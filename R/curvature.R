#' Calculate Curvature Scale Space (CSS) Descriptor
#'
#' Computes the Curvature Scale Space representation of a 2D closed contour.
#' This involves progressively smoothing the contour with Gaussian kernels and
#' identifying curvature zero-crossings (inflection points) at each scale.
#'
#' @param contour A matrix or data.frame with 'X' and 'Y' columns representing the contour.
#' @param n_points Integer. Number of points to resample the contour to (for uniform arc length).
#' @param sigma_max Numeric. Maximum smoothing scale (standard deviation of Gaussian).
#' @param sigma_step Numeric. Increment for the smoothing scale.
#' @return A list of class `css` containing:
#'   \itemize{
#'     \item \code{css_matrix}: A data.frame with columns `arc_length` and `sigma` for each zero-crossing.
#'     \item \code{n_points}: Number of points used in resampling.
#'     \item \code{contour_original}: The original contour.
#'   }
#' @export
#' @examples
#' # Example with a simple shape
#' theta <- seq(0, 2 * pi, length.out = 100)
#' x <- cos(theta) + 0.1 * cos(3 * theta)
#' y <- sin(theta) + 0.1 * sin(3 * theta)
#' contour <- data.frame(X = x, Y = y)
#' css_res <- calculate_css(contour)
#' # plot(css_res)
calculate_css <- function(contour, n_points = 200, sigma_max = 10, sigma_step = 0.5) {
    # Input validation
    if (is.null(contour)) stop("contour cannot be NULL")

    # Ensure matrix format
    if (is.data.frame(contour)) {
        coo <- as.matrix(contour[, c("X", "Y")])
    } else {
        coo <- as.matrix(contour)
    }

    # Remove duplicated last point if present (common in closed contours)
    if (nrow(coo) > 1 && all(coo[1, ] == coo[nrow(coo), ])) {
        coo <- coo[-nrow(coo), , drop = FALSE]
    }

    # 1. Resample contour for uniform arc length (Modern Improvement)
    # Ensure we have enough points for resampling
    if (requireNamespace("Momocs", quietly = TRUE)) {
        if (nrow(coo) < n_points) {
            coo <- Momocs::coo_interpolate(coo, n_points)
        } else {
            coo <- Momocs::coo_sample(coo, n_points)
        }
    } else {
        # Simple linear interpolation for resampling
        dists <- c(0, cumsum(sqrt(rowSums(diff(rbind(coo, coo[1, ]))^2))))
        total_len <- dists[length(dists)]
        target_dists <- seq(0, total_len, length.out = n_points + 1)[1:n_points]
        coo_x <- approx(dists, c(coo[, 1], coo[1, 1]), xout = target_dists)$y
        coo_y <- approx(dists, c(coo[, 2], coo[1, 2]), xout = target_dists)$y
        coo <- cbind(coo_x, coo_y)
    }

    # Prepare for smoothing
    N <- nrow(coo)
    sigmas <- seq(0, sigma_max, by = sigma_step)
    css_points <- list()

    # Manual convolution helper for circular padding
    circ_conv <- function(v, kernel) {
        k_len <- length(kernel)
        n_pad <- floor(k_len / 2)
        # Circular padding
        v_pad <- c(v[(N - n_pad + 1):N], v, v[1:n_pad])
        res <- rep(0, N)
        for (i in 1:N) {
            res[i] <- sum(v_pad[i:(i + k_len - 1)] * rev(kernel))
        }
        res
    }

    # Gaussian kernel function
    get_gaussian <- function(sigma) {
        if (sigma == 0) {
            return(1)
        }
        limit <- ceiling(4 * sigma)
        x <- seq(-limit, limit)
        g <- exp(-x^2 / (2 * sigma^2))
        g / sum(g)
    }

    # Gaussian derivative kernel
    get_gaussian_deriv <- function(sigma, order = 1) {
        limit <- ceiling(4 * sigma)
        x <- seq(-limit, limit)
        g <- exp(-x^2 / (2 * sigma^2))
        if (order == 1) {
            gd <- -(x / sigma^2) * g
        } else if (order == 2) {
            gd <- (x^2 / sigma^4 - 1 / sigma^2) * g
        }
        # For derivatives, we don't normalize to sum=1, but we keep scale consistent
        # The zero crossings are invariant to global scale anyway
        gd
    }

    for (sigma in sigmas) {
        if (sigma == 0) {
            # Compute curvature for raw data
            curv <- .compute_curvature_raw(coo)
        } else {
            # 2. Smooth and compute derivatives using Gaussian filters
            gd1 <- get_gaussian_deriv(sigma, 1)
            gd2 <- get_gaussian_deriv(sigma, 2)

            X_p <- circ_conv(coo[, 1], gd1)
            Y_p <- circ_conv(coo[, 2], gd1)
            X_pp <- circ_conv(coo[, 1], gd2)
            Y_pp <- circ_conv(coo[, 2], gd2)

            # 3. Curvature: K = (X'Y'' - Y'X'') / (X'^2 + Y'^2)^1.5
            num <- (X_p * Y_pp - Y_p * X_pp)
            den <- (X_p^2 + Y_p^2)^1.5
            curv <- num / ifelse(den == 0, 1e-10, den)
        }

        # 4. Find zero crossings (inflection points)
        # Use a small epsilon to avoid numerical noise
        curv_shifted <- c(curv[-1], curv[1])
        # A true zero crossing must cross a threshold or change sign significantly
        zeros <- which(curv * curv_shifted < 0)

        if (length(zeros) > 0) {
            # Interpolate for better precision
            zero_positions <- zeros + abs(curv[zeros]) / (abs(curv[zeros]) + abs(curv_shifted[zeros]) + 1e-12)
            # Normalize arc length to [0, 1]
            zero_positions <- (zero_positions - 1) / N
            css_points[[as.character(sigma)]] <- data.frame(arc_length = zero_positions, sigma = sigma)
        } else if (sigma > 1) {
            # If no zero crossings and we've smoothed enough, we can stop early
            break
        }
    }

    css_data <- do.call(rbind, css_points)
    rownames(css_data) <- NULL

    result <- list(
        css_matrix = css_data,
        n_points = n_points,
        contour_original = coo
    )
    class(result) <- "css"
    return(result)
}

#' Internal helper for raw curvature
#' @keywords internal
.compute_curvature_raw <- function(coo) {
    N <- nrow(coo)
    # Use simple finite differences for sigma=0
    i_prev <- c(N, 1:(N - 1))
    i_next <- c(2:N, 1)

    dx <- (coo[i_next, 1] - coo[i_prev, 1]) / 2
    dy <- (coo[i_next, 2] - coo[i_prev, 2]) / 2
    dxx <- coo[i_next, 1] - 2 * coo[, 1] + coo[i_prev, 1]
    dyy <- coo[i_next, 2] - 2 * coo[, 2] + coo[i_prev, 2]

    (dx * dyy - dy * dxx) / (dx^2 + dy^2)^1.5
}

#' Plot Curvature Scale Space Image
#'
#' Visualizes the CSS image, showing inflection points across different scales.
#'
#' @param x An object of class `css`.
#' @param ... Additional arguments passed to plot.
#' @export
plot.css <- function(x, ...) {
    if (is.null(x$css_matrix) || nrow(x$css_matrix) == 0) {
        plot(0, 0, type = "n", main = "No inflection points found", xlab = "Arc length", ylab = "Sigma")
        return(invisible(NULL))
    }

    plot(x$css_matrix$arc_length, x$css_matrix$sigma,
        pch = 20, cex = 0.5,
        xlab = "Arc Length (normalized)",
        ylab = "Scale (Sigma)",
        main = "Curvature Scale Space Image",
        ylim = c(0, max(x$css_matrix$sigma) * 1.1),
        xlim = c(0, 1),
        ...
    )
}
