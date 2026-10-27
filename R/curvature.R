#' Calculate Curvature Scale Space (CSS) Descriptor
#'
#' Computes the Curvature Scale Space representation of a 2D closed contour.
#' This involves progressively smoothing the contour with Gaussian kernels and
#' identifying curvature zero-crossings (inflection points) at each scale.
#' By default (\code{start_point = "auto"}), the contour is reordered so that
#' normalized arc length \code{0.0} starts at the anatomical rostrum (rightmost point),
#' matching the 4 AFORO zones used in \code{\link{calculate_distances}}.
#'
#' @param contour A matrix or data.frame with 'X' and 'Y' columns representing the contour.
#' @param n_points Integer, NULL, or \code{"max"}. Number of points to resample the contour to
#'                 for uniform arc length (default: 512, matching the AFORO standard).
#'                 If \code{NULL} or \code{"max"}, uses the full native resolution of the
#'                 extracted contour (\code{nrow(contour)} points).
#' @param sigma_max Numeric. Maximum smoothing scale (standard deviation of Gaussian).
#' @param sigma_step Numeric. Increment for the smoothing scale.
#' @param start_point Strategy or coordinates for selecting the contour starting point (\code{arc_length = 0.0}):
#'   \itemize{
#'     \item \code{"auto"} (default, or \code{NULL}): Automatically selects the rightmost point (\eqn{\max(X)}) as the rostrum.
#'     \item \code{"manual"}: Interactive mouse selection mode. Plots the contour with the automatic
#'           rostrum candidate (\code{o}) and prompts the user to click on the desired starting point
#'           with the mouse (snapping to the nearest contour vertex).
#'     \item Numeric coordinate pair \code{c(x, y)} (e.g. returned by \code{\link{select_starting_point}}),
#'           a 1x2 matrix, or an integer contour index (\code{1..nrow(contour)}).
#'   }
#' @return A list of class `css` containing:
#'   \itemize{
#'     \item \code{css_matrix}: A data.frame with columns `arc_length`, `sigma`, and `zone` for each zero-crossing.
#'     \item \code{n_points}: Number of points used in resampling.
#'     \item \code{contour_original}: The resampled contour (starting at the rostrum if aligned).
#'     \item \code{start_point_coords}: Numeric vector \code{c(X, Y)} of the starting point (\code{arc_length = 0}).
#'   }
#' @importFrom stats approx
#' @export
#' @examples
#' # Example with a simple shape
#' theta <- seq(0, 2 * pi, length.out = 100)
#' x <- cos(theta) + 0.1 * cos(3 * theta)
#' y <- sin(theta) + 0.1 * sin(3 * theta)
#' contour <- data.frame(X = x, Y = y)
#' css_res <- calculate_css(contour)
#' plot(css_res)
calculate_css <- function(contour, n_points = 512, sigma_max = 10, sigma_step = 0.5,
                          start_point = "auto") {
    if (!is.numeric(sigma_max) || length(sigma_max) != 1 || is.na(sigma_max) || sigma_max <= 0) {
        stop("sigma_max must be a positive number")
    }

    if (!is.numeric(sigma_step) || length(sigma_step) != 1 || is.na(sigma_step) || sigma_step <= 0) {
        stop("sigma_step must be a positive number")
    }

    # 0. Validate contour and align starting point (arc_length = 0.0) to rostrum
    resolved <- .resolve_contour_and_start(contour, start_point = start_point, remove_duplicate_end = TRUE)
    coo <- resolved$coords
    start_pt_coords <- resolved$start_point_coords

    if (is.null(n_points) || (is.character(n_points) && length(n_points) == 1 && tolower(n_points) == "max")) {
        n_points <- nrow(coo)
    } else if (!is.numeric(n_points) || length(n_points) != 1 || is.na(n_points) || n_points < 3) {
        stop("n_points must be an integer >= 3, NULL, or 'max'")
    } else {
        n_points <- as.integer(n_points)
    }

    # 1. Resample contour for uniform arc length (Modern Improvement)
    coo <- .coo_interpolate(coo, n_points)

    # Prepare for smoothing
    N <- nrow(coo)
    sigmas <- seq(0, sigma_max, by = sigma_step)
    css_points <- list()

    # Circular convolution using modular indexing
    circ_conv <- function(v, kernel) {
        k_len <- length(kernel)
        n_pad <- floor(k_len / 2)
        offsets <- (-n_pad):n_pad
        res <- numeric(N)
        for (i in 1:N) {
            idx <- ((i - 1 + offsets) %% N) + 1
            res[i] <- sum(v[idx] * rev(kernel))
        }
        res
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
        curv_shifted <- c(curv[-1], curv[1])
        zeros <- which(curv * curv_shifted < 0)

        if (length(zeros) > 0) {
            # Interpolate for better precision
            zero_positions <- zeros + abs(curv[zeros]) / (abs(curv[zeros]) + abs(curv_shifted[zeros]) + 1e-12)
            # Normalize arc length to [0, 1]
            zero_positions <- pmin(1, pmax(0, (zero_positions - 1) / N))
            zones_vec <- as.character(cut(
                zero_positions,
                breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                include.lowest = TRUE,
                labels = c("Zone 1", "Zone 2", "Zone 3", "Zone 4")
            ))
            css_points[[as.character(sigma)]] <- data.frame(
                arc_length = zero_positions,
                sigma = sigma,
                zone = zones_vec,
                stringsAsFactors = FALSE
            )
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
        contour_original = coo,
        start_point_coords = start_pt_coords
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

    num <- (dx * dyy - dy * dxx)
    den <- (dx^2 + dy^2)^1.5
    num / ifelse(den == 0, 1e-10, den)
}

#' Plot Curvature Scale Space Image
#'
#' Visualizes the CSS image, showing inflection points across different scales.
#' By default (\code{show_zones = TRUE}), points are colored according to the 4
#' AFORO contour sectors starting from the rostrum (\code{0.0}): Zone 1 (0-0.25),
#' Zone 2 (0.25-0.50), Zone 3 (0.50-0.75), and Zone 4 (0.75-1.0).
#'
#' @param x An object of class `css`.
#' @param show_zones Logical; if TRUE (default), colors points by the 4 AFORO
#'                   contour zones (unless \code{col} is explicitly supplied in \code{...})
#'                   and draws vertical zone boundary guides at 0.25, 0.50, and 0.75.
#' @param ... Additional arguments passed to plot.
#' @return An invisible NULL, called for its plotting side effect.
#' @importFrom utils modifyList
#' @importFrom graphics plot abline text mtext
#' @export
#' @examples
#' \dontrun{
#' css_res <- calculate_css(contour)
#' plot(css_res)
#' plot(css_res, main = "Custom Title", col = "blue")
#' }
plot.css <- function(x, show_zones = TRUE, ...) {
    if (!inherits(x, "css") || !is.list(x)) {
        stop("x must be an object of class 'css'")
    }

    user_args <- list(...)

    if (is.null(x$css_matrix) || nrow(x$css_matrix) == 0) {
        empty_defaults <- list(
            x = 0, y = 0, type = "n",
            main = "No inflection points found",
            xlab = "Arc length (0 = Rostrum)",
            ylab = "Sigma"
        )
        empty_args <- utils::modifyList(empty_defaults, user_args)
        do.call(graphics::plot, empty_args)
        return(invisible(NULL))
    }

    max_sigma <- max(x$css_matrix$sigma)
    ylim_val <- if (is.finite(max_sigma) && max_sigma > 0) c(0, max_sigma * 1.15) else c(0, 1)

    zone_palette <- c(
        "Zone 1" = "#00B8D4",
        "Zone 2" = "#E6B800",
        "Zone 3" = "#FF007F",
        "Zone 4" = "#00C853"
    )
    default_cols <- if (isTRUE(show_zones)) {
        z_idx <- as.character(cut(
            x$css_matrix$arc_length,
            breaks = c(0, 0.25, 0.50, 0.75, 1.0),
            include.lowest = TRUE,
            labels = c("Zone 1", "Zone 2", "Zone 3", "Zone 4")
        ))
        unname(zone_palette[z_idx])
    } else {
        "black"
    }

    defaults <- list(
        x = x$css_matrix$arc_length,
        y = x$css_matrix$sigma,
        pch = 20,
        cex = 0.65,
        col = default_cols,
        xlab = "Normalized Arc Length (0.0 = Rostrum)",
        ylab = "Scale (Sigma)",
        main = "Curvature Scale Space Image",
        ylim = ylim_val,
        xlim = c(0, 1)
    )

    plot_args <- utils::modifyList(defaults, user_args)
    do.call(graphics::plot, plot_args)

    if (isTRUE(show_zones)) {
        graphics::abline(v = c(0.25, 0.50, 0.75), lty = 3, col = "gray65", lwd = 1)
        y_top <- plot_args$ylim[2] * 0.96
        z_centers <- c(0.125, 0.375, 0.625, 0.875)
        z_names <- c("Zone 1", "Zone 2", "Zone 3", "Zone 4")
        for (k in 1:4) {
            graphics::text(z_centers[k], y_top, labels = z_names[k],
                           col = zone_palette[k], font = 2, cex = 0.75)
        }
    }

    invisible(NULL)
}
