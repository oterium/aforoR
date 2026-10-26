#' Calculate Adaptive Tri-Class Threshold
#'
#' Computes an adaptive threshold using a three-class Otsu variance maximization,
#' separating background noise from the translucent/shaded margins of the otolith.
#'
#' @param img An \code{Image} object or numeric matrix/vector of intensities in [0, 1].
#' @return A numeric scalar threshold value.
#' @keywords internal
calculate_adaptive_threshold <- function(img) {
  if (is.null(img)) return(0.2)

  if (inherits(img, "Image") && EBImage::colorMode(img) != EBImage::Grayscale) {
    img <- EBImage::channel(img, "grey")
  } else if (is.array(img) && length(dim(img)) >= 3) {
    if (dim(img)[3] >= 3) {
      img <- 0.2989 * img[, , 1] + 0.5870 * img[, , 2] + 0.1140 * img[, , 3]
    } else {
      img <- img[, , 1]
    }
  }

  if (is.numeric(img) && any(!is.na(img)) && max(img, na.rm = TRUE) > 1 && max(img, na.rm = TRUE) <= 255) {
    img <- img / 255
  }

  vals <- as.vector(img)
  vals <- vals[!is.na(vals) & vals >= 0 & vals <= 1]
  if (length(vals) < 100) return(0.2)

  th_std <- tryCatch(
    {
      t_raw <- EBImage::otsu(if (is.array(img)) img else as.matrix(img))
      if (length(t_raw) > 1) t_raw[1] else t_raw
    },
    error = function(e) 0.3
  )
  if (!is.numeric(th_std) || length(th_std) != 1 || is.na(th_std)) th_std <- 0.3

  h <- graphics::hist(vals, breaks = seq(0, 1, length.out = 257), plot = FALSE)
  counts <- h$counts
  total <- sum(counts)
  if (total == 0) return(th_std)

  p <- counts / total
  omega <- cumsum(p)
  mu <- cumsum(p * (1:256))
  mu_t <- mu[256]

  best_var <- -1
  best_t1 <- 1

  # Search range for T1 (background separator) and T2 (dense core separator)
  for (t1 in 2:80) {
    w0 <- omega[t1]
    if (w0 == 0) next
    u0 <- mu[t1] / w0

    for (t2 in (t1 + 5):200) {
      w1 <- omega[t2] - w0
      w2 <- 1 - omega[t2]
      if (w1 <= 0 || w2 <= 0) next

      u1 <- (mu[t2] - mu[t1]) / w1
      u2 <- (mu_t - mu[t2]) / w2

      sigma_b <- w0 * (u0 - mu_t)^2 + w1 * (u1 - mu_t)^2 + w2 * (u2 - mu_t)^2
      if (sigma_b > best_var) {
        best_var <- sigma_b
        best_t1 <- t1
      }
    }
  }

  if (best_var < 0) {
    return(th_std * 0.75)
  }

  th_adaptive <- h$mids[best_t1]

  # Safeguard: If adaptive threshold is unrealistically low (< 0.05) or
  # higher than standard Otsu, fall back to 0.75 * standard Otsu
  if (is.na(th_adaptive) || th_adaptive < 0.05 || th_adaptive >= th_std) {
    th_adaptive <- th_std * 0.75
  }

  return(th_adaptive)
}

#' Preprocess Image for Contour Analysis
#'
#' Prepares an image for shape analysis by applying grayscale conversion,
#' Gaussian smoothing to reduce noise, and binarization using the Otsu method
#' (default) or a fixed threshold.
#'
#' @param image_path A string specifying the path to the image file (e.g., .jpg,
#'                   .png, .tiff) OR an EBImage object.
#' @param threshold A numeric value (0-1) for binarization. If NULL, Otsu's method
#'                  is automatically applied.
#' @param blur_size Size of the Gaussian brush for smoothing in pixels. Can be an
#'                  odd integer (e.g., 15, 31), or \code{"auto"} / \code{NULL} (default)
#'                  to dynamically adapt the kernel size to ~1% of the image's minimum
#'                  dimension (\eqn{\min(\text{width}, \text{height}) \times 0.01}).
#'                  This preserves fine crenulations in small images while effectively
#'                  suppressing noise in high-resolution or 4K images.
#' @param blur_sigma Standard deviation of the Gaussian filter brush. Can be a numeric
#'                   value (e.g., 9), or \code{"auto"} / \code{NULL} (default) to
#'                   automatically scale with \eqn{\sigma \approx \text{blur\_size} / 3.5}.
#' @param downsample_threshold Logical; if TRUE, calculates Otsu threshold on a
#'                             downsampled thumbnail for faster preprocessing of
#'                             very large images (default: FALSE).
#' @param threshold_method Method for automatic binarization: "adaptive" (three-class Otsu
#'                         separating background, edge, and core), "otsu" (standard Otsu),
#'                         or "fixed". Default is "adaptive".
#' @param threshold_factor Multiplier applied to the calculated threshold (default: 1.0).
#' @return A \code{pixset} (EBImage) representing the binarized image.
#' @export
#' @examples
#' # Use the sample otolith image included in the package
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' processed_img <- preprocess_image(image_path)
preprocess_image <- function(image_path, threshold = NULL, blur_size = "auto",
                             blur_sigma = "auto", downsample_threshold = FALSE,
                             threshold_method = c("adaptive", "otsu", "fixed"),
                             threshold_factor = 1.0) {
  is_image <- inherits(image_path, "Image")

  if (!is_image && (!is.character(image_path) || length(image_path) != 1)) {
    stop("image_path must be a single character string or an Image object")
  }

  if (!is_image && !file.exists(image_path)) {
    stop(paste("Image file does not exist:", image_path))
  }

  if (!is.null(threshold) &&
    (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) || is.infinite(threshold) || threshold < 0 || threshold > 1)) {
    stop("threshold must be a single numeric value between 0 and 1, or NULL")
  }

  tryCatch(
    {
      if (is_image) {
        foto <- image_path
      } else {
        foto <- EBImage::readImage(image_path)
      }

      foto2 <- if (EBImage::colorMode(foto) != EBImage::Grayscale) {
        EBImage::channel(foto, "grey")
      } else {
        foto
      }

      # Compute adaptive blur_size if "auto" or NULL based on image resolution
      dims <- dim(foto2)
      min_dim <- min(dims[1], dims[2])

      if (!is.null(blur_size) && !identical(blur_size, "auto")) {
        if (!is.numeric(blur_size) || length(blur_size) != 1 || is.na(blur_size) || is.infinite(blur_size) || blur_size <= 0) {
          stop("blur_size must be a positive numeric value, 'auto', or NULL")
        }
      }

      if (!is.null(blur_sigma) && !identical(blur_sigma, "auto")) {
        if (!is.numeric(blur_sigma) || length(blur_sigma) != 1 || is.na(blur_sigma) || is.infinite(blur_sigma) || blur_sigma <= 0) {
          stop("blur_sigma must be a positive numeric value, 'auto', or NULL")
        }
      }

      if (is.null(blur_size) || identical(blur_size, "auto")) {
        calc_size <- as.integer(round(min_dim * 0.01))
        if (calc_size %% 2 == 0) calc_size <- calc_size + 1L
        blur_size <- max(5L, calc_size)
      } else {
        blur_size <- as.integer(round(blur_size))
        if (blur_size %% 2 == 0) blur_size <- blur_size + 1L
        blur_size <- max(3L, blur_size)
      }

      if (is.null(blur_sigma) || identical(blur_sigma, "auto")) {
        blur_sigma <- max(0.5, round(blur_size / 3.5, 2))
      } else {
        blur_sigma <- max(0.1, as.numeric(blur_sigma))
      }

      foto2 <- EBImage::filter2(
        foto2,
        EBImage::makeBrush(
          size = blur_size, shape = "gaussian",
          sigma = blur_sigma
        )
      )

      threshold_method <- match.arg(threshold_method)
      if (!is.null(threshold)) {
        threshold_method <- "fixed"
      }

      if (threshold_method == "fixed") {
        foto2 <- foto2 > threshold
      } else {
        if (isTRUE(downsample_threshold)) {
          w_small <- min(dim(foto2)[1], 800)
          foto2_th <- EBImage::resize(foto2, w = w_small)
        } else {
          foto2_th <- foto2
        }

        if (threshold_method == "adaptive") {
          th <- calculate_adaptive_threshold(foto2_th)
        } else {
          th <- EBImage::otsu(foto2_th)
        }

        if (!is.null(threshold_factor) && is.numeric(threshold_factor) && threshold_factor > 0) {
          th <- th * threshold_factor
        }

        foto2 <- foto2 > th
      }

      return(foto2)
    },
    error = function(e) {
      stop(paste("Error processing image:", e$message))
    }
  )
}

#' Extract Main Contour from Binarized Image
#'
#' Detects all contours in a binary image and identifies the largest one that
#' meets specific size criteria (points >= 200 and area > 5000), typically
#' representing the fish otolith.
#'
#' @param binary_image A binarized image produced by \code{preprocess_image()}.
#' @param min_points Minimum number of points required for a valid contour
#'                   (default: 200).
#' @param min_area Minimum area required for a valid contour (default: 5000).
#' @return A matrix of (X, Y) coordinates of the main contour, or NULL if no
#'         suitable contour is found.
#' @export
#' @examples
#' # Process sample image and extract its contour
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' binary_img <- preprocess_image(image_path)
#' contour <- extract_contour(binary_img)
#' head(contour)
extract_contour <- function(binary_image, min_points = 200, min_area = 5000) {
  if (is.null(binary_image)) {
    stop("binary_image cannot be NULL")
  }

  if (!is.numeric(min_points) || length(min_points) != 1 || is.na(min_points) ||
      is.infinite(min_points) || min_points < 3 || min_points != as.integer(min_points)) {
    stop("min_points must be an integer >= 3")
  }
  min_points <- as.integer(min_points)

  if (!is.numeric(min_area) || length(min_area) != 1 || is.na(min_area) || min_area < 0) {
    stop("min_area must be a numeric value >= 0")
  }

  tryCatch(
    {
      cont <- EBImage::ocontour(EBImage::bwlabel(binary_image))
      contorno <- NULL

      if (length(cont) == 0) {
        warning("No contours found in image")
        return(NULL)
      }

      max_area <- 0
      for (x in seq_along(cont)) {
        c_mat <- cont[[x]]
        if (is.null(c_mat) || !is.matrix(c_mat) || nrow(c_mat) < min_points) {
          next
        }
        current_area <- .coo_area(c_mat)
        if (current_area > min_area && current_area > max_area) {
          max_area <- current_area
          contorno <- c_mat
        }
      }

      if (is.null(contorno)) {
        warning(paste(
          "No contour found meeting size criteria (> ", min_points,
          " points and area > ", min_area, ")",
          sep = ""
        ))
      }

      return(contorno)
    },
    error = function(e) {
      stop(paste("Error extracting contour:", e$message))
    }
  )
}

#' Interactively Select Contour Starting Point (Rostrum)
#'
#' Select Starting Point Interactively
#'
#' Displays an otolith image with its extracted contour overlaid, prompting the user
#' to click on the anatomical rostrum or preferred starting point.
#'
#' @param image An \code{Image} object or file path.
#' @param contour A matrix or data.frame with (X, Y) contour coordinates.
#' @return A numeric vector of length 2 \code{c(x, y)} representing the clicked point coordinates,
#'         or NULL if cancelled.
#' @importFrom graphics lines title
#' @export
#' @examples
#' \dontrun{
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' pt <- select_starting_point(image_path)
#' }
select_starting_point <- function(image, contour = NULL) {
  if (is.character(image) && file.exists(image)) {
    img_obj <- EBImage::readImage(image)
  } else if (inherits(image, "Image")) {
    img_obj <- image
  } else {
    stop("image must be an Image object or a valid file path")
  }

  if (!interactive()) {
    message("Non-interactive session: using automatic rightmost point.")
    return(NULL)
  }

  plot(img_obj)
  if (!is.null(contour)) {
    lines(contour[, 1], contour[, 2], col = "#00E5FF", lwd = 2)
  }
  title(sub = "Click on the rostrum/starting point (or press Esc/right-click to use default)",
        col.sub = "red")

  pt <- tryCatch({
    graphics::locator(n = 1)
  }, error = function(e) NULL)

  if (is.null(pt) || length(pt$x) == 0) {
    message("No point clicked; using automatic rightmost point.")
    return(NULL)
  }

  return(c(x = pt$x[1], y = pt$y[1]))
}

#' Calculate Distance Measures from Contour
#'
#' Computes two types of distance measures from an otolith contour: polar
#' distances (radial distances from the centroid) and perimeter-based distances.
#' These measures are crucial for subsequent wavelet analysis.
#'
#' @param contour A matrix of coordinates from \code{extract_contour()}.
#' @param n_points Target number of points for resampling (default: 512).
#' @param start_point Optional numeric vector specifying the starting point (rostrum).
#'                    Can be either an integer index (1 to nrow(contour)) or coordinate
#'                    pair \code{c(x, y)}. If NULL (default), the automatic rightmost point
#'                    is used.
#' @return A list containing:
#'   \itemize{
#'     \item \code{polar}: Radii, coordinates, and normalized radii.
#'     \item \code{perimeter}: Distances, coordinates, and normalized distances.
#'     \item \code{reordered_coords}: The contour coordinates reordered from the
#'           starting point.
#'   }
#' @export
#' @examples
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' binary_img <- preprocess_image(image_path)
#' contour <- extract_contour(binary_img)
#' dists <- calculate_distances(contour, n_points = 512)
#' plot(dists$polar$normalized, type = "l")
calculate_distances <- function(contour, n_points = 512, start_point = NULL) {
  if (is.matrix(contour)) {
    if (ncol(contour) >= 2) {
      contour <- as.data.frame(contour[, 1:2, drop = FALSE])
      colnames(contour) <- c("X", "Y")
    } else {
      stop("contour must have at least 2 columns")
    }
  } else if (is.data.frame(contour)) {
    if (all(c("X", "Y") %in% colnames(contour))) {
      contour <- contour[, c("X", "Y")]
    } else if (all(c("x", "y") %in% colnames(contour))) {
      contour <- contour[, c("x", "y")]
      colnames(contour) <- c("X", "Y")
    } else if (ncol(contour) >= 2) {
      contour <- contour[, 1:2]
      colnames(contour) <- c("X", "Y")
    } else {
      stop("contour must have at least 2 columns")
    }
  } else {
    stop("contour must be a data.frame or matrix with X and Y columns")
  }

  if (nrow(contour) < 3) {
    stop("contour must contain at least 3 points")
  }

  if (!is.numeric(contour$X) || !is.numeric(contour$Y)) {
    stop("contour coordinates must be numeric")
  }

  if (any(is.na(contour$X)) || any(is.na(contour$Y))) {
    stop("contour coordinates cannot contain NA values")
  }

  if (any(is.infinite(contour$X)) || any(is.infinite(contour$Y))) {
    stop("contour coordinates cannot contain infinite values")
  }

  if (all(contour$X == contour$X[1]) && all(contour$Y == contour$Y[1])) {
    stop("contour is degenerate with identical coordinates")
  }

  if (!is.numeric(n_points) || length(n_points) != 1 || is.na(n_points) ||
      is.infinite(n_points) || n_points < 2 || n_points != as.integer(n_points)) {
    stop("n_points must be an integer >= 2")
  }

  # Determine starting point index
  if (!is.null(start_point)) {
    if (!is.numeric(start_point) || any(is.na(start_point)) || any(is.infinite(start_point))) {
      stop("start_point must be a valid numeric index or coordinate pair without NA/Inf")
    }
    if (length(start_point) == 1) {
      start_idx <- as.integer(start_point)
      if (start_idx < 1 || start_idx > nrow(contour)) {
        stop(paste("start_point index must be between 1 and", nrow(contour)))
      }
    } else if (length(start_point) == 2) {
      dists_sq <- (contour$X - start_point[1])^2 + (contour$Y - start_point[2])^2
      start_idx <- which.min(dists_sq)
    } else {
      stop("start_point must be an index (length 1) or coordinates c(x, y) (length 2)")
    }
  } else {
    # Find rightmost point automatically
    x_max_idx <- which.max(contour$X)
    x_max <- contour$X[x_max_idx]
    y_max <- contour$Y[x_max_idx]
    punto <- list(x = x_max, y = y_max)

    x_mean <- mean(contour$X)
    rightmost_candidates <- which(contour$X == punto$x & contour$X > x_mean)

    if (length(rightmost_candidates) == 0) {
      start_idx <- x_max_idx
    } else {
      start_idx <- rightmost_candidates[1]
    }
  }

  n_contour <- nrow(contour)
  if (start_idx == 1) {
    x <- contour$X
    y <- contour$Y
  } else {
    indices <- c(start_idx:n_contour, 1:(start_idx - 1))
    x <- contour$X[indices]
    y <- contour$Y[indices]
  }

  tryCatch(
    {
      dist <- regularradius(x, y, n = n_points)
      dist2 <- dper(x, y, n_points)
      max_r1 <- max(dist$radii)
      max_r2 <- max(dist2$dist)
      dist.norm <- if (max_r1 > 0) dist$radii / max_r1 else rep(0, length(dist$radii))
      dist.norm2 <- if (max_r2 > 0) dist2$dist / max_r2 else rep(0, length(dist2$dist))

      return(list(
        polar = list(
          radii = dist$radii, coord = dist$coord, normalized = dist.norm
        ),
        perimeter = list(
          dist = dist2$dist, coords = dist2$coords, normalized = dist.norm2
        ),
        reordered_coords = list(x = x, y = y)
      ))
    },
    error = function(e) {
      stop(paste("Error in distance calculation:", e$message))
    }
  )
}

#' Calculate Wavelet Analysis for Distance Data
#'
#' Performs multi-scale wavelet decomposition on normalized distance data
#' (polar or perimeter). This is useful for capturing shape variation at
#' different frequencies.
#'
#' @param distances A list from \code{calculate_distances()} containing normalized data.
#' @param n_scales Number of wavelet scales to extract (default: 9).
#' @param detail Logical; if FALSE (default), returns wavelet detail
#'               coefficients (high-to-medium frequency differences between
#'               scales, centered around zero). If TRUE, returns smoothed
#'               approximation coefficients.
#' @return A list with wavelet coefficients for 'polar' and 'perimeter' components.
#' @export
#' @examples
#' # Example using calculated distances from an otolith contour:
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' binary_img <- preprocess_image(image_path)
#' contour <- extract_contour(binary_img)
#' distances <- calculate_distances(contour)
#' wavelets <- calculate_wavelets_analysis(distances, n_scales = 5)
#' str(wavelets)
#' @seealso \code{\link{fwaveletspl_3}}
calculate_wavelets_analysis <- function(distances, n_scales = 9,
                                        detail = FALSE) {
  if (!is.list(distances)) {
    stop("distances must be a list")
  }

  if (!all(c("polar", "perimeter") %in% names(distances))) {
    stop("distances must contain 'polar' and 'perimeter' components")
  }

  if (!is.numeric(n_scales) || length(n_scales) != 1 || is.na(n_scales) || is.infinite(n_scales) || n_scales <= 0 || n_scales != as.integer(n_scales)) {
    stop("n_scales must be a positive integer")
  }

  if (!is.logical(detail) || length(detail) != 1 || is.na(detail)) {
    stop("detail must be a logical value (TRUE or FALSE)")
  }

  if (!("normalized" %in% names(distances$polar)) ||
    !("normalized" %in% names(distances$perimeter))) {
    stop("distances must contain normalized components")
  }

  tryCatch(
    {
      wavelets_polar <- fwaveletspl_3(
        distances$polar$normalized, n_scales,
        det = detail
      )
      wavelets_perimeter <- fwaveletspl_3(
        distances$perimeter$normalized, n_scales,
        det = detail
      )

      return(list(
        polar = wavelets_polar,
        perimeter = wavelets_perimeter
      ))
    },
    error = function(e) {
      stop(paste("Error in wavelet analysis:", e$message))
    }
  )
}

#' Save Visualization Images for Analysis
#'
#' Creates and saves visualization images showing contour analysis and wavelet results.
#' Overlays the extracted contour, centroid, rostrum, and polar analysis sectors
#' directly onto the original photograph in high contrast colors.
#'
#' @param image An \code{Image} object (original or binary) OR a string specifying
#'              the image file path. If NULL, draws on coordinate space (useful for GPA).
#' @param distances The distance measures from \code{\link{calculate_distances}}.
#' @param wavelets The wavelet analysis from \code{\link{calculate_wavelets_analysis}}.
#' @param image_name The base name of the image file.
#' @param output_dir The directory to save visualization images.
#' @param wavelet_scale The wavelet scale to visualize (default: 5).
#' @param plot_width Width in pixels for diagnostic JPEG images (default: 1200).
#' @param plot_height Height in pixels for diagnostic JPEG images (default: 800).
#' @param plot_res Resolution in dpi for diagnostic JPEG images (default: 150).
#' @param binary_image Deprecated alias for \code{image} for backwards compatibility.
#' @param scale_info Optional list with scale bar metrics (min_x, max_x, min_y, max_y, bar_length_px, scale_mm, px_per_mm).
#' @return An invisible NULL, called for its side effect of saving diagnostic visualization images.
#' @importFrom grDevices dev.off jpeg adjustcolor rgb
#' @importFrom graphics lines points polygon legend
#' @export
#' @examples
#' \dontrun{
#' # Assuming image, distances, and wavelets are already computed:
#' save_visualization(image, distances, wavelets, "sample_image", tempdir())
#' }
save_visualization <- function(image = NULL, distances, wavelets, image_name,
                               output_dir, wavelet_scale = 5,
                               plot_width = 1200, plot_height = 800, plot_res = 150,
                               binary_image = NULL, scale_info = NULL) {
  if (is.null(image) && !is.null(binary_image)) {
    image <- binary_image
  }

  if (!is.list(distances)) {
    stop("distances must be a list")
  }

  if (!is.null(wavelets) && !is.list(wavelets)) {
    stop("wavelets must be a list or NULL")
  }

  if (!is.character(image_name) || length(image_name) != 1) {
    stop("image_name must be a single character string")
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }

  if (!is.numeric(wavelet_scale) || length(wavelet_scale) != 1 || is.na(wavelet_scale) || wavelet_scale <= 0) {
    stop("wavelet_scale must be a positive integer")
  }

  if (!is.numeric(plot_width) || length(plot_width) != 1 || is.na(plot_width) || plot_width <= 0) {
    stop("plot_width must be a positive number")
  }

  if (!is.numeric(plot_height) || length(plot_height) != 1 || is.na(plot_height) || plot_height <= 0) {
    stop("plot_height must be a positive number")
  }

  if (!is.numeric(plot_res) || length(plot_res) != 1 || is.na(plot_res) || plot_res <= 0) {
    stop("plot_res must be a positive number")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cur_dev <- NULL
  on.exit({
    if (!is.null(cur_dev) && cur_dev > 1) {
      try(grDevices::dev.off(cur_dev), silent = TRUE)
      cur_dev <- NULL
    }
  }, add = TRUE)

  tryCatch(
    {
      zonal_file <- file.path(
        output_dir,
        paste0("zonal_", tools::file_path_sans_ext(image_name), ".jpg")
      )
      jpeg(filename = zonal_file, width = plot_width, height = plot_height,
           res = plot_res, quality = 90)
      cur_dev <- grDevices::dev.cur()

      has_img <- FALSE
      if (!is.null(image)) {
        if (is.character(image) && file.exists(image)) {
          img_to_plot <- EBImage::readImage(image)
          plot(img_to_plot)
          has_img <- TRUE
        } else if (inherits(image, "Image")) {
          plot(image)
          has_img <- TRUE
        }
      }

      cx <- mean(distances$reordered_coords$x)
      cy <- mean(distances$reordered_coords$y)

      if (!has_img) {
        xlims <- range(distances$polar$coord[, 1] + cx)
        ylims <- range(distances$polar$coord[, 2] + cy)
        xpad <- if (is.na(diff(xlims)) || diff(xlims) == 0) 1 else diff(xlims) * 0.1
        ypad <- if (is.na(diff(ylims)) || diff(ylims) == 0) 1 else diff(ylims) * 0.1
        plot(NA, xlim = c(xlims[1] - xpad, xlims[2] + xpad),
             ylim = c(ylims[1] - ypad, ylims[2] + ypad),
             xlab = "X", ylab = "Y", main = paste("Zonal Analysis -", image_name))
      }

      points(cx, cy, pch = 3, col = "red", lwd = 2.5, cex = 1.6)
      points(distances$reordered_coords$x[1], distances$reordered_coords$y[1],
             pch = 16, col = "red", cex = 1.6)

      n_points <- length(distances$polar$radii)
      b <- round(seq(0, n_points, length.out = 5))
      zones <- lapply(1:4, function(idx) if (b[idx + 1] > b[idx]) (b[idx] + 1L):b[idx + 1] else integer(0))
      colors <- c("#00E5FF", "#FFEA00", "#FF007F", "#00E676")

      for (i in seq_along(zones)) {
        if (length(zones[[i]]) > 0) {
          lines(
            distances$polar$coord[zones[[i]], 1] + cx,
            distances$polar$coord[zones[[i]], 2] + cy,
            col = colors[i], lwd = 2.5
          )
        }
      }

      if (!is.null(scale_info) && is.list(scale_info) && !is.null(scale_info$min_x)) {
        pad <- 15
        graphics::rect(scale_info$min_x - pad, scale_info$min_y - pad,
                       scale_info$max_x + pad, scale_info$max_y + pad,
                       border = "#00FF66", lwd = 2.5)
        scale_label <- sprintf("Scale: %g mm (%.0f px | %.1f px/mm)",
                               scale_info$scale_mm, scale_info$bar_length_px, scale_info$px_per_mm)
        img_w <- if (has_img) dim(image)[1] else max(scale_info$max_x + 100, 1000)
        if (mean(c(scale_info$min_x, scale_info$max_x)) > img_w * 0.5) {
          t_x <- min(img_w - 30, scale_info$max_x + pad)
          a_x <- 1
        } else {
          t_x <- max(30, scale_info$min_x - pad)
          a_x <- 0
        }
        t_y <- max(35, scale_info$min_y - pad - 20)
        graphics::text(t_x, t_y, labels = scale_label, adj = c(a_x, 0.5),
                       col = "#00FF66", font = 2, cex = 0.80)
      }

      leg_labels <- c("Zone 1", "Zone 2", "Zone 3", "Zone 4", "Centroid (+)", "Rostrum (o)")
      leg_cols   <- c(colors, "red", "red")
      leg_pch    <- c(rep(NA, 4), 3, 16)
      leg_lty    <- c(rep(1, 4), NA, NA)

      if (!is.null(scale_info) && is.list(scale_info) && !is.null(scale_info$min_x)) {
        leg_labels <- c(leg_labels, sprintf("Scale bar (%g mm)", scale_info$scale_mm))
        leg_cols   <- c(leg_cols, "#00FF66")
        leg_pch    <- c(leg_pch, NA)
        leg_lty    <- c(leg_lty, 1)
      }

      legend("topleft",
             legend = leg_labels,
             col = leg_cols, pch = leg_pch,
             lty = leg_lty, lwd = 2, bg = rgb(0, 0, 0, 0.65), text.col = "white",
             cex = 0.80, pt.cex = 1.1, inset = c(0.02, 0.03))
      if (!is.null(cur_dev) && cur_dev > 1) {
        try(grDevices::dev.off(cur_dev), silent = TRUE)
        cur_dev <- NULL
      }

      if (!is.null(wavelets) && !is.null(wavelets$polar) && nrow(wavelets$polar) > 0) {
        scale_to_plot <- min(as.integer(wavelet_scale), nrow(wavelets$polar))
        wav_file <- file.path(
          output_dir,
          paste0("wavelet_", tools::file_path_sans_ext(image_name), ".jpg")
        )
        jpeg(filename = wav_file, width = 800, height = 500, res = 100, quality = 85)
        cur_dev <- grDevices::dev.cur()
        plot(wavelets$polar[scale_to_plot, ],
          main = paste("Wavelet Scale", scale_to_plot, "-", image_name),
          xlab = "Points", ylab = "Coefficient", type = "l", col = "#00E5FF", lwd = 2
        )

        for (i in seq_along(zones)) {
          if (length(zones[[i]]) > 0) {
            polygon(
              x = c(min(zones[[i]]), zones[[i]], max(zones[[i]])),
              y = c(
                min(wavelets$polar[scale_to_plot, ]),
                wavelets$polar[scale_to_plot, zones[[i]]],
                min(wavelets$polar[scale_to_plot, ])
              ),
              col = adjustcolor(colors[i], alpha.f = 0.5),
              border = colors[i]
            )
          }
        }
        if (!is.null(cur_dev) && cur_dev > 1) {
          try(grDevices::dev.off(cur_dev), silent = TRUE)
          cur_dev <- NULL
        }
      }
    },
    error = function(e) {
      stop(paste("Error creating visualization:", e$message))
    }
  )
}

#' Save Analysis Results to CSV Files
#'
#' Exports analysis results including distances, wavelets, and coordinates to CSV files.
#'
#' @param analysis_results A list containing all analysis results for multiple images.
#' @param output_dir The directory to save CSV files.
#' @param result_type A string indicating the type of results ("polar" or "perimeter").
#' @return An invisible NULL, called for its side effect of writing analysis CSV files.
#' @importFrom utils write.table
#' @export
#' @examples
#' \dontrun{
#' # Assuming analysis_results is a list with output from calculate_distances etc.
#' save_analysis_results(analysis_results, tempdir(), "polar")
#' }
save_analysis_results <- function(analysis_results, output_dir,
                                  result_type = "polar") {
  if (!is.list(analysis_results)) {
    stop("analysis_results must be a list")
  }

  if (length(analysis_results) == 0) {
    warning("analysis_results is empty, no files will be saved")
    return(invisible(NULL))
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }

  if (!result_type %in% c("polar", "perimeter")) {
    stop("result_type must be either 'polar' or 'perimeter'")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  tryCatch(
    {
      analysis_results <- analysis_results[!sapply(analysis_results, is.null)]
      if (length(analysis_results) == 0) {
        warning("No valid analysis results to save")
        return(invisible(NULL))
      }

      image_names <- sapply(analysis_results, function(x) x$nombre)

      distances <- do.call(
        rbind, lapply(analysis_results, function(x) x$Distancia)
      )
      distances_norm <- do.call(
        rbind, lapply(analysis_results, function(x) x$Distancia_Norm)
      )
      coordinates <- do.call(
        rbind, lapply(analysis_results, function(x) x$coords)
      )

      write_analysis_csv(
        distances, file.path(output_dir, "DistanciaEN.csv"),
        image_names
      )
      write_analysis_csv(
        distances_norm,
        file.path(output_dir, "Distancia_NormEN.csv"),
        image_names
      )

      wavelet_names <- grep("Wavelet_", names(analysis_results[[1]]),
        value = TRUE
      )
      for (wavelet_name in wavelet_names) {
        wavelet_data <- do.call(
          rbind,
          lapply(analysis_results, function(x) x[[wavelet_name]])
        )
        if (!is.null(wavelet_data)) {
          write_analysis_csv(
            wavelet_data,
            file.path(
              output_dir,
              paste0(wavelet_name, "EN.csv")
            ),
            image_names
          )
        }
      }

      write.table(coordinates,
        file = file.path(output_dir, "Coords.csv"),
        dec = ".", sep = ";", col.names = FALSE
      )

      if ("elliptic" %in% names(analysis_results[[1]])) {
        elliptic_data <- do.call(
          rbind,
          lapply(analysis_results, function(x) x$elliptic)
        )
        n_harmonics <- ncol(elliptic_data) / 4
        colnames(elliptic_data) <- c(
          paste("A", 1:n_harmonics, sep = ""),
          paste("B", 1:n_harmonics, sep = ""),
          paste("C", 1:n_harmonics, sep = ""),
          paste("D", 1:n_harmonics, sep = "")
        )
        rownames(elliptic_data) <- image_names
        write.table(elliptic_data,
          file = file.path(output_dir, "EllipticCoeEN.csv"),
          dec = ".", sep = ";", col.names = TRUE
        )
      }
    },
    error = function(e) {
      stop(paste("Error saving analysis results:", e$message))
    }
  )
}

write_analysis_csv <- function(data, filename, row_names) {
  if (is.null(data)) stop("data cannot be NULL")
  if (!is.character(filename) || length(filename) != 1) stop("filename must be a single character string")
  if (!is.character(row_names)) stop("row_names must be a character vector")

  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }

  if (nrow(data) != length(row_names)) {
    stop("Number of rows in data must match length of row_names")
  }

  rownames(data) <- row_names
  write.table(data, file = filename, dec = ".", sep = ";", col.names = FALSE)
}

.setup_directories_and_files <- function(folder, coordinates = "both") {
  polar_dir <- file.path(folder, "Polar")
  cartesian_dir <- file.path(folder, "Cartesian")

  if (coordinates %in% c("both", "polar")) {
    if (!dir.exists(polar_dir)) dir.create(polar_dir, recursive = TRUE)
  }
  if (coordinates %in% c("both", "cartesian")) {
    if (!dir.exists(cartesian_dir)) dir.create(cartesian_dir, recursive = TRUE)
  }

  imag <- list.files(folder, pattern = "(?i)\\.(jpe?g|png|tiff?)$", full.names = TRUE)
  if (length(imag) == 0) {
    warning(paste("No supported image files (.jpg, .jpeg, .png, .tif, .tiff) found in folder:", folder))
    return(NULL)
  }

  list(polar_dir = polar_dir, cartesian_dir = cartesian_dir, imag = imag)
}

.process_single_image_worker <- function(image_file, threshold, blur_size, blur_sigma,
                                        detect_scale, pixels_per_mm, scale_mm,
                                        min_points, min_area, n_points, n_harmonics,
                                        wavelets, ef, testing, save,
                                        polar_dir, cartesian_dir,
                                        plot_width, plot_height, plot_res,
                                        downsample_threshold,
                                        coordinates = "both",
                                        start_point = NULL,
                                        threshold_method = "adaptive",
                                        threshold_factor = 1.0) {
  tryCatch({
    img_raw <- EBImage::readImage(image_file)

    current_pixels_per_mm <- pixels_per_mm
    scale_info <- NULL
    if (detect_scale) {
      img_grey <- if (EBImage::colorMode(img_raw) != EBImage::Grayscale) {
        EBImage::channel(img_raw, "grey")
      } else {
        img_raw
      }
      thresh_scale <- EBImage::otsu(img_grey)
      binary_scale <- img_grey > thresh_scale
      detected_px_per_mm <- find_scale_bar_length(binary_scale, scale_mm = scale_mm)
      if (!is.null(detected_px_per_mm)) {
        current_pixels_per_mm <- detected_px_per_mm
        scale_info <- attr(detected_px_per_mm, "scale_info")
      } else {
        warning(paste("Scale bar not found for image", basename(image_file), ". Using unscaled values."))
      }
      rm(img_grey, binary_scale)
    }

    binary_image <- preprocess_image(
      img_raw, threshold = threshold, blur_size = blur_size,
      blur_sigma = blur_sigma, downsample_threshold = downsample_threshold,
      threshold_method = threshold_method, threshold_factor = threshold_factor
    )

    contorno_raw <- extract_contour(binary_image, min_points = min_points, min_area = min_area)
    rm(binary_image)

    if (is.null(contorno_raw)) {
      warning(paste("No suitable contour found for image:", basename(image_file)))
      rm(img_raw)
      return(NULL)
    }

    coef_e <- if (ef) .efourier(contorno_raw, n_harmonics, norm = TRUE) else NULL
    distances <- calculate_distances(contorno_raw, n_points = n_points, start_point = start_point)
    wavelet_results <- if (wavelets) calculate_wavelets_analysis(distances, n_scales = floor(log2(n_points)), detail = FALSE) else NULL

    morpho <- calculate_morphometrics(contorno_raw, current_pixels_per_mm)
    morpho$Image <- basename(image_file)

    do_polar <- coordinates %in% c("both", "polar")
    do_cart  <- coordinates %in% c("both", "cartesian")

    if (testing) {
      if (do_polar) {
        save_visualization(
          image = img_raw, distances = distances, wavelets = wavelet_results,
          image_name = basename(image_file), output_dir = polar_dir,
          plot_width = plot_width, plot_height = plot_height, plot_res = plot_res,
          scale_info = scale_info
        )
      }
      if (do_cart) {
        save_visualization_perimeter(
          image = img_raw, distances = distances, wavelets = wavelet_results,
          image_name = basename(image_file), output_dir = cartesian_dir,
          plot_width = plot_width, plot_height = plot_height, plot_res = plot_res,
          scale_info = scale_info
        )
      }
    }
    rm(img_raw)

    res_polar <- if (save && do_polar) create_result_structure(distances$polar, wavelet_results$polar, basename(image_file), coef_e, "polar") else NULL
    res_cart  <- if (save && do_cart)  create_result_structure(distances$perimeter, wavelet_results$perimeter, basename(image_file), coef_e, "perimeter") else NULL

    list(
      image = basename(image_file),
      polar = res_polar,
      perimeter = res_cart,
      morpho = morpho,
      contour = contorno_raw
    )
  }, error = function(e) {
    warning(paste("Error processing image", basename(image_file), ":", e$message))
    return(NULL)
  })
}

.apply_gpa_alignment <- function(contours_raw, valid_indices, procrustes, n_points) {
  contours_to_process <- contours_raw
  if (procrustes && length(valid_indices) > 1) {
    message("\nPhase 1.5: Applying Generalized Procrustes Analysis (GPA)...")
    valid_contours <- lapply(contours_raw[valid_indices], function(c) .coo_interpolate(c, n_points))
    out_aligned <- tryCatch(
      .procrustes_align(valid_contours),
      error = function(e) {
        lapply(valid_contours, .coo_align)
      }
    )
    if (!is.null(out_aligned) && length(out_aligned) == length(valid_indices)) {
      for (j in seq_along(valid_indices)) {
        contours_to_process[[valid_indices[j]]] <- out_aligned[[j]]
      }
    }
  } else if (procrustes) {
    warning("Not enough valid contours found to apply Procrustes alignment.")
  }
  contours_to_process
}

.save_batch_csv_results <- function(result, result2, morpho_results, polar_dir, cartesian_dir, coordinates = "both") {
  morpho_results <- morpho_results[!sapply(morpho_results, is.null)]
  morpho_df <- if (length(morpho_results) > 0) {
    df <- do.call(rbind, lapply(morpho_results, as.data.frame))
    df[, c("Image", setdiff(names(df), "Image"))]
  } else NULL

  if (coordinates %in% c("both", "polar")) {
    result <- result[!sapply(result, is.null)]
    if (length(result) > 0) {
      save_analysis_results(result, polar_dir, "polar")
    }
    if (!is.null(morpho_df)) {
      utils::write.table(morpho_df, file = file.path(polar_dir, "MorphometricsEN.csv"), row.names = FALSE, sep = ";", dec = ".")
    }
  }
  if (coordinates %in% c("both", "cartesian")) {
    result2 <- result2[!sapply(result2, is.null)]
    if (length(result2) > 0) {
      save_analysis_results(result2, cartesian_dir, "perimeter")
    }
    if (!is.null(morpho_df) && coordinates == "cartesian") {
      utils::write.table(morpho_df, file = file.path(cartesian_dir, "MorphometricsEN.csv"), row.names = FALSE, sep = ";", dec = ".")
    }
  }
}

#' Process Images in Folder
#'
#' This function processes images in a specified folder, applying various transformations and analyses.
#' Supports multi-core parallel processing, recursive subfolder processing, manual/interactive
#' rostrum selection, and lightweight diagnostic visualization overlaid on real photographs.
#'
#' @param folder A string specifying the path to the folder containing images. Supported formats include `.jpg`, `.jpeg`, `.png`, `.tif`, and `.tiff`.
#' @param subfolder A logical value indicating whether to recursively process subfolders. Default is FALSE.
#' @param threshold A numeric value for binarization. If not provided, the Otsu method is used.
#' @param wavelets A logical value indicating whether to obtain wavelets. Default is TRUE.
#' @param ef A logical value indicating whether to obtain elliptic Fourier descriptors. Default is TRUE.
#' @param testing A logical value indicating whether to save diagnostic test images. Default is TRUE.
#' @param coordinates A string specifying which coordinate systems to compute and save:
#'                    \code{"both"} (default), \code{"polar"}, or \code{"cartesian"}.
#' @param save A logical value indicating whether to save the results as CSV files. Default is TRUE.
#' @param pixels_per_mm A numeric value specifying the scale (pixels per mm). Default is NULL.
#' @param detect_scale A logical value indicating whether to automatically detect the scale bar in the image. Default is FALSE.
#' @param blur_size Size of the Gaussian brush for smoothing in pixels. Can be an odd integer,
#'                  or \code{"auto"} (default) to dynamically calculate \eqn{\min(\text{width},\text{height}) \times 0.01}.
#' @param blur_sigma Standard deviation of the Gaussian filter brush. Can be a numeric value,
#'                   or \code{"auto"} (default) to calculate \eqn{\text{blur\_size} / 3.5}.
#' @param min_points Minimum number of points required for a valid contour (default: 200).
#' @param min_area Minimum area required for a valid contour (default: 5000).
#' @param n_harmonics Number of elliptic Fourier harmonics to calculate (default: 32).
#' @param n_points Target number of points for resampling (default: 512).
#' @param scale_mm Total length in mm of the scale bar for detection (default: 1).
#' @param procrustes A logical value indicating whether to apply Generalized Procrustes Analysis (GPA) to align all extracted contours before calculating shape descriptors. Default is FALSE.
#' @param cores Number of CPU cores to use for parallel processing. If NULL (default) or <= 1, processing is sequential.
#' @param plot_width Width in pixels for diagnostic JPEG images (default: 1200).
#' @param plot_height Height in pixels for diagnostic JPEG images (default: 800).
#' @param plot_res Resolution in dpi for diagnostic JPEG images (default: 150).
#' @param downsample_threshold Logical; if TRUE, calculates Otsu threshold on a downsampled thumbnail for faster preprocessing (default: FALSE).
#' @param threshold_method Method for automatic binarization: "adaptive", "otsu", or "fixed". Default is "adaptive".
#' @param threshold_factor Multiplier applied to the calculated threshold (default: 1.0).
#' @param interactive_rostrum Logical; if TRUE, sequentially prompts the user to click the rostrum on each image before dispatching the analysis (default: FALSE).
#' @param start_points Optional named list or data.frame with custom starting point coordinates \code{c(x, y)} or index for specific images.
#' @return Invisible NULL. The function is used for its side effects (creating files).
#' @export
#' @examples
#' \dontrun{
#' # Process all images in the "images" folder
#' process_images(folder = "path/to/images", pixels_per_mm = 100)
#'
#' # Process in parallel across 4 cores with polar coordinates only
#' process_images(folder = "path/to/images", cores = 4, coordinates = "polar")
#' }
process_images <- function(folder, subfolder = FALSE, threshold = NULL,
                           wavelets = TRUE, ef = TRUE, testing = TRUE,
                           coordinates = c("both", "polar", "cartesian"), save = TRUE,
                           pixels_per_mm = NULL, detect_scale = FALSE,
                           blur_size = "auto", blur_sigma = "auto", min_points = 200,
                           min_area = 5000, n_harmonics = 32, n_points = 512,
                           scale_mm = 1, procrustes = FALSE,
                           cores = NULL,
                           plot_width = 1200, plot_height = 800, plot_res = 150,
                           downsample_threshold = FALSE,
                           threshold_method = c("adaptive", "otsu", "fixed"),
                           threshold_factor = 1.0,
                           interactive_rostrum = FALSE,
                           start_points = NULL) {
  coordinates <- match.arg(coordinates)
  threshold_method <- match.arg(threshold_method)

  if (!is.numeric(threshold_factor) || length(threshold_factor) != 1 || is.na(threshold_factor) || threshold_factor <= 0) {
    stop("threshold_factor must be a positive number")
  }

  if (!is.character(folder) || length(folder) != 1) stop("folder must be a single character string")
  if (!dir.exists(folder)) stop(paste("Folder does not exist:", folder))
  if (!is.logical(wavelets) || !is.logical(ef) || !is.logical(testing) || !is.logical(save) || !is.logical(detect_scale) || !is.logical(procrustes) || !is.logical(subfolder) || !is.logical(downsample_threshold) || !is.logical(interactive_rostrum)) {
    stop("wavelets, ef, testing, save, detect_scale, procrustes, subfolder, downsample_threshold, and interactive_rostrum must be logical values")
  }
  if (!is.null(pixels_per_mm) && (!is.numeric(pixels_per_mm) || length(pixels_per_mm) != 1 || is.na(pixels_per_mm) || is.infinite(pixels_per_mm) || pixels_per_mm <= 0)) {
    stop("pixels_per_mm must be a positive numeric value")
  }
  if (!is.null(cores) && (!is.numeric(cores) || length(cores) != 1 || is.na(cores) || is.infinite(cores) || cores <= 0 || cores != as.integer(cores))) {
    stop("cores must be a positive integer or NULL")
  }
  if (!is.numeric(scale_mm) || length(scale_mm) != 1 || is.na(scale_mm) || is.infinite(scale_mm) || scale_mm <= 0) {
    stop("scale_mm must be a positive numeric value")
  }
  if (!is.numeric(n_points) || length(n_points) != 1 || is.na(n_points) || is.infinite(n_points) || n_points < 3 || n_points != as.integer(n_points)) {
    stop("n_points must be an integer >= 3")
  }
  if (!is.numeric(n_harmonics) || length(n_harmonics) != 1 || is.na(n_harmonics) || is.infinite(n_harmonics) || n_harmonics < 1 || n_harmonics != as.integer(n_harmonics)) {
    stop("n_harmonics must be a positive integer")
  }
  if (!is.numeric(min_points) || length(min_points) != 1 || is.na(min_points) || is.infinite(min_points) || min_points < 3 || min_points != as.integer(min_points)) {
    stop("min_points must be an integer >= 3")
  }
  if (!is.numeric(min_area) || length(min_area) != 1 || is.na(min_area) || is.infinite(min_area) || min_area < 0) {
    stop("min_area must be a numeric value >= 0")
  }

  if (subfolder) {
    all_dirs <- list.dirs(folder, full.names = TRUE, recursive = TRUE)
    all_dirs <- all_dirs[!grepl("(Polar|Cartesian)$", all_dirs)]
    dirs_with_images <- all_dirs[sapply(all_dirs, function(d) {
      length(list.files(d, pattern = "(?i)\\.(jpe?g|png|tiff?)$")) > 0
    })]
    if (length(dirs_with_images) == 0) {
      warning(paste("No subfolders containing supported images found in:", folder))
      return(invisible(NULL))
    }
    for (d in dirs_with_images) {
      message(paste("\n--- Processing subfolder:", d, "---"))
      process_images(
        folder = d, subfolder = FALSE, threshold = threshold,
        wavelets = wavelets, ef = ef, testing = testing,
        coordinates = coordinates, save = save,
        pixels_per_mm = pixels_per_mm, detect_scale = detect_scale,
        blur_size = blur_size, blur_sigma = blur_sigma,
        min_points = min_points, min_area = min_area,
        n_harmonics = n_harmonics, n_points = n_points,
        scale_mm = scale_mm, procrustes = procrustes,
        cores = cores,
        plot_width = plot_width, plot_height = plot_height, plot_res = plot_res,
        downsample_threshold = downsample_threshold,
        threshold_method = threshold_method,
        threshold_factor = threshold_factor,
        interactive_rostrum = interactive_rostrum,
        start_points = start_points
      )
    }
    return(invisible(NULL))
  }

  tryCatch({
    setup_data <- .setup_directories_and_files(folder, coordinates = coordinates)
    if (is.null(setup_data)) return(invisible(NULL))

    imag <- setup_data$imag
    n_imgs <- length(imag)

    # Interactive rostrum selection if requested
    chosen_rostrums <- vector("list", n_imgs)
    names(chosen_rostrums) <- basename(imag)

    if (isTRUE(interactive_rostrum)) {
      message("\nInteractive rostrum selection mode: Click on each otolith to set the starting point...")
      for (i in seq_along(imag)) {
        message(paste("Select rostrum for:", basename(imag[i])))
        img_preview <- EBImage::readImage(imag[i])
        bin_preview <- preprocess_image(img_preview, threshold = threshold, blur_size = blur_size,
                                        blur_sigma = blur_sigma, downsample_threshold = TRUE,
                                        threshold_method = threshold_method, threshold_factor = threshold_factor)
        cnt_preview <- extract_contour(bin_preview, min_points = min_points, min_area = min_area)
        chosen_rostrums[[i]] <- select_starting_point(img_preview, cnt_preview)
        rm(img_preview, bin_preview, cnt_preview)
      }
    }

    if (!is.null(start_points)) {
      if (is.data.frame(start_points) && all(c("Image", "X", "Y") %in% names(start_points))) {
        for (r in seq_len(nrow(start_points))) {
          img_name <- start_points$Image[r]
          chosen_rostrums[[img_name]] <- c(start_points$X[r], start_points$Y[r])
        }
      } else if (is.list(start_points)) {
        for (nm in names(start_points)) {
          chosen_rostrums[[nm]] <- start_points[[nm]]
        }
      }
    }

    # Determine parallelization
    use_parallel <- !is.null(cores) && cores > 1 && n_imgs > 1
    actual_cores <- if (use_parallel) min(as.integer(cores), n_imgs, parallel::detectCores()) else 1
    chk_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
    if (use_parallel && nzchar(chk_limit) && chk_limit != "false") {
      actual_cores <- min(actual_cores, 2L)
    }

    if (!procrustes) {
      message(paste("\nProcessing", n_imgs, "images in:", basename(folder),
                    if (use_parallel) paste("using", actual_cores, "cores") else "sequentially"))

      if (use_parallel) {
        cl <- parallel::makeCluster(actual_cores)
        on.exit({
          if (!is.null(cl)) {
            try(parallel::stopCluster(cl), silent = TRUE)
            cl <- NULL
          }
        }, add = TRUE)

        parallel::clusterEvalQ(cl, {
          suppressPackageStartupMessages(library(aforoR))
        })
        parallel::clusterExport(cl, c(".process_single_image_worker", "find_scale_bar_length", "calculate_adaptive_threshold"),
                                envir = environment(process_images))

        # Bundle image path with custom start point
        task_data <- lapply(imag, function(f) {
          list(file = f, start_point = chosen_rostrums[[basename(f)]])
        })

        processed_list <- parallel::parLapply(cl, task_data, function(item) {
          .process_single_image_worker(
            image_file = item$file, threshold = threshold, blur_size = blur_size,
            blur_sigma = blur_sigma, detect_scale = detect_scale,
            pixels_per_mm = pixels_per_mm, scale_mm = scale_mm,
            min_points = min_points, min_area = min_area,
            n_points = n_points, n_harmonics = n_harmonics,
            wavelets = wavelets, ef = ef, testing = testing, save = save,
            polar_dir = setup_data$polar_dir, cartesian_dir = setup_data$cartesian_dir,
            plot_width = plot_width, plot_height = plot_height, plot_res = plot_res,
            downsample_threshold = downsample_threshold,
            coordinates = coordinates,
            start_point = item$start_point,
            threshold_method = threshold_method,
            threshold_factor = threshold_factor
          )
        })
        try(parallel::stopCluster(cl), silent = TRUE)
        cl <- NULL
      } else {
        pb <- utils::txtProgressBar(max = n_imgs, style = 3)
        processed_list <- vector("list", n_imgs)
        for (i in seq_along(imag)) {
          utils::setTxtProgressBar(pb, i)
          processed_list[[i]] <- .process_single_image_worker(
            image_file = imag[i], threshold = threshold, blur_size = blur_size,
            blur_sigma = blur_sigma, detect_scale = detect_scale,
            pixels_per_mm = pixels_per_mm, scale_mm = scale_mm,
            min_points = min_points, min_area = min_area,
            n_points = n_points, n_harmonics = n_harmonics,
            wavelets = wavelets, ef = ef, testing = testing, save = save,
            polar_dir = setup_data$polar_dir, cartesian_dir = setup_data$cartesian_dir,
            plot_width = plot_width, plot_height = plot_height, plot_res = plot_res,
            downsample_threshold = downsample_threshold,
            coordinates = coordinates,
            start_point = chosen_rostrums[[basename(imag[i])]],
            threshold_method = threshold_method,
            threshold_factor = threshold_factor
          )
        }
        close(pb)
      }

      processed_list <- processed_list[!sapply(processed_list, is.null)]
      if (length(processed_list) == 0) {
        warning("No images were successfully processed")
        return(invisible(NULL))
      }

      result_polar <- lapply(processed_list, function(x) x$polar)
      result_cart  <- lapply(processed_list, function(x) x$perimeter)
      morpho_list  <- lapply(processed_list, function(x) x$morpho)

      if (save) {
        .save_batch_csv_results(result_polar, result_cart, morpho_list,
                                setup_data$polar_dir, setup_data$cartesian_dir,
                                coordinates = coordinates)
      }

    } else {
      # PROCRUSTES (GPA) PATH
      message(paste("\nPhase 1: Extracting contours for GPA alignment (", n_imgs, "images)..."))
      extract_worker <- function(item) {
        tryCatch({
          img_file <- item$file
          st_pt <- item$start_point
          img_raw <- EBImage::readImage(img_file)
          current_px <- pixels_per_mm
          if (detect_scale) {
            img_g <- if (EBImage::colorMode(img_raw) != EBImage::Grayscale) EBImage::channel(img_raw, "grey") else img_raw
            b_sc <- img_g > EBImage::otsu(img_g)
            d_px <- find_scale_bar_length(b_sc, scale_mm = scale_mm)
            if (!is.null(d_px)) current_px <- d_px
          }
          bin_img <- preprocess_image(
            img_raw, threshold = threshold, blur_size = blur_size,
            blur_sigma = blur_sigma, downsample_threshold = downsample_threshold,
            threshold_method = threshold_method, threshold_factor = threshold_factor
          )
          cnt <- extract_contour(bin_img, min_points, min_area)
          list(contour = cnt, px = current_px, start_point = st_pt)
        }, error = function(e) NULL)
      }

      task_data <- lapply(imag, function(f) {
        list(file = f, start_point = chosen_rostrums[[basename(f)]])
      })

      if (use_parallel) {
        cl <- parallel::makeCluster(actual_cores)
        on.exit({
          if (!is.null(cl)) {
            try(parallel::stopCluster(cl), silent = TRUE)
            cl <- NULL
          }
        }, add = TRUE)
        parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(library(aforoR)) })
        parallel::clusterExport(cl, c("find_scale_bar_length", "calculate_adaptive_threshold"), envir = environment(process_images))
        extr_list <- parallel::parLapply(cl, task_data, extract_worker)
        try(parallel::stopCluster(cl), silent = TRUE)
        cl <- NULL
      } else {
        pb <- utils::txtProgressBar(max = n_imgs, style = 3)
        extr_list <- vector("list", n_imgs)
        for (i in seq_along(imag)) {
          utils::setTxtProgressBar(pb, i)
          extr_list[[i]] <- extract_worker(task_data[[i]])
        }
        close(pb)
      }

      valid_indices <- which(!sapply(extr_list, is.null) & sapply(extr_list, function(x) !is.null(x$contour)))
      if (length(valid_indices) == 0) {
        warning("No valid contours found for Procrustes alignment.")
        return(invisible(NULL))
      }

      contours_raw <- lapply(extr_list, function(x) x$contour)
      px_list <- lapply(extr_list, function(x) x$px)
      st_list <- lapply(extr_list, function(x) x$start_point)

      contours_to_process <- .apply_gpa_alignment(contours_raw, valid_indices, procrustes = TRUE, n_points = n_points)

      message("\nPhase 2: Calculating shape descriptors on GPA-aligned contours...")
      result_polar <- list()
      result_cart  <- list()
      morpho_list  <- list()

      do_polar <- coordinates %in% c("both", "polar")
      do_cart  <- coordinates %in% c("both", "cartesian")

      pb2 <- utils::txtProgressBar(max = length(valid_indices), style = 3)
      for (idx_counter in seq_along(valid_indices)) {
        utils::setTxtProgressBar(pb2, idx_counter)
        i <- valid_indices[idx_counter]
        c_proc <- contours_to_process[[i]]
        c_raw  <- contours_raw[[i]]
        cur_px <- px_list[[i]]
        st_pt  <- st_list[[i]]

        coef_e <- if (ef) .efourier(c_proc, n_harmonics, norm = FALSE) else NULL
        distances <- calculate_distances(c_proc, n_points = n_points, start_point = st_pt)
        wavelet_results <- if (wavelets) calculate_wavelets_analysis(distances, n_scales = floor(log2(n_points)), detail = FALSE) else NULL
        morpho <- calculate_morphometrics(c_raw, cur_px)
        morpho$Image <- basename(imag[i])
        morpho_list[[length(morpho_list) + 1]] <- morpho

        if (save) {
          if (do_polar) result_polar[[length(result_polar) + 1]] <- create_result_structure(distances$polar, wavelet_results$polar, basename(imag[i]), coef_e, "polar")
          if (do_cart)  result_cart[[length(result_cart) + 1]]   <- create_result_structure(distances$perimeter, wavelet_results$perimeter, basename(imag[i]), coef_e, "perimeter")
        }

        if (testing) {
          if (do_polar) {
            save_visualization(
              image = NULL, distances = distances, wavelets = wavelet_results,
              image_name = basename(imag[i]), output_dir = setup_data$polar_dir,
              plot_width = plot_width, plot_height = plot_height, plot_res = plot_res
            )
          }
          if (do_cart) {
            save_visualization_perimeter(
              image = NULL, distances = distances, wavelets = wavelet_results,
              image_name = basename(imag[i]), output_dir = setup_data$cartesian_dir,
              plot_width = plot_width, plot_height = plot_height, plot_res = plot_res
            )
          }
        }
      }
      close(pb2)

      if (save) {
        .save_batch_csv_results(result_polar, result_cart, morpho_list,
                                setup_data$polar_dir, setup_data$cartesian_dir,
                                coordinates = coordinates)
      }
    }
  }, error = function(e) {
    stop(paste("Error in process_images:", e$message))
  })
}

find_scale_bar_length <- function(binary_image, scale_mm = 1) {
  if (!is.numeric(scale_mm) || length(scale_mm) != 1 || is.na(scale_mm) || is.infinite(scale_mm) || scale_mm <= 0) {
    stop("scale_mm must be a positive numeric value")
  }
  tryCatch(
    {
      labeled <- EBImage::bwlabel(binary_image)
      features <- EBImage::computeFeatures.shape(labeled)
      moments <- EBImage::computeFeatures.moment(labeled)

      if (is.null(features) || nrow(features) == 0) {
        return(NULL)
      }

      dims <- dim(binary_image)
      total_px <- prod(dims[1:2])
      max_area <- max(5000, min(30000, 0.005 * total_px))
      min_area <- max(50, 0.00005 * total_px)
      largest_idx <- which.max(features[, "s.area"])

      # Scale bar is a thin straight line and never the main specimen
      candidates <- which(
        seq_len(nrow(features)) != largest_idx &
        features[, "s.area"] > min_area &
        features[, "s.area"] < max_area &
        moments[, "m.eccentricity"] > 0.95
      )

      if (length(candidates) == 0) {
        candidates <- which(
          seq_len(nrow(features)) != largest_idx &
          features[, "s.area"] > min_area &
          features[, "s.area"] < max_area &
          moments[, "m.eccentricity"] > 0.88
        )
      }

      if (length(candidates) == 0) {
        return(NULL)
      }

      best_idx <- which.max(moments[candidates, "m.majoraxis"])
      best_candidate <- candidates[best_idx]

      val <- moments[best_candidate, "m.majoraxis"] / scale_mm
      bar_pts <- which(labeled == best_candidate, arr.ind = TRUE)
      if (length(bar_pts) > 0) {
        bbox <- list(
          min_x = min(bar_pts[, 1]),
          max_x = max(bar_pts[, 1]),
          min_y = min(bar_pts[, 2]),
          max_y = max(bar_pts[, 2]),
          bar_length_px = moments[best_candidate, "m.majoraxis"],
          scale_mm = scale_mm,
          px_per_mm = val
        )
        attr(val, "scale_info") <- bbox
      }
      return(val)
    },
    error = function(e) {
      warning(paste("Scale detection error:", e$message))
      invisible(NULL)
    }
  )
}

#' Create Standardized Analysis Result Structure
#'
#' Creates a standardized result structure for storing analysis data.
#'
#' @param distances The distance data (polar or perimeter).
#' @param wavelets The wavelet analysis results.
#' @param image_name The name of the processed image.
#' @param elliptic_coef The elliptic Fourier coefficients (optional).
#' @param result_type The type of result ("polar" or "perimeter").
#' @return A list containing the structured analysis results.
#' @export
#' @examples
#' dists <- list(radii = runif(512), dist = runif(512), normalized = runif(512))
#' res <- create_result_structure(dists, NULL, "sample", result_type = "polar")
#' str(res)
create_result_structure <- function(distances, wavelets, image_name,
                                    elliptic_coef = NULL,
                                    result_type = "polar") {
  if (!is.list(distances)) {
    stop("distances must be a list")
  }

  if (!is.character(image_name) || length(image_name) != 1) {
    stop("image_name must be a single character string")
  }

  if (!result_type %in% c("polar", "perimeter")) {
    stop("result_type must be either 'polar' or 'perimeter'")
  }

  result <- list(
    Distancia = if (result_type == "polar") distances$radii else distances$dist,
    Distancia_Norm = distances$normalized,
    nombre = image_name
  )

  if (!is.null(wavelets)) {
    if (is.vector(wavelets)) {
      wavelets <- matrix(wavelets, nrow = 1)
    } else if (!is.matrix(wavelets) && !is.data.frame(wavelets)) {
      stop("wavelets must be a matrix, data.frame, vector, or NULL")
    }
    for (i in seq_len(nrow(wavelets))) {
      result[[paste0("Wavelet_", i)]] <- wavelets[i, ]
    }
  }

  coord_data <- if (result_type == "polar") {
    if (!is.null(distances$coord)) distances$coord else distances$coords
  } else {
    if (!is.null(distances$coords)) distances$coords else distances$coord
  }
  if (!is.null(coord_data)) {
    result$coords <- cbind(image_name, coord_data)
  }

  if (!is.null(elliptic_coef)) {
    result$elliptic <- c(
      elliptic_coef$an, elliptic_coef$bn, elliptic_coef$cn, elliptic_coef$dn
    )
  }

  return(result)
}

#' Save Perimeter Visualization Images
#'
#' Creates and saves visualization images for perimeter-based analysis.
#' Overlays the extracted contour, centroid, rostrum, and perimeter sectors
#' directly onto the original photograph in high contrast colors.
#'
#' @param image An \code{Image} object (original or binary) OR a string specifying
#'              the image file path. If NULL, draws on coordinate space (useful for GPA).
#' @param distances The distance measures from \code{\link{calculate_distances}}.
#' @param wavelets The wavelet analysis from \code{\link{calculate_wavelets_analysis}}.
#' @param image_name The base name of the image file.
#' @param output_dir The directory to save visualization images.
#' @param wavelet_scale The wavelet scale to visualize (default: 5).
#' @param plot_width Width in pixels for diagnostic JPEG images (default: 1200).
#' @param plot_height Height in pixels for diagnostic JPEG images (default: 800).
#' @param plot_res Resolution in dpi for diagnostic JPEG images (default: 150).
#' @param binary_image Deprecated alias for \code{image} for backwards compatibility.
#' @param scale_info Optional list with scale bar metrics (min_x, max_x, min_y, max_y, bar_length_px, scale_mm, px_per_mm).
#' @return An invisible NULL, called for its side effect of saving diagnostic visualization images.
#' @importFrom grDevices dev.off jpeg adjustcolor rgb
#' @importFrom graphics lines points polygon legend
#' @export
#' @examples
#' \dontrun{
#' save_visualization_perimeter(image, distances, wavelets, "sample_image", tempdir())
#' }
save_visualization_perimeter <- function(image = NULL, distances, wavelets,
                                         image_name, output_dir,
                                         wavelet_scale = 5,
                                         plot_width = 1200, plot_height = 800, plot_res = 150,
                                         binary_image = NULL, scale_info = NULL) {
  if (is.null(image) && !is.null(binary_image)) {
    image <- binary_image
  }

  if (!is.list(distances)) {
    stop("distances must be a list")
  }

  if (!is.null(wavelets) && !is.list(wavelets)) {
    stop("wavelets must be a list or NULL")
  }

  if (!is.character(image_name) || length(image_name) != 1) {
    stop("image_name must be a single character string")
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }

  if (!is.numeric(wavelet_scale) || length(wavelet_scale) != 1 || is.na(wavelet_scale) || wavelet_scale <= 0) {
    stop("wavelet_scale must be a positive integer")
  }

  if (!is.numeric(plot_width) || length(plot_width) != 1 || is.na(plot_width) || plot_width <= 0) {
    stop("plot_width must be a positive number")
  }

  if (!is.numeric(plot_height) || length(plot_height) != 1 || is.na(plot_height) || plot_height <= 0) {
    stop("plot_height must be a positive number")
  }

  if (!is.numeric(plot_res) || length(plot_res) != 1 || is.na(plot_res) || plot_res <= 0) {
    stop("plot_res must be a positive number")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cur_dev <- NULL
  on.exit({
    if (!is.null(cur_dev) && cur_dev > 1) {
      try(grDevices::dev.off(cur_dev), silent = TRUE)
      cur_dev <- NULL
    }
  }, add = TRUE)

  tryCatch(
    {
      zonal_file <- file.path(
        output_dir,
        paste0("zonal_", tools::file_path_sans_ext(image_name), ".jpg")
      )
      jpeg(filename = zonal_file, width = plot_width, height = plot_height,
           res = plot_res, quality = 90)
      cur_dev <- grDevices::dev.cur()

      has_img <- FALSE
      if (!is.null(image)) {
        if (is.character(image) && file.exists(image)) {
          img_to_plot <- EBImage::readImage(image)
          plot(img_to_plot)
          has_img <- TRUE
        } else if (inherits(image, "Image")) {
          plot(image)
          has_img <- TRUE
        }
      }

      cx <- mean(distances$reordered_coords$x)
      cy <- mean(distances$reordered_coords$y)

      if (!has_img) {
        xlims <- range(distances$perimeter$coords[, 1])
        ylims <- range(distances$perimeter$coords[, 2])
        xpad <- if (is.na(diff(xlims)) || diff(xlims) == 0) 1 else diff(xlims) * 0.1
        ypad <- if (is.na(diff(ylims)) || diff(ylims) == 0) 1 else diff(ylims) * 0.1
        plot(NA, xlim = c(xlims[1] - xpad, xlims[2] + xpad),
             ylim = c(ylims[1] - ypad, ylims[2] + ypad),
             xlab = "X", ylab = "Y", main = paste("Perimeter Analysis -", image_name))
      }

      points(cx, cy, pch = 3, col = "red", lwd = 2.5, cex = 1.6)
      points(distances$reordered_coords$x[1], distances$reordered_coords$y[1],
             pch = 16, col = "red", cex = 1.6)

      n_points <- length(distances$perimeter$dist)
      b <- round(seq(0, n_points, length.out = 5))
      zones <- lapply(1:4, function(idx) if (b[idx + 1] > b[idx]) (b[idx] + 1L):b[idx + 1] else integer(0))
      colors <- c("#00E5FF", "#FFEA00", "#FF007F", "#00E676")

      for (i in seq_along(zones)) {
        if (length(zones[[i]]) > 0) {
          lines(
            distances$perimeter$coords[zones[[i]], 1],
            distances$perimeter$coords[zones[[i]], 2],
            col = colors[i], lwd = 2.5
          )
        }
      }

      if (!is.null(scale_info) && is.list(scale_info) && !is.null(scale_info$min_x)) {
        pad <- 15
        graphics::rect(scale_info$min_x - pad, scale_info$min_y - pad,
                       scale_info$max_x + pad, scale_info$max_y + pad,
                       border = "#00FF66", lwd = 2.5)
        scale_label <- sprintf("Scale: %g mm (%.0f px | %.1f px/mm)",
                               scale_info$scale_mm, scale_info$bar_length_px, scale_info$px_per_mm)
        img_w <- if (has_img) dim(image)[1] else max(scale_info$max_x + 100, 1000)
        if (mean(c(scale_info$min_x, scale_info$max_x)) > img_w * 0.5) {
          t_x <- min(img_w - 30, scale_info$max_x + pad)
          a_x <- 1
        } else {
          t_x <- max(30, scale_info$min_x - pad)
          a_x <- 0
        }
        t_y <- max(35, scale_info$min_y - pad - 20)
        graphics::text(t_x, t_y, labels = scale_label, adj = c(a_x, 0.5),
                       col = "#00FF66", font = 2, cex = 0.80)
      }

      leg_labels <- c("Zone 1", "Zone 2", "Zone 3", "Zone 4", "Centroid (+)", "Rostrum (o)")
      leg_cols   <- c(colors, "red", "red")
      leg_pch    <- c(rep(NA, 4), 3, 16)
      leg_lty    <- c(rep(1, 4), NA, NA)

      if (!is.null(scale_info) && is.list(scale_info) && !is.null(scale_info$min_x)) {
        leg_labels <- c(leg_labels, sprintf("Scale bar (%g mm)", scale_info$scale_mm))
        leg_cols   <- c(leg_cols, "#00FF66")
        leg_pch    <- c(leg_pch, NA)
        leg_lty    <- c(leg_lty, 1)
      }

      legend("topleft",
             legend = leg_labels,
             col = leg_cols, pch = leg_pch,
             lty = leg_lty, lwd = 2, bg = rgb(0, 0, 0, 0.65), text.col = "white",
             cex = 0.80, pt.cex = 1.1, inset = c(0.02, 0.03))
      if (!is.null(cur_dev) && cur_dev > 1) {
        try(grDevices::dev.off(cur_dev), silent = TRUE)
        cur_dev <- NULL
      }

      if (!is.null(wavelets) && !is.null(wavelets$perimeter) && nrow(wavelets$perimeter) > 0) {
        scale_to_plot <- min(as.integer(wavelet_scale), nrow(wavelets$perimeter))
        wav_file <- file.path(
          output_dir,
          paste0("wavelet_", tools::file_path_sans_ext(image_name), ".jpg")
        )
        jpeg(filename = wav_file, width = 800, height = 500, res = 100, quality = 85)
        cur_dev <- grDevices::dev.cur()
        plot(wavelets$perimeter[scale_to_plot, ],
          main = paste("Perimeter Wavelet Scale", scale_to_plot, "-", image_name),
          xlab = "Points", ylab = "Coefficient", type = "l", col = "#00E5FF", lwd = 2
        )

        for (i in seq_along(zones)) {
          if (length(zones[[i]]) > 0) {
            polygon(
              x = c(min(zones[[i]]), zones[[i]], max(zones[[i]])),
              y = c(
                min(wavelets$perimeter[scale_to_plot, ]),
                wavelets$perimeter[scale_to_plot, zones[[i]]],
                min(wavelets$perimeter[scale_to_plot, ])
              ),
              col = adjustcolor(colors[i], alpha.f = 0.5),
              border = colors[i]
            )
          }
        }
        if (!is.null(cur_dev) && cur_dev > 1) {
          try(grDevices::dev.off(cur_dev), silent = TRUE)
          cur_dev <- NULL
        }
      }
    },
    error = function(e) {
      stop(paste("Error creating perimeter visualization:", e$message))
    }
  )
}

