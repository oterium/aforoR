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

  if (inherits(img, "Image")) {
    if (length(dim(img)) == 3 && dim(img)[3] == 2) {
      img <- img[, , 1]
      EBImage::colorMode(img) <- 0L
    } else if (EBImage::colorMode(img) != EBImage::Grayscale) {
      img <- EBImage::channel(img, "grey")
    }
  } else if (is.array(img) && length(dim(img)) >= 3) {
    if (dim(img)[3] == 2) {
      img <- img[, , 1]
    } else if (dim(img)[3] >= 3) {
      img <- 0.2989 * img[, , 1] + 0.5870 * img[, , 2] + 0.1140 * img[, , 3]
    }
  }

  if (is.numeric(img) && any(!is.na(img)) && max(img, na.rm = TRUE) > 1.05 && max(img, na.rm = TRUE) <= 255) {
    img <- img / 255
  }

  # Strictly clamp intensities to [0, 1] to prevent histogram and Otsu range errors
  img[img < 0] <- 0
  img[img > 1] <- 1

  vals <- as.vector(img)
  vals <- vals[!is.na(vals) & vals >= 0 & vals <= 1]
  if (length(vals) < 100) return(0.2)

  th_std <- tryCatch(
    {
      img_for_otsu <- if (is.array(img)) img else as.matrix(img)
      img_for_otsu[img_for_otsu < 0] <- 0
      img_for_otsu[img_for_otsu > 1] <- 1
      t_raw <- EBImage::otsu(img_for_otsu)
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

  # Safeguard: If adaptive threshold is unrealistically low (< 0.05),
  # implausibly low compared to Otsu (< 0.65 * th_std) as in textured grey
  # backgrounds (e.g. caballa/PHOTOS), or higher than standard Otsu,
  # fall back to 0.75 * standard Otsu
  if (is.na(th_adaptive) || th_adaptive < 0.05 || th_adaptive < 0.65 * th_std || th_adaptive >= th_std) {
    th_adaptive <- th_std * 0.75
  }

  return(th_adaptive)
}

#' Read and Process Image with Multi-Channel Normalization
#'
#' Reads an image from disk or accepts an existing Image object, handles 2-channel
#' (grayscale + alpha) TIFF images by stripping the alpha channel and setting
#' Grayscale colormode, and returns a clean Image object.
#'
#' @param image A character string path to an image file or an \code{Image} object.
#' @return An \code{Image} object in Grayscale or Color mode, with any alpha channel removed.
#' @export
#' @examples
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' img <- read_and_process_image(image_path)
read_and_process_image <- function(image) {
  if (inherits(image, "Image")) {
    img <- image
  } else if (is.character(image) && length(image) == 1) {
    if (!file.exists(image)) {
      stop(paste("Image file does not exist:", image))
    }
    img <- EBImage::readImage(image)
  } else if (is.array(image)) {
    img <- EBImage::Image(image)
  } else {
    stop("image must be a single file path character string or an Image object")
  }

  if (length(dim(img)) == 3 && dim(img)[3] == 2) {
    img <- img[, , 1]
    EBImage::colorMode(img) <- 0L
  }

  return(img)
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

  if (is.null(threshold_factor) || !is.numeric(threshold_factor) || length(threshold_factor) != 1 || is.na(threshold_factor) || is.infinite(threshold_factor) || threshold_factor <= 0) {
    stop("threshold_factor must be a positive numeric value")
  }

  tryCatch(
    {
      foto <- read_and_process_image(image_path)

      if (length(dim(foto)) == 3 && dim(foto)[3] == 2) {
        foto <- foto[, , 1]
        EBImage::colorMode(foto) <- 0L
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

      # Clamp to [0, 1] to eliminate floating-point overshoot/undershoot from FFT convolution
      foto2[foto2 < 0] <- 0
      foto2[foto2 > 1] <- 1

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
          foto2_th[foto2_th < 0] <- 0
          foto2_th[foto2_th > 1] <- 1
        } else {
          foto2_th <- foto2
        }

        if (threshold_method == "adaptive") {
          th <- calculate_adaptive_threshold(foto2_th)
        } else {
          foto2_otsu <- if (is.array(foto2_th)) foto2_th else as.matrix(foto2_th)
          foto2_otsu[foto2_otsu < 0] <- 0
          foto2_otsu[foto2_otsu > 1] <- 1
          th <- EBImage::otsu(foto2_otsu)
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
#' @param min_area Minimum area required for a valid contour. Can be \code{"auto"} (default)
#'                 to dynamically adapt the threshold to 0.5% of total image area
#'                 (bounded between 500 and 5000 \eqn{\text{px}^2}), or a positive numeric value in pixels.
#' @return A matrix of (X, Y) coordinates of the main contour, or NULL if no
#'         suitable contour is found.
#' @export
#' @examples
#' # Process sample image and extract its contour
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' binary_img <- preprocess_image(image_path)
#' contour <- extract_contour(binary_img)
#' head(contour)
extract_contour <- function(binary_image, min_points = 200, min_area = "auto") {
  if (is.null(binary_image)) {
    stop("binary_image cannot be NULL")
  }

  if (!is.numeric(min_points) || length(min_points) != 1 || is.na(min_points) ||
      is.infinite(min_points) || min_points < 3 || min_points != as.integer(min_points)) {
    stop("min_points must be an integer >= 3")
  }
  min_points <- as.integer(min_points)

  dims <- dim(binary_image)
  if (is.null(dims) || length(dims) < 2) {
    stop("binary_image must be a 2D matrix or Image")
  }
  total_area <- as.numeric(dims[1]) * as.numeric(dims[2])
  max_allowed_area <- 0.85 * total_area

  is_auto_area <- identical(min_area, "auto") || is.null(min_area)
  if (is_auto_area) {
    calc_min_area <- min(5000, max(500, total_area * 0.005))
  } else if (is.numeric(min_area) && length(min_area) == 1 && !is.na(min_area) && !is.infinite(min_area) && min_area >= 0) {
    calc_min_area <- as.numeric(min_area)
  } else {
    stop("min_area must be 'auto', NULL, or a numeric value >= 0")
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
        # Discard any contour whose area is >= 85% of total image area (outer frame border)
        if (current_area > calc_min_area && current_area < max_allowed_area && current_area > max_area) {
          max_area <- current_area
          contorno <- c_mat
        }
      }

      if (is.null(contorno)) {
        warning(paste(
          "No contour found meeting size criteria (> ", min_points,
          " points and area > ", round(calc_min_area), ")",
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
