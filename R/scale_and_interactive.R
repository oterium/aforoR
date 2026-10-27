#' Interactively Select Contour Starting Point (Rostrum)
#'
#' Displays an otolith image with its extracted contour overlaid in cyan and the
#' automatic rightmost rostrum candidate highlighted with a yellow circle (\code{o}).
#' Prompts the user to click on the anatomical rostrum with the left mouse button
#' (snapping automatically to the nearest contour vertex and marking it in red),
#' or right-click / press \code{Esc} to keep the automatic rightmost point.
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
    img_obj <- read_and_process_image(image)
  } else if (inherits(image, "Image")) {
    img_obj <- image
  } else {
    stop("image must be an Image object or a valid file path")
  }

  if (!interactive() && is.null(getOption("aforoR.SimLocator"))) {
    message("Non-interactive session: using automatic rightmost point.")
    return(NULL)
  }

  if (is.null(contour)) {
    bin_auto <- tryCatch(preprocess_image(img_obj, downsample_threshold = TRUE), error = function(e) NULL)
    if (!is.null(bin_auto)) {
      contour <- tryCatch(extract_contour(bin_auto), error = function(e) NULL)
    }
  }

  plot(img_obj)
  if (!is.null(contour) && (is.matrix(contour) || is.data.frame(contour)) && nrow(contour) >= 3) {
    lines(contour[, 1], contour[, 2], col = "#00E5FF", lwd = 2)
    auto_idx <- which.max(contour[, 1])
    points(contour[auto_idx, 1], contour[auto_idx, 2], pch = 1, col = "#FFEA00", lwd = 2.5, cex = 1.8)
  }
  legend("topleft",
         legend = c("Click on Rostrum (starting point)", "Right-click / Esc = Keep Auto (o)"),
         col = c("red", "#FFEA00"), pch = c(16, 1),
         bg = rgb(0, 0, 0, 0.7), text.col = "white", cex = 0.85, inset = c(0.02, 0.02))

  sim_loc <- getOption("aforoR.SimLocator")
  pt <- tryCatch({
    if (is.function(sim_loc)) sim_loc(n = 1) else graphics::locator(n = 1)
  }, error = function(e) NULL)

  if (is.null(pt) || length(pt$x) == 0) {
    message("No point clicked; using automatic rightmost point.")
    return(NULL)
  }

  sel_x <- pt$x[1]
  sel_y <- pt$y[1]
  if (!is.null(contour) && (is.matrix(contour) || is.data.frame(contour)) && nrow(contour) >= 3) {
    d_sq <- (contour[, 1] - sel_x)^2 + (contour[, 2] - sel_y)^2
    snap_idx <- which.min(d_sq)
    sel_x <- unname(contour[snap_idx, 1])
    sel_y <- unname(contour[snap_idx, 2])
  }
  points(sel_x, sel_y, pch = 16, col = "red", cex = 1.6)

  return(c(x = sel_x, y = sel_y))
}

#' Interactively Select Calibration Scale Bar Endpoints
#'
#' Displays an otolith image and prompts the user to click two points with the mouse
#' corresponding to the start and end of the calibration scale bar (no contour snapping).
#' Computes the exact Euclidean pixel distance between the two clicked endpoints and
#' returns the calibration factor in pixels per millimeter (\code{px_per_mm}) along with
#' the \code{scale_info} attribute for diagnostic overlay rendering.
#'
#' @param image An \code{Image} object or file path.
#' @param scale_mm Numeric value specifying the physical length of the scale bar in mm (default: 1).
#' @return A numeric scalar representing \code{pixels_per_mm} with attribute \code{scale_info},
#'         or NULL if cancelled or in a non-interactive session.
#' @importFrom graphics lines points legend
#' @export
#' @examples
#' \dontrun{
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' px_mm <- select_scale_bar(image_path, scale_mm = 1)
#' }
select_scale_bar <- function(image, scale_mm = 1) {
  if (!is.numeric(scale_mm) || length(scale_mm) != 1 || is.na(scale_mm) || is.infinite(scale_mm) || scale_mm <= 0) {
    stop("scale_mm must be a positive numeric value")
  }

  if (is.character(image) && file.exists(image)) {
    img_obj <- read_and_process_image(image)
  } else if (inherits(image, "Image")) {
    img_obj <- image
  } else {
    stop("image must be an Image object or a valid file path")
  }

  if (!interactive() && is.null(getOption("aforoR.SimLocator"))) {
    message("Non-interactive session: skipping manual scale bar selection.")
    return(NULL)
  }

  plot(img_obj)
  legend("topleft",
         legend = c(sprintf("Click 2 points: start & end of %g mm scale bar", scale_mm),
                    "Right-click / Esc = Skip / Cancel"),
         col = c("#00FF66", "#FF0033"), pch = c(4, NA), lty = c(1, NA), lwd = 2,
         bg = rgb(0, 0, 0, 0.7), text.col = "white", cex = 0.85, inset = c(0.02, 0.02))

  sim_loc <- getOption("aforoR.SimLocator")
  pts <- tryCatch({
    if (is.function(sim_loc)) sim_loc(n = 2) else graphics::locator(n = 2)
  }, error = function(e) NULL)

  if (is.null(pts) || length(pts$x) < 2) {
    message("Fewer than 2 points clicked; scale bar selection cancelled.")
    return(NULL)
  }

  x1 <- as.integer(round(pts$x[1]))
  y1 <- as.integer(round(pts$y[1]))
  x2 <- as.integer(round(pts$x[2]))
  y2 <- as.integer(round(pts$y[2]))

  graphics::lines(c(x1, x2), c(y1, y2), col = "#00FF66", lwd = 3)
  graphics::points(c(x1, x2), c(y1, y2), pch = 4, col = "#FF0033", lwd = 3, cex = 1.6)
  graphics::points(c(x1, x2), c(y1, y2), pch = 1, col = "#00FF66", lwd = 2, cex = 1.8)

  bar_length_px <- sqrt((pts$x[2] - pts$x[1])^2 + (pts$y[2] - pts$y[1])^2)
  px_per_mm <- bar_length_px / scale_mm
  val <- px_per_mm
  scale_info <- list(
    min_x = min(x1, x2),
    max_x = max(x1, x2),
    min_y = min(y1, y2),
    max_y = max(y1, y2),
    x1 = x1,
    y1 = y1,
    x2 = x2,
    y2 = y2,
    endpoints = list(x1 = x1, y1 = y1, x2 = x2, y2 = y2),
    bar_length_px = bar_length_px,
    scale_mm = scale_mm,
    px_per_mm = px_per_mm,
    pixels_per_mm = px_per_mm,
    confidence = 1.0
  )
  attr(val, "scale_info") <- scale_info
  return(val)
}

#' Automatically Detect Calibration Scale Bar Length
#'
#' Detects horizontal, vertical, inverted plaque, or chromatic red calibration
#' scale bars in an image and calculates the calibration factor in pixels per mm.
#'
#' @param binary_image Binarized or raw image matrix/Image.
#' @param scale_mm Physical length of the scale bar in mm (default: 1).
#' @param bar_type Scale bar type (default: "auto").
#' @param raw_image Optional original color/grayscale Image object.
#' @return Numeric scalar \code{pixels_per_mm} with attribute \code{scale_info}, or NULL.
#' @importFrom stats median
#' @keywords internal
#' @noRd
find_scale_bar_length <- function(binary_image, scale_mm = 1, bar_type = "auto", raw_image = NULL) {
  if (!is.numeric(scale_mm) || length(scale_mm) != 1 || is.na(scale_mm) || is.infinite(scale_mm) || scale_mm <= 0) {
    stop("scale_mm must be a positive numeric value")
  }
  tryCatch(
    {
      # Handle 2-channel TIFFs (Grayscale + Alpha)
      if (length(dim(binary_image)) >= 3 && dim(binary_image)[3] == 2) {
        binary_image <- EBImage::Image(binary_image[, , 1], colormode = "Grayscale")
      }
      if (!is.null(raw_image) && length(dim(raw_image)) >= 3 && dim(raw_image)[3] == 2) {
        raw_image <- EBImage::Image(raw_image[, , 1], colormode = "Grayscale")
      }

      color_img <- NULL
      if (!is.null(raw_image) && (length(dim(raw_image)) >= 3 && dim(raw_image)[3] >= 3)) {
        color_img <- raw_image
      } else if ((length(dim(binary_image)) >= 3 && dim(binary_image)[3] >= 3) ||
                 (inherits(binary_image, "Image") && EBImage::colorMode(binary_image) == EBImage::Color && length(dim(binary_image)) >= 3 && dim(binary_image)[3] >= 3)) {
        color_img <- binary_image
      } else if (exists("img_raw", envir = parent.frame(), inherits = FALSE)) {
        parent_raw <- get("img_raw", envir = parent.frame())
        if ((length(dim(parent_raw)) >= 3 && dim(parent_raw)[3] >= 3) ||
            (inherits(parent_raw, "Image") && EBImage::colorMode(parent_raw) == EBImage::Color && length(dim(parent_raw)) >= 3 && dim(parent_raw)[3] >= 3)) {
          color_img <- parent_raw
        }
      }

      # 1. Chromatic red mask detection for thin red scale bars (e.g. Canaries_JosePajuelo)
      if (!is.null(color_img)) {
        r <- color_img[, , 1]
        g <- color_img[, , 2]
        b <- color_img[, , 3]
        red_mask <- (r > 0.45) & (r > 1.4 * pmax(g, b))
        if (any(red_mask)) {
          lab_red <- EBImage::bwlabel(red_mask)
          feat_red <- EBImage::computeFeatures.shape(lab_red)
          mom_red <- EBImage::computeFeatures.moment(lab_red)
          if (!is.null(feat_red) && nrow(feat_red) > 0) {
            # Isolated red scale bar: do NOT discard largest_idx
            red_cands <- which(feat_red[, "s.area"] > 15 & mom_red[, "m.eccentricity"] > 0.88)
            if (length(red_cands) > 0) {
              valid_cands <- list()
              for (k in red_cands) {
                pts <- which(lab_red == k, arr.ind = TRUE)
                dx <- max(pts[, 1]) - min(pts[, 1]) + 1
                dy <- max(pts[, 2]) - min(pts[, 2]) + 1
                span <- max(dx, dy)
                aspect <- span / max(1, min(dx, dy))
                if (span >= 40 && aspect >= 3) {
                  valid_cands[[length(valid_cands) + 1]] <- list(idx = k, span = span, pts = pts, dx = dx, dy = dy)
                }
              }
              if (length(valid_cands) > 0) {
                best_cand <- valid_cands[[which.max(sapply(valid_cands, function(vc) vc$span))]]
                bar_pts <- best_cand$pts
                x_coords <- bar_pts[, 1]
                y_coords <- bar_pts[, 2]

                min_x <- min(x_coords)
                max_x <- max(x_coords)
                min_y <- min(y_coords)
                max_y <- max(y_coords)
                delta_x <- max_x - min_x + 1
                delta_y <- max_y - min_y + 1

                if (delta_x >= delta_y) {
                  bar_length_px <- delta_x
                  x1 <- as.integer(min_x)
                  x2 <- as.integer(max_x)
                  left_pts <- y_coords[x_coords == min_x]
                  right_pts <- y_coords[x_coords == max_x]
                  y1 <- as.integer(round(median(if (length(left_pts) > 0) left_pts else y_coords)))
                  y2 <- as.integer(round(median(if (length(right_pts) > 0) right_pts else y_coords)))
                } else {
                  bar_length_px <- delta_y
                  y1 <- as.integer(min_y)
                  y2 <- as.integer(max_y)
                  top_pts <- x_coords[y_coords == min_y]
                  bottom_pts <- x_coords[y_coords == max_y]
                  x1 <- as.integer(round(median(if (length(top_pts) > 0) top_pts else x_coords)))
                  x2 <- as.integer(round(median(if (length(bottom_pts) > 0) bottom_pts else x_coords)))
                }

                if (abs(x2 - x1) > 10 && abs(y2 - y1) > 10) {
                  eucl_dist <- sqrt((x2 - x1)^2 + (y2 - y1)^2) + 1
                  if (eucl_dist > bar_length_px) {
                    bar_length_px <- as.integer(round(eucl_dist))
                  }
                }

                px_per_mm <- bar_length_px / scale_mm
                val <- px_per_mm
                scale_info <- list(
                  min_x = min_x,
                  max_x = max_x,
                  min_y = min_y,
                  max_y = max_y,
                  x1 = x1,
                  y1 = y1,
                  x2 = x2,
                  y2 = y2,
                  endpoints = list(x1 = x1, y1 = y1, x2 = x2, y2 = y2),
                  bar_length_px = bar_length_px,
                  scale_mm = scale_mm,
                  px_per_mm = px_per_mm,
                  pixels_per_mm = px_per_mm,
                  confidence = 1.0
                )
                attr(val, "scale_info") <- scale_info
                return(val)
              }
            }
          }
        }
      }

      # 2. Standard / Plaque / Vertical grayscale scale bar detection
      if (length(dim(binary_image)) >= 3 && dim(binary_image)[3] >= 3) {
        img_grey <- EBImage::channel(binary_image, "grey")
        binary_image <- img_grey > EBImage::otsu(img_grey)
      } else if (is.numeric(binary_image) && !is.logical(binary_image)) {
        if (max(binary_image, na.rm = TRUE) > 1.05 || (min(binary_image, na.rm = TRUE) >= 0 && max(binary_image, na.rm = TRUE) <= 1 && any(binary_image > 0 & binary_image < 1))) {
          binary_image <- binary_image > EBImage::otsu(binary_image)
        }
      }

      labeled <- EBImage::bwlabel(binary_image)
      features <- EBImage::computeFeatures.shape(labeled)
      moments <- EBImage::computeFeatures.moment(labeled)

      if (is.null(features) || nrow(features) == 0) {
        return(NULL)
      }

      dims <- dim(binary_image)
      total_px <- prod(dims[1:2])
      max_area <- max(5000, min(30000, 0.005 * total_px))
      min_area <- max(10, 0.00002 * total_px)
      largest_idx <- which.max(features[, "s.area"])

      # Scale bar is a thin straight line or plaque, never the main specimen
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

      # Fallback for plaque-style inverted bars if no candidate found in foreground
      if (length(candidates) == 0) {
        inv_bin <- !binary_image
        inv_lab <- EBImage::bwlabel(inv_bin)
        inv_feat <- EBImage::computeFeatures.shape(inv_lab)
        inv_mom <- EBImage::computeFeatures.moment(inv_lab)
        if (!is.null(inv_feat) && nrow(inv_feat) > 0) {
          inv_largest <- which.max(inv_feat[, "s.area"])
          inv_cands <- which(
            seq_len(nrow(inv_feat)) != inv_largest &
            inv_feat[, "s.area"] > min_area &
            inv_feat[, "s.area"] < max_area &
            inv_mom[, "m.eccentricity"] > 0.88
          )
          if (length(inv_cands) > 0) {
            labeled <- inv_lab
            features <- inv_feat
            moments <- inv_mom
            candidates <- inv_cands
          }
        }
      }

      if (length(candidates) == 0) {
        return(NULL)
      }

      best_idx <- which.max(moments[candidates, "m.majoraxis"])
      best_candidate <- candidates[best_idx]

      bar_pts <- which(labeled == best_candidate, arr.ind = TRUE)
      if (length(bar_pts) == 0) {
        return(NULL)
      }

      is_base_matrix <- is.matrix(binary_image) && !inherits(binary_image, "Image")
      if (is_base_matrix) {
        # Base R matrix: dim 1 is rows (Y), dim 2 is cols (X)
        x_coords <- bar_pts[, 2]
        y_coords <- bar_pts[, 1]
      } else {
        # EBImage Image: dim 1 is width (X), dim 2 is height (Y)
        x_coords <- bar_pts[, 1]
        y_coords <- bar_pts[, 2]
      }

      min_x <- min(x_coords)
      max_x <- max(x_coords)
      min_y <- min(y_coords)
      max_y <- max(y_coords)
      delta_x <- max_x - min_x + 1
      delta_y <- max_y - min_y + 1

      # Physical length and tips: horizontal vs vertical orientation
      if (delta_x >= delta_y) {
        bar_length_px <- delta_x
        x1 <- as.integer(min_x)
        x2 <- as.integer(max_x)
        left_pts <- y_coords[x_coords == min_x]
        right_pts <- y_coords[x_coords == max_x]
        y1 <- as.integer(round(median(if (length(left_pts) > 0) left_pts else y_coords)))
        y2 <- as.integer(round(median(if (length(right_pts) > 0) right_pts else y_coords)))
      } else {
        bar_length_px <- delta_y
        y1 <- as.integer(min_y)
        y2 <- as.integer(max_y)
        top_pts <- x_coords[y_coords == min_y]
        bottom_pts <- x_coords[y_coords == max_y]
        x1 <- as.integer(round(median(if (length(top_pts) > 0) top_pts else x_coords)))
        x2 <- as.integer(round(median(if (length(bottom_pts) > 0) bottom_pts else x_coords)))
      }

      if (abs(x2 - x1) > 10 && abs(y2 - y1) > 10) {
        eucl_dist <- sqrt((x2 - x1)^2 + (y2 - y1)^2) + 1
        if (eucl_dist > bar_length_px) {
          bar_length_px <- as.integer(round(eucl_dist))
        }
      }

      px_per_mm <- bar_length_px / scale_mm
      val <- px_per_mm
      scale_info <- list(
        min_x = min_x,
        max_x = max_x,
        min_y = min_y,
        max_y = max_y,
        x1 = x1,
        y1 = y1,
        x2 = x2,
        y2 = y2,
        endpoints = list(x1 = x1, y1 = y1, x2 = x2, y2 = y2),
        bar_length_px = bar_length_px,
        scale_mm = scale_mm,
        px_per_mm = px_per_mm,
        pixels_per_mm = px_per_mm,
        confidence = 1.0
      )
      attr(val, "scale_info") <- scale_info
      return(val)
    },
    error = function(e) {
      warning(paste("Scale detection error:", e$message))
      invisible(NULL)
    }
  )
}
