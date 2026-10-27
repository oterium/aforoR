#' Internal Shared Renderer for Polar and Perimeter Visualizations
#'
#' @param image An \code{Image} object, file path string, or NULL.
#' @param distances Distance list from \code{\link{calculate_distances}}.
#' @param wavelets Wavelet analysis list from \code{\link{calculate_wavelets_analysis}} or NULL.
#' @param image_name Base name of the image file.
#' @param output_dir Output directory path.
#' @param wavelet_scale Wavelet scale to plot.
#' @param plot_width Width in pixels.
#' @param plot_height Height in pixels.
#' @param plot_res Resolution in dpi.
#' @param binary_image Deprecated alias for \code{image}.
#' @param scale_info Optional scale calibration info list.
#' @param mode Either \code{"polar"} or \code{"perimeter"}.
#' @return Invisible NULL.
#' @keywords internal
#' @noRd
.render_zonal_and_wavelet_plots <- function(image = NULL, distances, wavelets,
                                            image_name, output_dir,
                                            wavelet_scale = 5,
                                            plot_width = 1200, plot_height = 800,
                                            plot_res = 150, binary_image = NULL,
                                            scale_info = NULL,
                                            mode = c("polar", "perimeter")) {
  mode <- match.arg(mode)

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

  if (!is.numeric(wavelet_scale) || length(wavelet_scale) != 1 || is.na(wavelet_scale) || is.infinite(wavelet_scale) || wavelet_scale <= 0) {
    stop("wavelet_scale must be a positive integer")
  }

  if (!is.numeric(plot_width) || length(plot_width) != 1 || is.na(plot_width) || is.infinite(plot_width) || plot_width <= 0) {
    stop("plot_width must be a positive number")
  }

  if (!is.numeric(plot_height) || length(plot_height) != 1 || is.na(plot_height) || is.infinite(plot_height) || plot_height <= 0) {
    stop("plot_height must be a positive number")
  }

  if (!is.numeric(plot_res) || length(plot_res) != 1 || is.na(plot_res) || is.infinite(plot_res) || plot_res <= 0) {
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

  err_prefix <- if (mode == "polar") "Error creating visualization:" else "Error creating perimeter visualization:"

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
      img_obj <- NULL
      if (!is.null(image)) {
        if (is.character(image) && length(image) == 1 && file.exists(image)) {
          img_obj <- read_and_process_image(image)
          plot(img_obj)
          has_img <- TRUE
        } else if (inherits(image, "Image")) {
          img_obj <- image
          plot(img_obj)
          has_img <- TRUE
        }
      }

      if (mode == "polar") {
        cx <- if (!is.null(distances$reordered_coords$x)) mean(distances$reordered_coords$x) else if (!is.null(distances$polar$centroid)) distances$polar$centroid[1] else mean(distances$polar$coord[, 1])
        cy <- if (!is.null(distances$reordered_coords$y)) mean(distances$reordered_coords$y) else if (!is.null(distances$polar$centroid)) distances$polar$centroid[2] else mean(distances$polar$coord[, 2])
        pts_x <- distances$polar$coord[, 1] + cx
        pts_y <- distances$polar$coord[, 2] + cy
        rx <- if (!is.null(distances$reordered_coords$x)) distances$reordered_coords$x[1] else if (!is.null(distances$polar$rostrum)) distances$polar$rostrum[1] else cx
        ry <- if (!is.null(distances$reordered_coords$y)) distances$reordered_coords$y[1] else if (!is.null(distances$polar$rostrum)) distances$polar$rostrum[2] else cy
        n_points <- length(distances$polar$radii)
        main_title <- paste("Zonal Analysis -", image_name)
        wav_mat <- if (!is.null(wavelets)) wavelets$polar else NULL
        wav_title_prefix <- "Wavelet Scale"
      } else {
        cx <- if (!is.null(distances$reordered_coords$x)) mean(distances$reordered_coords$x) else if (!is.null(distances$perimeter$centroid)) distances$perimeter$centroid[1] else mean(distances$perimeter$coords[, 1])
        cy <- if (!is.null(distances$reordered_coords$y)) mean(distances$reordered_coords$y) else if (!is.null(distances$perimeter$centroid)) distances$perimeter$centroid[2] else mean(distances$perimeter$coords[, 2])
        pts_x <- distances$perimeter$coords[, 1]
        pts_y <- distances$perimeter$coords[, 2]
        rx <- if (!is.null(distances$reordered_coords$x)) distances$reordered_coords$x[1] else if (!is.null(distances$perimeter$rostrum)) distances$perimeter$rostrum[1] else cx
        ry <- if (!is.null(distances$reordered_coords$y)) distances$reordered_coords$y[1] else if (!is.null(distances$perimeter$rostrum)) distances$perimeter$rostrum[2] else cy
        n_points <- length(distances$perimeter$dist)
        main_title <- paste("Perimeter Analysis -", image_name)
        wav_mat <- if (!is.null(wavelets)) wavelets$perimeter else NULL
        wav_title_prefix <- "Perimeter Wavelet Scale"
      }

      if (!has_img) {
        xlims <- range(pts_x)
        ylims <- range(pts_y)
        xpad <- if (is.na(diff(xlims)) || diff(xlims) == 0) 1 else diff(xlims) * 0.1
        ypad <- if (is.na(diff(ylims)) || diff(ylims) == 0) 1 else diff(ylims) * 0.1
        plot(NA, xlim = c(xlims[1] - xpad, xlims[2] + xpad),
             ylim = c(ylims[1] - ypad, ylims[2] + ypad),
             xlab = "X", ylab = "Y", main = main_title)
      }

      points(cx, cy, pch = 3, col = "red", lwd = 2.5, cex = 1.6)
      points(rx, ry, pch = 16, col = "red", cex = 1.6)

      b <- round(seq(0, n_points, length.out = 5))
      zones <- lapply(1:4, function(idx) if (b[idx + 1] > b[idx]) (b[idx] + 1L):b[idx + 1] else integer(0))
      colors <- c("#00E5FF", "#FFEA00", "#FF007F", "#00E676")

      for (i in seq_along(zones)) {
        if (length(zones[[i]]) > 0) {
          lines(
            pts_x[zones[[i]]],
            pts_y[zones[[i]]],
            col = colors[i], lwd = 2.5
          )
        }
      }

      if (!is.null(scale_info) && is.list(scale_info) && !is.null(scale_info$min_x)) {
        pad <- 15
        graphics::rect(scale_info$min_x - pad, scale_info$min_y - pad,
                       scale_info$max_x + pad, scale_info$max_y + pad,
                       border = "#00FF66", lwd = 2.5)

        # Visibly mark exact scale bar physical endpoints
        ep <- scale_info$endpoints
        ep_x1 <- if (!is.null(ep$x1)) ep$x1 else scale_info$x1
        ep_y1 <- if (!is.null(ep$y1)) ep$y1 else scale_info$y1
        ep_x2 <- if (!is.null(ep$x2)) ep$x2 else scale_info$x2
        ep_y2 <- if (!is.null(ep$y2)) ep$y2 else scale_info$y2
        if (!is.null(ep_x1) && !is.null(ep_y1) && !is.null(ep_x2) && !is.null(ep_y2)) {
          graphics::lines(c(ep_x1, ep_x2), c(ep_y1, ep_y2), col = "#00FF66", lwd = 3)
          graphics::points(c(ep_x1, ep_x2), c(ep_y1, ep_y2), pch = 4, col = "#FF0033", lwd = 3, cex = 1.6)
          graphics::points(c(ep_x1, ep_x2), c(ep_y1, ep_y2), pch = 1, col = "#00FF66", lwd = 2, cex = 1.8)
        }

        scale_label <- sprintf("Scale: %g mm (%.0f px | %.1f px/mm)",
                               scale_info$scale_mm, scale_info$bar_length_px, scale_info$px_per_mm)
        img_w <- if (has_img && !is.null(dim(img_obj))) dim(img_obj)[1] else max(scale_info$max_x + 100, 1000)
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
        leg_labels <- c(leg_labels, sprintf("Scale bar (%g mm)", scale_info$scale_mm), "Scale endpoints (x)")
        leg_cols   <- c(leg_cols, "#00FF66", "#FF0033")
        leg_pch    <- c(leg_pch, NA, 4)
        leg_lty    <- c(leg_lty, 1, NA)
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

      if (!is.null(wav_mat) && nrow(wav_mat) > 0) {
        scale_to_plot <- min(as.integer(wavelet_scale), nrow(wav_mat))
        wav_file <- file.path(
          output_dir,
          paste0("wavelet_", tools::file_path_sans_ext(image_name), ".jpg")
        )
        jpeg(filename = wav_file, width = 800, height = 500, res = 100, quality = 85)
        cur_dev <- grDevices::dev.cur()
        plot(wav_mat[scale_to_plot, ],
          main = paste(wav_title_prefix, scale_to_plot, "-", image_name),
          xlab = "Points", ylab = "Coefficient", type = "l", col = "#00E5FF", lwd = 2
        )

        for (i in seq_along(zones)) {
          if (length(zones[[i]]) > 0) {
            polygon(
              x = c(min(zones[[i]]), zones[[i]], max(zones[[i]])),
              y = c(
                min(wav_mat[scale_to_plot, ]),
                wav_mat[scale_to_plot, zones[[i]]],
                min(wav_mat[scale_to_plot, ])
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
      invisible(NULL)
    },
    error = function(e) {
      stop(paste(err_prefix, e$message))
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
  .render_zonal_and_wavelet_plots(
    image = image, distances = distances, wavelets = wavelets,
    image_name = image_name, output_dir = output_dir,
    wavelet_scale = wavelet_scale,
    plot_width = plot_width, plot_height = plot_height, plot_res = plot_res,
    binary_image = binary_image, scale_info = scale_info,
    mode = "polar"
  )
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
  .render_zonal_and_wavelet_plots(
    image = image, distances = distances, wavelets = wavelets,
    image_name = image_name, output_dir = output_dir,
    wavelet_scale = wavelet_scale,
    plot_width = plot_width, plot_height = plot_height, plot_res = plot_res,
    binary_image = binary_image, scale_info = scale_info,
    mode = "perimeter"
  )
}
