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

.extract_image_contour_and_scale <- function(image_file, threshold, blur_size, blur_sigma,
                                             detect_scale, pixels_per_mm, scale_mm,
                                             min_points, min_area, downsample_threshold,
                                             threshold_method = "adaptive",
                                             threshold_factor = 1.0,
                                             manual_scale = NULL,
                                             warn_on_missing_scale = TRUE,
                                             keep_raw_image = FALSE) {
  img_raw <- read_and_process_image(image_file)

  current_pixels_per_mm <- pixels_per_mm
  scale_info <- NULL
  is_auto_scale <- isTRUE(detect_scale) || (is.character(detect_scale) && length(detect_scale) == 1 && tolower(detect_scale) == "auto")
  if (!is.null(manual_scale)) {
    current_pixels_per_mm <- as.numeric(manual_scale)
    scale_info <- attr(manual_scale, "scale_info")
  } else if (is_auto_scale) {
    img_grey <- if (EBImage::colorMode(img_raw) != EBImage::Grayscale) {
      EBImage::channel(img_raw, "grey")
    } else {
      img_raw
    }
    thresh_scale <- EBImage::otsu(img_grey)
    binary_scale <- img_grey > thresh_scale
    detected_px_per_mm <- find_scale_bar_length(binary_scale, scale_mm = scale_mm, raw_image = img_raw)
    if (!is.null(detected_px_per_mm)) {
      current_pixels_per_mm <- detected_px_per_mm
      scale_info <- attr(detected_px_per_mm, "scale_info")
    } else if (warn_on_missing_scale) {
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

  if (!keep_raw_image) {
    rm(img_raw)
    img_raw <- NULL
  }

  list(
    img_raw = img_raw,
    contour = contorno_raw,
    pixels_per_mm = current_pixels_per_mm,
    scale_info = scale_info
  )
}

.run_batch_tasks <- function(task_data, worker_fn, use_parallel = FALSE, actual_cores = 1L) {
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
    parallel::clusterExport(
      cl,
      c(".process_single_image_worker", ".extract_image_contour_and_scale",
        ".render_zonal_and_wavelet_plots", ".resolve_contour_and_start",
        "find_scale_bar_length", "calculate_adaptive_threshold", "read_and_process_image"),
      envir = environment(process_images)
    )

    res <- parallel::parLapply(cl, task_data, worker_fn)
    try(parallel::stopCluster(cl), silent = TRUE)
    cl <- NULL
    return(res)
  } else {
    n_tasks <- length(task_data)
    pb <- utils::txtProgressBar(max = n_tasks, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
    res <- vector("list", n_tasks)
    for (i in seq_len(n_tasks)) {
      utils::setTxtProgressBar(pb, i)
      res[[i]] <- worker_fn(task_data[[i]])
    }
    return(res)
  }
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
                                        threshold_factor = 1.0,
                                        manual_scale = NULL) {
  tryCatch({
    extracted <- .extract_image_contour_and_scale(
      image_file = image_file, threshold = threshold, blur_size = blur_size,
      blur_sigma = blur_sigma, detect_scale = detect_scale,
      pixels_per_mm = pixels_per_mm, scale_mm = scale_mm,
      min_points = min_points, min_area = min_area,
      downsample_threshold = downsample_threshold,
      threshold_method = threshold_method, threshold_factor = threshold_factor,
      manual_scale = manual_scale, warn_on_missing_scale = TRUE,
      keep_raw_image = testing
    )

    contorno_raw <- extracted$contour
    current_pixels_per_mm <- extracted$pixels_per_mm
    scale_info <- extracted$scale_info
    img_raw <- extracted$img_raw
    rm(extracted)

    if (is.null(contorno_raw)) {
      warning(paste("No suitable contour found for image:", basename(image_file)))
      rm(img_raw)
      return(NULL)
    }

    coef_e <- if (ef) .efourier(contorno_raw, n_harmonics, norm = TRUE) else NULL
    distances <- calculate_distances(contorno_raw, n_points = n_points, start_point = start_point)
    wavelet_results <- if (wavelets) calculate_wavelets_analysis(distances, n_scales = floor(log2(n_points)), detail = TRUE) else NULL

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
#' @param detect_scale Strategy for calibrating the scale bar in each image. Accepts one of three options:
#'   \itemize{
#'     \item \code{"none"} (default): Scale bar detection is disabled (uses
#'           \code{pixels_per_mm} if provided, or outputs measurements in pixels).
#'     \item \code{"auto"}: Automatically detects the calibration scale bar
#'           in each image and measures its exact physical pixel span via \code{find_scale_bar_length()}.
#'     \item \code{"manual"}: Interactive mouse calibration mode. Sequentially displays each image
#'           so the user can click the start and end points of the calibration scale bar with the mouse
#'           via \code{\link{select_scale_bar}} (without snapping). Useful as a fallback when automatic
#'           detection fails on non-standard or low-contrast rulers.
#'   }
#'   (Logical \code{FALSE} and \code{TRUE} are also accepted as aliases for \code{"none"} and \code{"auto"}.)
#' @param blur_size Size of the Gaussian brush for smoothing in pixels. Can be an odd integer,
#'                  or \code{"auto"} (default) to dynamically calculate \eqn{\min(\text{width},\text{height}) \times 0.01}.
#' @param blur_sigma Standard deviation of the Gaussian filter brush. Can be a numeric value,
#'                   or \code{"auto"} (default) to calculate \eqn{\text{blur\_size} / 3.5}.
#' @param min_points Minimum number of points required for a valid contour (default: 200).
#' @param min_area Minimum area required for a valid contour. Can be \code{"auto"} (default)
#'                 to dynamically adapt the threshold to 0.5% of total image area
#'                 (bounded between 500 and 5000 \eqn{\text{px}^2}), or a positive numeric value in pixels.
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
#' @param start_points Strategy or coordinates for selecting the contour starting point (rostrum)
#'   on each image. Accepts one of three options:
#'   \itemize{
#'     \item \code{"auto"} (default): Automatically selects the rightmost contour
#'           point (\eqn{\max(X)}) as the anatomical rostrum on every image.
#'     \item \code{"manual"}: Interactive mouse selection mode. Sequentially displays
#'           each image with its extracted contour and the automatic rostrum candidate (\code{o})
#'           overlaid, allowing the user to click the desired starting point with the mouse via
#'           \code{\link{select_starting_point}} (or right-click / press \code{Esc} to keep \code{"auto"}).
#'     \item A 2-column numeric \code{matrix} (or \code{data.frame}) of \code{(X, Y)} starting point
#'           coordinates:
#'           \itemize{
#'             \item Must contain one row per image in \code{folder} (\code{nrow(start_points)} must
#'                   equal the total number of images in the folder).
#'             \item Cannot be combined with \code{subfolder = TRUE}; folders with custom starting point
#'                   matrices must be processed individually (\code{subfolder = FALSE}).
#'           }
#'   }
#' @details
#' \subsection{Scale Bar Calibration (\code{detect_scale})}{
#' Physical calibration (converting pixels to millimeters) is controlled via \code{detect_scale}
#' and \code{scale_mm}:
#' \enumerate{
#'   \item \strong{Disabled (\code{detect_scale = "none"}, default)}:
#'         Uses the fixed \code{pixels_per_mm} argument if supplied, or reports morphometrics in \code{px}.
#'   \item \strong{Automatic detection (\code{detect_scale = "auto"})}:
#'         Automatically segments the calibration bar (supporting white, dark-on-plaque, and pure red
#'         bars in horizontal or vertical orientation) and measures its exact tip-to-tip pixel span,
#'         marking both physical endpoints (\code{x}) on the diagnostic \code{zonal_*.jpg} image.
#'   \item \strong{Interactive manual selection (\code{detect_scale = "manual"})}:
#'         Opens each image sequentially and prompts the user to click the two endpoints (start and end)
#'         of the scale bar with the mouse via \code{\link{select_scale_bar}} (no contour snapping).
#'         Computes the exact Euclidean pixel distance between the two clicks and overlays the selected
#'         endpoints on the diagnostic plot.
#' }
#' }
#' \subsection{Starting Point (Rostrum) Selection (\code{start_points})}{
#' Contour descriptors (Wavelets, Polar/Perimeter distances, and CSS) require a consistent
#' anatomical starting point (\code{0.0}), which in the AFORO standard corresponds to the
#' \strong{rostrum}. The \code{start_points} argument controls how this starting point is
#' determined across the folder:
#' \enumerate{
#'   \item \strong{Automatic mode (\code{start_points = "auto"}, default)}:
#'         Automatically locates the rightmost contour vertex (\eqn{\max(X)}) on each image.
#'         Compatible with both single-folder (\code{subfolder = FALSE}) and recursive
#'         multi-folder (\code{subfolder = TRUE}) batch processing.
#'   \item \strong{Interactive mouse mode (\code{start_points = "manual"})}:
#'         Before running the batch analysis, sequentially opens each image in the R graphics
#'         device showing the extracted contour and the automatic rostrum candidate (\code{o}).
#'         Left-clicking sets and snaps the starting point to the nearest contour vertex;
#'         right-clicking or pressing \code{Esc} retains the automatic candidate for that image.
#'   \item \strong{Custom coordinate matrix (\code{start_points = matrix_or_df})}:
#'         Accepts a 2-column numeric \code{matrix} or \code{data.frame} of \code{(X, Y)} coordinates
#'         (or a \code{data.frame} with columns \code{Image}, \code{X}, \code{Y}). For each image,
#'         the contour vertex closest to the supplied \code{(X, Y)} coordinate is selected as
#'         the starting point. To prevent ambiguity across directories, \code{nrow(start_points)}
#'         must strictly equal the total number of images in \code{folder}, and recursive processing
#'         (\code{subfolder = TRUE}) is disabled when a custom matrix is supplied (process one
#'         folder at a time with \code{subfolder = FALSE}).
#' }
#' }
#' @return Invisible NULL. The function is used for its side effects (creating files).
#' @export
#' @examples
#' \dontrun{
#' # 1. Default automatic rostrum detection and automatic scale bar calibration
#' process_images(folder = "path/to/images", detect_scale = "auto", start_points = "auto")
#'
#' # 2. Interactive mouse selection of the scale bar and/or rostrum on each image
#' process_images(folder = "path/to/images", detect_scale = "manual", start_points = "manual")
#'
#' # 3. Custom N x 2 matrix of starting points (one row per image, subfolder = FALSE)
#' rostrum_mat <- matrix(c(1850, 920,
#'                         1790, 880,
#'                         1910, 950), ncol = 2, byrow = TRUE)
#' colnames(rostrum_mat) <- c("X", "Y")
#' process_images(folder = "path/to/3_images", subfolder = FALSE, start_points = rostrum_mat)
#'
#' # 4. Process in parallel across 4 cores with polar coordinates only
#' process_images(folder = "path/to/images", cores = 4, coordinates = "polar")
#' }
process_images <- function(folder, subfolder = FALSE, threshold = NULL,
                           wavelets = TRUE, ef = TRUE, testing = TRUE,
                           coordinates = c("both", "polar", "cartesian"), save = TRUE,
                           pixels_per_mm = NULL, detect_scale = "none",
                           blur_size = "auto", blur_sigma = "auto", min_points = 200,
                           min_area = "auto", n_harmonics = 32, n_points = 512,
                           scale_mm = 1, procrustes = FALSE,
                           cores = NULL,
                           plot_width = 1200, plot_height = 800, plot_res = 150,
                           downsample_threshold = FALSE,
                           threshold_method = c("adaptive", "otsu", "fixed"),
                           threshold_factor = 1.0,
                           start_points = "auto") {
  coordinates <- match.arg(coordinates)
  threshold_method <- match.arg(threshold_method)

  if (!is.numeric(threshold_factor) || length(threshold_factor) != 1 || is.na(threshold_factor) || is.infinite(threshold_factor) || threshold_factor <= 0) {
    stop("threshold_factor must be a positive number")
  }

  if (!is.character(folder) || length(folder) != 1) stop("folder must be a single character string")
  if (!dir.exists(folder)) stop(paste("Folder does not exist:", folder))
  is_scalar_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!is_scalar_flag(wavelets) || !is_scalar_flag(ef) || !is_scalar_flag(testing) ||
      !is_scalar_flag(save) || !is_scalar_flag(procrustes) || !is_scalar_flag(subfolder) ||
      !is_scalar_flag(downsample_threshold)) {
    stop("wavelets, ef, testing, save, procrustes, subfolder, and downsample_threshold must be logical values")
  }
  is_valid_ds_logical <- is.logical(detect_scale) && length(detect_scale) == 1 && !is.na(detect_scale)
  is_valid_ds_char <- is.character(detect_scale) && length(detect_scale) == 1 && tolower(detect_scale) %in% c("none", "auto", "manual", "interactive")
  if (!is_valid_ds_logical && !is_valid_ds_char) {
    stop("detect_scale must be 'none', 'auto', or 'manual'")
  }
  if (!is.null(start_points)) {
    is_valid_mode <- is.character(start_points) && length(start_points) == 1 && tolower(start_points) %in% c("auto", "manual", "interactive")
    is_valid_struct <- is.matrix(start_points) || is.data.frame(start_points)
    if (!is_valid_mode && !is_valid_struct) {
      stop("start_points must be 'auto', 'manual', or a 2-column matrix/data.frame of (X, Y) starting coordinates")
    }
    if (is_valid_struct) {
      if (subfolder) {
        stop("subfolder = TRUE is not allowed when providing a custom start_points matrix. Please process one folder at a time (subfolder = FALSE).")
      }
      if (ncol(start_points) < 2) {
        stop("start_points matrix or data.frame must have at least 2 columns (X, Y)")
      }
    }
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
  is_valid_area <- (identical(min_area, "auto") || is.null(min_area)) ||
                   (is.numeric(min_area) && length(min_area) == 1 && !is.na(min_area) && !is.infinite(min_area) && min_area >= 0)
  if (!is_valid_area) {
    stop("min_area must be 'auto', NULL, or a numeric value >= 0")
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

    # Resolve manual scale bar selection if detect_scale = "manual"
    chosen_scales <- vector("list", n_imgs)
    names(chosen_scales) <- basename(imag)
    is_manual_scale <- is.character(detect_scale) && length(detect_scale) == 1 && tolower(detect_scale) %in% c("manual", "interactive")
    is_auto_scale   <- isTRUE(detect_scale) || (is.character(detect_scale) && length(detect_scale) == 1 && tolower(detect_scale) == "auto")

    if (is_manual_scale) {
      message("\nManual scale bar calibration mode: Click start and end points of the scale bar on each image...")
      for (i in seq_along(imag)) {
        message(paste("Select scale bar endpoints for:", basename(imag[i])))
        img_preview <- read_and_process_image(imag[i])
        chosen_scales[[i]] <- select_scale_bar(img_preview, scale_mm = scale_mm)
        rm(img_preview)
      }
    }

    # Resolve starting points ("auto", "manual", or matrix/data.frame)
    chosen_rostrums <- vector("list", n_imgs)
    names(chosen_rostrums) <- basename(imag)

    is_manual_mode <- is.character(start_points) && length(start_points) == 1 && tolower(start_points) %in% c("manual", "interactive")

    if (is_manual_mode) {
      message("\nManual rostrum selection mode: Click on each otolith to set the starting point...")
      for (i in seq_along(imag)) {
        message(paste("Select rostrum for:", basename(imag[i])))
        img_preview <- read_and_process_image(imag[i])
        bin_preview <- preprocess_image(img_preview, threshold = threshold, blur_size = blur_size,
                                        blur_sigma = blur_sigma, downsample_threshold = TRUE,
                                        threshold_method = threshold_method, threshold_factor = threshold_factor)
        cnt_preview <- extract_contour(bin_preview, min_points = min_points, min_area = min_area)
        chosen_rostrums[[i]] <- select_starting_point(img_preview, cnt_preview)
        rm(img_preview, bin_preview, cnt_preview)
      }
    } else if (!is.null(start_points) && (is.matrix(start_points) || is.data.frame(start_points))) {
      if (nrow(start_points) != n_imgs) {
        stop(sprintf(
          "The number of rows in start_points (%d) must equal the total number of images in the folder (%d).",
          nrow(start_points), n_imgs
        ))
      }
      if (is.data.frame(start_points) && all(c("Image", "X", "Y") %in% names(start_points))) {
        for (r in seq_len(n_imgs)) {
          img_name <- as.character(start_points$Image[r])
          if (img_name %in% names(chosen_rostrums)) {
            chosen_rostrums[[img_name]] <- c(as.numeric(start_points$X[r]), as.numeric(start_points$Y[r]))
          } else {
            chosen_rostrums[[r]] <- c(as.numeric(start_points$X[r]), as.numeric(start_points$Y[r]))
          }
        }
      } else {
        rn <- rownames(start_points)
        if (!is.null(rn) && all(rn %in% names(chosen_rostrums))) {
          for (r in seq_len(n_imgs)) {
            chosen_rostrums[[rn[r]]] <- c(as.numeric(start_points[r, 1]), as.numeric(start_points[r, 2]))
          }
        } else {
          for (r in seq_len(n_imgs)) {
            chosen_rostrums[[r]] <- c(as.numeric(start_points[r, 1]), as.numeric(start_points[r, 2]))
          }
        }
      }
    }

    # Determine parallelization
    use_parallel <- !is.null(cores) && cores > 1 && n_imgs > 1
    actual_cores <- if (use_parallel) min(as.integer(cores), n_imgs, parallel::detectCores()) else 1L
    chk_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
    if (use_parallel && nzchar(chk_limit) && chk_limit != "false") {
      actual_cores <- min(actual_cores, 2L)
    }

    task_data <- lapply(imag, function(f) {
      list(
        file = f,
        start_point = chosen_rostrums[[basename(f)]],
        manual_scale = chosen_scales[[basename(f)]]
      )
    })

    if (!procrustes) {
      message(paste("\nProcessing", n_imgs, "images in:", basename(folder),
                    if (use_parallel) paste("using", actual_cores, "cores") else "sequentially"))

      processed_list <- .run_batch_tasks(
        task_data = task_data,
        worker_fn = function(item) {
          .process_single_image_worker(
            image_file = item$file, threshold = threshold, blur_size = blur_size,
            blur_sigma = blur_sigma, detect_scale = is_auto_scale,
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
            threshold_factor = threshold_factor,
            manual_scale = item$manual_scale
          )
        },
        use_parallel = use_parallel,
        actual_cores = actual_cores
      )

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
          extracted <- .extract_image_contour_and_scale(
            image_file = item$file, threshold = threshold, blur_size = blur_size,
            blur_sigma = blur_sigma, detect_scale = is_auto_scale,
            pixels_per_mm = pixels_per_mm, scale_mm = scale_mm,
            min_points = min_points, min_area = min_area,
            downsample_threshold = downsample_threshold,
            threshold_method = threshold_method, threshold_factor = threshold_factor,
            manual_scale = item$manual_scale, warn_on_missing_scale = FALSE,
            keep_raw_image = FALSE
          )
          list(contour = extracted$contour, px = extracted$pixels_per_mm, start_point = item$start_point)
        }, error = function(e) NULL)
      }

      extr_list <- .run_batch_tasks(
        task_data = task_data,
        worker_fn = extract_worker,
        use_parallel = use_parallel,
        actual_cores = actual_cores
      )

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
      on.exit(try(close(pb2), silent = TRUE), add = TRUE)
      for (idx_counter in seq_along(valid_indices)) {
        utils::setTxtProgressBar(pb2, idx_counter)
        i <- valid_indices[idx_counter]
        c_proc <- contours_to_process[[i]]
        c_raw  <- contours_raw[[i]]
        cur_px <- px_list[[i]]
        st_pt  <- st_list[[i]]

        coef_e <- if (ef) .efourier(c_proc, n_harmonics, norm = FALSE) else NULL
        distances <- calculate_distances(c_proc, n_points = n_points, start_point = st_pt)
        wavelet_results <- if (wavelets) calculate_wavelets_analysis(distances, n_scales = floor(log2(n_points)), detail = TRUE) else NULL
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
