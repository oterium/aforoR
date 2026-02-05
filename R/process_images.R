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
#' @param blur_size Size of the Gaussian brush for smoothing (default: 31).
#' @param blur_sigma Sigma of the Gaussian brush for smoothing (default: 9).
#' @return A \code{pixset} (EBImage) representing the binarized image.
#' @export
#' @examples
#' # Use the sample otolith image included in the package
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' processed_img <- preprocess_image(image_path)
#' # display(processed_img) # if using EBImage
preprocess_image <- function(image_path, threshold = NULL, blur_size = 31,
                             blur_sigma = 9) {
  # Input validation
  is_image <- inherits(image_path, "Image")

  if (!is_image && (!is.character(image_path) || length(image_path) != 1)) {
    stop(paste(
      "image_path must be a single character string",
      "or an Image object"
    ))
  }

  if (!is_image && !file.exists(image_path)) {
    stop(paste("Image file does not exist:", image_path))
  }

  if (!is.null(threshold) &&
    (!is.numeric(threshold) || length(threshold) != 1)) {
    stop("threshold must be a single numeric value or NULL")
  }

  tryCatch(
    {
      # Read and process image
      if (is_image) {
        foto <- image_path
      } else {
        foto <- EBImage::readImage(image_path)
      }

      foto2 <- EBImage::channel(foto, "grey")
      foto2 <- EBImage::filter2(
        foto2,
        EBImage::makeBrush(
          size = blur_size, shape = "gaussian",
          sigma = blur_sigma
        )
      )

      if (is.null(threshold)) {
        foto2 <- foto2 > EBImage::otsu(foto2)
      } else {
        foto2 <- foto2 > threshold
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
#' meets specific size criteria (points > 200 and area > 5000), typically
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
  # Input validation
  if (is.null(binary_image)) {
    stop("binary_image cannot be NULL")
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
        current_area <- Momocs::coo_area(cont[[x]])
        if (length(cont[[x]])[1] > min_points && current_area > min_area) {
          if (current_area > max_area) {
            max_area <- current_area
            contorno <- cont[[x]]
          }
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

#' Calculate Distance Measures from Contour
#'
#' Computes two types of distance measures from an otolith contour: polar
#' distances (radial distances from the centroid) and perimeter-based distances.
#' These measures are crucial for subsequent wavelet analysis.
#'
#' @param contour A matrix of coordinates from \code{extract_contour()}.
#' @param n_points Target number of points for resampling (default: 512).
#' @return A list containing:
#'   \itemize{
#'     \item \code{polar}: Radii, coordinates, and normalized radii.
#'     \item \code{perimeter}: Distances, coordinates, and normalized distances.
#'     \item \code{reordered_coords}: The contour coordinates reordered from the
#'           rightmost point.
#'   }
#' @export
#' @examples
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' binary_img <- preprocess_image(image_path)
#' contour <- extract_contour(binary_img)
#' dists <- calculate_distances(contour, n_points = 512)
#' plot(dists$polar$normalized, type = "l")
calculate_distances <- function(contour, n_points = 512) {
  # Input validation
  if (is.matrix(contour)) {
    contour <- as.data.frame(contour)
    if (ncol(contour) == 2 &&
      (is.null(colnames(contour)) ||
        !all(c("X", "Y") %in% colnames(contour)))) {
      colnames(contour) <- c("X", "Y")
    }
  }

  if (is.null(contour) || !is.data.frame(contour)) {
    stop("contour must be a data.frame or matrix with X and Y columns")
  }

  if (!"X" %in% names(contour) || !"Y" %in% names(contour)) {
    stop("contour must have columns named 'X' and 'Y'")
  }

  if (n_points <= 0 || !is.numeric(n_points)) {
    stop("n_points must be a positive integer")
  }

  # Find rightmost point more safely
  x_max_idx <- which.max(contour$X)
  x_max <- contour$X[x_max_idx]
  y_max <- contour$Y[x_max_idx]
  punto <- list(x = x_max, y = y_max)

  # More robust contour reordering
  # Find points where X coordinate matches the maximum and X > mean(X)
  x_mean <- mean(contour$X)
  rightmost_candidates <- which(contour$X == punto$x & contour$X > x_mean)

  if (length(rightmost_candidates) == 0) {
    # Fallback: use the point with maximum X coordinate
    start_idx <- x_max_idx
  } else {
    start_idx <- rightmost_candidates[1]
  }

  # Safely reorder the contour
  n_contour <- nrow(contour)
  if (start_idx == 1) {
    # No reordering needed
    x <- contour$X
    y <- contour$Y
  } else {
    # Reorder starting from start_idx
    indices <- c(start_idx:n_contour, 1:(start_idx - 1))
    x <- contour$X[indices]
    y <- contour$Y[indices]
  }

  # Calculate distances with error handling
  tryCatch(
    {
      dist <- regularradius(x, y, n = n_points)
      dist2 <- dper(x, y, n_points)
      dist.norm <- dist$radii / max(dist$radii)
      dist.norm2 <- dist2$dist / max(dist2$dist)

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
#' @param detail Logical; if TRUE, returns detail coefficients instead of
#'               approximations.
#' @return A list with wavelet coefficients for 'polar' and 'perimeter' components.
#' @export
#' @examples
#' # Example using existing data for wavelet analysis
#' data(Aphanopus_W5)
#' # This data frame already contains coefficients at scale 5.
#' @seealso \code{\link{fwaveletspl_3}}
calculate_wavelets_analysis <- function(distances, n_scales = 9,
                                        detail = FALSE) {
  # Input validation
  if (!is.list(distances)) {
    stop("distances must be a list")
  }

  if (!all(c("polar", "perimeter") %in% names(distances))) {
    stop("distances must contain 'polar' and 'perimeter' components")
  }

  if (!is.numeric(n_scales) || length(n_scales) != 1 || n_scales <= 0) {
    stop("n_scales must be a positive integer")
  }

  if (!is.logical(detail) || length(detail) != 1) {
    stop("detail must be a logical value")
  }

  # Check if normalized distances exist
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
#'
#' @param binary_image The binarized image from preprocessing.
#' @param distances The distance measures from calculate_distances().
#' @param wavelets The wavelet analysis from calculate_wavelets_analysis().
#' @param image_name The base name of the image file.
#' @param output_dir The directory to save visualization images.
#' @param wavelet_scale The wavelet scale to visualize (default: 5).
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics lines points polygon
#' @export
save_visualization <- function(binary_image, distances, wavelets, image_name,
                               output_dir, wavelet_scale = 5) {
  # Input validation
  if (is.null(binary_image)) {
    stop("binary_image cannot be NULL")
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

  if (!is.numeric(wavelet_scale) || length(wavelet_scale) != 1 || wavelet_scale <= 0) {
    stop("wavelet_scale must be a positive integer")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  tryCatch(
    {
      # Save zonal analysis plot
      jpeg(
        filename = file.path(output_dir, paste("zonal_", tools::file_path_sans_ext(image_name), ".jpg", sep = "")),
        res = 600, width = 600, height = 400, units = "mm"
      )
      plot(binary_image)
      points(mean(distances$reordered_coords$x),
        mean(distances$reordered_coords$y),
        pch = 16
      )
      points(distances$reordered_coords$x[1], distances$reordered_coords$y[1],
        col = 2, pch = 16
      )

      # Draw colored zones dynamically based on number of points
      n_points <- length(distances$polar$radii)
      z1 <- 1:floor(n_points * 0.25)
      z2 <- (max(z1) + 1):floor(n_points * 0.5)
      z3 <- (max(z2) + 1):floor(n_points * 0.75)
      z4 <- (max(z3) + 1):n_points

      zones <- list(z1, z2, z3, z4)
      colors <- c(4, 5, 6, 7)

      for (i in seq_along(zones)) {
        lines(
          distances$polar$coord[zones[[i]], 1] +
            mean(distances$reordered_coords$x),
          distances$polar$coord[zones[[i]], 2] +
            mean(distances$reordered_coords$y),
          col = colors[i], pch = 16, lwd = 2, type = "b"
        )
      }
      dev.off()

      # Save wavelet plot if available
      if (!is.null(wavelets)) {
        jpeg(filename = file.path(
          output_dir,
          paste("wavelet_", tools::file_path_sans_ext(image_name), ".jpg", sep = "")
        ))
        plot(wavelets$polar[wavelet_scale, ],
          main = paste("Wavelet", wavelet_scale, image_name),
          xlab = "", ylab = "", type = "l", col = 3, lwd = 2
        )

        # Add colored polygons for zones
        for (i in seq_along(zones)) {
          polygon(
            x = c(min(zones[[i]]), zones[[i]], max(zones[[i]])),
            y = c(
              min(wavelets$polar[wavelet_scale, ]),
              wavelets$polar[wavelet_scale, zones[[i]]],
              min(wavelets$polar[wavelet_scale, ])
            ),
            col = colors[i]
          )
        }
        dev.off()
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
#' @importFrom utils write.table
#' @export
save_analysis_results <- function(analysis_results, output_dir,
                                  result_type = "polar") {
  # Input validation
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
      # Filter out NULL results
      analysis_results <- analysis_results[!sapply(analysis_results, is.null)]
      if (length(analysis_results) == 0) {
        warning("No valid analysis results to save")
        return(invisible(NULL))
      }

      # Extract data for each type of measurement
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

      # Save distance data
      write_analysis_csv(
        distances, file.path(output_dir, "DistanciaEN.csv"),
        image_names
      )
      write_analysis_csv(
        distances_norm,
        file.path(output_dir, "Distancia_NormEN.csv"),
        image_names
      )

      # Save wavelet data (detect number of scales from result)
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

      # Save coordinates
      write.table(coordinates,
        file = file.path(output_dir, "Coords.csv"),
        dec = ".", sep = ";", col.names = FALSE
      )

      # Save elliptic coefficients if available
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

#' Helper function to write analysis data to CSV
#'
#' Internal helper function to write analysis data matrices to CSV files with proper formatting.
#'
#' @param data The data matrix or vector to write.
#' @param filename The name of the CSV file.
#' @param row_names The row names for the data.
#' @noRd
write_analysis_csv <- function(data, filename, row_names) {
  # Input validation
  if (is.null(data)) {
    stop("data cannot be NULL")
  }

  if (!is.character(filename) || length(filename) != 1) {
    stop("filename must be a single character string")
  }

  if (!is.character(row_names)) {
    stop("row_names must be a character vector")
  }

  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }

  if (nrow(data) != length(row_names)) {
    stop("Number of rows in data must match length of row_names")
  }

  rownames(data) <- row_names
  write.table(data, file = filename, dec = ".", sep = ";", col.names = FALSE)
}

#' Process Images in Folder
#'
#' This function processes images in a specified folder, applying various transformations and analyses.
#'
#' @param folder A string specifying the path to the folder containing images.
#' @param subfolder A logical value indicating whether the folder has subfolders.
#' @param threshold A numeric value for binarization. If not provided, the Otsu method is used.
#' @param wavelets A logical value indicating whether to obtain wavelets. Default is TRUE.
#' @param ef A logical value indicating whether to obtain elliptic Fourier descriptors. Default is TRUE.
#' @param testing A logical value indicating whether to save test images. Default is TRUE.
#' @param pseudolandmarks A string specifying the type of pseudolandmarks. Options are 'curvilinear', 'polar', 'both'. Default is 'both'.
#' @param save A logical value indicating whether to save the results as CSV files. Default is TRUE.
#' @param pixels_per_mm A numeric value specifying the scale (pixels per mm). Default is NULL.
#' @param detect_scale A logical value indicating whether to automatically detect the scale bar in the image. Default is FALSE.
#' @param blur_size Size of the Gaussian brush for smoothing (default: 31).
#' @param blur_sigma Sigma of the Gaussian brush for smoothing (default: 9).
#' @param min_points Minimum number of points required for a valid contour (default: 200).
#' @param min_area Minimum area required for a valid contour (default: 5000).
#' @param n_harmonics Number of elliptic Fourier harmonics to calculate (default: 32).
#' @param n_points Target number of points for resampling (default: 512).
#' @param scale_mm Total length in mm of the scale bar for detection (default: 1).
#' @return Invisible NULL. The function is used for its side effects (creating files).
#' @details
#' The function generates the following files in the `Polar` (and `Cartesian`) directories:
#' \itemize{
#'   \item \code{MorphometricsEN.csv}: Table with shape indices (Area, Perimeter, Circularity, etc.).
#'   \item \code{DistanciaEN.csv} and \code{Distancia_NormEN.csv}: Raw and normalized distances (radii or perimeter).
#'   \item \code{Wavelet_xEN.csv}: Wavelet coefficients for scales 1 to 9.
#'   \item \code{EllipticCoeEN.csv}: Elliptic Fourier coefficients (if \code{ef = TRUE}).
#'   \item \code{Coords.csv}: Contour coordinates.
#' }
#'
#' If \code{testing = TRUE}, diagnostic images are also saved:
#' \itemize{
#'   \item \code{zonal_[ImageName].jpg}: Visualization of the contour and analysis sectors.
#'   \item \code{wavelet_[ImageName].jpg}: Visualization of the wavelet transform.
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Process all images in the "images" folder
#' process_images(folder = "path/to/images", pixels_per_mm = 100)
#' }
process_images <- function(folder, subfolder = FALSE, threshold = NULL,
                           wavelets = TRUE, ef = TRUE, testing = TRUE,
                           pseudolandmarks = "both", save = TRUE,
                           pixels_per_mm = NULL, detect_scale = FALSE,
                           blur_size = 31, blur_sigma = 9, min_points = 200,
                           min_area = 5000, n_harmonics = 32, n_points = 512,
                           scale_mm = 1) {
  # Input validation
  if (!is.character(folder) || length(folder) != 1) {
    stop("folder must be a single character string")
  }

  if (!dir.exists(folder)) {
    stop(paste("Folder does not exist:", folder))
  }

  if (!is.logical(wavelets) || !is.logical(ef) || !is.logical(testing) ||
    !is.logical(save) || !is.logical(detect_scale)) {
    stop(paste(
      "wavelets, ef, testing, save, and detect_scale",
      "must be logical values"
    ))
  }

  if (!pseudolandmarks %in% c("curvilinear", "polar", "both")) {
    stop("pseudolandmarks must be one of: 'curvilinear', 'polar', 'both'")
  }

  if (!is.null(pixels_per_mm) && (!is.numeric(pixels_per_mm) || pixels_per_mm <= 0)) {
    stop("pixels_per_mm must be a positive numeric value")
  }

  tryCatch(
    {
      # Create processing directories if they don't exist
      polar_dir <- file.path(folder, "Polar")
      cartesian_dir <- file.path(folder, "Cartesian")

      if (!dir.exists(polar_dir)) {
        dir.create(polar_dir, recursive = TRUE)
      }
      if (!dir.exists(cartesian_dir)) {
        dir.create(cartesian_dir, recursive = TRUE)
      }

      # Get list of image files (JPEG, PNG, TIFF - case-insensitive)
      # EBImage::readImage supports these formats natively
      imag <- list.files(folder,
        pattern = "(?i)\\.(jpe?g|png|tiff?)$",
        full.names = TRUE
      )

      if (length(imag) == 0) {
        warning(paste(
          "No supported image files (.jpg, .jpeg, .png, .tif, .tiff)",
          "found in folder:", folder
        ))
        return(invisible(NULL))
      }

      pb <- utils::txtProgressBar(max = length(imag), style = 3)

      # Initialize results
      result <- list()
      result2 <- list()
      morpho_results <- list()

      for (i in seq_along(imag)) {
        utils::setTxtProgressBar(pb, i)

        tryCatch(
          {
            # Step 1: Read image ONCE
            img_raw <- EBImage::readImage(imag[i])

            # Step 2: Preprocess image
            binary_image <- preprocess_image(img_raw, threshold, blur_size = blur_size, blur_sigma = blur_sigma)

            # Scale detection logic
            current_pixels_per_mm <- pixels_per_mm
            if (detect_scale) {
              # Use simple Otsu without heavy smoothing for scale detection
              img_grey <- if (EBImage::colorMode(img_raw) != EBImage::Grayscale) {
                EBImage::channel(img_raw, "grey")
              } else {
                img_raw
              }
              thresh_scale <- EBImage::otsu(img_grey)
              binary_scale <- img_grey > thresh_scale

              detected_px_per_mm <- find_scale_bar_length(binary_scale,
                scale_mm = scale_mm
              )
              if (!is.null(detected_px_per_mm)) {
                current_pixels_per_mm <- detected_px_per_mm
              } else {
                warning(paste(
                  "Scale bar not found for image", basename(imag[i]),
                  ". Using unscaled values."
                ))
              }
            }

            # Step 3: Extract contour (with max area logic and params)
            contorno <- extract_contour(binary_image, min_points = min_points, min_area = min_area)

            if (!is.null(contorno)) {
              # Step 4: Calculate elliptic Fourier descriptors if requested
              coef_e <- NULL
              if (ef) {
                coef_e <- Momocs::efourier(contorno, n_harmonics)
              }

              # Step 5: Calculate distance measures
              distances <- calculate_distances(contorno, n_points = n_points)

              # Step 6: Calculate wavelets if requested
              wavelet_results <- NULL
              if (wavelets) {
                # log2(n_points) is the max scales
                n_scales_calc <- floor(log2(n_points))
                wavelet_results <- calculate_wavelets_analysis(
                  distances,
                  n_scales = n_scales_calc, detail = FALSE
                )
              }

              # Step 7: Calculate Morphometrics
              morpho <- calculate_morphometrics(contorno, current_pixels_per_mm)
              # Add image name to morpho results
              morpho$Image <- basename(imag[i])
              morpho_results[[i]] <- morpho

              if (save) {
                result[[i]] <- create_result_structure(
                  distances$polar, wavelet_results$polar, basename(imag[i]),
                  coef_e, "polar"
                )
                result2[[i]] <- create_result_structure(
                  distances$perimeter, wavelet_results$perimeter,
                  basename(imag[i]), coef_e, "perimeter"
                )
              }

              # Step 9: Save visualizations if requested
              if (testing) {
                save_visualization(
                  binary_image, distances, wavelet_results,
                  basename(imag[i]), polar_dir
                )
                save_visualization_perimeter(
                  binary_image, distances,
                  wavelet_results,
                  basename(imag[i]), cartesian_dir
                )
              }
            } else {
              warning(paste("No suitable contour found for image:", basename(imag[i])))
            }
          },
          error = function(e) {
            warning(paste(
              "Error processing image", basename(imag[i]), ":", e$message
            ))
          }
        )
      }

      # Step 10: Save results to CSV files
      if (save) {
        if (length(result) > 0) {
          save_analysis_results(result, polar_dir, "polar")
          save_analysis_results(result2, cartesian_dir, "perimeter")
        }

        if (length(morpho_results) > 0) {
          # Combine morphometrics into a data frame
          morpho_results <- morpho_results[!sapply(morpho_results, is.null)]
          if (length(morpho_results) > 0) {
            morpho_df <- do.call(rbind, lapply(morpho_results, as.data.frame))
            # Reorder to put Image first
            morpho_df <- morpho_df[, c("Image", setdiff(names(morpho_df), "Image"))]

            write.table(morpho_df,
              file = file.path(polar_dir, "MorphometricsEN.csv"),
              row.names = FALSE, sep = ";", dec = "."
            )
          }
        }
      }

      close(pb)
    },
    error = function(e) {
      stop(paste("Error in process_images:", e$message))
    }
  )
}

#' @param binary_image A binarized image.
#' @param scale_mm The actual length of the scale bar in millimeters (default: 1).
#' @return The number of pixels per mm, or NULL if not found.
#' @noRd
find_scale_bar_length <- function(binary_image, scale_mm = 1) {
  tryCatch(
    {
      labeled <- EBImage::bwlabel(binary_image)
      features <- EBImage::computeFeatures.shape(labeled)
      moments <- EBImage::computeFeatures.moment(labeled)

      if (is.null(features)) {
        return(NULL)
      }

      # Filter for scale bar candidates
      candidates <- which(features[, "s.area"] > 100 &
        features[, "s.area"] < 5000 &
        moments[, "m.eccentricity"] > 0.8)

      if (length(candidates) == 0) {
        return(NULL)
      }

      # For now, picking the candidate with the highest eccentricity
      best_idx <- which.max(moments[candidates, "m.eccentricity"])
      best_candidate <- candidates[best_idx]

      # Return pixels per mm (Major Axis length / scale_mm)
      return(moments[best_candidate, "m.majoraxis"] / scale_mm)
    },
    error = function(e) {
      warning(paste("Scale detection error:", e$message))
      invisible(NULL)
    }
  )
}

#' Create Result Structure for Analysis Data
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
create_result_structure <- function(distances, wavelets, image_name,
                                    elliptic_coef = NULL,
                                    result_type = "polar") {
  # Input validation
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

  # Add wavelet data if available
  if (!is.null(wavelets)) {
    for (i in seq_len(nrow(wavelets))) {
      result[[paste0("Wavelet_", i)]] <- wavelets[i, ]
    }
  }

  # Add coordinates
  if (result_type == "polar") {
    result$coords <- cbind(image_name, distances$coord)
  } else {
    result$coords <- cbind(image_name, distances$coords)
  }

  # Add elliptic coefficients if available
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
#'
#' @param binary_image The binarized image from preprocessing.
#' @param distances The distance measures from calculate_distances().
#' @param wavelets The wavelet analysis from calculate_wavelets_analysis().
#' @param image_name The base name of the image file.
#' @param output_dir The directory to save visualization images.
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics legend lines points polygon
#' @export
save_visualization_perimeter <- function(binary_image, distances, wavelets,
                                         image_name, output_dir) {
  # Input validation
  if (is.null(binary_image)) {
    stop("binary_image cannot be NULL")
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

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  tryCatch(
    {
      # Save zonal analysis plot for perimeter
      jpeg(
        filename = file.path(output_dir, paste("zonal_", tools::file_path_sans_ext(image_name), ".jpg", sep = "")),
        res = 600, width = 600, height = 400, units = "mm"
      )
      plot(binary_image)
      points(mean(distances$reordered_coords$x),
        mean(distances$reordered_coords$y),
        pch = 16
      )
      points(distances$reordered_coords$x[1], distances$reordered_coords$y[1],
        col = 3, pch = 16
      )

      # Draw colored zones for perimeter coordinates dynamically
      n_points <- length(distances$perimeter$dist)
      z1 <- 1:floor(n_points * 0.25)
      z2 <- (max(z1) + 1):floor(n_points * 0.5)
      z3 <- (max(z2) + 1):floor(n_points * 0.75)
      z4 <- (max(z3) + 1):n_points

      zones <- list(z1, z2, z3, z4)
      colors <- c(4, 5, 6, 7)

      for (i in seq_along(zones)) {
        lines(distances$perimeter$coords[zones[[i]], 1],
          distances$perimeter$coords[zones[[i]], 2],
          col = colors[i], pch = 16, lwd = 2, type = "b"
        )
      }
      dev.off()

      # Save wavelet plot for perimeter if available
      if (!is.null(wavelets)) {
        jpeg(filename = file.path(
          output_dir,
          paste("wavelet_", tools::file_path_sans_ext(image_name), ".jpg", sep = "")
        ))
        plot(wavelets$perimeter[5, ],
          main = paste("Wavelet 5", image_name),
          xlab = "", ylab = "", type = "l", col = 3, lwd = 2
        )

        # Add colored polygons for zones
        for (i in seq_along(zones)) {
          polygon(
            x = c(min(zones[[i]]), zones[[i]], max(zones[[i]])),
            y = c(
              min(wavelets$perimeter[5, ]),
              wavelets$perimeter[5, zones[[i]]],
              min(wavelets$perimeter[5, ])
            ),
            col = colors[i]
          )
        }

        legend("topleft",
          legend = paste(
            "Zona",
            sapply(zones, function(z) paste0(min(z), ":", max(z)))
          ),
          pch = c(15, 15, 15, 15), col = colors, border = NULL, bg = "white"
        )
        dev.off()
      }
    },
    error = function(e) {
      stop(paste("Error creating perimeter visualization:", e$message))
    }
  )
}
