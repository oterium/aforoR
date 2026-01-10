#' Create Zones for Visualization
#' 
#' @param n_points Total number of points
#' @param n_zones Number of zones to create
#' @return A list of integer vectors representing zone indices
#' @noRd
create_zones <- function(n_points, n_zones = 4) {
  zone_size <- n_points / n_zones
  lapply(1:n_zones, function(i) {
    start <- floor((i-1) * zone_size) + 1
    end <- if (i == n_zones) n_points else floor(i * zone_size)
    start:end
  })
}

#' Preprocess Image for Contour Analysis
#'
#' Applies preprocessing steps to an image including grayscale conversion,
#' Gaussian filtering, and binarization.
#'
#' @param image_path A string specifying the path to the image file.
#' @param threshold A numeric value for binarization. If NULL, uses Otsu method.
#' @param gaussian_size Integer specifying the size of the Gaussian filter brush (default: 31).
#' @param gaussian_sigma Numeric value specifying the sigma parameter for Gaussian filtering (default: 9).
#' @return A binarized image ready for contour detection.
#' @export
#' @examples
#' \dontrun{
#' # Example usage:
#' processed_img <- preprocess_image("path/to/image.jpg")
#' }
preprocess_image <- function(image_path, threshold = NULL, gaussian_size = 31, gaussian_sigma = 9) {
  # Input validation
  if (!is.character(image_path) || length(image_path) != 1) {
    stop("image_path must be a single character string")
  }

  if (!file.exists(image_path)) {
    stop(paste("Image file does not exist:", image_path))
  }

  if (!is.null(threshold) && (!is.numeric(threshold) || length(threshold) != 1)) {
    stop("threshold must be a single numeric value or NULL")
  }

  if (!is.null(threshold) && (threshold < 0 || threshold > 1)) {
    stop("threshold must be between 0 and 1")
  }

  tryCatch({
    # Read and process image
    foto <- EBImage::readImage(image_path)
    foto2 <- EBImage::channel(foto, "grey")
    foto2 <- EBImage::filter2(foto2, EBImage::makeBrush(size = gaussian_size, shape = 'gaussian', sigma = gaussian_sigma))

    if (is.null(threshold)) {
      foto2 <- foto2 > EBImage::otsu(foto2)
    } else {
      foto2 <- foto2 > threshold
    }

    return(foto2)
  }, error = function(e) {
    stop(paste("Error processing image:", e$message))
  })
}

#' Extract Main Contour from Binarized Image
#'
#' Detects and extracts the largest contour from a binarized image that meets
#' size criteria for otolith analysis.
#'
#' @param binary_image A binarized image from preprocess_image().
#' @param min_points Integer specifying the minimum number of points for a valid contour (default: 1000).
#' @param min_area Numeric value specifying the minimum area for a valid contour (default: 10^5).
#' @return A matrix containing the contour coordinates, or NULL if no suitable contour found.
#' @export
#' @examples
#' \dontrun{
#' # Example usage:
#' binary_img <- preprocess_image("path/to/image.jpg")
#' contour <- extract_contour(binary_img)
#' }
extract_contour <- function(binary_image, min_points = 1000, min_area = 10^5) {
  # Input validation
  if (is.null(binary_image)) {
    stop("binary_image cannot be NULL")
  }

  tryCatch({
    cont <- EBImage::ocontour(EBImage::bwlabel(binary_image))
    contorno <- NULL

    if (length(cont) == 0) {
      warning("No contours found in image")
      return(NULL)
    }

    for (x in 1:length(cont)) {
      if (length(cont[[x]])[1] > min_points && Momocs::coo_area(cont[[x]]) > min_area) {
        contorno <- cont[[x]]
        break
      }
    }

    if (is.null(contorno)) {
      warning(sprintf("No contour found meeting size criteria (>%d points and area >%g)", min_points, min_area))
    }

    return(contorno)
  }, error = function(e) {
    stop(paste("Error extracting contour:", e$message))
  })
}

#' Calculate Distance Measures from Contour
#'
#' Computes both polar and perimeter distance measures from an otolith contour.
#'
#' @param contour A matrix containing contour coordinates from extract_contour().
#' @param n_points An integer specifying the number of points to sample (default: 512).
#' @return A list containing:
#'   \describe{
#'     \item{polar}{List with components:
#'       \itemize{
#'         \item radii: numeric vector of distances from centroid
#'         \item coord: matrix of polar coordinates
#'         \item normalized: normalized radii values
#'       }
#'     }
#'     \item{perimeter}{List with components:
#'       \itemize{
#'         \item dist: numeric vector of perimeter distances
#'         \item coords: matrix of perimeter coordinates
#'         \item normalized: normalized distance values
#'       }
#'     }
#'     \item{reordered_coords}{List with x and y coordinate vectors}
#'   }
#' @export
#' @examples
#' \dontrun{
#' # Example usage:
#' distances <- calculate_distances(contour, n_points = 512)
#' }
calculate_distances <- function(contour, n_points = 512) {
  # Input validation
  if (is.null(contour) || !is.data.frame(contour) && !is.matrix(contour)) {
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
    indices <- c(start_idx:n_contour, 1:(start_idx-1))
    x <- contour$X[indices]
    y <- contour$Y[indices]
  }

  # Calculate distances with error handling
  tryCatch({
    dist <- regularradius(x, y, n = n_points)
    dist2 <- dper(x, y, n_points)
    dist.norm <- dist$radii / max(dist$radii)
    dist.norm2 <- dist2$dist / max(dist2$dist)

    return(list(
      polar = list(radii = dist$radii, coord = dist$coord, normalized = dist.norm),
      perimeter = list(dist = dist2$dist, coords = dist2$coords, normalized = dist.norm2),
      reordered_coords = list(x = x, y = y)
    ))
  }, error = function(e) {
    stop(paste("Error in distance calculation:", e$message))
  })
}

#' Calculate Wavelet Analysis for Distance Data
#'
#' Computes wavelet transforms for both polar and perimeter distance measures.
#'
#' @param distances A list from calculate_distances() containing distance measures.
#' @param n_scales An integer specifying the number of wavelet scales (default: 9).
#' @param detail A logical indicating whether to return detail coefficients (default: FALSE).
#' @return A list containing wavelet transforms for both polar and perimeter distances.
#' @export
#' @examples
#' \dontrun{
#' # Example usage:
#' wavelets <- calculate_wavelets_analysis(distances, n_scales = 9)
#' }
calculate_wavelets_analysis <- function(distances, n_scales = 9, detail = FALSE) {
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
  if (!("normalized" %in% names(distances$polar)) || !("normalized" %in% names(distances$perimeter))) {
    stop("distances must contain normalized components")
  }

  tryCatch({
    wavelets_polar <- fwaveletspl_3(distances$polar$normalized, n_scales, det = detail)
    wavelets_perimeter <- fwaveletspl_3(distances$perimeter$normalized, n_scales, det = detail)

    return(list(
      polar = wavelets_polar,
      perimeter = wavelets_perimeter
    ))
  }, error = function(e) {
    stop(paste("Error in wavelet analysis:", e$message))
  })
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
#' @export
save_visualization <- function(binary_image, distances, wavelets, image_name, output_dir, wavelet_scale = 5) {
  # Input validation
  if (is.null(binary_image)) {
    stop("binary_image cannot be NULL")
  }

  if (!is.list(distances) || !is.list(wavelets)) {
    stop("distances and wavelets must be lists")
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

  original_dir <- getwd()

  tryCatch({
    setwd(output_dir)

    # Save zonal analysis plot
    jpeg(file = paste("zonal ", image_name, '.jpg', sep = ''),
         res = 600, width = 600, height = 400, units = 'mm')
    plot(binary_image)
    points(mean(distances$reordered_coords$x), mean(distances$reordered_coords$y), pch = 16)
    points(distances$reordered_coords$x[1], distances$reordered_coords$y[1], col = 2, pch = 16)

    # Draw colored zones
    zones <- create_zones(n_points = 512, n_zones = 4)
    colors <- c(4, 5, 6, 7)

    for (i in 1:length(zones)) {
      lines(distances$polar$coord[zones[[i]], 1] + mean(distances$reordered_coords$x),
            distances$polar$coord[zones[[i]], 2] + mean(distances$reordered_coords$y),
            col = colors[i], pch = 16, lwd = 2, type = 'b')
    }
    dev.off()

    # Save wavelet plot
    jpeg(file = paste("wavelet ", image_name, '.jpg', sep = ''))
    plot(wavelets$polar[wavelet_scale, ],
         main = paste("Wavelet", wavelet_scale, image_name),
         xlab = "", ylab = "", type = "l", col = 3, lwd = 2)

    # Add colored polygons for zones
    for (i in 1:length(zones)) {
      polygon(x = c(min(zones[[i]]), zones[[i]], max(zones[[i]])),
              y = c(min(wavelets$polar[wavelet_scale, ]),
                   wavelets$polar[wavelet_scale, zones[[i]]],
                   min(wavelets$polar[wavelet_scale, ])),
              col = colors[i])
    }
    dev.off()
  }, error = function(e) {
    stop(paste("Error creating visualization:", e$message))
  }, finally = {
    setwd(original_dir)
  })
}

#' Save Analysis Results to CSV Files
#'
#' Exports analysis results including distances, wavelets, and coordinates to CSV files.
#'
#' @param analysis_results A list containing all analysis results for multiple images.
#' @param output_dir The directory to save CSV files.
#' @param result_type A string indicating the type of results ("polar" or "perimeter").
#' @export
save_analysis_results <- function(analysis_results, output_dir, result_type = "polar") {
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

  original_dir <- getwd()

  tryCatch({
    setwd(output_dir)

    # Extract data for each type of measurement
    distances <- sapply(analysis_results, function(x) x$Distancia)
    distances_norm <- sapply(analysis_results, function(x) x$Distancia_Norm)
    coordinates <- do.call(rbind, lapply(analysis_results, function(x) x$coords))
    image_names <- sapply(analysis_results, function(x) x$nombre)

    # Save distance data
    write_analysis_csv(distances, "DistanciaEN.csv", image_names)
    write_analysis_csv(distances_norm, "Distancia_NormEN.csv", image_names)

    # Save wavelet data
    for (i in 1:9) {
      wavelet_data <- sapply(analysis_results, function(x) x[[paste0("Wavelet_", i)]])
      write_analysis_csv(wavelet_data, paste0("Wavelet_", i, "EN.csv"), image_names)
    }

    # Save coordinates
    write.table(coordinates, file = "Coords.csv", dec = ".", sep = ";", col.names = FALSE)

    # Save elliptic coefficients if available
    if ("elliptic" %in% names(analysis_results[[1]])) {
      elliptic_data <- do.call(rbind, lapply(analysis_results, function(x) x$elliptic))
      colnames(elliptic_data) <- c(paste("A", 1:32, sep = ""), paste("B", 1:32, sep = ""),
                                   paste("C", 1:32, sep = ""), paste("D", 1:32, sep = ""))
      rownames(elliptic_data) <- image_names
      write.table(elliptic_data, file = "EllipticCoeEN.csv", dec = ".", sep = ";", col.names = TRUE)
    }
  }, error = function(e) {
    stop(paste("Error saving analysis results:", e$message))
  }, finally = {
    setwd(original_dir)
  })
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

  if (!is.matrix(data) && !is.data.frame(data)) {
    if (is.vector(data)) {
      data <- matrix(data, nrow = 1)
    } else {
      stop("data must be a matrix, data frame, or vector")
    }
  }

  if (is.null(nrow(data)) || nrow(data) != length(row_names)) {
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
#' @export
#'
process_images <- function(folder, subfolder = FALSE, threshold = NULL, wavelets = TRUE, ef = TRUE, testing = TRUE, pseudolandmarks = "both", save = TRUE) {

  # Input validation
  if (!is.character(folder) || length(folder) != 1) {
    stop("folder must be a single character string")
  }

  if (!dir.exists(folder)) {
    stop(paste("Folder does not exist:", folder))
  }

  if (!is.logical(wavelets) || !is.logical(ef) || !is.logical(testing) || !is.logical(save)) {
    stop("wavelets, ef, testing, and save must be logical values")
  }

  if (!pseudolandmarks %in% c('curvilinear', 'polar', 'both')) {
    stop("pseudolandmarks must be one of: 'curvilinear', 'polar', 'both'")
  }

  # Store original directory
  original_dir <- getwd()

  tryCatch({
    # Change to working directory safely
    tryCatch({
      setwd(folder)
    }, error = function(e) {
      stop(paste("Cannot change to folder:", folder, "-", e$message))
    })

    # Create processing directories if they don't exist
    if (!dir.exists("Procesamiento")) {
      dir.create("Procesamiento", recursive = TRUE)
    }
    if (!dir.exists("Procesamiento2")) {
      dir.create("Procesamiento2", recursive = TRUE)
    }

    # Get list of image files
    imag <- list.files(pattern = "\\.jpg$", full.names = TRUE)

    if (length(imag) == 0) {
      warning("No .jpg files found in the specified folder")
      return(invisible(NULL))
    }

    pb <- utils::txtProgressBar(max = length(imag), style = 3)

    # Initialize results
    result <- list()
    result2 <- list()

    for (i in 1:length(imag)) {
      utils::setTxtProgressBar(pb, i)

      tryCatch({
        # Step 1: Preprocess image
        binary_image <- preprocess_image(imag[i], threshold)

        # Step 2: Extract contour
        contorno <- extract_contour(binary_image)

        if (!is.null(contorno)) {
          # Step 3: Calculate elliptic Fourier descriptors if requested
          coefE <- NULL
          if (ef) {
            coefE <- Momocs::efourier(contorno, 32)
          }

          # Step 4: Calculate distance measures
          distances <- calculate_distances(contorno, n_points = 512)

          # Step 5: Calculate wavelets if requested
          wavelet_results <- NULL
          if (wavelets) {
            wavelet_results <- calculate_wavelets_analysis(distances, n_scales = 9, detail = FALSE)
          }

          # Step 6: Save results
          if (save) {
            result[[i]] <- create_result_structure(distances$polar, wavelet_results$polar, basename(imag[i]), coefE, "polar")
            result2[[i]] <- create_result_structure(distances$perimeter, wavelet_results$perimeter, basename(imag[i]), coefE, "perimeter")
          }

          # Step 7: Save visualizations if requested
          if (testing) {
            save_visualization(binary_image, distances, wavelet_results, basename(imag[i]), "./Procesamiento/")
            save_visualization_perimeter(binary_image, distances, wavelet_results, basename(imag[i]), "./Procesamiento2/")
          }
        } else {
          warning(paste("No suitable contour found for image:", basename(imag[i])))
        }
      }, error = function(e) {
        warning(paste("Error processing image", basename(imag[i]), ":", e$message))
      })
    }

    # Step 8: Save results to CSV files
    if (save && length(result) > 0) {
      save_analysis_results(result, "./Procesamiento/", "polar")
      save_analysis_results(result2, "./Procesamiento2/", "perimeter")
    }

    close(pb)

  }, error = function(e) {
    stop(paste("Error in process_images:", e$message))
  }, finally = {
    # Always return to original directory
    setwd(original_dir)
  })
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
create_result_structure <- function(distances, wavelets, image_name, elliptic_coef = NULL, result_type = "polar") {
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
    for (i in 1:9) {
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
    result$elliptic <- c(elliptic_coef$an, elliptic_coef$bn, elliptic_coef$cn, elliptic_coef$dn)
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
#' @export
save_visualization_perimeter <- function(binary_image, distances, wavelets, image_name, output_dir) {
  # Input validation
  if (is.null(binary_image)) {
    stop("binary_image cannot be NULL")
  }

  if (!is.list(distances) || !is.list(wavelets)) {
    stop("distances and wavelets must be lists")
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

  original_dir <- getwd()

  tryCatch({
    setwd(output_dir)

    # Save zonal analysis plot for perimeter
    jpeg(file = paste("zonal ", image_name, '.jpg', sep = ''),
         res = 600, width = 600, height = 400, units = 'mm')
    plot(binary_image)
    points(mean(distances$reordered_coords$x), mean(distances$reordered_coords$y), pch = 16)
    points(distances$reordered_coords$x[1], distances$reordered_coords$y[1], col = 3, pch = 16)

    # Draw colored zones for perimeter coordinates
    zones <- create_zones(n_points = 512, n_zones = 4)
    colors <- c(4, 5, 6, 7)

    for (i in 1:length(zones)) {
      lines(distances$perimeter$coords[zones[[i]], 1],
            distances$perimeter$coords[zones[[i]], 2],
            col = colors[i], pch = 16, lwd = 2, type = 'b')
    }
    dev.off()

    # Save wavelet plot for perimeter
    jpeg(file = paste("wavelet ", image_name, '.jpg', sep = ''))
    plot(wavelets$perimeter[5, ],
         main = paste("Wavelet 5", image_name),
         xlab = "", ylab = "", type = "l", col = 3, lwd = 2)

    # Add colored polygons for zones
    for (i in 1:length(zones)) {
      polygon(x = c(min(zones[[i]]), zones[[i]], max(zones[[i]])),
              y = c(min(wavelets$perimeter[5, ]),
                   wavelets$perimeter[5, zones[[i]]],
                   min(wavelets$perimeter[5, ])),
              col = colors[i])
    }

    legend("topleft", legend = c("Zona 1:126", "Zona 127:256", "Zona 257:374", "Zona 375:512"),
           pch = c(15, 15, 15, 15), col = colors, border = NULL, bg = "white")
    dev.off()
  }, error = function(e) {
    stop(paste("Error creating perimeter visualization:", e$message))
  }, finally = {
    setwd(original_dir)
  })
}


