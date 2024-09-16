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

  # Set working directory
  setwd(folder)

  # Create processing directories if they don't exist
  if (!dir.exists("Procesamiento")) {
    dir.create("Procesamiento")
  }
  if (!dir.exists("Procesamiento2")) {
    dir.create("Procesamiento2")
  }

  # Get list of image files
  imag <- list.files(pattern = "\\.jpg$", full.names = TRUE)
  pb <- utils::txtProgressBar(max = length(imag), style = 3)

  # Initialize results
  result <- list()
  result2 <- list()

  for (i in 1:length(imag)) {
    setTxtProgressBar(pb, i)

    # Read and process image
    foto <- EBImage::readImage(imag[i])
    foto2 <- EBImage::channel(foto, "grey")
    foto2 <- filter2(foto2, EBImage::makeBrush(size = 31, shape = 'gaussian', sigma = 9))

    if (is.null(threshold)) {
      foto2 <- foto2 > EBImage::otsu(foto2)
    } else {
      foto2 <- foto2 > threshold
    }

    cont <- EBImage::ocontour(EBImage::bwlabel(foto2))
    contorno <- NULL
    for (x in 1:length(cont)) {
      if (length(cont[[x]])[1] > 1000 && Momocs::coo_area(cont[[x]]) > 10^5) {
        contorno <- cont[[x]]
      }
    }

    if (!is.null(contorno)) {
      coefE <- Momocs::efourier(contorno, 32)
      coords <- Momocs::efourier_i(coefE, nb.pts = dim(contorno)[1])

      # Reorder contour
      x_max <- max(contorno$X)
      y_max <- contorno$Y[which.max(contorno$X)]
      punto <- list(x = x_max, y = y_max)

      x <- c(contorno$X[which(contorno$X == (punto$x) & contorno$X > mean(contorno$X))[1]:1],
             contorno$X[length(contorno$X):(which(contorno$Y == (punto$y) & contorno$X > mean(contorno$X))[1] + 1)])
      y <- c(contorno$Y[which(contorno$X == (punto$x) & contorno$X > mean(contorno$X))[1]:1],
             contorno$Y[length(contorno$X):(which(contorno$Y == (punto$y) & contorno$X > mean(contorno$X))[1] + 1)])

      # Calculate distances
      dist <- regularradius(x, y, n = 512)
      dist2 <- dper(x, y, 512)
      dist.norm <- dist$radii / max(dist$radii)
      dist.norm2 <- dist2$dist / max(dist2$dist)

      # Calculate wavelets
      if (wavelets) {
        wavelets <- fwaveletspl_3(dist.norm, 9, det = FALSE)
        wavelets2 <- fwaveletspl_3(dist.norm2, 9, det = FALSE)
      }

      # Save results
      if (save) {
        result[[i]] <- list(
          Distancia = dist$radii,
          Distancia_Norm = dist.norm,
          Wavelet_1 = wavelets[1, ],
          Wavelet_2 = wavelets[2, ],
          Wavelet_3 = wavelets[3, ],
          Wavelet_4 = wavelets[4, ],
          Wavelet_5 = wavelets[5, ],
          Wavelet_6 = wavelets[6, ],
          Wavelet_7 = wavelets[7, ],
          Wavelet_8 = wavelets[8, ],
          Wavelet_9 = wavelets[9, ],
          nombre = basename(imag[i]),
          coords = cbind(basename(imag[i]), dist$coord)
        )

        result2[[i]] <- list(
          Distancia = dist2$dist,
          Distancia_Norm = dist.norm2,
          elliptic = c(coefE$an, coefE$bn, coefE$cn, coefE$dn),
          Wavelet_1 = wavelets2[1, ],
          Wavelet_2 = wavelets2[2, ],
          Wavelet_3 = wavelets2[3, ],
          Wavelet_4 = wavelets2[4, ],
          Wavelet_5 = wavelets2[5, ],
          Wavelet_6 = wavelets2[6, ],
          Wavelet_7 = wavelets2[7, ],
          Wavelet_8 = wavelets2[8, ],
          Wavelet_9 = wavelets2[9, ],
          nombre = basename(imag[i]),
          coords = cbind(basename(imag[i]), dist2$coords)
        )
      }

      # Save test images
      if (testing) {
        setwd("./Procesamiento/")
        jpeg(file = paste("zonal ", basename(imag[i]), '.jpg', sep = ''), res = 600, width = 600, height = 400, units = 'mm')
        plot(foto2)
        points(mean(x), mean(y), pch = 16)
        points(x[1], y[1], col = 2, pch = 16)
        lines(dist$coord[1:126, 1] + mean(x), dist$coord[1:126, 2] + mean(y), col = 4, pch = 16, lwd = 2, type = 'b')
        lines(dist$coord[127:256, 1] + mean(x), dist$coord[127:256, 2] + mean(y), col = 5, pch = 16, lwd = 2, type = 'b')
        lines(dist$coord[257:374, 1] + mean(x), dist$coord[257:374, 2] + mean(y), col = 6, pch = 16, lwd = 2, type = 'b')
        lines(dist$coord[375:512, 1] + mean(x), dist$coord[375:512, 2] + mean(y), col = 7, pch = 16, lwd = 2, type = 'b')
        dev.off()

        jpeg(file = paste("wavelet ", basename(imag[i]), '.jpg', sep = ''))
        plot(wavelets[5, ], main = paste("Wavelet 5", basename(imag[i])), xlab = "", ylab = "", type = "l", col = 3, lwd = 2)
        polygon(x = c(1, 1:126, 126), y = c(min(wavelets[5, ]), wavelets[5, 1:126], min(wavelets[5, ])), col = 4)
        polygon(x = c(127, 127:256, 256), y = c(min(wavelets[5, ]), wavelets[5, 127:256], min(wavelets[5, ])), col = 5)
        polygon(x = c(257, 257:374, 374), y = c(min(wavelets[5, ]), wavelets[5, 257:374], min(wavelets[5, ])), col = 6)
        polygon(x = c(375, 375:512, 512), y = c(min(wavelets[5, ]), wavelets[5, 375:512], min(wavelets[5, ])), col = 7)
        dev.off()

        setwd("./Procesamiento2/")
        jpeg(file = paste("zonal ", basename(imag[i]), '.jpg', sep = ''), res = 600, width = 600, height = 400, units = 'mm')
        plot(foto2)
        points(mean(x), mean(y), pch = 16)
        points(x[1], y[1], col = 3, pch = 16)
        lines(dist2$coords[1:126, 1], dist2$coords[1:126, 2], col = 4, pch = 16, lwd = 2, type = 'b')
        lines(dist2$coords[127:256, 1], dist2$coords[127:256, 2], col = 5, pch = 16, lwd = 2, type = 'b')
        lines(dist2$coords[257:374, 1], dist2$coords[257:374, 2], col = 6, pch = 16, lwd = 2, type = 'b')
        lines(dist2$coords[375:512, 1], dist2$coords[375:512, 2], col = 7, pch = 16, lwd = 2, type = 'b')
        dev.off()

        jpeg(file = paste("wavelet ", basename(imag[i]), '.jpg', sep = ''))
        plot(wavelets2[5, ], main = paste("Wavelet 5", basename(imag[i])), xlab = "", ylab = "", type = "l", col = 3, lwd = 2)
        polygon(x = c(1, 1:126, 126), y = c(min(wavelets2[5, ]), wavelets2[5, 1:126], min(wavelets2[5, ])), col = 4)
        polygon(x = c(127, 127:256, 256), y = c(min(wavelets2[5, ]), wavelets2[5, 127:256], min(wavelets2[5, ])), col = 5)
        polygon(x = c(257, 257:374, 374), y = c(min(wavelets2[5, ]), wavelets2[5, 257:374], min(wavelets2[5, ])), col = 6)
        polygon(x = c(375, 375:512, 512), y = c(min(wavelets2[5, ]), wavelets2[5, 375:512], min(wavelets2[5, ])), col = 7)
        legend("topleft", legend = c("Zona 1:126", "Zona 127:256", "Zona 257:374", "Zona 375:512"),
               pch = c(15, 15, 15, 15), col = c(4, 5, 6, 7), border = NULL, bg = "white")
        dev.off()
      }
    }
  }

  # Save results to CSV files
  if (save) {
    save_results_to_csv(result, "Procesamiento")
    save_results_to_csv(result2, "Procesamiento2")
  }

  close(pb)
}

# Helper function to save results to CSV files
save_results_to_csv <- function(result, directory) {
  setwd(directory)

  write_csv <- function(data, file_name) {
    z <- data
    rownames(z) <- result$nombre
    write.table(z, file = file_name, dec = ".", sep = ";", col.names = FALSE)
  }

  write_csv(result$Distancia, "DistanciaEN.csv")
  write_csv(result$Distancia_Norm, "Distancia_NormEN.csv")
  write_csv(result$Wavelet_1, "Wavelet_1EN.csv")
  write_csv(result$Wavelet_2, "Wavelet_2EN.csv")
  write_csv(result$Wavelet_3, "Wavelet_3EN.csv")
  write_csv(result$Wavelet_4, "Wavelet_4EN.csv")
  write_csv(result$Wavelet_5, "Wavelet_5EN.csv")
  write_csv(result$Wavelet_6, "Wavelet_6EN.csv")
  write_csv(result$Wavelet_7, "Wavelet_7EN.csv")
  write_csv(result$Wavelet_8, "Wavelet_8EN.csv")
  write_csv(result$Wavelet_9, "Wavelet_9EN.csv")
  write_csv(result$coords, "Coords.csv")

  if ("elliptic" %in% names(result)) {
    z <- result$elliptic
    colnames(z) <- c(paste("A", 1:32, sep = ""), paste("B", 1:32, sep = ""),
                     paste("C", 1:32, sep = ""), paste("D", 1:32, sep = ""))
    rownames(z) <- result$nombre
    write.table(z, file = "EllipticCoeEN.csv", dec = ".", sep = ";", col.names = TRUE)
  }
}
