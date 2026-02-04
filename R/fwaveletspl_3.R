#' B-Spline Wavelet Transform
#'
#' Computes the stationary wavelet transform of a signal using B-spline filters.
#' This specific implementation is optimized for 1D signals like otolith radial
#' distances.
#'
#' @param ra0 A numeric vector representing the input signal (must be a power of
#'            2 length for some applications, though here it uses the signal
#'            length).
#' @param sjm An integer specifying the number of decomposition scales.
#' @param det A logical value. If FALSE (default), it returns the approximation
#'            coefficients (low-frequency). If TRUE, it returns the detail
#'            coefficients (high-frequency).
#' @return A matrix of dimensions (sjm x length(ra0)) containing the
#'         coefficients for each scale.
#' @export
#' @examples
#' # Decompose a signal into 5 scales
#' sig <- rnorm(512)
#' sjm <- 9 # Number of scales
#' wavelets <- fwaveletspl_3(sig, sjm)
#' # plot(wavelets[1,], type="l") # Plot first scale
fwaveletspl_3 <- function(ra0, sjm, det = FALSE) {
  # Input validation
  if (!is.numeric(ra0)) {
    stop("ra0 must be a numeric vector")
  }

  if (length(ra0) == 0) {
    stop("ra0 cannot be empty")
  }

  if (!is.numeric(sjm) || length(sjm) != 1 || sjm <= 0) {
    stop("sjm must be a positive integer")
  }

  if (!is.logical(det) || length(det) != 1) {
    stop("det must be a logical value")
  }

  if (any(is.na(ra0))) {
    stop("ra0 cannot contain NA values")
  }

  # Filters
  rh_ <- c(0.0625, 0.25, 0.375, 0.25, 0.0625)
  snh_min <- -2
  snh_max <- 2

  # Number of samples
  s_n <- length(ra0)

  # Number of scales
  # This is used for safety, recalculates sjm based on ra_0.
  sj_max <- ceiling(log(s_n, 2))

  if (sjm > sj_max) {
    stop(paste(
      "sjm cannot exceed the max scales for this signal length:",
      sj_max
    ))
  }

  raj_1_anterior <- ra0

  wave <- matrix(nrow = sjm, ncol = s_n)
  for (sj in 1:sjm) {
    raj_1 <- fatrous1d(raj_1_anterior, rh_, snh_min, snh_max, sj)
    rd_j <- raj_1_anterior - raj_1 # Approximation
    raj_1_anterior <- raj_1
    if (det == FALSE) {
      wave[sj, ] <- rd_j
    } else {
      wave[sj, ] <- raj_1
    }
  }
  return(wave)
}
