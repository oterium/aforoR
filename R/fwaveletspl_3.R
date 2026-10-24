#' B-Spline Wavelet Transform
#'
#' Computes the stationary wavelet transform of a signal using B-spline filters.
#' This specific implementation is optimized for 1D signals like otolith radial
#' distances.
#'
#' @param ra0 A numeric vector representing the input signal (must be a power of
#'            2 length for some applications, though here it uses the signal
#'            length).
#' @param sjm A positive integer specifying the number of decomposition scales (e.g., 9).
#' @param det A logical value. If FALSE (default), it returns the wavelet
#'            detail coefficients (high-to-medium frequency differences between
#'            consecutive scales, centered around zero, capturing localized
#'            features such as lobation). If TRUE, it returns the smoothed
#'            approximation coefficients (low-frequency smoothed outline).
#' @return A matrix of dimensions (sjm x length(ra0)) containing the
#'         coefficients for each scale.
#' @export
#' @examples
#' # Decompose a signal into 9 scales
#' sig <- rnorm(512)
#' sjm <- 9
#' wavelets <- fwaveletspl_3(sig, sjm)
#' dim(wavelets)
fwaveletspl_3 <- function(ra0, sjm, det = FALSE) {
  # Input validation
  if (!is.numeric(ra0)) {
    stop("ra0 must be a numeric vector")
  }

  if (length(ra0) < 2) {
    stop("ra0 must contain at least 2 points")
  }

  if (any(is.na(ra0))) {
    stop("ra0 cannot contain NA values")
  }

  if (any(is.infinite(ra0))) {
    stop("ra0 cannot contain infinite values")
  }

  if (!is.numeric(sjm) || length(sjm) != 1 || is.na(sjm) || is.infinite(sjm) || sjm <= 0 || sjm != as.integer(sjm)) {
    stop("sjm must be a positive integer")
  }
  sjm <- as.integer(sjm)

  if (!is.logical(det) || length(det) != 1 || is.na(det)) {
    stop("det must be a logical value (TRUE or FALSE)")
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
    rd_j <- raj_1_anterior - raj_1 # Detail (wavelet plane: difference between consecutive approximations)
    raj_1_anterior <- raj_1
    if (det == FALSE) {
      wave[sj, ] <- rd_j
    } else {
      wave[sj, ] <- raj_1
    }
  }
  return(wave)
}
