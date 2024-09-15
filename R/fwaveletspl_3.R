#' Wavelet Transform
#'
#' Computes the wavelet transform of a signal using the specified filters and scales.
#'
#' @param ra0 A numeric vector representing the input signal.
#' @param sjm An integer specifying the number of scales.
#' @param det A logical value indicating whether to return the detail coefficients (TRUE) or the approximation coefficients (FALSE). Default is FALSE.
#' @return A matrix containing the wavelet coefficients.
#' @export
#' @examples
#' # Example usage:
#' ra0 <- rnorm(1024)  # Example signal
#' sjm <- 5  # Number of scales
#' wavelets <- fwaveletspl_3(ra0, sjm)
#' print(wavelets)
fwaveletspl_3 <- function(ra0, sjm, det = F) {
  # Filters
  rh_ <- c(0.0625, 0.25, 0.375, 0.25, 0.0625)
  snh_min <- -2
  snh_max <- 2

  # Number of samples
  sN <- length(ra0)

  # Number of scales
  sj <- ceiling(log(sN, 2))  # This is used for safety, recalculates sjm based on ra0.

  if (sj != sjm) stop("Revisa los parametros")

  raj_1Anterior <- ra0

  wave <- matrix(nrow = sjm, ncol = sN)
  for (sj in 1:sjm) {
    raj_1 <- fatrous1d(raj_1Anterior, rh_, snh_min, snh_max, sj)
    rdj <- raj_1Anterior - raj_1  # Approximation
    raj_1Anterior <- raj_1
    if (det == F) {
      wave[sj, ] <- rdj
    } else {
      wave[sj, ] <- raj_1
    }
  }
  return(wave)
}
