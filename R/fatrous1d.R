#' Atrous Convolution
#'
#' Computes the atrous convolution of a signal.
#'
#' @param rx A numeric vector representing the input signal.
#' @param rh A numeric vector representing the filter coefficients.
#' @param snmin The minimum scale index.
#' @param snmax The maximum scale index.
#' @param sj The scale factor (must be >= 0).
#' @return A matrix representing the atrous convolution result.
#' @export
#' @examples
#' # Example usage:
#' ra1 <- fatrous1d(ra0, rh_, snh_min, snh_max, 1)
#' rd1 <- fatrous1d(ra0, rg_, sng_min, sng_max, 1)

fatrous1d <- function(rx, rh, snmin, snmax, sj) {
  # sj: scale
  # Computes the atrous convolution:
  #
  # y = x[n] * h?[n];
  #              *: Convolution product
  #              h? = h with 2^j-1 zeros between each sample
  #
  # Example:
  #  [ra0, rt] = facb3transient(100, 10, 10e-3, 6e-3, 2, 5120);
  #  rh_ = sqrt(2) * c(0.125, 0.375, 0.375, 0.125);
  #  snh_min = -2; snh_max = 1;
  #  rg_ = sqrt(2) * c(0.5, -0.5);
  #  sng_min = -1; sng_max = 0;
  #  ra1 = fatrous1d(ra0, rh_, snh_min, snh_max, 1);
  #  rd1 = fatrous1d(ra0, rg_, sng_min, sng_max, 1);

  if (sj < 0) stop("sj tiene que ser >=0")

  sj = sj - 1
  rg = rh[seq(snmax - snmin + 1, 1, -1)]
  smmin = -snmax
  smmax = -snmin
  slrx = length(rx)
  ry = matrix(0, nrow = 1, ncol = slrx)

  for (sn in 0:(slrx - 1)) {
    for (sm in seq(smmin, smmax, 1)) {
      sp = rg[sm - smmin + 1] * rx[(1 + (sn + sm * (2^sj)) %% slrx)] ## %% Modulus
      ry[, sn + 1] = ry[, sn + 1] + sp
    }
  }

  return(ry)
}
