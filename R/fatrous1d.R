#' Discrete Wavelet Transform (à trous algorithm)
#'
#' Computes the 1D discrete wavelet transform using the "à trous" (with holes)
#' algorithm. This is a non-decimated transform, meaning the signal length
#' remains constant across scales.
#'
#' @param rx A numeric vector representing the input signal (e.g., radial distances).
#' @param rh A numeric vector representing the filter coefficients.
#' @param snmin The minimum scale index for the filter.
#' @param snmax The maximum scale index for the filter.
#' @param sj The current scale index (j >= 1).
#' @return A numeric matrix (1 x length(rx)) representing the convolution result.
#' @export
#' @examples
#' # Internal use example:
#' sig <- rnorm(512)
#' filter <- c(0.0625, 0.25, 0.375, 0.25, 0.0625)
#' res <- fatrous1d(sig, filter, -2, 2, 1)
fatrous1d <- function(rx, rh, snmin, snmax, sj) {
  # Input validation
  if (!is.numeric(rx) || !is.numeric(rh)) {
    stop("rx and rh must be numeric vectors")
  }

  if (length(rx) == 0 || length(rh) == 0) {
    stop("rx and rh cannot be empty")
  }

  if (!is.numeric(snmin) || !is.numeric(snmax) || !is.numeric(sj)) {
    stop("snmin, snmax, and sj must be numeric values")
  }

  if (length(snmin) != 1 || length(snmax) != 1 || length(sj) != 1) {
    stop("snmin, snmax, and sj must be single values")
  }

  if (sj < 0) stop("sj must be >= 0")

  if (snmax < snmin) {
    stop("snmax must be greater than or equal to snmin")
  }

  sj <- sj - 1
  rg <- rh[seq(snmax - snmin + 1, 1, -1)]
  smmin <- -snmax
  smmax <- -snmin
  slrx <- length(rx)
  ry <- matrix(0, nrow = 1, ncol = slrx)

  sn_indices <- 0:(slrx - 1)
  sm_values <- smmin:smmax
  
  # Outer matrix of 1-based indices (dimensions: length(sm) x slrx)
  idx_matrix <- outer(sm_values, sn_indices, function(sm, sn) {
    1 + (sn + sm * (2^sj)) %% slrx
  })

  # Map values from rx and multiply by weights in rg
  rx_mapped <- matrix(rx[idx_matrix], nrow = length(sm_values))
  ry[1, ] <- colSums(rg * rx_mapped)

  return(ry)
}
