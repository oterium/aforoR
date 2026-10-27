#' Discrete Wavelet Transform (a trous algorithm)
#'
#' Computes the 1D discrete wavelet transform using the a trous (with holes)
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

  if (any(is.na(rx)) || any(is.na(rh))) {
    stop("rx and rh cannot contain NA values")
  }

  if (any(is.infinite(rx)) || any(is.infinite(rh))) {
    stop("rx and rh cannot contain infinite values")
  }

  if (!is.numeric(snmin) || !is.numeric(snmax) || !is.numeric(sj)) {
    stop("snmin, snmax, and sj must be numeric values")
  }

  if (length(snmin) != 1 || length(snmax) != 1 || length(sj) != 1) {
    stop("snmin, snmax, and sj must be single values")
  }

  if (is.na(snmin) || is.na(snmax) || is.na(sj) ||
      is.infinite(snmin) || is.infinite(snmax) || is.infinite(sj)) {
    stop("snmin, snmax, and sj must be finite values")
  }

  if (snmin != as.integer(snmin) || snmax != as.integer(snmax)) {
    stop("snmin and snmax must be integers")
  }
  snmin <- as.integer(snmin)
  snmax <- as.integer(snmax)

  if (sj < 1 || sj != as.integer(sj)) {
    stop("sj must be an integer >= 1")
  }
  sj <- as.integer(sj)

  if (snmax < snmin) {
    stop("snmax must be greater than or equal to snmin")
  }

  if (length(rh) != (snmax - snmin + 1)) {
    stop("rh must have length equal to snmax - snmin + 1")
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
#' @param det A logical value. If TRUE (default), it returns the wavelet
#'            detail coefficients (high-to-medium frequency differences between
#'            consecutive scales, centered around zero, capturing localized
#'            features such as lobation). If FALSE, it returns the smoothed
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
fwaveletspl_3 <- function(ra0, sjm, det = TRUE) {
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
    if (det) {
      wave[sj, ] <- rd_j
    } else {
      wave[sj, ] <- raj_1
    }
  }
  return(wave)
}

#' Calculate Wavelet Analysis for Distance Data
#'
#' Performs multi-scale wavelet decomposition on normalized distance data
#' (polar or perimeter). This is useful for capturing shape variation at
#' different frequencies.
#'
#' @param distances A list from \code{calculate_distances()} containing normalized data.
#' @param n_scales Number of wavelet scales to extract (default: 9).
#' @param detail Logical; if TRUE (default), returns wavelet detail
#'               coefficients (high-to-medium frequency differences between
#'               scales, centered around zero). If FALSE, returns smoothed
#'               approximation coefficients.
#' @return A list with wavelet coefficients for 'polar' and 'perimeter' components.
#' @export
#' @examples
#' # Example using calculated distances from an otolith contour:
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' binary_img <- preprocess_image(image_path)
#' contour <- extract_contour(binary_img)
#' distances <- calculate_distances(contour)
#' wavelets <- calculate_wavelets_analysis(distances, n_scales = 5)
#' str(wavelets)
#' @seealso \code{\link{fwaveletspl_3}}
calculate_wavelets_analysis <- function(distances, n_scales = 9,
                                        detail = TRUE) {
  if (!is.list(distances)) {
    stop("distances must be a list")
  }

  if (!all(c("polar", "perimeter") %in% names(distances))) {
    stop("distances must contain 'polar' and 'perimeter' components")
  }

  if (!is.numeric(n_scales) || length(n_scales) != 1 || is.na(n_scales) || is.infinite(n_scales) || n_scales <= 0 || n_scales != as.integer(n_scales)) {
    stop("n_scales must be a positive integer")
  }

  if (!is.logical(detail) || length(detail) != 1 || is.na(detail)) {
    stop("detail must be a logical value (TRUE or FALSE)")
  }

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
