#' 1D Discrete Wavelet Convolution (A Trous Algorithm)
#'
#' Computes the 1D stationary (non-decimated) wavelet convolution using the
#' à trous ("with holes") algorithm across dyadic scales.
#'
#' @param signal A numeric vector representing the 1D input signal (e.g., radial distances).
#' @param filter_weights A numeric vector representing the convolution filter coefficients.
#' @param scale_min Integer; the minimum index for the filter window.
#' @param scale_max Integer; the maximum index for the filter window.
#' @param scale_j Integer >= 1; the current decomposition scale index.
#' @return A numeric matrix of dimensions (1 x length(signal)) containing the convolution result.
#' @export
#' @examples
#' sig <- rnorm(512)
#' bspline_filter <- c(0.0625, 0.25, 0.375, 0.25, 0.0625)
#' res <- wavelet_atrous_1d(sig, bspline_filter, -2, 2, 1)
wavelet_atrous_1d <- function(signal, filter_weights, scale_min, scale_max, scale_j) {
  # Input validation
  if (!is.numeric(signal) || !is.numeric(filter_weights)) {
    stop("signal and filter_weights must be numeric vectors")
  }

  if (length(signal) == 0 || length(filter_weights) == 0) {
    stop("signal and filter_weights cannot be empty")
  }

  if (any(is.na(signal)) || any(is.na(filter_weights))) {
    stop("signal and filter_weights cannot contain NA values")
  }

  if (any(is.infinite(signal)) || any(is.infinite(filter_weights))) {
    stop("signal and filter_weights cannot contain infinite values")
  }

  if (!is.numeric(scale_min) || !is.numeric(scale_max) || !is.numeric(scale_j)) {
    stop("scale_min, scale_max, and scale_j must be numeric values")
  }

  if (length(scale_min) != 1 || length(scale_max) != 1 || length(scale_j) != 1) {
    stop("scale_min, scale_max, and scale_j must be single values")
  }

  if (is.na(scale_min) || is.na(scale_max) || is.na(scale_j) ||
      is.infinite(scale_min) || is.infinite(scale_max) || is.infinite(scale_j)) {
    stop("scale_min, scale_max, and scale_j must be finite values")
  }

  if (scale_min != as.integer(scale_min) || scale_max != as.integer(scale_max)) {
    stop("scale_min and scale_max must be integers")
  }
  scale_min <- as.integer(scale_min)
  scale_max <- as.integer(scale_max)

  if (scale_j < 1 || scale_j != as.integer(scale_j)) {
    stop("scale_j must be an integer >= 1")
  }
  scale_j <- as.integer(scale_j)

  if (scale_max < scale_min) {
    stop("scale_max must be greater than or equal to scale_min")
  }

  if (length(filter_weights) != (scale_max - scale_min + 1)) {
    stop("filter_weights must have length equal to scale_max - scale_min + 1")
  }

  dilation_power <- scale_j - 1L
  step_size <- 2L^dilation_power
  rev_weights <- rev(filter_weights)
  shift_range <- (-scale_max):(-scale_min)
  n_pts <- length(signal)

  sample_indices <- 0L:(n_pts - 1L)
  index_grid <- outer(shift_range, sample_indices, function(m, n) {
    1L + (n + m * step_size) %% n_pts
  })

  mapped_samples <- matrix(signal[index_grid], nrow = length(shift_range))
  smoothed_output <- matrix(colSums(rev_weights * mapped_samples), nrow = 1L, ncol = n_pts)

  return(smoothed_output)
}

#' Multi-Scale Spline Wavelet Transform
#'
#' Computes the multi-scale stationary wavelet decomposition of a 1D closed signal
#' using cubic B-spline scaling filters.
#'
#' @param signal A numeric vector representing the input signal (e.g., radial distances).
#' @param n_scales A positive integer specifying the number of decomposition scales (e.g., 9).
#' @param detail A logical value. If TRUE (default), it returns the wavelet
#'               detail coefficients (high-to-medium frequency differences between
#'               consecutive scales, centered around zero, capturing localized
#'               features such as lobation). If FALSE, it returns the smoothed
#'               approximation coefficients (low-frequency smoothed outline).
#' @return A matrix of dimensions (n_scales x length(signal)) containing the
#'         coefficients for each scale.
#' @export
#' @examples
#' # Decompose a signal into 9 scales
#' sig <- rnorm(512)
#' n_scales <- 9
#' wavelets <- wavelet_spline_transform(sig, n_scales)
#' dim(wavelets)
wavelet_spline_transform <- function(signal, n_scales, detail = TRUE) {
  # Input validation
  if (!is.numeric(signal)) {
    stop("signal must be a numeric vector")
  }

  if (length(signal) < 2) {
    stop("signal must contain at least 2 points")
  }

  if (any(is.na(signal))) {
    stop("signal cannot contain NA values")
  }

  if (any(is.infinite(signal))) {
    stop("signal cannot contain infinite values")
  }

  if (!is.numeric(n_scales) || length(n_scales) != 1 || is.na(n_scales) || is.infinite(n_scales) || n_scales <= 0 || n_scales != as.integer(n_scales)) {
    stop("n_scales must be a positive integer")
  }
  n_scales <- as.integer(n_scales)

  if (!is.logical(detail) || length(detail) != 1 || is.na(detail)) {
    stop("detail must be a logical value (TRUE or FALSE)")
  }

  # Cubic B-spline filter kernel
  spline_kernel <- c(0.0625, 0.25, 0.375, 0.25, 0.0625)
  n_pts <- length(signal)
  max_allowable_scales <- ceiling(log2(n_pts))

  if (n_scales > max_allowable_scales) {
    stop(paste(
      "n_scales cannot exceed the max scales for this signal length:",
      max_allowable_scales
    ))
  }

  current_approx <- signal
  coeffs_matrix <- matrix(numeric(n_scales * n_pts), nrow = n_scales, ncol = n_pts)

  for (s_idx in 1:n_scales) {
    next_approx <- as.vector(wavelet_atrous_1d(current_approx, spline_kernel, -2, 2, s_idx))
    diff_detail <- current_approx - next_approx
    current_approx <- next_approx
    if (detail) {
      coeffs_matrix[s_idx, ] <- diff_detail
    } else {
      coeffs_matrix[s_idx, ] <- next_approx
    }
  }

  return(coeffs_matrix)
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
#' @seealso \code{\link{wavelet_spline_transform}}
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
      wavelets_polar <- wavelet_spline_transform(
        distances$polar$normalized, n_scales,
        detail = detail
      )
      wavelets_perimeter <- wavelet_spline_transform(
        distances$perimeter$normalized, n_scales,
        detail = detail
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
