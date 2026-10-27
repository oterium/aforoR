#' Calculate Euclidean Distance Between Two Points
#'
#' Calculates the Euclidean distance between two points in 2D space.
#'
#' @param p1 A numeric vector of length 2 representing the first point (x, y).
#' @param p2 A numeric vector of length 2 representing the second point (x, y).
#' @return A numeric value representing the Euclidean distance between the two points.
#' @export
#' @examples
#' # Example usage:
#' point1 <- c(0, 0)
#' point2 <- c(3, 4)
#' distance <- ild(point1, point2)
#' print(distance) # Should be 5
ild <- function(p1, p2) {
  # Input validation
  if (!is.numeric(p1) || !is.numeric(p2)) {
    stop("Both p1 and p2 must be numeric vectors")
  }

  if (length(p1) != 2 || length(p2) != 2) {
    stop("Both p1 and p2 must be vectors of length 2")
  }

  if (any(is.na(p1)) || any(is.na(p2))) {
    stop("Points cannot contain NA values")
  }

  if (any(is.infinite(p1)) || any(is.infinite(p2))) {
    stop("Points cannot contain infinite values")
  }

  sqrt(sum((p1 - p2)^2))
}

#' Compute Distances to Centroid (Polar Resampling)
#'
#' Calculates radial distances from the centroid to a specified number of points
#' sampled along the contour using polar coordinates. This method ensures regular
#' angular sampling, which is a prerequisite for standard AFORO wavelet analysis.
#'
#' @param Rx A numeric vector of x-coordinates.
#' @param Ry A numeric vector of y-coordinates.
#' @param n Number of points to sample (e.g., 512).
#' @return A list containing:
#'   \itemize{
#'     \item \code{pixindices}: Indices of the sampled points in the original contour.
#'     \item \code{radii}: Distances from the centroid to each sampled point.
#'     \item \code{coord}: Coordinates of the sampled points relative to the centroid.
#'   }
#' @export
#' @examples
#' # Example using sample data
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' contour <- extract_contour(preprocess_image(image_path))
#' polar_dists <- regularradius(contour[, 1], contour[, 2], n = 512)
#' plot(polar_dists$radii, type = "l", main = "Radial Distances")
regularradius <- function(Rx, Ry, n) {
  # Input validation
  if (!is.numeric(Rx) || !is.numeric(Ry)) {
    stop("Rx and Ry must be numeric vectors")
  }

  if (length(Rx) != length(Ry)) {
    stop("Rx and Ry must have the same length")
  }

  if (!is.numeric(n) || length(n) != 1 || is.na(n) || is.infinite(n) || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer")
  }

  if (any(is.na(Rx)) || any(is.na(Ry))) {
    stop("Rx and Ry cannot contain NA values")
  }

  if (any(is.infinite(Rx)) || any(is.infinite(Ry))) {
    stop("Rx and Ry cannot contain infinite values")
  }

  if (length(Rx) < 3) {
    stop("Rx and Ry must contain at least 3 points")
  }

  if (all(Rx == Rx[1]) && all(Ry == Ry[1])) {
    stop("Contour is degenerate with identical coordinates")
  }

  n <- as.integer(n)
  le <- length(Rx)

  M1 <- matrix(c(Rx - mean(Rx), Ry - mean(Ry)), le, 2)
  V1 <- complex(real = M1[, 1], imaginary = M1[, 2])
  M2 <- matrix(c(Arg(V1), Mod(V1)), le, 2)
  V2 <- integer(n)
  for (i in 0:(n - 1)) {
    V2[i + 1] <- which.max((cos(M2[, 1] - 2 * i * pi / n)))
  }
  V2 <- sort(V2)
  list("pixindices" = V2, "radii" = M2[V2, 2], "coord" = M1[V2, , drop = FALSE])
}

#' Compute Perimeter Distances (Curvilinear Resampling)
#'
#' Calculates distances from the centroid to points sampled at regular intervals
#' along the perimeter of the shape. This provides an alternative to polar sampling
#' for shapes with complex bounds.
#'
#' @param x A numeric vector of x-coordinates.
#' @param y A numeric vector of y-coordinates.
#' @param n Number of points to sample along the perimeter.
#' @return A list containing:
#'   \itemize{
#'     \item \code{dist}: Radial distances from the centroid to the sampled points.
#'     \item \code{coords}: Coordinates of the sampled points.
#'   }
#' @export
#' @examples
#' # Example using sample data
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' contour <- extract_contour(preprocess_image(image_path))
#' perim_dists <- dper(contour[, 1], contour[, 2], n = 512)
#' plot(perim_dists$dist, type = "l", main = "Perimeter-based Distances")
dper <- function(x, y, n) {
  # Input validation
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors")
  }

  if (length(x) != length(y)) {
    stop("x and y must have the same length")
  }

  if (!is.numeric(n) || length(n) != 1 || is.na(n) || is.infinite(n) || n < 2 || n != as.integer(n)) {
    stop("n must be an integer >= 2")
  }

  if (any(is.na(x)) || any(is.na(y))) {
    stop("x and y cannot contain NA values")
  }

  if (any(is.infinite(x)) || any(is.infinite(y))) {
    stop("x and y cannot contain infinite values")
  }

  if (length(x) < 3) {
    stop("Need at least 3 points to calculate perimeter distances")
  }

  if (all(x == x[1]) && all(y == y[1])) {
    stop("Contour is degenerate with identical coordinates")
  }

  # Sample n points along the perimeter at strictly equidistant metric arc lengths
  coo_sampled <- .coo_interpolate(cbind(x, y), n = as.integer(n))
  centroid <- c(mean(coo_sampled[, 1]), mean(coo_sampled[, 2]))
  distances <- sqrt((coo_sampled[, 1] - centroid[1])^2 + (coo_sampled[, 2] - centroid[2])^2)

  return(list(dist = distances, coords = coo_sampled))
}

#' Calculate Distance Measures from Contour
#'
#' Computes two types of distance measures from an otolith contour: polar
#' distances (radial distances from the centroid) and perimeter-based distances.
#' These measures are crucial for subsequent wavelet analysis.
#'
#' @param contour A matrix of coordinates from \code{extract_contour()}.
#' @param n_points Target number of points for resampling (default: 512).
#' @param start_point Strategy or coordinates for selecting the contour starting point (rostrum):
#'   \itemize{
#'     \item \code{"auto"} (default, or \code{NULL}): Automatically selects the rightmost point (\eqn{\max(X)}).
#'     \item \code{"manual"}: Interactive mouse selection mode. Plots the contour with the automatic
#'           rostrum candidate (\code{o}) and prompts the user to click on the desired starting point
#'           with the mouse (snapping to the nearest contour vertex).
#'     \item Numeric coordinate pair \code{c(x, y)} (or 1x2 matrix) or integer index \code{1..nrow(contour)}:
#'           Starts at the contour point closest to \code{(x, y)} or at the specified row index.
#'   }
#' @return A list containing:
#'   \itemize{
#'     \item \code{polar}: Radii, coordinates, and normalized radii.
#'     \item \code{perimeter}: Distances, coordinates, and normalized distances.
#'     \item \code{reordered_coords}: The contour coordinates reordered from the
#'           starting point.
#'   }
#' @export
#' @examples
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' binary_img <- preprocess_image(image_path)
#' contour <- extract_contour(binary_img)
#' dists <- calculate_distances(contour, n_points = 512)
#' plot(dists$polar$normalized, type = "l")
calculate_distances <- function(contour, n_points = 512, start_point = "auto") {
  if (!is.numeric(n_points) || length(n_points) != 1 || is.na(n_points) ||
      is.infinite(n_points) || n_points < 2 || n_points != as.integer(n_points)) {
    stop("n_points must be an integer >= 2")
  }

  resolved <- .resolve_contour_and_start(contour, start_point = start_point, remove_duplicate_end = FALSE)
  x <- unname(resolved$coords[, 1])
  y <- unname(resolved$coords[, 2])

  tryCatch(
    {
      dist <- regularradius(x, y, n = n_points)
      dist2 <- dper(x, y, n_points)
      max_r1 <- max(dist$radii)
      max_r2 <- max(dist2$dist)
      dist.norm <- if (max_r1 > 0) dist$radii / max_r1 else rep(0, length(dist$radii))
      dist.norm2 <- if (max_r2 > 0) dist2$dist / max_r2 else rep(0, length(dist2$dist))

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
