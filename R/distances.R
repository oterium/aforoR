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

  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("n must be a positive integer")
  }

  if (any(is.na(Rx)) || any(is.na(Ry))) {
    stop("Rx and Ry cannot contain NA values")
  }

  le <- length(Rx)

  M1 <- matrix(c(Rx - mean(Rx), Ry - mean(Ry)), le, 2)
  V1 <- complex(real = M1[, 1], imaginary = M1[, 2])
  M2 <- matrix(c(Arg(V1), Mod(V1)), le, 2)
  V2 <- NA
  for (i in 0:(n - 1)) {
    V2[i + 1] <- which.max((cos(M2[, 1] - 2 * i * pi / n)))
  }
  V2 <- sort(V2)
  list("pixindices" = V2, "radii" = M2[V2, 2], "coord" = M1[V2, ])
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

  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("n must be a positive integer")
  }

  if (any(is.na(x)) || any(is.na(y))) {
    stop("x and y cannot contain NA values")
  }

  if (length(x) < 2) {
    stop("Need at least 2 points to calculate perimeter distances")
  }

  # x values ordered
  # y values ordered
  # n values of interest (e.g., 512)

  xr <- (x[round(seq(1, length(x), length.out = n + 1))])[-1]
  yr <- (y[round(seq(1, length(y), length.out = n + 1))])[-1]

  distancias <- NULL
  for (i in 1:n) {
    distancias[i] <- ild(c(mean(x), mean(y)), c(xr[i], yr[i]))
    # Alternative: distancias[i] <- ild(c(xr[1], yr[1]), c(xr[i], yr[i]))
  }
  return(list(dist = distancias, coords = cbind(xr, yr)))
}
