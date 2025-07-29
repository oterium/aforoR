#' Calculate Euclidean Distance Between Two Points
#'
#' Calculates the Euclidean distance between two points in 2D space.
#'
#' @param p1 A numeric vector of length 2 representing the first point (x, y).
#' @param p2 A numeric vector of length 2 representing the second point (x, y).
#' @return A numeric value representing the Euclidean distance between the two points.
#' @examples
#' # Example usage:
#' point1 <- c(0, 0)
#' point2 <- c(3, 4)
#' distance <- ild(point1, point2)
#' print(distance)  # Should be 5
ild <- function(p1, p2) {
  sqrt(sum((p1 - p2)^2))
}

#' Compute Distances to Centroid
#'
#' Calculates the distances from points to the centroid and returns the indices, radii, and coordinates of the points.
#'
#' @param Rx A numeric vector representing the x-coordinates of the points.
#' @param Ry A numeric vector representing the y-coordinates of the points.
#' @param n An integer specifying the number of points to sample.
#' @return A list containing:
#'   - `pixindices`: A numeric vector of indices of the sampled points.
#'   - `radii`: A numeric vector of radii from the centroid to the sampled points.
#'   - `coord`: A matrix of coordinates of the sampled points.
#' @export
#' @examples
#' # Example usage:
#' Rx <- rnorm(100)
#' Ry <- rnorm(100)
#' n <- 50
#' result <- regularradius(Rx, Ry, n)
#' print(result)
regularradius <- function(Rx, Ry, n) {
  le <- length(Rx)
  M <- matrix(c(Rx, Ry), le, 2)
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

#' Compute Perimeter Distances
#'
#' Calculates the distances from the centroid to points along the perimeter of a shape.
#'
#' @param x A numeric vector representing the x-coordinates of the perimeter points.
#' @param y A numeric vector representing the y-coordinates of the perimeter points.
#' @param n An integer specifying the number of points to sample along the perimeter.
#' @return A list containing:
#'   - `dist`: A numeric vector of distances from the centroid to the sampled points.
#'   - `coords`: A matrix of coordinates of the sampled points.
#' @export
#' @examples
#' # Example usage:
#' x <- rnorm(100)
#' y <- rnorm(100)
#' n <- 50
#' result <- dper(x, y, n)
#' print(result)
dper <- function(x, y, n) {
  # x values ordered
  # y values ordered
  # n values of interest (e.g., 512)

  xr <- (x[round(seq(1, length(x), length.out = n + 1))])[-1]
  yr <- (y[round(seq(1, length(y), length.out = n + 1))])[-1]

  distancias <- NULL
  for (i in 1:n) {
    distancias[i] <- ild(c(mean(x), mean(y)), c(xr[i], yr[i]))
    # distancias[i] <- ild(c(xr[1], yr[1], c(xr[i], yr[i]))
  }
  return(list(dist = distancias, coords = cbind(xr, yr)))
}
