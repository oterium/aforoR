# Native Geometric and Contour Shape Utilities (replacing external dependencies)

#' Calculate Polygon Area using Gauss's Area Formula (Shoelace Formula)
#'
#' @param coo Matrix or data.frame with X and Y columns.
#' @return Numeric scalar representing the polygon area.
#' @keywords internal
#' @noRd
.coo_area <- function(coo) {
  coo <- as.matrix(coo)
  n <- nrow(coo)
  if (n < 3) return(0)
  x <- coo[, 1]
  y <- coo[, 2]
  j <- c(2:n, 1)
  0.5 * abs(sum(x * y[j] - y * x[j]))
}

#' Calculate Polygon Perimeter
#'
#' Computes the sum of euclidean distances between consecutive vertices.
#'
#' @param coo Matrix or data.frame with X and Y columns.
#' @return Numeric scalar representing the perimeter.
#' @keywords internal
#' @noRd
.coo_perim <- function(coo) {
  coo <- as.matrix(coo)
  n <- nrow(coo)
  if (n < 2) return(0)
  dx <- diff(coo[, 1])
  dy <- diff(coo[, 2])
  sum(sqrt(dx^2 + dy^2))
}

#' Align Contour along Principal Axes using SVD
#'
#' Centers the contour and rotates it to align with its principal components.
#'
#' @param coo Matrix or data.frame with X and Y columns.
#' @return Matrix of aligned coordinates.
#' @keywords internal
#' @noRd
.coo_align <- function(coo) {
  coo <- as.matrix(coo)
  m <- colMeans(coo)
  cnt <- sweep(coo, 2, m, "-")
  s <- svd(stats::var(cnt))
  cnt %*% s$u
}

#' Calculate Major Axis Length (Length)
#'
#' @param coo Matrix or data.frame with X and Y columns.
#' @return Numeric scalar of length along major principal axis.
#' @keywords internal
#' @noRd
.coo_length <- function(coo) {
  aligned <- .coo_align(coo)
  r <- apply(aligned, 2, range)
  abs(r[2, 1] - r[1, 1])
}

#' Calculate Minor Axis Width (Width)
#'
#' @param coo Matrix or data.frame with X and Y columns.
#' @return Numeric scalar of width along minor principal axis.
#' @keywords internal
#' @noRd
.coo_width <- function(coo) {
  aligned <- .coo_align(coo)
  r <- apply(aligned, 2, range)
  abs(r[2, 2] - r[1, 2])
}

#' Maximum Feret Diameter (Calliper Length)
#'
#' @param coo Matrix or data.frame with X and Y columns.
#' @return Numeric scalar of maximum pairwise distance.
#' @keywords internal
#' @noRd
.coo_calliper <- function(coo) {
  coo <- as.matrix(coo)
  max(stats::dist(coo))
}

#' Resample Contour to N Points along Cumulative Arc Length
#'
#' Uniformly interpolates closed or open contours to exactly \code{n} points.
#'
#' @param coo Matrix or data.frame with X and Y columns.
#' @param n Target number of points.
#' @return Matrix of size \code{n x 2}.
#' @keywords internal
#' @noRd
.coo_interpolate <- function(coo, n) {
  coo <- as.matrix(coo)
  nr <- nrow(coo)
  if (nr == 0) return(matrix(numeric(0), nrow = 0, ncol = 2))
  if (nr == 1) return(matrix(rep(coo[1, ], n), nrow = n, byrow = TRUE))

  # Ensure closed loop for perimeter calculation
  if (!isTRUE(all.equal(coo[1, ], coo[nr, ]))) {
    coo_c <- rbind(coo, coo[1, ])
  } else {
    coo_c <- coo
  }

  dx <- diff(coo_c[, 1])
  dy <- diff(coo_c[, 2])
  dists <- sqrt(dx^2 + dy^2)
  orig <- c(0, cumsum(dists))
  total_perim <- orig[length(orig)]

  if (total_perim < 1e-12) {
    return(matrix(rep(coo[1, ], n), nrow = n, byrow = TRUE))
  }

  targ <- seq(0, total_perim, length.out = n + 1)[-(n + 1)]
  res <- matrix(NA_real_, nrow = n, ncol = 2)
  res[1, ] <- coo[1, ]

  for (i in 2:n) {
    k <- max(which(orig <= targ[i]))
    if (k >= nrow(coo_c)) {
      res[i, ] <- coo_c[nrow(coo_c), ]
    } else {
      seg_len <- orig[k + 1] - orig[k]
      r <- if (seg_len > 1e-12) (targ[i] - orig[k]) / seg_len else 0
      res[i, ] <- coo_c[k, ] + r * (coo_c[k + 1, ] - coo_c[k, ])
    }
  }
  colnames(res) <- c("X", "Y")
  res
}

#' Elliptic Fourier Analysis (Kuhl & Giardina 1982)
#'
#' Decomposes closed contour coordinates into harmonic elliptic Fourier coefficients.
#'
#' @param x Matrix or data.frame of contour coordinates.
#' @param nb.h Number of harmonics.
#' @param norm Logical; if TRUE, normalizes coefficients for size, rotation, and starting point.
#' @return List with harmonic coefficients \code{an}, \code{bn}, \code{cn}, \code{dn}.
#' @keywords internal
#' @noRd
.efourier <- function(x, nb.h = 12, norm = TRUE) {
  coo <- as.matrix(x)
  nr <- nrow(coo)
  if (isTRUE(all.equal(coo[1, ], coo[nr, ]))) {
    coo <- coo[-nr, , drop = FALSE]
    nr <- nrow(coo)
  }

  if (nb.h * 2 > nr) {
    nb.h <- max(1, floor(nr / 2))
  }

  Dx <- coo[, 1] - coo[, 1][c(nr, (1:(nr - 1)))]
  Dy <- coo[, 2] - coo[, 2][c(nr, (1:(nr - 1)))]
  Dt <- sqrt(Dx^2 + Dy^2)
  Dt[Dt < 1e-10] <- 1e-10
  t1 <- cumsum(Dt)
  t1m1 <- c(0, t1[-nr])
  T_total <- sum(Dt)

  an <- bn <- cn <- dn <- numeric(nb.h)
  for (i in 1:nb.h) {
    Ti <- (T_total / (2 * pi^2 * i^2))
    r <- 2 * i * pi
    an[i] <- Ti * sum((Dx / Dt) * (cos(r * t1 / T_total) - cos(r * t1m1 / T_total)))
    bn[i] <- Ti * sum((Dx / Dt) * (sin(r * t1 / T_total) - sin(r * t1m1 / T_total)))
    cn[i] <- Ti * sum((Dy / Dt) * (cos(r * t1 / T_total) - cos(r * t1m1 / T_total)))
    dn[i] <- Ti * sum((Dy / Dt) * (sin(r * t1 / T_total) - sin(r * t1m1 / T_total)))
  }
  ao <- 2 * sum(coo[, 1] * Dt / T_total)
  co <- 2 * sum(coo[, 2] * Dt / T_total)
  raw_res <- list(an = an, bn = bn, cn = cn, dn = dn, ao = ao, co = co)

  if (!norm) {
    return(raw_res)
  }

  # Normalization for size, rotation, and starting phase
  A1 <- an[1]
  B1 <- bn[1]
  C1 <- cn[1]
  D1 <- dn[1]

  theta <- (0.5 * atan(2 * (A1 * B1 + C1 * D1) / (A1^2 + C1^2 - B1^2 - D1^2))) %% pi
  phaseshift <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  M2 <- matrix(c(A1, C1, B1, D1), 2, 2) %*% phaseshift
  v <- colSums(M2^2)
  if (v[1] < v[2]) {
    theta <- theta + pi / 2
  }
  theta <- (theta + pi / 2) %% pi - pi / 2
  Aa <- A1 * cos(theta) + B1 * sin(theta)
  Cc <- C1 * cos(theta) + D1 * sin(theta)
  scale <- sqrt(Aa^2 + Cc^2)
  psi <- atan(Cc / Aa) %% pi
  if (Aa < 0) {
    psi <- psi + pi
  }
  size <- if (scale > 1e-12) 1 / scale else 1
  rotation <- matrix(c(cos(psi), -sin(psi), sin(psi), cos(psi)), 2, 2)

  A <- B <- C <- D <- numeric(nb.h)
  for (i in 1:nb.h) {
    mat <- size * rotation %*% matrix(c(an[i], cn[i], bn[i], dn[i]), 2, 2) %*%
           matrix(c(cos(i * theta), sin(i * theta), -sin(i * theta), cos(i * theta)), 2, 2)
    A[i] <- mat[1, 1]
    B[i] <- mat[1, 2]
    C[i] <- mat[2, 1]
    D[i] <- mat[2, 2]
  }

  list(an = A, bn = B, cn = C, dn = D, ao = ao, co = co,
       A = A, B = B, C = C, D = D, size = scale, theta = theta, psi = psi)
}

#' Generalized Procrustes Analysis (GPA) for Shape Alignment
#'
#' Registers a list of contours to a consensus shape via translation, scaling, and rotation.
#'
#' @param shapes List of 2D coordinate matrices.
#' @param max_iter Maximum alignment iterations.
#' @param tol Convergence tolerance.
#' @return List of aligned 2D coordinate matrices.
#' @keywords internal
#' @noRd
.procrustes_align <- function(shapes, max_iter = 30, tol = 1e-5) {
  n <- length(shapes)
  if (n <= 1) return(shapes)

  # Center and scale to centroid size = 1
  aligned <- lapply(shapes, function(s) {
    m <- colMeans(s)
    cnt <- sweep(s, 2, m, "-")
    cs <- sqrt(sum(cnt^2))
    if (cs < 1e-12) cnt else cnt / cs
  })

  mean_shape <- aligned[[1]]

  for (iter in seq_len(max_iter)) {
    diff_acc <- 0
    new_aligned <- vector("list", n)
    for (i in seq_len(n)) {
      si <- aligned[[i]]
      H <- t(mean_shape) %*% si
      sv <- svd(H)
      d <- sign(det(sv$v %*% t(sv$u)))
      V <- sv$v
      V[, 2] <- V[, 2] * d
      rot_mat <- V %*% t(sv$u)
      rot_si <- si %*% rot_mat
      diff_acc <- diff_acc + sum((rot_si - si)^2)
      new_aligned[[i]] <- rot_si
    }
    mean_mat <- Reduce("+", new_aligned) / n
    cs_m <- sqrt(sum(mean_mat^2))
    mean_shape <- if (cs_m < 1e-12) mean_mat else mean_mat / cs_m
    aligned <- new_aligned
    if (diff_acc < tol) break
  }
  aligned
}
