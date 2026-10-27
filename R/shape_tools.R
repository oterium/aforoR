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

#' Validate Contour Matrix and Resolve Rostrum Starting Point
#'
#' Validates that a contour has at least 3 finite, non-NA, non-degenerate numeric
#' points and reorders the contour vertices so row 1 corresponds to the resolved
#' starting point (\code{"auto"}, \code{"manual"}, index, or \code{c(x, y)}).
#'
#' @param contour A matrix or data.frame with X and Y columns.
#' @param start_point Starting point mode or coordinates (\code{"auto"}, \code{"manual"}, index, or \code{c(x, y)}).
#' @param remove_duplicate_end Logical; if TRUE, drops a duplicated closing vertex when row 1 equals row N.
#' @return A list with \code{coords} (reordered N x 2 numeric matrix with colnames X, Y),
#'         \code{start_idx} (integer index in original contour), and
#'         \code{start_point_coords} (named numeric vector \code{c(X, Y)}).
#' @keywords internal
#' @noRd
.resolve_contour_and_start <- function(contour, start_point = "auto", remove_duplicate_end = FALSE) {
  if (is.null(contour)) stop("contour cannot be NULL")

  if (is.data.frame(contour)) {
    if (all(c("X", "Y") %in% colnames(contour))) {
      coo <- as.matrix(contour[, c("X", "Y")])
    } else if (all(c("x", "y") %in% colnames(contour))) {
      coo <- as.matrix(contour[, c("x", "y")])
      colnames(coo) <- c("X", "Y")
    } else if (ncol(contour) >= 2) {
      coo <- as.matrix(contour[, 1:2])
      colnames(coo) <- c("X", "Y")
    } else {
      stop("contour must have at least 2 columns")
    }
  } else if (is.matrix(contour)) {
    if (ncol(contour) >= 2) {
      coo <- contour[, 1:2, drop = FALSE]
      colnames(coo) <- c("X", "Y")
    } else {
      stop("contour must have at least 2 columns")
    }
  } else {
    stop("contour must be a data.frame or matrix with X and Y columns")
  }

  if (nrow(coo) < 3) {
    stop("contour must contain at least 3 points")
  }

  if (!is.numeric(coo)) {
    stop("contour coordinates must be numeric")
  }

  if (any(is.na(coo))) {
    stop("contour coordinates cannot contain NA values")
  }

  if (any(is.infinite(coo))) {
    stop("contour coordinates cannot contain infinite values")
  }

  n_orig <- nrow(coo)
  had_dup_end <- FALSE
  if (remove_duplicate_end && nrow(coo) > 1 && isTRUE(all(coo[1, ] == coo[nrow(coo), ]))) {
    coo <- coo[-nrow(coo), , drop = FALSE]
    had_dup_end <- TRUE
  }

  if (nrow(coo) < 3) {
    stop("contour must contain at least 3 points")
  }

  if (all(coo[, 1] == coo[1, 1]) && all(coo[, 2] == coo[1, 2])) {
    stop("contour is degenerate with identical coordinates")
  }

  n_raw <- nrow(coo)
  if (is.matrix(start_point) || is.data.frame(start_point)) {
    if (nrow(start_point) != 1L || ncol(start_point) < 1L) {
      stop("start_point matrix or data.frame must have exactly 1 row and 1 or 2 columns")
    }
    start_point <- as.numeric(start_point[1, ])
  }

  if (is.null(start_point) || isTRUE(start_point) ||
      (is.character(start_point) && length(start_point) == 1 && tolower(start_point) == "auto")) {
    x_max_idx <- which.max(coo[, 1])
    x_max <- coo[x_max_idx, 1]
    x_mean <- mean(coo[, 1])
    rightmost_candidates <- which(coo[, 1] == x_max & coo[, 1] > x_mean)
    start_idx <- if (length(rightmost_candidates) == 0) x_max_idx else rightmost_candidates[1]
  } else if (is.character(start_point) && length(start_point) == 1 && tolower(start_point) %in% c("manual", "interactive")) {
    auto_idx <- which.max(coo[, 1])
    if (!interactive() && is.null(getOption("aforoR.SimLocator"))) {
      message("Non-interactive session: using automatic rightmost point.")
      start_idx <- auto_idx
    } else {
      graphics::plot(coo[, 1], coo[, 2], type = "l", col = "#00E5FF", lwd = 2, asp = 1,
                     xlab = "X", ylab = "Y", main = "Select Starting Point (Rostrum)")
      graphics::points(coo[auto_idx, 1], coo[auto_idx, 2], pch = 1, col = "#FFEA00", lwd = 2.5, cex = 1.8)
      graphics::legend("topleft",
                       legend = c("Click on Rostrum (starting point)", "Right-click / Esc = Keep Auto (o)"),
                       col = c("red", "#FFEA00"), pch = c(16, 1),
                       bg = grDevices::rgb(0, 0, 0, 0.7), text.col = "white", cex = 0.85, inset = c(0.02, 0.02))
      sim_loc <- getOption("aforoR.SimLocator")
      pt <- tryCatch(if (is.function(sim_loc)) sim_loc(n = 1) else graphics::locator(n = 1), error = function(e) NULL)
      if (is.null(pt) || length(pt$x) == 0) {
        start_idx <- auto_idx
      } else {
        d_sq <- (coo[, 1] - pt$x[1])^2 + (coo[, 2] - pt$y[1])^2
        start_idx <- which.min(d_sq)
        graphics::points(coo[start_idx, 1], coo[start_idx, 2], pch = 16, col = "red", cex = 1.6)
      }
    }
  } else {
    if (!is.numeric(start_point) || any(is.na(start_point)) || any(is.infinite(start_point))) {
      stop("start_point must be 'auto', 'manual', or a valid numeric index or coordinate pair c(x, y) without NA/Inf")
    }
    if (length(start_point) == 1) {
      if (start_point != as.integer(start_point)) {
        stop(paste("start_point index must be an integer between 1 and", n_raw))
      }
      start_idx <- as.integer(start_point)
      if (had_dup_end && start_idx == n_orig) {
        start_idx <- 1L
      }
      if (start_idx < 1L || start_idx > n_raw) {
        stop(paste("start_point index must be between 1 and", n_raw))
      }
    } else if (length(start_point) == 2) {
      dists_sq <- (coo[, 1] - start_point[1])^2 + (coo[, 2] - start_point[2])^2
      start_idx <- which.min(dists_sq)
    } else {
      stop("start_point must be an index (length 1) or coordinates c(x, y) (length 2)")
    }
  }

  if (start_idx > 1L) {
    indices <- c(start_idx:n_raw, 1:(start_idx - 1L))
    coo <- coo[indices, , drop = FALSE]
  }
  start_pt_coords <- c(X = unname(coo[1, 1]), Y = unname(coo[1, 2]))

  list(
    coords = coo,
    start_idx = start_idx,
    start_point_coords = start_pt_coords
  )
}
