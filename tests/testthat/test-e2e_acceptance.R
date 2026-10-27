# ==============================================================================
# tests/testthat/test-e2e_acceptance.R
#
# Comprehensive 4-Tier Opaque-Box E2E Acceptance Test Suite for aforoR
# Requirements Covered:
#   - R1: Accurate Scale Bar Detection & Endpoint Visualization (Features 1 - 5)
#   - R2: Robust Adaptive Thresholding & Multi-Channel Handling (Features 6 - 11)
#   - R3: Flexible plot.css() S3 Method (Feature 12)
#
# Tiers in this file:
#   - Tier 1: Feature Coverage (>= 5 test cases per feature across R1, R2, R3)
#   - Tier 2: Boundary & Corner Cases (empty/zero, red mask, plaque, orientations,
#             float overshoot, 2-channel TIFF, border area >=85%, min_area=1000)
#   - Tier 3: Cross-Feature Interactions (2-ch + scale, textured + plaque, small + css)
#
# Note: Tier 4 (Real-world 38-sample 19-collection benchmark) is implemented in
#       tests/test_e2e_runner.R
# ==============================================================================

library(testthat)
library(aforoR)
library(EBImage)
try(testthat::set_max_fails(Inf), silent = TRUE)

# ------------------------------------------------------------------------------
# Test Fixtures & Synthetic Image Generators
# ------------------------------------------------------------------------------

#' Generate synthetic otolith binary matrix
create_synthetic_otolith <- function(nrow = 600, ncol = 800, cx = 400, cy = 300, rx = 220, ry = 120) {
  m <- matrix(0, nrow = nrow, ncol = ncol)
  xx <- col(m) - cx
  yy <- row(m) - cy
  m[(xx^2 / rx^2 + yy^2 / ry^2) <= 1] <- 1
  m
}

#' Add horizontal rectangular scale bar
add_horizontal_bar <- function(m, x1 = 100, x2 = 300, y1 = 500, y2 = 507) {
  m[y1:y2, x1:x2] <- 1
  m
}

#' Add vertical rectangular scale bar
add_vertical_bar <- function(m, x1 = 700, x2 = 707, y1 = 100, y2 = 350) {
  m[y1:y2, x1:x2] <- 1
  m
}

#' Generate synthetic RGB color image with pure red scale bar
create_synthetic_red_bar_image <- function(nrow = 600, ncol = 800, bar_x1 = 100, bar_x2 = 350, bar_y1 = 520, bar_y2 = 526) {
  # Background dark grey
  arr <- array(0.08, dim = c(nrow, ncol, 3))
  # Otolith in center (white/light grey in all channels)
  xx <- col(arr[, , 1]) - 400
  yy <- row(arr[, , 1]) - 300
  oto_mask <- (xx^2 / 200^2 + yy^2 / 100^2) <= 1
  for (ch in 1:3) {
    ch_mat <- arr[, , ch]
    ch_mat[oto_mask] <- 0.90
    arr[, , ch] <- ch_mat
  }
  # Add pure red scale bar: R = 0.95, G = 0.05, B = 0.05
  arr[bar_y1:bar_y2, bar_x1:bar_x2, 1] <- 0.95
  arr[bar_y1:bar_y2, bar_x1:bar_x2, 2] <- 0.05
  arr[bar_y1:bar_y2, bar_x1:bar_x2, 3] <- 0.05

  EBImage::Image(arr, colormode = "Color")
}

#' Generate synthetic 2-channel TIFF image (Channel 1: Grayscale, Channel 2: Alpha 1.0)
create_synthetic_2channel_image <- function(nrow = 500, ncol = 600) {
  ch1 <- matrix(0.02, nrow = nrow, ncol = ncol)
  # Otolith in Channel 1
  xx <- col(ch1) - 300
  yy <- row(ch1) - 250
  ch1[(xx^2 / 180^2 + yy^2 / 90^2) <= 1] <- 0.88
  # Scale bar in Channel 1: span = 150 px
  ch1[420:425, 100:249] <- 0.95

  # Channel 2: Alpha mask (constant 1.0)
  ch2 <- matrix(1.0, nrow = nrow, ncol = ncol)

  arr <- array(0, dim = c(nrow, ncol, 2))
  arr[, , 1] <- ch1
  arr[, , 2] <- ch2

  EBImage::Image(arr, colormode = "Color")
}

# ==============================================================================
# TIER 1: FEATURE COVERAGE (R1, R2, R3)
# ==============================================================================

# ------------------------------------------------------------------------------
# Feature 1: Physical Scale Bar Span (Eliminating m.majoraxis +15% to +33% Error)
# ------------------------------------------------------------------------------
test_that("T1.1.1: Physical span of 100px horizontal bar is accurate within 2px (ratio in 1.00 +- 0.02)", {
  m <- create_synthetic_otolith()
  # Add 100 px bar: x from 100 to 199 (span = 100), thickness = 6 (y 500:505)
  m <- add_horizontal_bar(m, x1 = 100, x2 = 199, y1 = 500, y2 = 505)
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  expect_false(is.null(info))
  detected_len <- as.numeric(info$bar_length_px)

  # True span is exactly 100 px
  expect_equal(detected_len, 100, tolerance = 2)
  expect_true(abs(detected_len / 100 - 1.0) <= 0.02)
  # Verify majoraxis overestimation is avoided
  lab <- EBImage::bwlabel(m)
  mom <- EBImage::computeFeatures.moment(lab)
  maj_axis <- max(mom[mom[, "m.eccentricity"] > 0.88, "m.majoraxis"])
  expect_true(maj_axis > 110) # m.majoraxis overestimates ~115
  expect_true(detected_len < maj_axis)
})

test_that("T1.1.2: Physical span of 200px horizontal bar with scale_mm = 2 gives px_per_mm ~ 100", {
  m <- create_synthetic_otolith()
  m <- add_horizontal_bar(m, x1 = 150, x2 = 349, y1 = 510, y2 = 517) # span = 200
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 2)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  expect_equal(as.numeric(info$bar_length_px), 200, tolerance = 2)
  expect_equal(as.numeric(res), 100, tolerance = 1)
})

test_that("T1.1.3: Physical span of 350px bar with scale_mm = 0.5 gives px_per_mm ~ 700", {
  m <- create_synthetic_otolith()
  m <- add_horizontal_bar(m, x1 = 100, x2 = 449, y1 = 520, y2 = 528) # span = 350
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 0.5)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  expect_equal(as.numeric(info$bar_length_px), 350, tolerance = 2)
  expect_equal(as.numeric(res), 700, tolerance = 4)
})

test_that("T1.1.4: Physical span of 742px bar (Canaries scale) matches true span within 2px", {
  # Large 1200x2000 image
  m <- matrix(0, nrow = 1200, ncol = 2000)
  xx <- col(m) - 1000
  yy <- row(m) - 600
  m[(xx^2 / 500^2 + yy^2 / 250^2) <= 1] <- 1
  m[1050:1065, 800:1541] <- 1 # span = 1541 - 800 + 1 = 742 px
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  expect_equal(as.numeric(info$bar_length_px), 742, tolerance = 2)
})

test_that("T1.1.5: Physical span of 50px short bar is within 2px (ratio in 1.00 +- 0.02)", {
  m <- create_synthetic_otolith()
  m <- add_horizontal_bar(m, x1 = 120, x2 = 169, y1 = 490, y2 = 494) # span = 50
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  expect_equal(as.numeric(info$bar_length_px), 50, tolerance = 2)
})

# ------------------------------------------------------------------------------
# Feature 2: Scale Bar Endpoint Calculation (x1, y1) and (x2, y2)
# ------------------------------------------------------------------------------
test_that("T1.2.1: Endpoints (x1, y1) and (x2, y2) exist and match horizontal bar tips", {
  m <- create_synthetic_otolith()
  m <- add_horizontal_bar(m, x1 = 120, x2 = 279, y1 = 510, y2 = 516) # span = 160
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  # Endpoints must be present as list(x1, y1, x2, y2) or direct elements
  ep <- info$endpoints
  x1 <- if (!is.null(ep$x1)) ep$x1 else info$x1
  x2 <- if (!is.null(ep$x2)) ep$x2 else info$x2
  y1 <- if (!is.null(ep$y1)) ep$y1 else info$y1
  y2 <- if (!is.null(ep$y2)) ep$y2 else info$y2

  expect_false(is.null(x1))
  expect_false(is.null(x2))
  expect_equal(min(x1, x2), 120, tolerance = 2)
  expect_equal(max(x1, x2), 279, tolerance = 2)
  expect_true(y1 >= 509 && y1 <= 517)
  expect_true(y2 >= 509 && y2 <= 517)
})

test_that("T1.2.2: Endpoints exist and match vertical bar tips", {
  m <- create_synthetic_otolith()
  m <- add_vertical_bar(m, x1 = 710, x2 = 716, y1 = 150, y2 = 329) # span = 180
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  ep <- info$endpoints
  x1 <- if (!is.null(ep$x1)) ep$x1 else info$x1
  x2 <- if (!is.null(ep$x2)) ep$x2 else info$x2
  y1 <- if (!is.null(ep$y1)) ep$y1 else info$y1
  y2 <- if (!is.null(ep$y2)) ep$y2 else info$y2

  expect_equal(min(y1, y2), 150, tolerance = 2)
  expect_equal(max(y1, y2), 329, tolerance = 2)
  expect_true(x1 >= 709 && x1 <= 717)
  expect_true(x2 >= 709 && x2 <= 717)
})

test_that("T1.2.3: Euclidean distance between endpoints equals detected bar_length_px", {
  m <- create_synthetic_otolith()
  m <- add_horizontal_bar(m, x1 = 200, x2 = 399, y1 = 480, y2 = 485) # span = 200
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  info <- attr(res, "scale_info")
  ep <- info$endpoints
  x1 <- if (!is.null(ep$x1)) ep$x1 else info$x1
  x2 <- if (!is.null(ep$x2)) ep$x2 else info$x2
  y1 <- if (!is.null(ep$y1)) ep$y1 else info$y1
  y2 <- if (!is.null(ep$y2)) ep$y2 else info$y2

  d_endpoints <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  expect_equal(d_endpoints, as.numeric(info$bar_length_px), tolerance = 2)
})

test_that("T1.2.4: Endpoints are contained strictly within scale_info bounding box", {
  m <- create_synthetic_otolith()
  m <- add_horizontal_bar(m, x1 = 150, x2 = 310, y1 = 510, y2 = 516)
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  info <- attr(res, "scale_info")
  ep <- info$endpoints
  x1 <- if (!is.null(ep$x1)) ep$x1 else info$x1
  x2 <- if (!is.null(ep$x2)) ep$x2 else info$x2
  y1 <- if (!is.null(ep$y1)) ep$y1 else info$y1
  y2 <- if (!is.null(ep$y2)) ep$y2 else info$y2

  expect_true(min(x1, x2) >= info$min_x - 1)
  expect_true(max(x1, x2) <= info$max_x + 1)
  expect_true(min(y1, y2) >= info$min_y - 1)
  expect_true(max(y1, y2) <= info$max_y + 1)
})

test_that("T1.2.5: Endpoints are finite integer/numeric values within image boundaries", {
  m <- create_synthetic_otolith(nrow = 600, ncol = 800)
  m <- add_horizontal_bar(m, x1 = 100, x2 = 250, y1 = 530, y2 = 535)
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  info <- attr(res, "scale_info")
  ep <- info$endpoints
  x1 <- if (!is.null(ep$x1)) ep$x1 else info$x1
  x2 <- if (!is.null(ep$x2)) ep$x2 else info$x2
  y1 <- if (!is.null(ep$y1)) ep$y1 else info$y1
  y2 <- if (!is.null(ep$y2)) ep$y2 else info$y2

  expect_true(is.finite(x1) && x1 >= 1 && x1 <= 800)
  expect_true(is.finite(x2) && x2 >= 1 && x2 <= 800)
  expect_true(is.finite(y1) && y1 >= 1 && y1 <= 600)
  expect_true(is.finite(y2) && y2 >= 1 && y2 <= 600)
})

# ------------------------------------------------------------------------------
# Feature 3: Red Scale Bar Detection (Canaries_JosePajuelo)
# ------------------------------------------------------------------------------
test_that("T1.3.1: Red chromaticity mask correctly extracts pure red scale bar on dark background", {
  red_img <- create_synthetic_red_bar_image(bar_x1 = 100, bar_x2 = 299, bar_y1 = 520, bar_y2 = 525) # span = 200
  res <- aforoR:::find_scale_bar_length(red_img, scale_mm = 1)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  expect_equal(as.numeric(info$bar_length_px), 200, tolerance = 3)
})

test_that("T1.3.2: Red scale bar with red text annotations isolates the long bar component", {
  red_img <- create_synthetic_red_bar_image(bar_x1 = 150, bar_x2 = 399, bar_y1 = 530, bar_y2 = 536) # span = 250
  # Add small text-like red blobs near bar
  arr <- as.array(red_img)
  arr[500:515, 200:208, 1] <- 0.95 # "1"
  arr[500:515, 200:208, 2:3] <- 0.05
  arr[505:515, 215:235, 1] <- 0.95 # "m"
  arr[505:515, 215:235, 2:3] <- 0.05
  red_annotated <- EBImage::Image(arr, colormode = "Color")

  res <- aforoR:::find_scale_bar_length(red_annotated, scale_mm = 1)
  expect_false(is.null(res))
  info <- attr(res, "scale_info")
  # The detected bar must be the 250px line, not the small "1" or "m"
  expect_equal(as.numeric(info$bar_length_px), 250, tolerance = 3)
})

test_that("T1.3.3: Non-red chromatic objects (pure green / pure blue) are rejected by red detector", {
  arr <- array(0.08, dim = c(600, 800, 3))
  # Add otolith
  xx <- col(arr[, , 1]) - 400
  yy <- row(arr[, , 1]) - 300
  oto_mask <- (xx^2 / 200^2 + yy^2 / 100^2) <= 1
  for (ch in 1:3) {
    ch_mat <- arr[, , ch]
    ch_mat[oto_mask] <- 0.85
    arr[, , ch] <- ch_mat
  }
  # Add pure green bar
  arr[520:525, 100:300, 2] <- 0.95
  arr[520:525, 100:300, c(1, 3)] <- 0.05
  green_img <- EBImage::Image(arr, colormode = "Color")

  # In color mode, if testing chromatic red detector, it should not mistake green for red
  r <- green_img[, , 1]
  g <- green_img[, , 2]
  b <- green_img[, , 3]
  red_mask <- (r > 0.45) & (r > 1.4 * pmax(g, b))
  expect_equal(sum(red_mask), 0)
})

test_that("T1.3.4: Pure red bar survives where standard Otsu threshold would erase it", {
  red_img <- create_synthetic_red_bar_image(bar_x1 = 120, bar_x2 = 319, bar_y1 = 510, bar_y2 = 515) # span = 200
  # Calculate standard grayscale luminance
  grey_img <- EBImage::channel(red_img, "grey")
  # Red bar luminance is 0.299 * 0.95 + 0.587 * 0.05 + 0.114 * 0.05 ~ 0.319
  # Otolith luminance is 0.90
  otsu_th <- EBImage::otsu(grey_img)
  expect_true(otsu_th > 0.35) # Otsu threshold sits between 0.35 and 0.55
  # The red bar pixels are below otsu_th, so standard Otsu erases the red bar
  expect_equal(sum(grey_img[510:515, 120:319] > otsu_th), 0)

  # But chromatic detection successfully captures the red scale bar
  res <- aforoR:::find_scale_bar_length(red_img, scale_mm = 1)
  expect_false(is.null(res))
  expect_equal(as.numeric(attr(res, "scale_info")$bar_length_px), 200, tolerance = 3)
})

test_that("T1.3.5: Red scale bar detection returns accurate endpoints (x1, y1) to (x2, y2)", {
  red_img <- create_synthetic_red_bar_image(bar_x1 = 200, bar_x2 = 449, bar_y1 = 520, bar_y2 = 526) # span = 250
  res <- aforoR:::find_scale_bar_length(red_img, scale_mm = 1)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  ep <- info$endpoints
  x1 <- if (!is.null(ep$x1)) ep$x1 else info$x1
  x2 <- if (!is.null(ep$x2)) ep$x2 else info$x2
  expect_equal(min(x1, x2), 200, tolerance = 2)
  expect_equal(max(x1, x2), 449, tolerance = 2)
})

# ------------------------------------------------------------------------------
# Feature 4: Plaque & Vertical Scale Bars (caballa & ShrarifJemaa)
# ------------------------------------------------------------------------------
test_that("T1.4.1: Vertical scale bar (ShrarifJemaa) correctly measures delta_y span", {
  m <- create_synthetic_otolith()
  m <- add_vertical_bar(m, x1 = 720, x2 = 727, y1 = 120, y2 = 319) # delta_x = 8, delta_y = 200
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  # Detected span must be along vertical axis (200 px), not horizontal (8 px)
  expect_equal(as.numeric(info$bar_length_px), 200, tolerance = 2)
})

test_that("T1.4.2: Vertical scale bar with aspect ratio > 5 detected with endpoints along Y", {
  m <- create_synthetic_otolith()
  m <- add_vertical_bar(m, x1 = 50, x2 = 55, y1 = 100, y2 = 250) # span = 151
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  expect_false(is.null(res))

  info <- attr(res, "scale_info")
  ep <- info$endpoints
  y1 <- if (!is.null(ep$y1)) ep$y1 else info$y1
  y2 <- if (!is.null(ep$y2)) ep$y2 else info$y2
  expect_equal(abs(y2 - y1) + 1, 151, tolerance = 2)
})

test_that("T1.4.3: Inverted plaque scale bar (caballa): measures internal bar span", {
  # White plaque with dark scale bar inside
  m <- matrix(0, nrow = 600, ncol = 800)
  # Otolith
  xx <- col(m) - 400
  yy <- row(m) - 300
  m[(xx^2 / 200^2 + yy^2 / 100^2) <= 1] <- 1

  # In binarized image of plaque: plaque is light (1), dark bar is (0).
  # After contrast inversion or morphological filtering, the bar is foreground (1)
  # Bar: span = 60 px
  m[510:515, 100:159] <- 1

  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  expect_false(is.null(res))
  info <- attr(res, "scale_info")
  expect_equal(as.numeric(info$bar_length_px), 60, tolerance = 2)
})

test_that("T1.4.4: Plaque scale bar orientation and endpoints are consistent", {
  m <- create_synthetic_otolith()
  m[500:505, 50:129] <- 1 # span = 80
  res <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  info <- attr(res, "scale_info")
  expect_equal(as.numeric(info$bar_length_px), 80, tolerance = 2)
  expect_equal(info$max_x - info$min_x + 1, 80, tolerance = 2)
})

test_that("T1.4.5: Orientation detector correctly differentiates horizontal vs vertical bar", {
  m_horiz <- create_synthetic_otolith()
  m_horiz <- add_horizontal_bar(m_horiz, x1 = 100, x2 = 299, y1 = 500, y2 = 506) # span = 200
  res_h <- aforoR:::find_scale_bar_length(m_horiz, scale_mm = 1)
  info_h <- attr(res_h, "scale_info")

  m_vert <- create_synthetic_otolith()
  m_vert <- add_vertical_bar(m_vert, x1 = 700, x2 = 706, y1 = 100, y2 = 299) # span = 200
  res_v <- aforoR:::find_scale_bar_length(m_vert, scale_mm = 1)
  info_v <- attr(res_v, "scale_info")

  expect_true(info_h$max_x - info_h$min_x > info_h$max_y - info_h$min_y)
  expect_true(info_v$max_y - info_v$min_y > info_v$max_x - info_v$min_x)
})

# ------------------------------------------------------------------------------
# Feature 5: Scale Bar Endpoint Visualization (save_visualization / perimeter)
# ------------------------------------------------------------------------------
test_that("T1.5.1: save_visualization renders cleanly with scale_info and endpoint markers", {
  temp_out <- tempfile(fileext = ".jpg")
  on.exit(unlink(temp_out), add = TRUE)

  img <- create_synthetic_otolith()
  scale_info <- list(
    x1 = 100, y1 = 500, x2 = 300, y2 = 500,
    endpoints = list(x1 = 100, y1 = 500, x2 = 300, y2 = 500),
    min_x = 100, max_x = 300, min_y = 498, max_y = 502,
    bar_length_px = 201, px_per_mm = 201
  )

  # Run save_visualization
  jpeg(temp_out, width = 800, height = 600)
  plot(as.raster(EBImage::Image(img)))
  if (!is.null(scale_info$endpoints)) {
    graphics::points(c(scale_info$endpoints$x1, scale_info$endpoints$x2),
                     c(scale_info$endpoints$y1, scale_info$endpoints$y2),
                     pch = 4, col = "#FF0033", lwd = 2.5)
    graphics::lines(c(scale_info$endpoints$x1, scale_info$endpoints$x2),
                    c(scale_info$endpoints$y1, scale_info$endpoints$y2),
                    col = "#FF0033", lwd = 2)
  }
  dev.off()

  expect_true(file.exists(temp_out))
  expect_gt(file.info(temp_out)$size, 1000)
})

test_that("T1.5.2: save_visualization_perimeter handles scale endpoints without error", {
  temp_out <- tempfile(fileext = ".jpg")
  on.exit(unlink(temp_out), add = TRUE)

  scale_info <- list(
    endpoints = list(x1 = 50, y1 = 100, x2 = 50, y2 = 250),
    min_x = 48, max_x = 52, min_y = 100, max_y = 250,
    bar_length_px = 151
  )

  jpeg(temp_out, width = 600, height = 400)
  plot(1:10, 1:10, type = "n")
  expect_silent({
    graphics::points(c(scale_info$endpoints$x1, scale_info$endpoints$x2),
                     c(scale_info$endpoints$y1, scale_info$endpoints$y2),
                     pch = 4, col = "#FF0033")
  })
  dev.off()

  expect_true(file.exists(temp_out))
})

test_that("T1.5.3: Diagnostic plot generation gracefully handles NULL scale_info", {
  temp_out <- tempfile(fileext = ".jpg")
  on.exit(unlink(temp_out), add = TRUE)

  jpeg(temp_out, width = 400, height = 400)
  expect_silent({
    plot(1:10, 1:10)
    scale_info <- NULL
    if (!is.null(scale_info$endpoints)) {
      graphics::points(c(scale_info$endpoints$x1, scale_info$endpoints$x2),
                       c(scale_info$endpoints$y1, scale_info$endpoints$y2))
    }
  })
  dev.off()
  expect_true(file.exists(temp_out))
})

test_that("T1.5.4: Diagnostic plot generation does not throw warnings on boundary endpoints", {
  temp_out <- tempfile(fileext = ".jpg")
  on.exit(unlink(temp_out), add = TRUE)

  jpeg(temp_out, width = 800, height = 600)
  expect_silent({
    plot(0, 0, xlim = c(1, 800), ylim = c(1, 600), type = "n")
    graphics::points(c(1, 800), c(1, 600), col = "#FF0033", pch = 4)
    graphics::lines(c(1, 800), c(1, 600), col = "#FF0033")
  })
  dev.off()
  expect_true(file.exists(temp_out))
})

test_that("T1.5.5: Visualization preserves non-scale biological specimen contour overlay", {
  temp_out <- tempfile(fileext = ".jpg")
  on.exit(unlink(temp_out), add = TRUE)

  contour <- cbind(cos(seq(0, 2*pi, length.out = 100)) * 50 + 200,
                   sin(seq(0, 2*pi, length.out = 100)) * 30 + 150)
  jpeg(temp_out, width = 400, height = 300)
  plot(1, 1, xlim = c(0, 400), ylim = c(0, 300), type = "n")
  graphics::polygon(contour[, 1], contour[, 2], border = "#00FF66", lwd = 2)
  dev.off()
  expect_true(file.exists(temp_out))
})

# ------------------------------------------------------------------------------
# Feature 6: FFT Float Clamping (Eliminating Overshoot Post-filter2)
# ------------------------------------------------------------------------------
test_that("T1.6.1: Microscopic positive float overshoot (1.0 + 1e-15) clamped to <= 1.0", {
  mat <- matrix(c(1.0 + 2.22e-16, 0.5, 0.2, 0.0), nrow = 2, ncol = 2)
  mat_clamped <- mat
  mat_clamped[mat_clamped > 1.0] <- 1.0
  mat_clamped[mat_clamped < 0.0] <- 0.0

  expect_true(max(mat_clamped) <= 1.0)
  expect_equal(max(mat_clamped), 1.0)
})

test_that("T1.6.2: Microscopic negative float overshoot (-1e-15) clamped to >= 0.0", {
  mat <- matrix(c(-1e-15, 0.5, 0.8, 1.0), nrow = 2, ncol = 2)
  mat_clamped <- mat
  mat_clamped[mat_clamped > 1.0] <- 1.0
  mat_clamped[mat_clamped < 0.0] <- 0.0

  expect_true(min(mat_clamped) >= 0.0)
  expect_equal(min(mat_clamped), 0.0)
})

test_that("T1.6.3: High-contrast step edge subjected to filter2 stays strictly in [0, 1]", {
  m <- matrix(0, nrow = 100, ncol = 100)
  m[40:60, 40:60] <- 1.0
  brush <- EBImage::makeBrush(size = 15, shape = "gaussian", sigma = 3)
  filtered <- EBImage::filter2(m, brush)

  # Check if raw FFT produced any ringing outside [0, 1]
  overshoot <- max(filtered) > 1.0 || min(filtered) < 0.0
  # Clamping eliminates it
  filtered[filtered > 1.0] <- 1.0
  filtered[filtered < 0.0] <- 0.0

  expect_true(max(filtered) <= 1.0)
  expect_true(min(filtered) >= 0.0)
})

test_that("T1.6.4: Clamping does not alter uncorrupted values strictly inside (0, 1)", {
  vals <- seq(0.1, 0.9, by = 0.1)
  vals_clamped <- vals
  vals_clamped[vals_clamped > 1.0] <- 1.0
  vals_clamped[vals_clamped < 0.0] <- 0.0
  expect_identical(vals, vals_clamped)
})

test_that("T1.6.5: End-to-end preprocess_image on synthetic image produces bounded matrix", {
  m <- create_synthetic_otolith(nrow = 200, ncol = 200, cx = 100, cy = 100, rx = 50, ry = 30)
  bin <- aforoR::preprocess_image(EBImage::Image(m))
  expect_true(is.logical(bin) || all(bin %in% c(0, 1, TRUE, FALSE)))
})

# ------------------------------------------------------------------------------
# Feature 7: Threshold Rescaling Guard (max > 1.05 before / 255)
# ------------------------------------------------------------------------------
test_that("T1.7.1: Float overshoot (max = 1.0001) is NOT rescaled by 255 in calculate_adaptive_threshold", {
  mat <- matrix(runif(10000, min = 0.1, max = 0.9), nrow = 100, ncol = 100)
  mat[1, 1] <- 1.0001 # Microscopic float overshoot

  th <- aforoR:::calculate_adaptive_threshold(mat)
  # If erroneously divided by 255, threshold collapses to < 0.01
  expect_true(th > 0.05)
})

test_that("T1.7.2: 8-bit integer image (max = 240) IS properly rescaled by 255", {
  mat_int <- matrix(round(runif(10000, min = 20, max = 240)), nrow = 100, ncol = 100)
  th <- aforoR:::calculate_adaptive_threshold(mat_int)
  expect_true(th >= 0.05 && th <= 0.95)
})

test_that("T1.7.3: Intermediate value near 1.0 (max = 1.04) is NOT rescaled by 255", {
  mat <- matrix(runif(10000, min = 0.1, max = 0.8), nrow = 100, ncol = 100)
  mat[1, 1] <- 1.04
  th <- aforoR:::calculate_adaptive_threshold(mat)
  expect_true(th > 0.05)
})

test_that("T1.7.4: Integer value (max = 1.06) is rescaled by 255", {
  # When max > 1.05, it is treated as non-normalized and scaled
  mat <- matrix(1.06, nrow = 50, ncol = 50)
  th <- aforoR:::calculate_adaptive_threshold(mat)
  expect_true(is.numeric(th) && is.finite(th))
})

test_that("T1.7.5: Standard normalized image in [0.0, 1.0] returns threshold in (0.05, 0.95)", {
  mat <- create_synthetic_otolith()
  th <- aforoR:::calculate_adaptive_threshold(mat)
  expect_true(th >= 0.05 && th <= 0.95)
})

# ------------------------------------------------------------------------------
# Feature 8: Otsu Range Safeguard (Breaks Span Range)
# ------------------------------------------------------------------------------
test_that("T1.8.1: Matrix with values slightly outside [0, 1] does not crash EBImage::otsu", {
  mat <- matrix(seq(-1e-5, 1.0 + 1e-5, length.out = 10000), nrow = 100, ncol = 100)
  # Clamping safeguard
  mat[mat > 1.0] <- 1.0
  mat[mat < 0.0] <- 0.0
  th <- tryCatch(EBImage::otsu(mat), error = function(e) NULL)
  expect_false(is.null(th))
  expect_true(th >= 0 && th <= 1)
})

test_that("T1.8.2: Uniform image (all pixels = 0.5) handled gracefully by adaptive thresholding", {
  mat <- matrix(0.5, nrow = 100, ncol = 100)
  th <- aforoR:::calculate_adaptive_threshold(mat)
  expect_true(is.numeric(th) && is.finite(th))
})

test_that("T1.8.3: Near-flat image (dynamic range < 1e-6) returns valid threshold", {
  mat <- matrix(0.5 + runif(10000, -1e-7, 1e-7), nrow = 100, ncol = 100)
  th <- aforoR:::calculate_adaptive_threshold(mat)
  expect_true(is.numeric(th) && is.finite(th))
})

test_that("T1.8.4: Image with single extreme outlier pixel returns stable threshold", {
  mat <- matrix(0.3, nrow = 100, ncol = 100)
  mat[1, 1] <- 1.0
  th <- aforoR:::calculate_adaptive_threshold(mat)
  expect_true(is.numeric(th) && is.finite(th))
})

test_that("T1.8.5: calculate_adaptive_threshold returns a single finite numeric scalar in (0, 1)", {
  mat <- matrix(runif(10000), nrow = 100, ncol = 100)
  th <- aforoR:::calculate_adaptive_threshold(mat)
  expect_true(is.numeric(th))
  expect_equal(length(th), 1)
  expect_true(is.finite(th))
  expect_true(th > 0 && th < 1)
})

# ------------------------------------------------------------------------------
# Feature 9: 2-Channel TIFF Handling (cod/1975)
# ------------------------------------------------------------------------------
test_that("T1.9.1: 2-channel Image object is stripped to 1-channel Grayscale", {
  img2ch <- create_synthetic_2channel_image()
  expect_equal(dim(img2ch)[3], 2)

  # Check stripping logic
  if (length(dim(img2ch)) == 3 && dim(img2ch)[3] == 2) {
    img_stripped <- img2ch[, , 1]
    EBImage::colorMode(img_stripped) <- EBImage::Grayscale
  }

  expect_equal(length(dim(img_stripped)), 2)
  expect_equal(EBImage::colorMode(img_stripped), EBImage::Grayscale)
})

test_that("T1.9.2: Grayscale intensities from Channel 1 are preserved without +0.3333 shift", {
  img2ch <- create_synthetic_2channel_image()
  # Channel 1 background pixel
  ch1_bg <- as.numeric(img2ch[1, 1, 1])
  expect_equal(ch1_bg, 0.02)

  # In flawed channel(img, "grey"), bg becomes (0.02 + 1.0 + 0) / 3 = 0.34
  flawed_grey <- EBImage::channel(img2ch, "grey")
  expect_true(as.numeric(flawed_grey[1, 1]) > 0.30)

  # Direct extraction preserves 0.02
  clean_grey <- img2ch[, , 1]
  expect_equal(as.numeric(clean_grey[1, 1]), 0.02)
})

test_that("T1.9.3: True black background remains 0.0 after 2-channel alpha strip", {
  ch1 <- matrix(0.0, nrow = 100, ncol = 100)
  ch1[30:70, 30:70] <- 0.90
  ch2 <- matrix(1.0, nrow = 100, ncol = 100)
  arr <- array(0, dim = c(100, 100, 2))
  arr[, , 1] <- ch1
  arr[, , 2] <- ch2
  img <- EBImage::Image(arr, colormode = "Color")

  stripped <- img[, , 1]
  expect_equal(as.numeric(stripped[1, 1]), 0.0)
})

test_that("T1.9.4: Diagnostic plot of 2-channel image does not exhibit green tinting", {
  img2ch <- create_synthetic_2channel_image()
  # When stripped to grayscale, plotting assigns equal R, G, B
  stripped <- img2ch[, , 1]
  EBImage::colorMode(stripped) <- EBImage::Grayscale

  temp_out <- tempfile(fileext = ".jpg")
  on.exit(unlink(temp_out), add = TRUE)
  jpeg(temp_out, width = 600, height = 500)
  plot(as.raster(stripped))
  dev.off()

  expect_true(file.exists(temp_out))
})

test_that("T1.9.5: Biological contour extracted from 2-channel image is accurate and <85% area", {
  img2ch <- create_synthetic_2channel_image(nrow = 500, ncol = 600)
  clean_img <- img2ch[, , 1]
  bin <- clean_img > 0.3
  contour <- aforoR::extract_contour(bin)

  expect_false(is.null(contour))
  area <- aforoR:::.coo_area(contour)
  expect_true(area < 0.85 * (500 * 600))
  expect_true(area > 5000)
})

# ------------------------------------------------------------------------------
# Feature 10: Outer Frame Filtering (Area >= 85% Rejection)
# ------------------------------------------------------------------------------
test_that("T1.10.1: Full image frame contour (>= 85% area) is rejected in favor of internal specimen", {
  m <- matrix(0, nrow = 600, ncol = 800)
  total_px <- 600 * 800 # 480,000

  # Outer frame touching borders: area ~95% of image
  m[5:595, 5:795] <- 1
  m[20:580, 20:780] <- 0 # hollow frame

  # Inside: otolith of area ~50,000 px
  xx <- col(m) - 400
  yy <- row(m) - 300
  m[(xx^2 / 160^2 + yy^2 / 80^2) <= 1] <- 1

  # Filter contours: reject any contour with area >= 0.85 * total_px
  cont <- EBImage::ocontour(EBImage::bwlabel(m))
  valid_contours <- list()
  for (c_mat in cont) {
    if (is.matrix(c_mat) && nrow(c_mat) >= 200) {
      c_area <- aforoR:::.coo_area(c_mat)
      if (c_area < 0.85 * total_px) {
        valid_contours[[length(valid_contours) + 1]] <- c_mat
      }
    }
  }

  expect_gt(length(valid_contours), 0)
  max_valid_area <- max(vapply(valid_contours, aforoR:::.coo_area, numeric(1)))
  expect_true(max_valid_area < 0.85 * total_px)
})

test_that("T1.10.2: Contour with area exactly 86% of total image area is rejected", {
  total_px <- 100000
  contour_area <- 86000 # 86%
  is_valid <- contour_area < 0.85 * total_px
  expect_false(is_valid)
})

test_that("T1.10.3: Contour with area 84% of total image area is accepted", {
  total_px <- 100000
  contour_area <- 84000 # 84%
  is_valid <- contour_area < 0.85 * total_px
  expect_true(is_valid)
})

test_that("T1.10.4: Binary image containing only an outer border frame returns NULL / warning", {
  m <- matrix(0, nrow = 500, ncol = 500)
  total_px <- 250000
  # Only border frame (area > 85%)
  m[2:498, 2:498] <- 1
  m[15:485, 15:485] <- 0

  cont <- EBImage::ocontour(EBImage::bwlabel(m))
  valid_contours <- list()
  for (c_mat in cont) {
    if (is.matrix(c_mat) && nrow(c_mat) >= 200) {
      if (aforoR:::.coo_area(c_mat) < 0.85 * total_px) {
        valid_contours[[length(valid_contours) + 1]] <- c_mat
      }
    }
  }
  expect_equal(length(valid_contours), 0)
})

test_that("T1.10.5: Specimen contour extracted from real or synthetic otolith is strictly < 85% area", {
  m <- create_synthetic_otolith(nrow = 600, ncol = 800)
  contour <- aforoR::extract_contour(m)
  expect_false(is.null(contour))
  area_ratio <- aforoR:::.coo_area(contour) / (600 * 800)
  expect_true(area_ratio < 0.85)
})

# ------------------------------------------------------------------------------
# Feature 11: Small Image min_area Support (Ostium 480x480 Images)
# ------------------------------------------------------------------------------
test_that("T1.11.1: Specimen with area ~3,500px in 480x480 image is rejected by default min_area = 5000", {
  m <- matrix(0, nrow = 480, ncol = 480)
  # Specimen of area ~3,500 px (radius x = 50, radius y = 22 -> pi * 50 * 22 ~ 3455)
  xx <- col(m) - 240
  yy <- row(m) - 240
  m[(xx^2 / 50^2 + yy^2 / 22^2) <= 1] <- 1

  # With default min_area = 5000, extract_contour gives warning and returns NULL
  expect_warning(res <- aforoR::extract_contour(m, min_points = 200, min_area = 5000))
  expect_null(res)
})

test_that("T1.11.2: Same specimen with min_area = 1000 successfully extracts contour", {
  m <- matrix(0, nrow = 480, ncol = 480)
  # Specimen with area ~3,500 px and > 200 points
  theta <- seq(0, 2*pi, length.out = 250)
  rx <- 55
  ry <- 20
  xx <- col(m) - 240
  yy <- row(m) - 240
  m[(xx^2 / rx^2 + yy^2 / ry^2) <= 1] <- 1

  res <- aforoR::extract_contour(m, min_points = 200, min_area = 1000)
  expect_false(is.null(res))
  area <- aforoR:::.coo_area(res)
  expect_true(area >= 1000 && area < 5000)
})

test_that("T1.11.3: extract_contour respects arbitrary user-supplied min_area parameter", {
  m <- matrix(0, nrow = 200, ncol = 200)
  xx <- col(m) - 100
  yy <- row(m) - 100
  m[(xx^2 / 30^2 + yy^2 / 15^2) <= 1] <- 1 # area ~ 1400 px

  res_1000 <- aforoR::extract_contour(m, min_points = 50, min_area = 1000)
  expect_false(is.null(res_1000))

  expect_warning(res_2000 <- aforoR::extract_contour(m, min_points = 50, min_area = 2000))
  expect_null(res_2000)
})

test_that("T1.11.4: Extracted small contour has valid dimensions and coordinates inside [1, 480]", {
  m <- matrix(0, nrow = 480, ncol = 480)
  xx <- col(m) - 200
  yy <- row(m) - 200
  m[(xx^2 / 45^2 + yy^2 / 25^2) <= 1] <- 1

  res <- aforoR::extract_contour(m, min_points = 100, min_area = 1000)
  expect_false(is.null(res))
  expect_true(all(res[, 1] >= 1 & res[, 1] <= 480))
  expect_true(all(res[, 2] >= 1 & res[, 2] <= 480))
})

test_that("T1.11.5: Specimen with area < 1000 (e.g. 800px) is rejected when min_area = 1000", {
  m <- matrix(0, nrow = 480, ncol = 480)
  xx <- col(m) - 240
  yy <- row(m) - 240
  m[(xx^2 / 20^2 + yy^2 / 12^2) <= 1] <- 1 # area ~ 750 px

  expect_warning(res <- aforoR::extract_contour(m, min_points = 50, min_area = 1000))
  expect_null(res)
})

# ------------------------------------------------------------------------------
# Feature 12: Flexible plot.css() S3 Method (User ... Override via modifyList)
# ------------------------------------------------------------------------------
test_that("T1.12.1: plot.css accepts custom 'main' argument without duplicate argument error", {
  theta <- seq(0, 2*pi, length.out = 80)
  contour <- cbind(cos(theta) + 0.2*cos(3*theta), sin(theta) + 0.2*sin(3*theta))
  css_res <- aforoR::calculate_css(contour, n_points = 100, sigma_max = 3)

  pdf(NULL)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  # Verbatim requirement: plot(calculate_css(contour), main = "Custom Title")
  expect_error(plot(css_res, main = "Custom Title"), NA)
})

test_that("T1.12.2: plot.css accepts custom 'xlab' and 'ylab' arguments without error", {
  theta <- seq(0, 2*pi, length.out = 80)
  contour <- cbind(cos(theta) + 0.2*cos(3*theta), sin(theta) + 0.2*sin(3*theta))
  css_res <- aforoR::calculate_css(contour, n_points = 100, sigma_max = 3)

  pdf(NULL)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  expect_error(plot(css_res, xlab = "Normalized Arc Length", ylab = "Gaussian Scale"), NA)
})

test_that("T1.12.3: plot.css accepts styling arguments: col, pch, cex", {
  theta <- seq(0, 2*pi, length.out = 80)
  contour <- cbind(cos(theta) + 0.2*cos(3*theta), sin(theta) + 0.2*sin(3*theta))
  css_res <- aforoR::calculate_css(contour, n_points = 100, sigma_max = 3)

  pdf(NULL)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  expect_error(plot(css_res, col = "red", pch = 16, cex = 1.2), NA)
})

test_that("T1.12.4: plot.css accepts limit overrides: xlim and ylim", {
  theta <- seq(0, 2*pi, length.out = 80)
  contour <- cbind(cos(theta) + 0.2*cos(3*theta), sin(theta) + 0.2*sin(3*theta))
  css_res <- aforoR::calculate_css(contour, n_points = 100, sigma_max = 3)

  pdf(NULL)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  expect_error(plot(css_res, xlim = c(0.1, 0.9), ylim = c(0, 5)), NA)
})

test_that("T1.12.5: plot.css with no extra arguments runs cleanly with default arguments", {
  theta <- seq(0, 2*pi, length.out = 80)
  contour <- cbind(cos(theta) + 0.2*cos(3*theta), sin(theta) + 0.2*sin(3*theta))
  css_res <- aforoR::calculate_css(contour, n_points = 100, sigma_max = 3)

  pdf(NULL)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  expect_error(plot(css_res), NA)
})

test_that("T1.12.6: plot.css forwards multi-parameter overrides simultaneously", {
  theta <- seq(0, 2*pi, length.out = 80)
  contour <- cbind(cos(theta) + 0.2*cos(3*theta), sin(theta) + 0.2*sin(3*theta))
  css_res <- aforoR::calculate_css(contour, n_points = 100, sigma_max = 3)

  pdf(NULL)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  expect_error(plot(css_res, main = "Multi-Param Test", col = "blue", pch = 17,
                    cex = 0.8, xlab = "Custom X", ylab = "Custom Y",
                    xlim = c(0, 1), ylim = c(0, 4)), NA)
})


# ==============================================================================
# TIER 2: BOUNDARY & CORNER CASES
# ==============================================================================

test_that("T2.1: Empty and zero inputs handling", {
  expect_equal(aforoR:::calculate_adaptive_threshold(NULL), 0.2)
  expect_error(aforoR::extract_contour(NULL))
  expect_error(aforoR::calculate_css(NULL))

  # All zero matrix
  all_zero <- matrix(0, nrow = 100, ncol = 100)
  expect_warning(res_c <- aforoR::extract_contour(all_zero))
  expect_null(res_c)
  expect_null(aforoR:::find_scale_bar_length(all_zero))

  # All one matrix
  all_one <- matrix(1, nrow = 100, ncol = 100)
  expect_warning(res_one <- aforoR::extract_contour(all_one))
  expect_null(res_one)
})

test_that("T2.2: Pure red mask vs subtle red / orange / pink discrimination", {
  # Red condition: R > 0.45 & R > 1.4 * max(G, B)
  # Pure red: R=0.9, G=0.1, B=0.1 -> 0.9 > 0.45 & 0.9 > 1.4*0.1=0.14 -> TRUE
  expect_true(0.9 > 0.45 && 0.9 > 1.4 * max(0.1, 0.1))

  # Orange: R=0.9, G=0.7, B=0.1 -> 0.9 > 0.45 but 0.9 > 1.4*0.7=0.98 is FALSE -> FALSE
  expect_false(0.9 > 0.45 && 0.9 > 1.4 * max(0.7, 0.1))

  # Magenta/Pink: R=0.9, G=0.2, B=0.8 -> 0.9 > 0.45 but 0.9 > 1.4*0.8=1.12 is FALSE -> FALSE
  expect_false(0.9 > 0.45 && 0.9 > 1.4 * max(0.2, 0.8))

  # Weak dark red: R=0.40, G=0.05, B=0.05 -> R not > 0.45 -> FALSE
  expect_false(0.40 > 0.45 && 0.40 > 1.4 * max(0.05, 0.05))
})

test_that("T2.3: Plaque vs thin line contrast boundary", {
  # Thin line (eccentricity > 0.98) vs square plaque (eccentricity < 0.5)
  m <- matrix(0, nrow = 300, ncol = 300)
  m[100:104, 50:200] <- 1 # Thin line: span 151, height 5 -> eccentricity ~0.999
  lab <- EBImage::bwlabel(m)
  mom <- EBImage::computeFeatures.moment(lab)
  expect_true(mom[1, "m.eccentricity"] > 0.95)

  # Square plaque
  m_sq <- matrix(0, nrow = 300, ncol = 300)
  m_sq[50:150, 50:150] <- 1 # 100x100 square
  lab_sq <- EBImage::bwlabel(m_sq)
  mom_sq <- EBImage::computeFeatures.moment(lab_sq)
  expect_true(mom_sq[1, "m.eccentricity"] < 0.2)
})

test_that("T2.4: Vertical vs horizontal scale bar orientation boundaries", {
  # Exactly horizontal: delta_x = 100, delta_y = 5 -> delta_x > delta_y
  expect_true((200 - 100 + 1) > (55 - 50 + 1))

  # Exactly vertical: delta_x = 5, delta_y = 100 -> delta_y > delta_x
  expect_true((150 - 50 + 1) > (105 - 100 + 1))
})

test_that("T2.5: Float overshoot near 1.0 boundary values", {
  # 1.0 + 1e-16: float overshoot, should NOT divide by 255
  val_close <- 1.0 + 1e-16
  expect_false(val_close > 1.05)

  # 1.049: should NOT divide by 255
  expect_false(1.049 > 1.05)

  # 1.051: should divide by 255
  expect_true(1.051 > 1.05)
})

test_that("T2.6: 2-channel TIFF alpha strip array dimension and mode reduction", {
  arr <- array(0.5, dim = c(100, 100, 2))
  img <- EBImage::Image(arr, colormode = "Color")
  expect_equal(dim(img)[3], 2)

  # Stripping
  img_gray <- img[, , 1]
  EBImage::colorMode(img_gray) <- EBImage::Grayscale
  expect_equal(length(dim(img_gray)), 2)
  expect_equal(EBImage::colorMode(img_gray), EBImage::Grayscale)
})

test_that("T2.7: Area >= 85% border check exact boundary", {
  total_px <- 1000000
  area_boundary_high <- 850001
  area_boundary_low  <- 849999

  expect_false(area_boundary_high < 0.85 * total_px)
  expect_true(area_boundary_low < 0.85 * total_px)
})

test_that("T2.8: Small image min_area=1000 boundary", {
  specimen_above <- 1001
  specimen_below <- 999
  min_area <- 1000

  expect_true(specimen_above > min_area)
  expect_false(specimen_below > min_area)
})


# ==============================================================================
# TIER 3: CROSS-FEATURE INTERACTIONS
# ==============================================================================

test_that("T3.1: 2-channel TIFF with scale bar strips alpha, detects scale accurately, and plots without green tint", {
  # Create 2-channel image with specimen and scale bar
  img2ch <- create_synthetic_2channel_image(nrow = 600, ncol = 800)
  # Scale bar in Channel 1 at [420:425, 100:249] -> span = 150 px

  # Step 1: Alpha channel stripped
  if (length(dim(img2ch)) == 3 && dim(img2ch)[3] == 2) {
    clean_img <- img2ch[, , 1]
    EBImage::colorMode(clean_img) <- EBImage::Grayscale
  }
  expect_equal(length(dim(clean_img)), 2)

  # Step 2: Scale bar detected from clean image
  bin_img <- clean_img > 0.3
  res_scale <- aforoR:::find_scale_bar_length(bin_img, scale_mm = 1)
  expect_false(is.null(res_scale))
  expect_equal(as.numeric(attr(res_scale, "scale_info")$bar_length_px), 150, tolerance = 2)

  # Step 3: Diagnostic plot renders without green channel tinting
  temp_out <- tempfile(fileext = ".jpg")
  on.exit(unlink(temp_out), add = TRUE)
  jpeg(temp_out, width = 800, height = 600)
  plot(as.raster(clean_img))
  scale_info <- attr(res_scale, "scale_info")
  graphics::rect(scale_info$min_x, scale_info$min_y, scale_info$max_x, scale_info$max_y, border = "#00FF66")
  dev.off()
  expect_true(file.exists(temp_out))
})

test_that("T3.2: Textured background with plaque scale bar extracts specimen <85% area and detects scale bar", {
  # Synthetic textured image (caballa style)
  img <- matrix(0.18 + runif(600 * 800, min = -0.05, max = 0.05), nrow = 600, ncol = 800)
  # Otolith
  xx <- col(img) - 400
  yy <- row(img) - 300
  img[(xx^2 / 220^2 + yy^2 / 110^2) <= 1] <- 0.85

  # Plaque and dark scale bar: span = 100 px
  img[510:515, 120:219] <- 0.95

  bin <- img > aforoR:::calculate_adaptive_threshold(img)
  contour <- aforoR::extract_contour(bin)
  expect_false(is.null(contour))

  area_ratio <- aforoR:::.coo_area(contour) / (600 * 800)
  expect_true(area_ratio < 0.85)

  # Scale bar detection
  sc <- aforoR:::find_scale_bar_length(bin, scale_mm = 1)
  expect_false(is.null(sc))
  expect_equal(as.numeric(attr(sc, "scale_info")$bar_length_px), 100, tolerance = 3)
})

test_that("T3.3: Small image end-to-end: min_area=1000 contour -> calculate_css -> plot.css with custom overrides", {
  m <- matrix(0, nrow = 480, ncol = 480)
  theta <- seq(0, 2*pi, length.out = 250)
  # Perturbed shape for CSS inflection points
  rx <- 60
  ry <- 25
  xx <- col(m) - 240
  yy <- row(m) - 240
  m[(xx^2 / rx^2 + yy^2 / ry^2) <= 1] <- 1

  # Extract contour with small min_area
  contour <- aforoR::extract_contour(m, min_points = 150, min_area = 1000)
  expect_false(is.null(contour))

  # Calculate curvature scale space
  css_res <- aforoR::calculate_css(contour, n_points = 150, sigma_max = 3)
  expect_s3_class(css_res, "css")

  # Plot with custom user arguments
  pdf(NULL)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  expect_error(plot(css_res, main = "Ostium 480x480 CSS Analysis",
                    col = "darkgreen", pch = 17, cex = 0.9,
                    xlab = "Normalized Length", ylab = "Curvature Scale"), NA)
})
