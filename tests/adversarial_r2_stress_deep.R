#!/usr/bin/env Rscript
# ==============================================================================
# tests/adversarial_r2_stress_deep.R
#
# Deep Adversarial Stress Testing: Extreme Edge Cases & Hostile Inputs for R2
# - Constant/uniform images (all 0, all 1, all 0.5)
# - Degenerate geometries (extreme aspect ratios: 10x1000, 1000x10)
# - Missing/NaN/Inf values in image arrays
# - Donut/hollow shapes and multi-component contours
# - Boundary touching contours (x=1, y=1, x=W, y=H)
# - Invalid argument attacks against threshold_factor, blur_size, blur_sigma
# ==============================================================================

suppressPackageStartupMessages({
  library(EBImage)
  library(aforoR)
})

cat("================================================================================\n")
cat("      DEEP ADVERSARIAL STRESS TESTING: EXTREME EDGE CASES & HOSTILE INPUTS     \n")
cat("================================================================================\n\n")

calc_adaptive_th <- aforoR:::calculate_adaptive_threshold

total_cases <- 0
passed_cases <- 0
failed_cases <- 0

check_case <- function(name, passed, detail = "") {
  total_cases <<- total_cases + 1
  if (passed) {
    passed_cases <<- passed_cases + 1
    cat(sprintf("  [PASS] %s\n", name))
  } else {
    failed_cases <<- failed_cases + 1
    cat(sprintf("  [FAIL] %s - Detail: %s\n", name, detail))
  }
}

# ------------------------------------------------------------------------------
# 1. Constant & Zero-Variance Inputs
# ------------------------------------------------------------------------------
cat("\n--- 1. Constant & Zero-Variance Inputs ---\n")

# Pitch black (all 0)
m_black <- matrix(0, nrow = 100, ncol = 100)
th_black <- calc_adaptive_th(m_black)
check_case(
  "1.1 All-zero pitch black matrix does not crash and returns fallback threshold",
  is.numeric(th_black) && length(th_black) == 1 && !is.na(th_black) && th_black > 0,
  sprintf("th = %s", as.character(th_black))
)

# Pure white (all 1)
m_white <- matrix(1, nrow = 100, ncol = 100)
th_white <- calc_adaptive_th(m_white)
check_case(
  "1.2 All-one pure white matrix does not crash and returns fallback threshold",
  is.numeric(th_white) && length(th_white) == 1 && !is.na(th_white) && th_white > 0,
  sprintf("th = %s", as.character(th_white))
)

# Mid-gray (all 0.5)
m_gray <- matrix(0.5, nrow = 100, ncol = 100)
th_gray <- calc_adaptive_th(m_gray)
check_case(
  "1.3 All 0.5 mid-gray matrix does not crash and returns fallback threshold",
  is.numeric(th_gray) && length(th_gray) == 1 && !is.na(th_gray) && th_gray > 0,
  sprintf("th = %s", as.character(th_gray))
)

# ------------------------------------------------------------------------------
# 2. Extreme Aspect Ratios
# ------------------------------------------------------------------------------
cat("\n--- 2. Extreme Aspect Ratios ---\n")

# Tall thin image (1000 x 20)
m_tall <- matrix(0.05, nrow = 1000, ncol = 20)
m_tall[200:800, 5:15] <- 0.9
bin_tall <- tryCatch(preprocess_image(EBImage::Image(m_tall)), error = function(e) NULL)
check_case(
  "2.1 Extreme tall aspect ratio (1000x20) preprocessed cleanly",
  !is.null(bin_tall) && inherits(bin_tall, "Image"),
  "Failed on tall aspect ratio"
)

# Wide thin image (20 x 1000)
m_wide <- matrix(0.05, nrow = 20, ncol = 1000)
m_wide[5:15, 200:800] <- 0.9
bin_wide <- tryCatch(preprocess_image(EBImage::Image(m_wide)), error = function(e) NULL)
check_case(
  "2.2 Extreme wide aspect ratio (20x1000) preprocessed cleanly",
  !is.null(bin_wide) && inherits(bin_wide, "Image"),
  "Failed on wide aspect ratio"
)

# ------------------------------------------------------------------------------
# 3. Missing & Corrupted Float Inputs (NA, NaN, Inf)
# ------------------------------------------------------------------------------
cat("\n--- 3. Missing & Corrupted Values ---\n")

m_corrupt <- matrix(runif(10000, 0.1, 0.8), nrow = 100, ncol = 100)
m_corrupt[10:15, 10:15] <- NA_real_
m_corrupt[20:25, 20:25] <- NaN
m_corrupt[1, 1] <- Inf
m_corrupt[1, 2] <- -Inf

th_corrupt <- tryCatch(calc_adaptive_th(m_corrupt), error = function(e) NULL)
check_case(
  "3.1 Matrix with NA, NaN, Inf handled gracefully without crash",
  !is.null(th_corrupt) && is.numeric(th_corrupt) && th_corrupt >= 0 && th_corrupt <= 1,
  sprintf("th = %s", as.character(th_corrupt))
)

# All NA matrix
m_all_na <- matrix(NA_real_, nrow = 100, ncol = 100)
th_all_na <- tryCatch(calc_adaptive_th(m_all_na), error = function(e) NULL)
check_case(
  "3.2 All NA matrix returns default fallback threshold",
  !is.null(th_all_na) && th_all_na == 0.2,
  sprintf("th = %s", as.character(th_all_na))
)

# ------------------------------------------------------------------------------
# 4. Hollow / Donut Geometry & Edge-Touching Contours
# ------------------------------------------------------------------------------
cat("\n--- 4. Hollow / Donut Geometry & Edge-Touching Contours ---\n")

# Donut shape: outer radius 100, inner hole radius 50
m_donut <- matrix(0, nrow = 400, ncol = 400)
xx <- col(m_donut) - 200
yy <- row(m_donut) - 200
r2 <- xx^2 + yy^2
m_donut[r2 <= 100^2 & r2 >= 50^2] <- 1
cnt_donut <- extract_contour(m_donut, min_points = 100, min_area = 5000)
check_case(
  "4.1 Hollow donut otolith contour extracted without error",
  !is.null(cnt_donut) && is.matrix(cnt_donut) && nrow(cnt_donut) >= 100,
  "Failed on donut shape"
)

# Edge-touching specimen (reaches x = 1, y = 1, x = 600, y = 800)
m_touch <- matrix(0, nrow = 600, ncol = 800)
m_touch[1:500, 1:400] <- 1 # covers top-left quadrant touching edge
cnt_touch <- extract_contour(m_touch, min_points = 150, min_area = 5000)
check_case(
  "4.2 Edge-touching contour handled correctly",
  !is.null(cnt_touch) && (aforoR:::.coo_area(cnt_touch) / (600 * 800)) < 0.85,
  "Failed on edge touching contour"
)

# ------------------------------------------------------------------------------
# 5. Parameter Validation & Defensive Error Paths
# ------------------------------------------------------------------------------
cat("\n--- 5. Parameter Validation & Defensive Error Paths ---\n")

ex_img <- EBImage::Image(matrix(runif(10000), 100, 100))

# Invalid threshold_factor
err_tf1 <- inherits(tryCatch(preprocess_image(ex_img, threshold_factor = -1), error = function(e) e), "error")
err_tf2 <- inherits(tryCatch(preprocess_image(ex_img, threshold_factor = 0), error = function(e) e), "error")
err_tf3 <- inherits(tryCatch(preprocess_image(ex_img, threshold_factor = "invalid"), error = function(e) e), "error")
check_case(
  "5.1 Invalid threshold_factor (<= 0 or non-numeric) rejected with clear error",
  err_tf1 && err_tf2 && err_tf3,
  sprintf("tf1: %s, tf2: %s, tf3: %s", err_tf1, err_tf2, err_tf3)
)

# Invalid min_points and min_area in extract_contour
err_mp1 <- inherits(tryCatch(extract_contour(m_donut, min_points = 2), error = function(e) e), "error")
err_mp2 <- inherits(tryCatch(extract_contour(m_donut, min_points = -5), error = function(e) e), "error")
err_ma1 <- inherits(tryCatch(extract_contour(m_donut, min_area = -10), error = function(e) e), "error")
check_case(
  "5.2 Invalid min_points (< 3) or min_area (< 0) rejected with clear error",
  err_mp1 && err_mp2 && err_ma1,
  sprintf("mp1: %s, mp2: %s, ma1: %s", err_mp1, err_mp2, err_ma1)
)

# NULL binary_image in extract_contour
err_null_bin <- inherits(tryCatch(extract_contour(NULL), error = function(e) e), "error")
check_case(
  "5.3 NULL binary_image in extract_contour rejected with error",
  err_null_bin,
  "Did not reject NULL binary_image"
)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
cat("\n================================================================================\n")
cat("                  DEEP ADVERSARIAL STRESS TEST SUMMARY                          \n")
cat("================================================================================\n")
cat(sprintf("Total Cases Evaluated : %d\n", total_cases))
cat(sprintf("Passed                : %d (%.1f%%)\n", passed_cases, (passed_cases / total_cases) * 100))
cat(sprintf("Failed                : %d (%.1f%%)\n", failed_cases, (failed_cases / total_cases) * 100))

if (failed_cases > 0) {
  quit(status = 1, save = "no")
} else {
  cat("\n[SUCCESS] All deep adversarial edge cases passed cleanly!\n")
  quit(status = 0, save = "no")
}
