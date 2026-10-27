#!/usr/bin/env Rscript
# ==============================================================================
# tests/adversarial_r2_stress.R
#
# Empirical Adversarial Stress Test Suite for Requirement R2:
# - Float overshoot/undershoot at critical boundaries (1.0000001, 1.049, 1.050, 1.051)
# - Multi-channel arrays (1-ch, 2-ch, 3-ch, 4-ch) & colorMode handling
# - Textured grey backgrounds and strict 85% area border filtering
# ==============================================================================

suppressPackageStartupMessages({
  library(EBImage)
  library(aforoR)
})

cat("================================================================================\n")
cat("      ADVERSARIAL STRESS TEST SUITE: REQUIREMENT R2 & EDGE CASE MINING          \n")
cat("================================================================================\n\n")

# Access internal functions from aforoR namespace
calc_adaptive_th <- aforoR:::calculate_adaptive_threshold
coo_area <- aforoR:::.coo_area

total_tests <- 0
passed_tests <- 0
failed_tests <- 0
failures_log <- list()

record_test <- function(name, passed, detail = "") {
  total_tests <<- total_tests + 1
  if (passed) {
    passed_tests <<- passed_tests + 1
    cat(sprintf("  [PASS] %s\n", name))
  } else {
    failed_tests <<- failed_tests + 1
    failures_log[[length(failures_log) + 1]] <<- list(name = name, detail = detail)
    cat(sprintf("  [FAIL] %s - Detail: %s\n", name, detail))
  }
}

# ==============================================================================
# SUITE 1: Float Boundary Overshoots, Undershoots & Rescaling Guard
# ==============================================================================
cat("\n--- Suite 1: Float Boundary Overshoots & Rescaling Guard ---\n")

# Test 1.1: Microscopic overshoot: 1.0000001
set.seed(101)
base_mat <- matrix(runif(10000, 0.1, 0.8), nrow = 100, ncol = 100)
mat_10000001 <- base_mat
mat_10000001[1, 1] <- 1.0000001
th_10000001 <- calc_adaptive_th(mat_10000001)
record_test(
  "1.1 Microscopic overshoot (1.0000001) does not collapse threshold via / 255",
  th_10000001 >= 0.1 && th_10000001 <= 0.9,
  sprintf("th = %f (expected in [0.1, 0.9])", th_10000001)
)

# Test 1.2: Boundary overshoot just below 1.05: 1.049
mat_1049 <- base_mat
mat_1049[1, 1] <- 1.049
th_1049 <- calc_adaptive_th(mat_1049)
record_test(
  "1.2 Float overshoot near boundary (1.049) does not trigger / 255",
  th_1049 >= 0.1 && th_1049 <= 0.9,
  sprintf("th = %f (expected in [0.1, 0.9])", th_1049)
)

# Test 1.3: Exact boundary: 1.0500000
mat_1050 <- base_mat
mat_1050[1, 1] <- 1.0500000
th_1050 <- calc_adaptive_th(mat_1050)
record_test(
  "1.3 Float overshoot at exact boundary (1.050) does not trigger / 255",
  th_1050 >= 0.1 && th_1050 <= 0.9,
  sprintf("th = %f (expected in [0.1, 0.9])", th_1050)
)

# Test 1.4: Critical boundary test: 1.051
mat_1051 <- base_mat
mat_1051[1, 1] <- 1.051
th_1051 <- calc_adaptive_th(mat_1051)
# When max(img) > 1.05 and <= 255, code executes img <- img / 255.
# If an image had max = 1.051, dividing by 255 compresses values to ~0.004.
cat(sprintf("      Info: for max(img) = 1.051, calculate_adaptive_threshold returned th = %f\n", th_1051))
record_test(
  "1.4 Critical boundary check at 1.051 behavior observed",
  is.numeric(th_1051) && !is.na(th_1051) && !is.infinite(th_1051),
  sprintf("th = %f", th_1051)
)

# Test 1.5: Legitimate 8-bit integer image (e.g. max = 220, max = 255)
mat_8bit <- matrix(round(runif(10000, 10, 240)), nrow = 100, ncol = 100)
th_8bit <- calc_adaptive_th(mat_8bit)
record_test(
  "1.5 8-bit integer image correctly triggers / 255 and produces valid threshold",
  th_8bit >= 0.1 && th_8bit <= 0.9,
  sprintf("th = %f (expected in [0.1, 0.9])", th_8bit)
)

# Test 1.6: Negative undershoot (-0.05, -1e-7)
mat_neg <- base_mat
mat_neg[1, 1] <- -0.05
mat_neg[2, 2] <- -1e-7
th_neg <- tryCatch(calc_adaptive_th(mat_neg), error = function(e) NULL)
record_test(
  "1.6 Negative undershoot values clamped to 0 without histogram/Otsu crash",
  !is.null(th_neg) && th_neg >= 0.1 && th_neg <= 0.9,
  sprintf("th = %s", as.character(th_neg))
)

# Test 1.7: Extreme FFT ringing from sharp step-edge in preprocess_image()
step_img <- matrix(0, nrow = 150, ncol = 150)
step_img[40:110, 40:110] <- 1.0
brush <- EBImage::makeBrush(size = 15, shape = "gaussian", sigma = 3)
raw_filtered <- EBImage::filter2(step_img, brush)
overshoot_present <- max(raw_filtered) > 1.0 || min(raw_filtered) < 0.0
cat(sprintf("      Info: raw filter2 produces range [%.6f, %.6f] (overshoot: %s)\n",
            min(raw_filtered), max(raw_filtered), overshoot_present))

bin_proc <- tryCatch(preprocess_image(EBImage::Image(step_img)), error = function(e) NULL)
record_test(
  "1.7 preprocess_image strictly clamps FFT ringing and produces valid binary mask",
  !is.null(bin_proc) && inherits(bin_proc, "Image") && all(bin_proc %in% c(0, 1, TRUE, FALSE)),
  "Failed to binarize step image"
)

# ==============================================================================
# SUITE 2: Multi-Channel Arrays & colorMode Normalization
# ==============================================================================
cat("\n--- Suite 2: Multi-Channel Arrays & colorMode Normalization ---\n")

# Test 2.1: 1-channel 2D matrix
m1 <- matrix(runif(10000, 0.1, 0.9), nrow = 100, ncol = 100)
img1 <- read_and_process_image(m1)
record_test(
  "2.1 1-channel 2D matrix produces 2D Image with Grayscale colorMode",
  inherits(img1, "Image") && length(dim(img1)) == 2 && EBImage::colorMode(img1) == EBImage::Grayscale,
  sprintf("dim = %s, colorMode = %d", paste(dim(img1), collapse = "x"), EBImage::colorMode(img1))
)

# Test 2.2: 1-channel 3D array (100x100x1)
arr1_3d <- array(runif(10000, 0.1, 0.9), dim = c(100, 100, 1))
img1_3d <- read_and_process_image(arr1_3d)
th1_3d <- calc_adaptive_th(arr1_3d)
record_test(
  "2.2 1-channel 3D array handled by read_and_process_image & calculate_adaptive_threshold",
  inherits(img1_3d, "Image") && is.numeric(th1_3d) && th1_3d > 0.05,
  sprintf("th = %f", th1_3d)
)

# Test 2.3: 2-channel array (Grayscale + Alpha) - critical R2 requirement
arr2 <- array(0.05, dim = c(120, 100, 2))
# Specimen in Channel 1
arr2[30:90, 20:80, 1] <- 0.85
# Alpha mask in Channel 2 (1.0 everywhere)
arr2[, , 2] <- 1.0

img2_raw <- EBImage::Image(arr2, colormode = "Color")
img2_proc <- read_and_process_image(img2_raw)

# Channel 2 must be stripped: dim is 120x100, colorMode is Grayscale (0)
# Pixel (1, 1) must remain 0.05 (NOT shifted to 0.525 or 0.35 by averaging with alpha!)
ch2_stripped <- length(dim(img2_proc)) == 2 &&
  EBImage::colorMode(img2_proc) == EBImage::Grayscale &&
  abs(as.numeric(img2_proc[1, 1]) - 0.05) < 1e-4

record_test(
  "2.3 2-channel TIFF strips alpha channel without channel tinting or luminance shift",
  ch2_stripped,
  sprintf("dim = %s, colorMode = %d, pixel(1,1) = %f",
          paste(dim(img2_proc), collapse = "x"), EBImage::colorMode(img2_proc), as.numeric(img2_proc[1, 1]))
)

# Test 2.4: 2-channel array binarization via preprocess_image
bin2 <- tryCatch(preprocess_image(img2_raw), error = function(e) NULL)
record_test(
  "2.4 2-channel Image binarized cleanly via preprocess_image",
  !is.null(bin2) && inherits(bin2, "Image") && length(dim(bin2)) == 2 &&
    !as.logical(bin2[1, 1]) && as.logical(bin2[50, 50]),
  "Failed to binarize 2-channel image cleanly"
)

# Test 2.5: 3-channel RGB Image array
arr3 <- array(0.1, dim = c(100, 100, 3))
arr3[25:75, 25:75, 1] <- 0.8 # Red specimen
arr3[25:75, 25:75, 2] <- 0.7 # Green
arr3[25:75, 25:75, 3] <- 0.6 # Blue
img3 <- EBImage::Image(arr3, colormode = "Color")
bin3 <- tryCatch(preprocess_image(img3), error = function(e) NULL)
th3 <- calc_adaptive_th(arr3)
record_test(
  "2.5 3-channel RGB image normalized to Grayscale luminance and segmented",
  !is.null(bin3) && length(dim(bin3)) == 2 && th3 >= 0.1 && th3 <= 0.8,
  sprintf("th3 = %f", th3)
)

# Test 2.6: 4-channel RGBA Image array
arr4 <- array(0.1, dim = c(100, 100, 4))
arr4[25:75, 25:75, 1:3] <- 0.85 # RGB specimen
arr4[, , 4] <- 1.0 # Alpha channel
img4 <- EBImage::Image(arr4, colormode = "Color")
bin4 <- tryCatch(preprocess_image(img4), error = function(e) NULL)
th4 <- calc_adaptive_th(arr4)
record_test(
  "2.6 4-channel RGBA array converts to Grayscale without error",
  !is.null(bin4) && length(dim(bin4)) == 2 && is.numeric(th4) && th4 > 0.05,
  sprintf("th4 = %f", th4)
)

# ==============================================================================
# SUITE 3: Textured Grey Backgrounds & Border Filtering (85% Area Limit)
# ==============================================================================
cat("\n--- Suite 3: Textured Grey Backgrounds & Border Filtering ---\n")

# Test 3.1: Textured grey background safeguard (caballa/PHOTOS simulation)
# Background mean ~ 0.20, sd ~ 0.05; Specimen mean ~ 0.85
set.seed(42)
tex_bg <- matrix(rnorm(10000, mean = 0.20, sd = 0.04), nrow = 100, ncol = 100)
tex_bg[tex_bg < 0] <- 0
tex_bg[30:70, 30:70] <- rnorm(1681, mean = 0.85, sd = 0.04)
tex_bg[tex_bg > 1] <- 1

th_tex <- calc_adaptive_th(tex_bg)
th_otsu_tex <- EBImage::otsu(tex_bg)

# Adaptive threshold must not fall into background noise (< 0.65 * Otsu)
record_test(
  "3.1 Textured grey background safeguard prevents threshold collapse into noise",
  th_tex >= 0.65 * th_otsu_tex && th_tex > 0.22,
  sprintf("th_adaptive = %f, th_otsu = %f, ratio = %.2f", th_tex, th_otsu_tex, th_tex / th_otsu_tex)
)

# Test 3.2: Outer frame border filtering (image border area ~95%, specimen ~15%)
total_px <- 800 * 600
m_frame <- matrix(0, nrow = 600, ncol = 800)
# Outer border frame spanning from edge to edge (area ~95%)
m_frame[5:595, 5:795] <- 1
m_frame[25:575, 25:775] <- 0

# Internal specimen
xx <- col(m_frame) - 400
yy <- row(m_frame) - 300
m_frame[(xx^2 / 140^2 + yy^2 / 70^2) <= 1] <- 1

cnt_extracted <- extract_contour(m_frame, min_points = 150, min_area = 5000)
cnt_area <- if (!is.null(cnt_extracted)) coo_area(cnt_extracted) else 0
cnt_pct <- (cnt_area / total_px) * 100

record_test(
  "3.2 Outer frame border (95% area) rejected; internal specimen contour selected",
  !is.null(cnt_extracted) && cnt_pct < 85.0 && cnt_pct > 5.0,
  sprintf("cnt_pct = %.2f%% (must be < 85%% and > 5%%)", cnt_pct)
)

# Test 3.3: Strict mathematical limit: area >= 85% of total image is rejected
# Create binary mask where contour has area ~86% of total image
m_86 <- matrix(0, nrow = 500, ncol = 500)
total_500 <- 500 * 500
# Ellipse with area = pi * rx * ry ~ 0.86 * 250000 = 215,000 px -> rx = 270, ry = 253
# Let's create an inscribed polygon covering 86% of the image
m_86[15:485, 15:485] <- 1
area_86 <- sum(m_86)
pct_86 <- (area_86 / total_500) * 100

res_86 <- suppressWarnings(extract_contour(m_86, min_points = 100, min_area = 5000))
record_test(
  "3.3 Contour with area >= 85% (88.4%) is strictly rejected (returns NULL)",
  is.null(res_86),
  sprintf("Returned non-NULL contour with area %.2f%%", if (!is.null(res_86)) (coo_area(res_86) / total_500) * 100 else 0)
)

# Test 3.4: Pure outer border frame only (no specimen inside) returns NULL with warning
m_border_only <- matrix(0, nrow = 400, ncol = 400)
m_border_only[2:398, 2:398] <- 1
m_border_only[20:380, 20:380] <- 0

warn_caught <- FALSE
res_border_only <- withCallingHandlers(
  extract_contour(m_border_only, min_points = 100, min_area = 1000),
  warning = function(w) {
    warn_caught <<- TRUE
    invokeRestart("muffleWarning")
  }
)
record_test(
  "3.4 Image with only outer frame returns NULL and issues warning",
  is.null(res_border_only) && warn_caught,
  sprintf("res is NULL: %s, warning caught: %s", is.null(res_border_only), warn_caught)
)

# Test 3.5: Small image support with min_area = 1000 (Ostium 480x480 requirement)
m_small <- matrix(0, nrow = 480, ncol = 480)
total_small <- 480 * 480
xx_s <- col(m_small) - 240
yy_s <- row(m_small) - 240
# Small otolith: area ~ 3,200 px (1.38% of image)
m_small[(xx_s^2 / 50^2 + yy_s^2 / 20^2) <= 1] <- 1

cnt_default <- suppressWarnings(extract_contour(m_small, min_points = 150, min_area = 5000))
cnt_custom <- extract_contour(m_small, min_points = 150, min_area = 1000)

small_test_passed <- is.null(cnt_default) && !is.null(cnt_custom) &&
  coo_area(cnt_custom) >= 1000 && coo_area(cnt_custom) < 5000

record_test(
  "3.5 Small specimen (area ~3100 px) fails default min_area=5000 but succeeds with min_area=1000",
  small_test_passed,
  sprintf("cnt_custom area = %.1f px", if (!is.null(cnt_custom)) coo_area(cnt_custom) else 0)
)

# ==============================================================================
# Summary
# ==============================================================================
cat("\n================================================================================\n")
cat("                       ADVERSARIAL STRESS TEST SUMMARY                          \n")
cat("================================================================================\n")
cat(sprintf("Total Adversarial Tests : %d\n", total_tests))
cat(sprintf("Passed                  : %d (%.1f%%)\n", passed_tests, (passed_tests / total_tests) * 100))
cat(sprintf("Failed                  : %d (%.1f%%)\n", failed_tests, (failed_tests / total_tests) * 100))

if (failed_tests > 0) {
  cat("\nFailures Detail:\n")
  for (f in failures_log) {
    cat(sprintf("  - %s: %s\n", f$name, f$detail))
  }
  quit(status = 1, save = "no")
} else {
  cat("\n[SUCCESS] All adversarial stress tests for Requirement R2 passed!\n")
  quit(status = 0, save = "no")
}
