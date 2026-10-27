# ==============================================================================
# Adversarial Empirical Verification Suite for aforoR Requirements R1 & R3
# Author: teamwork_preview_challenger_1
# ==============================================================================

suppressMessages(library(EBImage))
suppressMessages(devtools::load_all(".", quiet = TRUE))

cat("==============================================================================\n")
cat("ADVERSARIAL VERIFICATION SUITE: REQUIREMENTS R1 & R3\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# 1. BUG REPRODUCTION: Vertical Chromatic Red Scale Bar Endpoint Transposition
# ------------------------------------------------------------------------------
cat("--- 1. EMPIRICAL TEST: Vertical Chromatic Red Scale Bar Endpoints ---\n")
w <- 1200; h <- 800
arr <- array(0, dim = c(w, h, 3))
xx <- row(arr[, , 1]) - 600; yy <- col(arr[, , 1]) - 400
arr[, , 1][(xx^2 / 300^2 + yy^2 / 180^2) <= 1] <- 0.9
arr[, , 2][(xx^2 / 300^2 + yy^2 / 180^2) <= 1] <- 0.9
arr[, , 3][(xx^2 / 300^2 + yy^2 / 180^2) <= 1] <- 0.9

# Physical vertical red scale bar at X in 100:107 (8 px wide), Y in 200:399 (200 px tall)
arr[100:107, 200:399, 1] <- 1.0
arr[100:107, 200:399, 2] <- 0.05
arr[100:107, 200:399, 3] <- 0.05

img_vert_red <- Image(arr, colormode = "Color")
res_vert_red <- aforoR:::find_scale_bar_length(img_vert_red, scale_mm = 1)
info_vert_red <- attr(res_vert_red, "scale_info")

cat("Ground Truth Physical Geometry:\n")
cat("  X range: [100, 107] (width 8 px), Y range: [200, 399] (height 200 px)\n")
cat("  Physical Endpoints: (x1=104, y1=200) and (x2=104, y2=399)\n")
cat("Reported by find_scale_bar_length():\n")
cat(sprintf("  min_x = %d, max_x = %d, min_y = %d, max_y = %d\n",
            info_vert_red$min_x, info_vert_red$max_x, info_vert_red$min_y, info_vert_red$max_y))
cat(sprintf("  endpoints: (x1=%d, y1=%d) to (x2=%d, y2=%d)\n",
            info_vert_red$endpoints$x1, info_vert_red$endpoints$y1,
            info_vert_red$endpoints$x2, info_vert_red$endpoints$y2))

bug_1_confirmed <- (info_vert_red$endpoints$x1 == 200 && info_vert_red$endpoints$y1 == 104)
cat(sprintf("=> BUG CONFIRMED (Coordinates Transposed X<->Y): %s\n\n", bug_1_confirmed))


# ------------------------------------------------------------------------------
# 2. EMPIRICAL TEST: Slanted / Rotated Scale Bar True Span Underestimation
# ------------------------------------------------------------------------------
cat("--- 2. EMPIRICAL TEST: Rotated / Slanted Scale Bars ---\n")
angles <- c(0, 5, 10, 15, 30, 45)
bar_len <- 200
hl <- bar_len / 2; ht <- 4
for (a in angles) {
  m <- matrix(0, 1000, 1000)
  xx <- col(m) - 500; yy <- row(m) - 400
  m[(xx^2 / 200^2 + yy^2 / 120^2) <= 1] <- 1
  rad <- a * pi / 180
  grid_x <- col(m) - 750; grid_y <- row(m) - 750
  rot_u <- grid_x * cos(-rad) - grid_y * sin(-rad)
  rot_v <- grid_x * sin(-rad) + grid_y * cos(-rad)
  m[(abs(rot_u) <= hl) & (abs(rot_v) <= ht)] <- 1
  
  res_rot <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  info_rot <- attr(res_rot, "scale_info")
  det_len <- info_rot$bar_length_px
  err <- abs(det_len - bar_len)
  eucl_span <- sqrt((info_rot$endpoints$x2 - info_rot$endpoints$x1)^2 + 
                    (info_rot$endpoints$y2 - info_rot$endpoints$y1)^2)
  cat(sprintf("  Angle %2d deg: max(dx,dy)=%.0f px (err=%.0f px) | Euclidean span=%.1f px\n",
              a, det_len, err, eucl_span))
}
cat("\n")


# ------------------------------------------------------------------------------
# 3. EMPIRICAL TEST: R1 Scale Bar Thickness & Noise Robustness
# ------------------------------------------------------------------------------
cat("--- 3. EMPIRICAL TEST: Line Thicknesses & Red Channel Robustness ---\n")
# Thickness 1 to 20 px
th_pass <- TRUE
for (th in c(1, 2, 5, 10, 20)) {
  m <- matrix(0, 1000, 1000)
  xx <- col(m) - 500; yy <- row(m) - 400
  m[(xx^2 / 200^2 + yy^2 / 120^2) <= 1] <- 1
  m[750:(750+th-1), 400:599] <- 1
  res_th <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
  info_th <- attr(res_th, "scale_info")
  if (is.null(info_th) || abs(info_th$bar_length_px - 200) > 2) th_pass <- FALSE
}
cat(sprintf("  Varying line thicknesses (1 to 20 px) within <= 2 px: %s\n", th_pass))

# Red channel text discrimination
arr_text <- arr
arr_text[880:920, 670:685, 1] <- 1.0 # simulated "1 mm" text block
res_text <- aforoR:::find_scale_bar_length(Image(arr_text, colormode="Color"), scale_mm = 1)
info_text <- attr(res_text, "scale_info")
cat(sprintf("  Red bar selected over red text label (length = %d px): %s\n\n",
            if(!is.null(info_text)) info_text$bar_length_px else 0,
            !is.null(info_text) && info_text$bar_length_px == 200))


# ------------------------------------------------------------------------------
# 4. EMPIRICAL TEST: R3 plot.css() Parameter Flexibility & Empty Matrix Handling
# ------------------------------------------------------------------------------
cat("--- 4. EMPIRICAL TEST: plot.css() S3 Method Flexibility ---\n")
theta <- seq(0, 2 * pi, length.out = 80)
css_obj <- calculate_css(cbind(cos(theta) + 0.2*cos(4*theta), sin(theta) + 0.2*sin(4*theta)), n_points = 80, sigma_max = 3)

pdf(NULL)
on.exit(try(dev.off(), silent = TRUE), add = TRUE)

r3_tests <- list(
  "Custom title, col, pch, cex" = function() plot(css_obj, main = "Title", col = "red", pch = 16, cex = 1.2),
  "Custom xlim, ylim, xlab, ylab" = function() plot(css_obj, xlim = c(0, 0.5), ylim = c(0, 2), xlab = "S", ylab = "Sigma"),
  "Unusual graphical params (bty, las, xaxt, yaxt)" = function() plot(css_obj, bty = "n", las = 2, xaxt = "n", yaxt = "n"),
  "Vector col and pch" = function() plot(css_obj, col = 1:nrow(css_obj$css_matrix), pch = 1:nrow(css_obj$css_matrix)),
  "panel.first hook" = function() plot(css_obj, panel.first = quote(graphics::grid())),
  "Inverted zoom limits" = function() plot(css_obj, xlim = c(1, 0), ylim = c(3, 0)),
  "NULL css_matrix with custom title" = function() plot(structure(list(css_matrix = NULL), class = "css"), main = "Empty"),
  "0-row data.frame css_matrix" = function() plot(structure(list(css_matrix = data.frame(arc_length=numeric(0), sigma=numeric(0))), class = "css"), col = "blue")
)

r3_all_pass <- TRUE
for (nm in names(r3_tests)) {
  err <- NULL
  tryCatch(r3_tests[[nm]](), error = function(e) err <<- conditionMessage(e))
  if (is.null(err)) {
    cat(sprintf("  [PASS] %s\n", nm))
  } else {
    cat(sprintf("  [FAIL] %s: %s\n", nm, err))
    r3_all_pass <- FALSE
  }
}
cat(sprintf("\n=> R3 Overall Stress Test Status: %s\n", if(r3_all_pass) "ALL PASSED" else "FAILURES DETECTED"))

cat("==============================================================================\n")
