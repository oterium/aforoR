test_that("calculate_distances handles degenerate contours and boundary cases", {
    # Less than 3 points
    expect_error(calculate_distances(data.frame(X = numeric(0), Y = numeric(0))), "contour must contain at least 3 points")
    expect_error(calculate_distances(data.frame(X = c(1, 2), Y = c(1, 2))), "contour must contain at least 3 points")

    # Identical coordinates (degenerate boundary)
    expect_error(calculate_distances(data.frame(X = rep(10, 5), Y = rep(10, 5))), "contour is degenerate")

    # NA coordinates
    expect_error(calculate_distances(data.frame(X = c(1, 2, NA), Y = c(1, 2, 3))), "contour coordinates cannot contain NA values")

    # Invalid start_point
    cnt <- data.frame(X = c(0, 5, 5, 0), Y = c(0, 0, 5, 5))
    expect_error(calculate_distances(cnt, start_point = c(NA, 1)), "start_point must be 'auto', 'manual', or a valid numeric")
    expect_error(calculate_distances(cnt, start_point = c(Inf, 0)), "start_point must be 'auto', 'manual', or a valid numeric")
    expect_error(calculate_distances(cnt, start_point = 0), "start_point index must be between 1 and")
    expect_error(calculate_distances(cnt, start_point = 10), "start_point index must be between 1 and")

    # Negative coordinates
    cnt_neg <- data.frame(X = c(-10, -5, -5, -10), Y = c(-10, -10, -5, -5))
    res_neg <- calculate_distances(cnt_neg, n_points = 16)
    expect_type(res_neg, "list")
    expect_length(res_neg$polar$radii, 16)
    expect_true(all(is.finite(res_neg$polar$normalized)))
    expect_true(all(is.finite(res_neg$perimeter$normalized)))
})

test_that("regularradius, dper, and ild handle edge cases cleanly", {
    # ild with infinite values
    expect_error(ild(c(Inf, 0), c(0, 0)), "cannot contain infinite values")
    expect_error(ild(c(0, 0), c(0, -Inf)), "cannot contain infinite values")

    # regularradius degenerate inputs
    expect_error(regularradius(c(1, 2), c(1, 2), n = 10), "must contain at least 3 points")
    expect_error(regularradius(rep(3, 10), rep(4, 10), n = 10), "degenerate")
    expect_error(regularradius(1:5, 1:5, n = 0), "positive integer")

    # regularradius with n = 1 retains 2D matrix structure in coord
    res_n1 <- regularradius(c(0, 5, 2), c(0, 0, 4), n = 1)
    expect_true(is.matrix(res_n1$coord))
    expect_equal(nrow(res_n1$coord), 1)

    # dper degenerate inputs
    expect_error(dper(c(1, 2), c(1, 2), n = 10), "at least 3 points")
    expect_error(dper(rep(2, 5), rep(2, 5), n = 10), "degenerate")
    expect_error(dper(1:5, 1:5, n = 1), "integer >= 2")
})

test_that("calculate_morphometrics and relative indices handle edge cases", {
    # calculate_morphometrics validation
    expect_error(calculate_morphometrics(data.frame(X = numeric(0), Y = numeric(0))), "at least 3 points")
    expect_error(calculate_morphometrics(data.frame(X = c(1, 2), Y = c(1, 2))), "at least 3 points")
    expect_error(calculate_morphometrics(data.frame(X = rep(1, 5), Y = rep(1, 5))), "degenerate")
    expect_error(calculate_morphometrics(data.frame(X = c(1, 2, NA), Y = c(1, 2, 3))), "cannot contain NA")

    # calculate_rol
    expect_error(calculate_rol(numeric(0), 10), "cannot be empty")
    expect_error(calculate_rol(5, numeric(0)), "cannot be empty")
    expect_error(calculate_rol(-1, 10), "cannot be negative")
    expect_error(calculate_rol(5, 0), "must be strictly positive")
    expect_error(calculate_rol(5, -10), "must be strictly positive")
    expect_error(calculate_rol(c(1, 2), c(10, 20, 30)), "Lengths of otolith_length and fish_length must match")

    # calculate_ros
    expect_error(calculate_ros(numeric(0), 10), "cannot be empty")
    expect_error(calculate_ros(-5, 10), "cannot be negative")
    expect_error(calculate_ros(5, 0), "must be strictly positive")
    expect_error(calculate_ros(c(1, 2), c(10, 20, 30)), "Lengths of otolith_area and fish_length must match")

    # calculate_relative_indices
    expect_error(calculate_relative_indices(c(1, 2), c(1, 2, 3), 10), "Lengths of otolith_length and otolith_area must match")
    rel <- calculate_relative_indices(5, 15, 50)
    expect_s3_class(rel, "data.frame")
    expect_equal(rel$ROL, 10)
})

test_that("calculate_css handles boundary conditions", {
    expect_error(calculate_css(data.frame(X = numeric(0), Y = numeric(0))), "at least 3 points")
    expect_error(calculate_css(data.frame(X = c(1, 2), Y = c(1, 2))), "at least 3 points")
    expect_error(calculate_css(data.frame(X = rep(5, 10), Y = rep(5, 10))), "degenerate")
    expect_error(calculate_css(data.frame(X = c(1, 2, NA), Y = c(1, 2, 3))), "cannot contain NA values")
    expect_error(calculate_css(data.frame(X = 1:5, Y = 1:5), n_points = 2), "integer >= 3")
})

test_that("fatrous1d validates filter length", {
    # rh length must match snmax - snmin + 1 (here 2 - (-2) + 1 = 5)
    expect_error(fatrous1d(1:10, c(0.1, 0.2), -2, 2, 1), "rh must have length equal to snmax", fixed = FALSE)
    # Correct length succeeds
    flt <- c(0.0625, 0.25, 0.375, 0.25, 0.0625)
    res <- fatrous1d(1:16, flt, -2, 2, 1)
    expect_true(is.matrix(res))
    expect_equal(ncol(res), 16)
})

test_that("extract_contour validates min_points and min_area", {
    dummy_img <- matrix(0, 10, 10)
    expect_error(extract_contour(dummy_img, min_points = 1), "min_points must be an integer >= 3")
    expect_error(extract_contour(dummy_img, min_area = -5), "min_area must be 'auto', NULL, or a numeric value >= 0")
    expect_error(extract_contour(dummy_img, min_area = "invalid"), "min_area must be 'auto', NULL, or a numeric value >= 0")
})

test_that("calculate_adaptive_threshold works on numeric vectors and handles extreme values", {
    # Numeric vector (array conversion for EBImage::otsu)
    vals <- runif(200, 0.1, 0.9)
    th <- aforoR:::calculate_adaptive_threshold(vals)
    expect_true(is.numeric(th) && length(th) == 1)
    expect_true(th > 0 && th < 1)

    # Constant all 0s
    th_zero <- aforoR:::calculate_adaptive_threshold(matrix(0, 20, 20))
    expect_true(is.numeric(th_zero) && length(th_zero) == 1)

    # Constant all 1s
    th_one <- aforoR:::calculate_adaptive_threshold(matrix(1, 20, 20))
    expect_true(is.numeric(th_one) && length(th_one) == 1)

    # High range values
    th_high <- aforoR:::calculate_adaptive_threshold(runif(500, 0.85, 0.95))
    expect_true(is.numeric(th_high) && length(th_high) == 1)
})

test_that("save_visualization and save_visualization_perimeter clamp scale gracefully", {
    skip_if_not_installed("EBImage")

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    bin <- preprocess_image(example_image)
    cnt <- extract_contour(bin)
    dists <- calculate_distances(cnt, n_points = 256)
    # Only 3 scales computed
    wavs <- calculate_wavelets_analysis(dists, n_scales = 3)

    temp_out <- file.path(tempdir(), "test_audit_viz")
    if (dir.exists(temp_out)) unlink(temp_out, recursive = TRUE)
    dir.create(temp_out)

    # Requesting scale 5 when only 3 exist: should not crash with out-of-bounds error
    expect_silent(save_visualization(
        image = NULL, distances = dists, wavelets = wavs,
        image_name = "scale_clamp_test.jpg", output_dir = temp_out,
        wavelet_scale = 5, plot_width = 400, plot_height = 300, plot_res = 72
    ))
    expect_true(file.exists(file.path(temp_out, "wavelet_scale_clamp_test.jpg")))

    expect_silent(save_visualization_perimeter(
        image = NULL, distances = dists, wavelets = wavs,
        image_name = "scale_clamp_perim.jpg", output_dir = temp_out,
        wavelet_scale = 5, plot_width = 400, plot_height = 300, plot_res = 72
    ))
    expect_true(file.exists(file.path(temp_out, "wavelet_scale_clamp_perim.jpg")))

    unlink(temp_out, recursive = TRUE)
})

test_that("select_points_hall input validation", {
    # Mismatched rows and grouping length
    expect_error(select_points_hall(matrix(1:10, 5, 2), factor(c("A", "B"))), "nrow\\(data\\) must match length\\(grouping\\)")

    # Insufficient observations
    expect_error(select_points_hall(matrix(1:2, 1, 2), factor("A")), "at least 2 observations")

    # Single class
    expect_error(select_points_hall(matrix(1:10, 5, 2), factor(rep("A", 5))), "at least 2 distinct classes")

    # Non-numeric data
    expect_error(select_points_hall(matrix(c("a", "b", "c", "d"), 2, 2), factor(c("A", "B"))), "data must be numeric")

    # NA in data or grouping
    expect_error(select_points_hall(matrix(c(1, NA, 3, 4), 2, 2), factor(c("A", "B"))), "data cannot contain NA values")
    expect_error(select_points_hall(matrix(1:4, 2, 2), factor(c("A", NA))), "grouping cannot contain NA values")

    # Invalid parameters
    expect_error(select_points_hall(matrix(1:4, 2, 2), factor(c("A", "B")), cv = "yes"), "cv must be a logical value")
    expect_error(select_points_hall(matrix(1:4, 2, 2), factor(c("A", "B")), parallel = 1), "parallel must be a logical value")
    expect_error(select_points_hall(matrix(1:4, 2, 2), factor(c("A", "B")), p = -0.5), "p must be a positive numeric value")
    expect_error(select_points_hall(matrix(1:4, 2, 2), factor(c("A", "B")), delta_t = 0), "delta_t must be a numeric value between 0 and 1")
    expect_error(select_points_hall(matrix(1:4, 2, 2), factor(c("A", "B")), method = "unsupported"), "'arg' should be one of")
})

test_that("select_points_hall runs with parallel = TRUE across cores without export leaks", {
    skip_if_not_installed("MASS")
    set.seed(123)
    data <- matrix(rnorm(30 * 8), 30, 8)
    data[1:15, 2] <- data[1:15, 2] + 4
    data[16:30, 5] <- data[16:30, 5] + 4
    groups <- factor(rep(c("Class1", "Class2"), each = 15))

    res_par <- select_points_hall(data, groups, method = "lda", cv = FALSE, parallel = TRUE, n_cores = 2)
    expect_type(res_par, "list")
    expect_true(length(res_par$selected_indices) >= 1)
    expect_true(all(res_par$selected_indices %in% 1:8))
    expect_true(is.table(res_par$confusion_matrix))
})

test_that("calculate_adaptive_threshold handles 3D arrays, 0..255 scale, and NULL input", {
    # NULL input
    expect_equal(aforoR:::calculate_adaptive_threshold(NULL), 0.2)

    # 3D array (RGB format: width x height x channels)
    set.seed(42)
    rgb_arr <- array(runif(300, 0, 1), dim = c(10, 10, 3))
    th_rgb <- aforoR:::calculate_adaptive_threshold(rgb_arr)
    expect_true(is.numeric(th_rgb) && length(th_rgb) == 1)
    expect_false(is.na(th_rgb))
    expect_true(th_rgb > 0 && th_rgb < 1)

    # 8-bit scale inputs (0..255)
    mat_255 <- matrix(sample(0:255, 400, replace = TRUE), 20, 20)
    th_255 <- aforoR:::calculate_adaptive_threshold(mat_255)
    expect_true(is.numeric(th_255) && length(th_255) == 1)
    expect_true(th_255 > 0 && th_255 <= 1)
})

test_that("fatrous1d validates scale index sj >= 1 and rejects invalid/NA inputs", {
    flt <- c(0.0625, 0.25, 0.375, 0.25, 0.0625)
    expect_error(fatrous1d(1:10, flt, -2, 2, sj = 0), "sj must be an integer >= 1")
    expect_error(fatrous1d(1:10, flt, -2, 2, sj = -1), "sj must be an integer >= 1")
    expect_error(fatrous1d(1:10, flt, -2, 2, sj = NA), "must be numeric values")
    expect_error(fatrous1d(1:10, flt, -2, 2, sj = NA_real_), "must be finite values")
    expect_error(fatrous1d(1:10, flt, -2, 2, sj = Inf), "must be finite values")
    expect_error(fatrous1d(c(1, NA, 3), flt, -2, 2, sj = 1), "rx and rh cannot contain NA values")
    expect_error(fatrous1d(1:10, c(flt[1:4], NA), -2, 2, sj = 1), "rx and rh cannot contain NA values")
    expect_error(fatrous1d(c(1, Inf, 3), flt, -2, 2, sj = 1), "rx and rh cannot contain infinite values")
})

test_that("calculate_css handles 2-column non-standard data.frames and validates sigma arguments", {
    theta <- seq(0, 2 * pi, length.out = 50)
    # Lowercase or custom column names
    cnt_df <- data.frame(pos_x = cos(theta) + 0.1 * cos(3 * theta), pos_y = sin(theta) + 0.1 * sin(3 * theta))
    res_css <- calculate_css(cnt_df, n_points = 50, sigma_max = 5, sigma_step = 1)
    expect_s3_class(res_css, "css")

    # Invalid sigma
    expect_error(calculate_css(cnt_df, sigma_max = -1), "sigma_max must be a positive number")
    expect_error(calculate_css(cnt_df, sigma_max = 0), "sigma_max must be a positive number")
    expect_error(calculate_css(cnt_df, sigma_step = -0.5), "sigma_step must be a positive number")
    expect_error(calculate_css(cnt_df, sigma_step = 0), "sigma_step must be a positive number")

    # Non-numeric or infinite coordinates
    expect_error(calculate_css(data.frame(x = c("a", "b", "c", "d"), y = 1:4)), "coordinates must be numeric")
    expect_error(calculate_css(data.frame(x = c(1, Inf, 3, 4), y = 1:4)), "cannot contain infinite values")
})

test_that("calculate_morphometrics validates numeric coordinates and infinite values", {
    expect_error(calculate_morphometrics(data.frame(X = c("a", "b", "c"), Y = 1:3)), "coordinates must be numeric")
    expect_error(calculate_morphometrics(data.frame(X = c(1, Inf, 3), Y = 1:3)), "cannot contain infinite values")
})

test_that("calculate_distances strictly validates n_points", {
    cnt <- data.frame(X = c(0, 5, 5, 0), Y = c(0, 0, 5, 5))
    expect_error(calculate_distances(cnt, n_points = c(10, 20)), "n_points must be an integer >= 2")
    expect_error(calculate_distances(cnt, n_points = 1), "n_points must be an integer >= 2")
    expect_error(calculate_distances(cnt, n_points = NA), "n_points must be an integer >= 2")
})

test_that("select_starting_point gracefully returns NULL in non-interactive sessions", {
    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    # In test runner (non-interactive), should return NULL without plotting to device
    res_pt <- select_starting_point(example_image)
    expect_null(res_pt)
})

test_that("circ_conv handles n_pad >= N without negative indexing errors", {
    # High-frequency contour with N = 20 points and sigma = 6 (n_pad = 24 >= 20)
    theta <- seq(0, 2 * pi, length.out = 20)
    x <- cos(theta) + 0.3 * cos(5 * theta)
    y <- sin(theta) + 0.3 * sin(5 * theta)
    cnt_df <- data.frame(x = x, y = y)

    expect_no_error(res <- calculate_css(cnt_df, n_points = 20, sigma_max = 6, sigma_step = 6))
    expect_s3_class(res, "css")

    # plot.css validation and return value
    expect_error(aforoR:::plot.css(NULL), "x must be an object of class 'css'")
    expect_error(aforoR:::plot.css(list(a = 1)), "x must be an object of class 'css'")

    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    plot_ret <- plot(res)
    expect_null(plot_ret)
})

test_that("fwaveletspl_3 validates NA, infinite, non-integer, and small signals", {
    # sjm validation
    expect_error(fwaveletspl_3(1:8, sjm = NA_real_), "sjm must be a positive integer")
    expect_error(fwaveletspl_3(1:8, sjm = Inf), "sjm must be a positive integer")
    expect_error(fwaveletspl_3(1:8, sjm = 2.5), "sjm must be a positive integer")
    expect_error(fwaveletspl_3(1:8, sjm = -1), "sjm must be a positive integer")

    # det validation
    expect_error(fwaveletspl_3(1:8, sjm = 2, det = NA), "det must be a logical value")

    # ra0 validation
    expect_error(fwaveletspl_3(c(1, Inf, 3, 4), sjm = 2), "ra0 cannot contain infinite values")
    expect_error(fwaveletspl_3(1, sjm = 1), "ra0 must contain at least 2 points")
})

test_that("fatrous1d strictly rejects non-integer filter scale parameters", {
    flt <- c(0.0625, 0.25, 0.375, 0.25, 0.0625)
    expect_error(fatrous1d(1:10, flt, -2.5, 2.5, 1), "snmin and snmax must be integers")
    expect_error(fatrous1d(1:10, flt, -2, 2.5, 1), "snmin and snmax must be integers")
    expect_error(fatrous1d(1:10, flt, -2, 2, sj = 1.5), "sj must be an integer >= 1")
})

test_that("calculate_distances supports 2-column lowercase frames and validates coordinates", {
    # Lowercase column names x, y
    cnt_lower <- data.frame(x = c(0, 5, 5, 0), y = c(0, 0, 5, 5))
    res_lower <- calculate_distances(cnt_lower, n_points = 16)
    expect_type(res_lower, "list")
    expect_named(res_lower, c("polar", "perimeter", "reordered_coords"))

    # Arbitrary 2 columns
    cnt_arb <- data.frame(col1 = c(0, 5, 5, 0), col2 = c(0, 0, 5, 5))
    res_arb <- calculate_distances(cnt_arb, n_points = 16)
    expect_type(res_arb, "list")

    # Non-numeric or infinite coordinates
    expect_error(calculate_distances(data.frame(X = c("a", "b", "c"), Y = 1:3)), "coordinates must be numeric")
    expect_error(calculate_distances(data.frame(X = c(1, Inf, 3), Y = 1:3)), "cannot contain infinite values")
    expect_error(calculate_distances(cnt_lower, n_points = 2.5), "n_points must be an integer >= 2")
})

test_that("calculate_morphometrics supports 2-column lowercase frames and validates infinite values", {
    cnt_lower <- data.frame(x = c(0, 5, 5, 0), y = c(0, 0, 5, 5))
    res_morpho <- calculate_morphometrics(cnt_lower)
    expect_type(res_morpho, "list")
    expect_true(res_morpho$Area > 0)

    # ROL and ROS with infinite inputs
    expect_error(calculate_rol(c(1, Inf), c(10, 20)), "cannot contain infinite values")
    expect_error(calculate_rol(c(1, 2), c(10, Inf)), "cannot contain infinite values")
    expect_error(calculate_ros(c(1, Inf), c(10, 20)), "cannot contain infinite values")
    expect_error(calculate_ros(c(1, 2), c(10, Inf)), "cannot contain infinite values")
})

test_that("regularradius and dper strictly reject infinite coordinates and non-integers", {
    expect_error(regularradius(c(1, Inf, 3), c(1, 2, 3), n = 10), "cannot contain infinite values")
    expect_error(regularradius(c(1, 2, 3), c(1, -Inf, 3), n = 10), "cannot contain infinite values")
    expect_error(regularradius(1:5, 1:5, n = 5.5), "n must be a positive integer")
    expect_error(regularradius(1:5, 1:5, n = Inf), "n must be a positive integer")

    expect_error(dper(c(1, Inf, 3), c(1, 2, 3), n = 10), "cannot contain infinite values")
    expect_error(dper(c(1, 2, 3), c(1, -Inf, 3), n = 10), "cannot contain infinite values")
    expect_error(dper(1:5, 1:5, n = 5.5), "n must be an integer >= 2")
    expect_error(dper(1:5, 1:5, n = Inf), "n must be an integer >= 2")
})

test_that("preprocess_image and find_scale_bar_length validate edge cases", {
    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    # Invalid threshold values
    expect_error(preprocess_image(example_image, threshold = NA_real_), "between 0 and 1")
    expect_error(preprocess_image(example_image, threshold = -0.1), "between 0 and 1")
    expect_error(preprocess_image(example_image, threshold = 1.5), "between 0 and 1")

    # Invalid blur arguments
    expect_error(preprocess_image(example_image, blur_size = NA), "positive numeric value")
    expect_error(preprocess_image(example_image, blur_size = -5), "positive numeric value")
    expect_error(preprocess_image(example_image, blur_sigma = NA), "positive numeric value")
    expect_error(preprocess_image(example_image, blur_sigma = -1), "positive numeric value")

    # extract_contour non-integer min_points
    dummy_img <- matrix(0, 10, 10)
    expect_error(extract_contour(dummy_img, min_points = 3.5), "min_points must be an integer >= 3")

    # find_scale_bar_length invalid scale_mm
    expect_error(aforoR:::find_scale_bar_length(dummy_img, scale_mm = 0), "positive numeric value")
    expect_error(aforoR:::find_scale_bar_length(dummy_img, scale_mm = -1), "positive numeric value")
})

test_that("save_visualization validates plot dimensions and handles flat bounding box", {
    dists <- list(
        polar = list(radii = 1:10, coord = matrix(rep(1, 20), 10, 2), normalized = runif(10)),
        perimeter = list(dist = 1:10, coords = matrix(rep(1, 20), 10, 2), normalized = runif(10)),
        reordered_coords = list(x = rep(1, 10), y = rep(1, 10))
    )
    tmp_out <- file.path(tempdir(), "test_viz_bounds")
    dir.create(tmp_out, showWarnings = FALSE)
    on.exit(unlink(tmp_out, recursive = TRUE), add = TRUE)

    # Invalid plot dimensions
    expect_error(save_visualization(distances = dists, wavelets = NULL, image_name = "test.jpg",
                                    output_dir = tmp_out, plot_width = -10), "plot_width must be a positive number")
    expect_error(save_visualization(distances = dists, wavelets = NULL, image_name = "test.jpg",
                                    output_dir = tmp_out, plot_height = 0), "plot_height must be a positive number")
    expect_error(save_visualization(distances = dists, wavelets = NULL, image_name = "test.jpg",
                                    output_dir = tmp_out, plot_res = NA), "plot_res must be a positive number")

    # Flat bounding box (diff(xlims) == 0): should not fail
    expect_silent(save_visualization(distances = dists, wavelets = NULL, image_name = "flat_box.jpg",
                                     output_dir = tmp_out, plot_width = 400, plot_height = 300, plot_res = 72))
    expect_silent(save_visualization_perimeter(distances = dists, wavelets = NULL, image_name = "flat_box_perim.jpg",
                                               output_dir = tmp_out, plot_width = 400, plot_height = 300, plot_res = 72))
})

test_that("create_result_structure handles vector wavelets safely", {
    dists <- list(radii = runif(10), dist = runif(10), normalized = runif(10), coord = matrix(1:20, 10, 2))
    res_vec <- create_result_structure(dists, 1:10, "test_img", result_type = "polar")
    expect_true("Wavelet_1" %in% names(res_vec))
    expect_equal(res_vec$Wavelet_1, 1:10)
})

test_that("process_images runs GPA procrustes alignment without unbound variable errors", {
    skip_if_not_installed("EBImage")

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    tmp_dir <- file.path(tempdir(), "test_gpa_audit")
    if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    # Copy two rotated versions of the test image
    bin <- preprocess_image(example_image)
    cnt <- extract_contour(bin)
    file.copy(example_image, file.path(tmp_dir, "otolith_1.jpg"))
    file.copy(example_image, file.path(tmp_dir, "otolith_2.jpg"))

    # Direct test of .apply_gpa_alignment with instantiated Out object
    c1 <- aforoR:::.coo_interpolate(cnt, 128)
    th <- 0.05; c2 <- c1 %*% matrix(c(cos(th), -sin(th), sin(th), cos(th)), 2, 2)
    aligned_res <- aforoR:::.apply_gpa_alignment(list(c1, c2), 1:2, procrustes = TRUE, n_points = 128)
    expect_length(aligned_res, 2)
    expect_true(is.matrix(aligned_res[[1]]))
    expect_equal(nrow(aligned_res[[1]]), 128)
})

test_that("save_visualization and save_visualization_perimeter handle small n_points (2 and 3) without subscript out of bounds", {
    tmp_viz <- file.path(tempdir(), "test_small_n_points_viz")
    if (dir.exists(tmp_viz)) unlink(tmp_viz, recursive = TRUE)
    dir.create(tmp_viz)
    on.exit(unlink(tmp_viz, recursive = TRUE), add = TRUE)

    cnt <- data.frame(X = c(1, 4, 2), Y = c(1, 5, 3))

    # Test n_points = 2
    dists_2 <- calculate_distances(cnt, n_points = 2)
    expect_silent(save_visualization(NULL, dists_2, NULL, "test_n2.jpg", tmp_viz, plot_width = 300, plot_height = 200, plot_res = 72))
    expect_silent(save_visualization_perimeter(NULL, dists_2, NULL, "test_n2_perim.jpg", tmp_viz, plot_width = 300, plot_height = 200, plot_res = 72))

    # Test n_points = 3
    dists_3 <- calculate_distances(cnt, n_points = 3)
    expect_silent(save_visualization(NULL, dists_3, NULL, "test_n3.jpg", tmp_viz, plot_width = 300, plot_height = 200, plot_res = 72))
    expect_silent(save_visualization_perimeter(NULL, dists_3, NULL, "test_n3_perim.jpg", tmp_viz, plot_width = 300, plot_height = 200, plot_res = 72))

    # Test with wavelets on small n_points
    wavs_3 <- calculate_wavelets_analysis(dists_3, n_scales = 1)
    expect_silent(save_visualization(NULL, dists_3, wavs_3, "test_n3_wav.jpg", tmp_viz, wavelet_scale = 1, plot_width = 300, plot_height = 200, plot_res = 72))
    expect_silent(save_visualization_perimeter(NULL, dists_3, wavs_3, "test_n3_wav_perim.jpg", tmp_viz, wavelet_scale = 1, plot_width = 300, plot_height = 200, plot_res = 72))
})

test_that(".compute_curvature_raw and calculate_css avoid division by zero and NaN on zero-velocity points", {
    # Contour with zero velocity between i_prev and i_next
    coo_degenerate <- matrix(c(1, 1.666667, 1.666667, 1, 1.666667, 1.666667), 3, 2)
    raw_curv <- aforoR:::.compute_curvature_raw(coo_degenerate)
    expect_false(any(is.nan(raw_curv)))
    expect_true(all(is.finite(raw_curv)))

    # calculate_css with small closed shape
    res_css <- calculate_css(matrix(c(1, 2, 3, 1, 1, 4, 2, 1), 4, 2), n_points = 3, sigma_max = 1)
    expect_s3_class(res_css, "css")
    expect_true(is.null(res_css$css_matrix) || !any(is.nan(res_css$css_matrix$arc_length)))
})

test_that("calculate_morphometrics strictly validates pixels_per_mm", {
    coo <- matrix(1:6, 3, 2)
    expect_error(calculate_morphometrics(coo, pixels_per_mm = NA_real_), "pixels_per_mm must be a positive numeric value")
    expect_error(calculate_morphometrics(coo, pixels_per_mm = c(10, 20)), "pixels_per_mm must be a positive numeric value")
    expect_error(calculate_morphometrics(coo, pixels_per_mm = Inf), "pixels_per_mm must be a positive numeric value")
    expect_error(calculate_morphometrics(coo, pixels_per_mm = -5), "pixels_per_mm must be a positive numeric value")
    expect_error(calculate_morphometrics(coo, pixels_per_mm = 0), "pixels_per_mm must be a positive numeric value")
})

test_that("process_images validates cores, pixels_per_mm, scale_mm, and n_points", {
    expect_error(process_images(tempdir(), cores = NA_real_), "cores must be a positive integer or NULL")
    expect_error(process_images(tempdir(), cores = c(2, 4)), "cores must be a positive integer or NULL")
    expect_error(process_images(tempdir(), cores = 2.5), "cores must be a positive integer or NULL")
    expect_error(process_images(tempdir(), pixels_per_mm = NA_real_), "pixels_per_mm must be a positive numeric value")
    expect_error(process_images(tempdir(), scale_mm = -1), "scale_mm must be a positive numeric value")
    expect_error(process_images(tempdir(), n_points = 2), "n_points must be an integer >= 3")
})

test_that("process_images saves MorphometricsEN.csv when coordinates = 'cartesian'", {
    tmp_cart <- file.path(tempdir(), "test_cart_morpho_save")
    if (dir.exists(tmp_cart)) unlink(tmp_cart, recursive = TRUE)
    dir.create(tmp_cart)
    on.exit(unlink(tmp_cart, recursive = TRUE), add = TRUE)

    # Test .save_batch_csv_results directly with coordinates = "cartesian"
    dummy_morpho <- list(list(Image = "img1.jpg", Area = 100, Perimeter = 40, Length = 15, Width = 8, Feret_Max = 16, Feret_Min = 7, PCA_Angle = 45, Units = "px"))
    dummy_res <- list(list(Distancia = runif(10), Distancia_Norm = runif(10), nombre = "img1.jpg", coords = matrix(1:20, 10, 2)))
    cart_dir <- file.path(tmp_cart, "Cartesian")
    dir.create(cart_dir)

    aforoR:::.save_batch_csv_results(result = list(), result2 = dummy_res, morpho_results = dummy_morpho,
                                     polar_dir = file.path(tmp_cart, "Polar"), cartesian_dir = cart_dir,
                                     coordinates = "cartesian")

    expect_true(file.exists(file.path(cart_dir, "MorphometricsEN.csv")))
    expect_true(file.exists(file.path(cart_dir, "DistanciaEN.csv")))
})

test_that("create_result_structure handles missing coord/coords cleanly", {
    dists <- list(radii = runif(10), dist = runif(10), normalized = runif(10))
    res <- create_result_structure(dists, NULL, "sample_no_coord", result_type = "polar")
    expect_null(res$coords)
})



test_that("process_images enforces subfolder = FALSE and matching row count for custom start_points matrix", {
    expect_error(
        process_images(tempdir(), subfolder = TRUE, start_points = matrix(1:4, nrow = 2, ncol = 2)),
        "subfolder = TRUE is not allowed when providing a custom start_points matrix"
    )
    expect_error(
        process_images(tempdir(), subfolder = FALSE, start_points = matrix(1:2, nrow = 2, ncol = 1)),
        "at least 2 columns"
    )

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    tmp_sp <- file.path(tempdir(), "test_sp_matrix")
    if (dir.exists(tmp_sp)) unlink(tmp_sp, recursive = TRUE)
    dir.create(tmp_sp)
    on.exit(unlink(tmp_sp, recursive = TRUE), add = TRUE)

    file.copy(example_image, file.path(tmp_sp, "oto1.jpg"))
    file.copy(example_image, file.path(tmp_sp, "oto2.jpg"))

    # Mismatched row count (3 rows for 2 images) should error
    expect_error(
        process_images(tmp_sp, subfolder = FALSE, start_points = matrix(1:6, nrow = 3, ncol = 2),
                       wavelets = FALSE, ef = FALSE, testing = FALSE, save = FALSE),
        "must equal the total number of images in the folder"
    )
})

test_that("select_scale_bar and detect_scale modes validate inputs and handle non-interactive sessions", {
    expect_error(select_scale_bar("nonexistent.jpg", scale_mm = 0), "scale_mm must be a positive numeric value")
    expect_error(select_scale_bar("nonexistent.jpg", scale_mm = -1), "scale_mm must be a positive numeric value")
    expect_error(select_scale_bar("nonexistent.jpg", scale_mm = 1), "image must be an Image object or a valid file path")
    expect_error(process_images(tempdir(), detect_scale = "invalid_mode"), "detect_scale must be 'none', 'auto', or 'manual'")

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    # In non-interactive session, select_scale_bar returns NULL cleanly
    expect_null(select_scale_bar(example_image, scale_mm = 1))
})

test_that("fwaveletspl_3 and calculate_wavelets_analysis return detail when det/detail = TRUE (default) and approximation when FALSE", {
    theta <- seq(0, 2 * pi, length.out = 64)
    sig <- 10 + cos(theta) + 0.3 * sin(4 * theta)

    w_def <- fwaveletspl_3(sig, sjm = 4)
    w_det <- fwaveletspl_3(sig, sjm = 4, det = TRUE)
    w_app <- fwaveletspl_3(sig, sjm = 4, det = FALSE)

    expect_equal(w_def, w_det)
    # Detail coefficients are centered around 0; approximation retains signal mean (~10)
    expect_lt(abs(mean(w_det[1, ])), 1e-6)
    expect_equal(mean(w_app[1, ]), 10, tolerance = 0.2)
    # Reconstruction check: sig = w_app[4, ] + colSums(w_det[1:4, ])
    expect_equal(as.numeric(w_app[4, ] + colSums(w_det[1:4, ])), sig, tolerance = 1e-10)
})

test_that("select_starting_point and select_scale_bar work end-to-end with simulated mouse clicks in manual mode", {
    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    pdf(NULL)
    on.exit({
        options(aforoR.SimLocator = NULL)
        try(dev.off(), silent = TRUE)
    }, add = TRUE)

    # 1. Test select_starting_point clicking near the leftmost tip (x = 50, y = 300)
    options(aforoR.SimLocator = function(n = 1) list(x = 50, y = 300))
    pt_snapped <- select_starting_point(example_image)
    expect_true(is.numeric(pt_snapped))
    expect_length(pt_snapped, 2)
    expect_named(pt_snapped, c("x", "y"))
    # Snapped point must be on the left half of the otolith, not at x = 50
    expect_gt(pt_snapped["x"], 50)

    # 2. Test select_scale_bar clicking 2 points 200 px apart (no snapping)
    options(aforoR.SimLocator = function(n = 2) list(x = c(100, 300), y = c(500, 500)))
    px_mm <- select_scale_bar(example_image, scale_mm = 2)
    expect_equal(as.numeric(px_mm), 100, tolerance = 1e-6)
    sc_info <- attr(px_mm, "scale_info")
    expect_equal(sc_info$bar_length_px, 200, tolerance = 1e-6)
    expect_equal(sc_info$x1, 100L)
    expect_equal(sc_info$x2, 300L)

    # 3. Test process_images with detect_scale = "manual" and start_points = "manual"
    tmp_man <- file.path(tempdir(), "test_manual_pipeline")
    if (dir.exists(tmp_man)) unlink(tmp_man, recursive = TRUE)
    dir.create(tmp_man)
    on.exit(unlink(tmp_man, recursive = TRUE), add = TRUE)
    file.copy(example_image, file.path(tmp_man, "oto_manual.jpg"))

    options(aforoR.SimLocator = function(n = 1) {
        if (n == 2) list(x = c(100, 250), y = c(550, 550)) else list(x = 80, y = 300)
    })
    suppressMessages(
        process_images(tmp_man, detect_scale = "manual", scale_mm = 1.5,
                       start_points = "manual", wavelets = TRUE, ef = FALSE,
                       testing = TRUE, save = TRUE)
    )
    morpho_df <- read.table(file.path(tmp_man, "Polar", "MorphometricsEN.csv"), header = TRUE, sep = ";")
    expect_equal(morpho_df$Units[1], "mm")
    expect_true(file.exists(file.path(tmp_man, "Polar", "zonal_oto_manual.jpg")))
    options(aforoR.SimLocator = NULL)
})

test_that("refactored internal helpers (.resolve_contour_and_start, .render_zonal_and_wavelet_plots, .run_batch_tasks) handle edge cases", {
    # 1. .resolve_contour_and_start with remove_duplicate_end = TRUE vs FALSE
    closed_cnt <- rbind(
        c(0, 0),
        c(10, 0),
        c(5, 5),
        c(0, 0)
    )
    res_keep <- aforoR:::.resolve_contour_and_start(closed_cnt, start_point = "auto", remove_duplicate_end = FALSE)
    res_drop <- aforoR:::.resolve_contour_and_start(closed_cnt, start_point = "auto", remove_duplicate_end = TRUE)
    expect_equal(nrow(res_keep$coords), 4L)
    expect_equal(nrow(res_drop$coords), 3L)
    expect_equal(unname(res_drop$start_point_coords), c(10, 0))

    # Closing vertex index (nrow(closed_cnt) == 4) maps cleanly to vertex 1 when remove_duplicate_end = TRUE
    res_end_idx <- aforoR:::.resolve_contour_and_start(closed_cnt, start_point = 4L, remove_duplicate_end = TRUE)
    expect_equal(res_end_idx$start_idx, 1L)
    expect_equal(unname(res_end_idx$start_point_coords), c(0, 0))

    # Non-numeric 3-row closed data.frame reports numeric validation error even when remove_duplicate_end = TRUE
    char_df <- data.frame(X = c("a", "b", "a"), Y = c("c", "d", "c"), stringsAsFactors = FALSE)
    expect_error(
        aforoR:::.resolve_contour_and_start(char_df, remove_duplicate_end = TRUE),
        "contour coordinates must be numeric"
    )

    # Reject invalid start_point length, empty matrix, 3-col matrix, and non-integer index
    expect_error(
        aforoR:::.resolve_contour_and_start(closed_cnt, start_point = c(1, 2, 3)),
        "start_point must be an index \\(length 1\\) or coordinates"
    )
    expect_error(
        aforoR:::.resolve_contour_and_start(closed_cnt, start_point = matrix(numeric(0), 0, 2)),
        "start_point matrix or data.frame must have exactly 1 row"
    )
    expect_error(
        aforoR:::.resolve_contour_and_start(closed_cnt, start_point = matrix(1:3, 1, 3)),
        "start_point must be an index \\(length 1\\) or coordinates"
    )
    expect_error(
        aforoR:::.resolve_contour_and_start(closed_cnt, start_point = 2.5),
        "start_point index must be an integer"
    )

    # 2. save_visualization / save_visualization_perimeter with file path string + scale_info
    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    if (file.exists(example_image)) {
        theta <- seq(0, 2 * pi, length.out = 65)[-65]
        cnt <- cbind(X = 300 + 100 * cos(theta), Y = 300 + 60 * sin(theta))
        dists <- calculate_distances(cnt, n_points = 64)
        wavs <- calculate_wavelets_analysis(dists, n_scales = 5)
        sc_info <- list(
            min_x = 50, max_x = 150, min_y = 50, max_y = 55,
            x1 = 50L, y1 = 52L, x2 = 150L, y2 = 52L,
            endpoints = list(x1 = 50L, y1 = 52L, x2 = 150L, y2 = 52L),
            bar_length_px = 100, scale_mm = 1, px_per_mm = 100
        )
        tmp_vis <- file.path(tempdir(), "test_vis_filepath_scale")
        on.exit(unlink(tmp_vis, recursive = TRUE), add = TRUE)
        expect_silent(
            save_visualization(
                image = example_image, distances = dists, wavelets = wavs,
                image_name = "oto_path.jpg", output_dir = tmp_vis, scale_info = sc_info
            )
        )
        expect_true(file.exists(file.path(tmp_vis, "zonal_oto_path.jpg")))
        expect_silent(
            save_visualization_perimeter(
                image = example_image, distances = dists, wavelets = wavs,
                image_name = "oto_perim_path.jpg", output_dir = tmp_vis, scale_info = sc_info
            )
        )
        expect_true(file.exists(file.path(tmp_vis, "zonal_oto_perim_path.jpg")))
        expect_error(
            save_visualization(
                image = NULL, distances = dists, wavelets = wavs,
                image_name = "oto_inf.jpg", output_dir = tmp_vis, wavelet_scale = Inf
            ),
            "wavelet_scale must be a positive integer"
        )
    }

    # 3. process_images rejects NA logical flags and Inf threshold_factor
    expect_error(
        process_images(tempdir(), wavelets = NA),
        "must be logical values"
    )
    expect_error(
        process_images(tempdir(), threshold_factor = Inf),
        "threshold_factor must be a positive number"
    )

    # 4. .run_batch_tasks sequential and parallel execution
    items <- list(list(x = 2), list(x = 5))
    seq_res <- aforoR:::.run_batch_tasks(items, function(it) it$x * 10, use_parallel = FALSE)
    expect_equal(seq_res, list(20, 50))

    par_res <- aforoR:::.run_batch_tasks(items, function(it) it$x * 10, use_parallel = TRUE, actual_cores = 2L)
    expect_equal(par_res, list(20, 50))
})
