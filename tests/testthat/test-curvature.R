library(testthat)
library(aforoR)

test_that("calculate_css works with a basic shape", {
    # Create a petal-like shape (with known inflection points)
    theta <- seq(0, 2 * pi, length.out = 100)
    # Increase the perturbation to ensure inflection points survive resampling
    x <- cos(theta) + 0.2 * cos(3 * theta)
    y <- sin(theta) + 0.2 * sin(3 * theta)
    contour <- data.frame(X = x, Y = y)

    res <- calculate_css(contour, n_points = 200, sigma_max = 5)

    expect_s3_class(res, "css")
    expect_type(res$css_matrix, "list")
    expect_true(nrow(res$css_matrix) > 0)
    expect_true(all(res$css_matrix$arc_length >= 0 & res$css_matrix$arc_length <= 1))
    expect_true(all(res$css_matrix$sigma >= 0))
})

test_that("calculate_css handles a circle (no inflection points)", {
    theta <- seq(0, 2 * pi, length.out = 100)
    x <- cos(theta)
    y <- sin(theta)
    contour <- data.frame(X = x, Y = y)

    # A perfect circle should have 0 inflection points at scales above noise
    res <- calculate_css(contour, sigma_max = 2)
    # Check only for sigmas > 0.1 to avoid raw numerical jitter
    inflections_above_noise <- subset(res$css_matrix, sigma > 0.1)
    expect_true(is.null(inflections_above_noise) || nrow(inflections_above_noise) == 0)
})

test_that("calculate_css handles invalid inputs", {
    expect_error(calculate_css(NULL))
})

test_that("plot.css works without error", {
    theta <- seq(0, 2 * pi, length.out = 50)
    x <- cos(theta) + 0.2 * cos(4 * theta)
    y <- sin(theta) + 0.2 * sin(4 * theta)
    res <- calculate_css(cbind(x, y), n_points = 50, sigma_max = 2)

    # Prevent Rplots.pdf generation
    pdf(NULL)
    on.exit(try(dev.off(), silent = TRUE), add = TRUE)

    # Just verify it doesn't crash
    expect_error(plot(res), NA)
    expect_null(plot(res))
})

test_that("plot.css accepts user-supplied graphical parameters and overrides defaults cleanly", {
    theta <- seq(0, 2 * pi, length.out = 50)
    x <- cos(theta) + 0.2 * cos(4 * theta)
    y <- sin(theta) + 0.2 * sin(4 * theta)
    contour <- cbind(X = x, Y = y)
    res <- calculate_css(contour, n_points = 50, sigma_max = 2)

    pdf(NULL)
    on.exit(try(dev.off(), silent = TRUE), add = TRUE)

    # 1. Exact dispatch requirement: override main, col, pch, cex
    expect_error(plot(calculate_css(contour), main = "Custom Title", col = "blue", pch = 16, cex = 1.2), NA)

    # 2. Overriding all default formal arguments (which previously caused duplicate argument errors)
    expect_error(plot(res,
        main = "Custom Title",
        xlab = "Normalized Arc Length",
        ylab = "Gaussian Scale Sigma",
        col = "darkgreen",
        pch = 17,
        cex = 1.5,
        xlim = c(0, 0.5),
        ylim = c(0, 3)
    ), NA)

    # 3. Return value contract: returns invisible NULL
    ret <- plot(res, main = "Custom Title", col = "blue")
    expect_null(ret)
})

test_that("plot.css handles empty inflection matrix with custom parameters", {
    pdf(NULL)
    on.exit(try(dev.off(), silent = TRUE), add = TRUE)

    # Empty data.frame css_matrix
    empty_css_df <- structure(list(
        css_matrix = data.frame(arc_length = numeric(0), sigma = numeric(0)),
        n_points = 50,
        contour_original = matrix(0, 10, 2)
    ), class = "css")

    expect_error(plot(empty_css_df, main = "No inflections found", col = "red", xlab = "Normalized s"), NA)
    expect_null(plot(empty_css_df, main = "No inflections found"))

    # NULL css_matrix
    empty_css_null <- structure(list(
        css_matrix = NULL,
        n_points = 50,
        contour_original = matrix(0, 10, 2)
    ), class = "css")

    expect_error(plot(empty_css_null, main = "NULL inflections", ylab = "Scale"), NA)
    expect_null(plot(empty_css_null))
})

test_that("plot.css handles line and type arguments and secondary graphical calls", {
    theta <- seq(0, 2 * pi, length.out = 50)
    x <- cos(theta) + 0.2 * cos(4 * theta)
    y <- sin(theta) + 0.2 * sin(4 * theta)
    res <- calculate_css(cbind(x, y), n_points = 50, sigma_max = 2)

    pdf(NULL)
    on.exit(try(dev.off(), silent = TRUE), add = TRUE)

    # Line-oriented type arguments
    expect_error(plot(res, type = "b", lty = 2, lwd = 1.5, col = "purple"), NA)
    expect_error(plot(res, type = "l", lty = 1, lwd = 2, col = "darkblue"), NA)

    # Secondary graphical calls on the active plot device execute cleanly
    expect_error({
        plot(res, main = "Overlay Test", col = "navy")
        graphics::lines(c(0, 1), c(0.5, 0.5), col = "red", lty = 2, lwd = 2)
        graphics::points(c(0.2, 0.8), c(0.5, 0.5), pch = 4, col = "black")
    }, NA)
})

test_that("plot.css validates class of input object", {
    expect_error(aforoR:::plot.css(NULL), "x must be an object of class 'css'")
    expect_error(aforoR:::plot.css(list(a = 1)), "x must be an object of class 'css'")
    expect_error(aforoR:::plot.css(matrix(1:4, 2, 2)), "x must be an object of class 'css'")
})


test_that("calculate_css aligns starting point (0.0) to rostrum and assigns 4 AFORO zones", {
    theta <- seq(0, 2 * pi, length.out = 120)[-120]
    theta_shifted <- theta + pi * 0.7
    x <- 2 * cos(theta_shifted) + 0.25 * cos(4 * theta_shifted)
    y <- sin(theta_shifted) + 0.25 * sin(4 * theta_shifted)
    contour <- cbind(X = x, Y = y)

    res_aligned <- calculate_css(contour, n_points = 128, sigma_max = 4, start_point = "auto")
    max_x_idx <- which.max(contour[, 1])
    expect_equal(unname(res_aligned$start_point_coords["X"]), unname(contour[max_x_idx, 1]), tolerance = 1e-6)
    expect_true("zone" %in% colnames(res_aligned$css_matrix))
    expect_true(all(res_aligned$css_matrix$zone %in% c("Zone 1", "Zone 2", "Zone 3", "Zone 4")))

    expect_error(calculate_css(contour, n_points = 128, sigma_max = 4, start_point = "none"), "start_point must be 'auto'")

    pdf(NULL)
    options(aforoR.SimLocator = function(n = 1) list(x = contour[10, 1], y = contour[10, 2]))
    res_manual <- calculate_css(contour, n_points = 128, sigma_max = 4, start_point = "manual")
    options(aforoR.SimLocator = NULL)
    try(dev.off(), silent = TRUE)
    expect_equal(unname(res_manual$start_point_coords), unname(contour[10, ]), tolerance = 1e-6)

    res_custom <- calculate_css(contour, n_points = 128, sigma_max = 4, start_point = 10)
    expect_equal(unname(res_custom$start_point_coords), unname(contour[10, ]), tolerance = 1e-6)

    res_mat <- calculate_css(contour, n_points = 128, sigma_max = 4, start_point = matrix(contour[15, ], nrow = 1))
    expect_equal(unname(res_mat$start_point_coords), unname(contour[15, ]), tolerance = 1e-6)

    pdf(NULL)
    on.exit(try(dev.off(), silent = TRUE), add = TRUE)
    expect_error(plot(res_aligned, show_zones = TRUE), NA)
    expect_error(plot(res_aligned, show_zones = FALSE), NA)
})

test_that("calculate_css defaults to n_points = 512 and supports n_points = NULL / 'max'", {
    theta <- seq(0, 2 * pi, length.out = 301)[-301]
    x <- 2 * cos(theta) + 0.25 * cos(4 * theta)
    y <- sin(theta) + 0.25 * sin(4 * theta)
    contour <- cbind(X = x, Y = y)

    res_def <- calculate_css(contour, sigma_max = 2)
    expect_equal(res_def$n_points, 512L)
    expect_equal(nrow(res_def$contour_original), 512L)

    res_null <- calculate_css(contour, n_points = NULL, sigma_max = 2)
    expect_equal(res_null$n_points, 300L)
    expect_equal(nrow(res_null$contour_original), 300L)

    res_max <- calculate_css(contour, n_points = "max", sigma_max = 2)
    expect_equal(res_max$n_points, 300L)
    expect_equal(nrow(res_max$contour_original), 300L)
})
