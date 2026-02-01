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

    # Just verify it doesn't crash
    expect_error(plot(res), NA)
})
