library(testthat)
library(aforoR)

test_that("select_points_hall basic functionality works", {
    # Create dummy functional-like data for classification
    set.seed(42)
    n <- 40
    p <- 100
    groups <- factor(rep(c("A", "B"), each = 20))

    # Group A has higher values at indices 10 and 50
    # Group B has higher values at indices 20 and 80
    data <- matrix(rnorm(n * p), n, p)
    data[1:20, 10] <- data[1:20, 10] + 5
    data[1:20, 50] <- data[1:20, 50] + 5
    data[21:40, 20] <- data[21:40, 20] + 5
    data[21:40, 80] <- data[21:40, 80] + 5

    # Run selection
    # Note: MASS must be installed for this test to pass if using LDA
    if (requireNamespace("MASS", quietly = TRUE)) {
        res <- select_points_hall(data, groups, method = "lda", cv = FALSE, p = 0.05)

        expect_type(res, "list")
        expect_true(length(res$selected_indices) >= 1)
        expect_true(all(res$selected_indices %in% 1:p))
        expect_true(res$error_path[length(res$error_path)] <= res$error_path[1])
    } else {
        skip("MASS package not available")
    }
})

test_that("select_points_hall handles invalid inputs", {
    data <- matrix(rnorm(100), 10, 10)
    groups <- rep(c("A", "B"), each = 5)

    expect_error(select_points_hall(data, groups, method = "invalid_method"))
    # If we try a method whose package is missing (not likely in dev environment but good to check)
})

test_that("select_points_hall stop criterion works", {
    set.seed(42)
    data <- matrix(rnorm(200), 20, 10)
    groups <- factor(rep(c("A", "B"), each = 10))

    # If data is pure noise, it should stop early or improvements should be small
    if (requireNamespace("MASS", quietly = TRUE)) {
        res <- select_points_hall(data, groups, method = "lda", p = 0.5)
        # High p should mean it stops very quickly
        expect_true(length(res$selected_indices) <= 2)
    }
})
