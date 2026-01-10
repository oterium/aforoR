test_that("ild calculates Euclidean distance correctly", {
  expect_equal(ild(c(0, 0), c(3, 4)), 5)
  expect_equal(ild(c(1, 1), c(1, 1)), 0)
  expect_error(ild(c(1), c(1, 2)), "must be vectors of length 2")
  expect_error(ild(c(1, NA), c(1, 2)), "cannot contain NA values")
})

test_that("regularradius validates inputs properly", {
  expect_error(regularradius(c(1, 2), c(1, 2, 3), 2), "must have the same length")
  expect_error(regularradius(c(1, 2), c(1, 2), -1), "must be a positive integer")
  expect_error(regularradius(c(1, 2), c(1, 2), 5), "cannot exceed the number")
})

test_that("dper calculates perimeter distances", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  result <- dper(x, y, 3)
  
  expect_type(result, "list")
  expect_named(result, c("dist", "coords"))
  expect_length(result$dist, 3)
  expect_equal(nrow(result$coords), 3)
})
