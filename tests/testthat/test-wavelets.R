test_that("fwaveletspl_3 validates inputs", {
  expect_error(fwaveletspl_3(c(), 5), "cannot be empty")
  expect_error(fwaveletspl_3(1:10, -1), "must be a positive integer")
  expect_error(fwaveletspl_3(1:10, 5, "not logical"), "must be a logical value")
})

test_that("fwaveletspl_3 produces correct output dimensions", {
  ra0 <- rnorm(512)
  sjm <- 9
  result <- fwaveletspl_3(ra0, sjm)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), sjm)
  expect_equal(ncol(result), length(ra0))
})

test_that("fatrous1d validates parameters", {
  expect_error(fatrous1d(c(), c(1, 2), 1, 2, 1), "cannot be empty")
  expect_error(fatrous1d(1:10, 1:5, 1, 2, -1), "sj must be >= 0")
  expect_error(fatrous1d(1:10, 1:5, 3, 2, 1), "snmax must be greater than or equal to snmin")
})
