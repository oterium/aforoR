test_that("calculate_morphometrics works correctly", {
    skip_if_not_installed("Momocs")

    # Create a simple square contour 10x10
    # Area = 100, Perimeter = 40
    contour <- data.frame(
        X = c(0, 10, 10, 0, 0),
        Y = c(0, 0, 10, 10, 0)
    )

    # Tested without scale
    metrics <- calculate_morphometrics(contour)

    expect_equal(metrics$Units, "px")
    expect_equal(metrics$Area, 100, tolerance = 0.01)
    expect_equal(metrics$Perimeter, 40, tolerance = 0.01)
    expect_equal(metrics$Length, sqrt(200), tolerance = 0.01) # Momocs measures max dimension (diagonal for square)
    expect_gt(metrics$Length, 10)

    # Shape indices for a square:
    # Rectangularity = Area / (L * W)
    # Circularity = P^2 / Area = 1600 / 100 = 16
    expect_equal(metrics$Circularity, 16, tolerance = 0.1)

    # Test with scale: 10 pixels = 1 mm
    # Area should be 100 / 100 = 1 mm2
    # Perimeter should be 40 / 10 = 4 mm
    metrics_mm <- calculate_morphometrics(contour, pixels_per_mm = 10)

    expect_equal(metrics_mm$Units, "mm")
    expect_equal(metrics_mm$Area, 1, tolerance = 0.01)
    expect_equal(metrics_mm$Perimeter, 4, tolerance = 0.01)
})

test_that("process_images integrates morphometrics", {
    skip_if_not_installed("EBImage")
    skip_if_not_installed("Momocs")

    # Setup temp dir
    temp_dir <- file.path(tempdir(), "test_morpho")
    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir)

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") {
        example_image <- file.path("../../inst/extdata", "otolith.jpg")
    }
    file.copy(example_image, file.path(temp_dir, "otolith.jpg"))

    # Run process_images with scale
    suppressWarnings(process_images(temp_dir, wavelets = FALSE, ef = FALSE, save = TRUE, pixels_per_mm = 100))

    # Check for CSV
    results_file <- file.path(temp_dir, "Polar", "MorphometricsEN.csv")
    expect_true(file.exists(results_file))

    # Read CSV and check columns
    data <- read.table(results_file, header = TRUE, sep = ";", dec = ".")
    expect_true("Area" %in% names(data))
    expect_true("Roundness" %in% names(data))
    expect_equal(data$Units[1], "mm")

    unlink(temp_dir, recursive = TRUE)
})
