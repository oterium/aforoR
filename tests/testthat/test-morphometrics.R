test_that("calculate_morphometrics works correctly", {

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
    expect_equal(metrics$Length, sqrt(200), tolerance = 0.01) # Measures max dimension (diagonal for square)
    expect_gt(metrics$Length, 10)
    expect_equal(metrics$Feret_Max, sqrt(200), tolerance = 0.01) # Exact max Feret is diagonal
    expect_equal(metrics$Feret_Min, 10, tolerance = 0.02)       # Min Feret is the side length
    expect_equal(metrics$PCA_Angle, 45, tolerance = 0.01)       # Diagonal symmetry rotates 45 deg

    # Test with scale: 10 pixels = 1 mm
    # Area should be 100 / 100 = 1 mm2
    # Perimeter should be 40 / 10 = 4 mm
    metrics_mm <- calculate_morphometrics(contour, pixels_per_mm = 10)

    expect_equal(metrics_mm$Units, "mm")
    expect_equal(metrics_mm$Area, 1, tolerance = 0.01)
    expect_equal(metrics_mm$Perimeter, 4, tolerance = 0.01)
    expect_equal(metrics_mm$Length, sqrt(200) / 10, tolerance = 0.01)
    expect_equal(metrics_mm$Feret_Max, sqrt(200) / 10, tolerance = 0.01)
    expect_equal(metrics_mm$Feret_Min, 1, tolerance = 0.02)
    expect_equal(metrics_mm$PCA_Angle, 45, tolerance = 0.01) # PCA angle is scale-invariant
})

test_that("process_images integrates morphometrics", {
    skip_if_not_installed("EBImage")

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
    expect_equal(data$Units[1], "mm")

    unlink(temp_dir, recursive = TRUE)
})

test_that("calculate_rol and calculate_ros work as expected", {
    # ROL = (OL * 100) / FL
    expect_equal(calculate_rol(5, 50), 10)
    expect_equal(calculate_rol(c(5, 10), c(50, 100)), c(10, 10))

    # ROS = (OA * 1000) / (FL^3)
    # OA = 10, FL = 10 -> (10 * 1000) / 1000 = 10
    expect_equal(calculate_ros(10, 10), 10)

    # Relative indices table
    df <- calculate_relative_indices(otolith_length = 5, otolith_area = 10, fish_length = 50)
    expect_s3_class(df, "data.frame")
    expect_equal(names(df), c("ROL", "ROS"))
    expect_equal(df$ROL, 10)
    expect_equal(df$ROS, (10 * 1000) / (50^3))

    # Error handling
    expect_error(calculate_rol("5", 50))
    expect_error(calculate_ros(10, -5))
})

