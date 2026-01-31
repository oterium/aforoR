test_that("process_images works with example image", {
    # Skip if dependencies are not installed
    skip_if_not_installed("EBImage")
    skip_if_not_installed("Momocs")

    # Create a temporary directory for the test
    temp_dir <- file.path(tempdir(), "test_images")
    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir)

    # Copy the example image to the temp directory
    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")

    # Fallback if the package is not installed but we are running tests locally
    if (example_image == "") {
        example_image <- file.path("../../inst/extdata", "otolith.jpg")
    }

    expect_true(file.exists(example_image), info = "Example image should exist")

    file.copy(example_image, file.path(temp_dir, "otolith.jpg"))

    # Run process_images
    # We suppress warnings because Momocs might produce some warnings about contours
    suppressWarnings(process_images(temp_dir, wavelets = TRUE, ef = TRUE, save = TRUE, testing = FALSE))

    # Check if results were created
    # process_images creates a 'Polar' folder
    results_dir <- file.path(temp_dir, "Polar")
    expect_true(dir.exists(results_dir), info = "Polar directory should be created")

    # Check for CSV files
    expect_true(file.exists(file.path(results_dir, "Distancia_NormEN.csv")))
    expect_true(file.exists(file.path(results_dir, "Wavelet_1EN.csv")))
    expect_true(file.exists(file.path(results_dir, "Coords.csv")))
    expect_true(file.exists(file.path(results_dir, "EllipticCoeEN.csv")))

    # Check for perimeter results
    results_dir2 <- file.path(temp_dir, "Cartesian")
    expect_true(dir.exists(results_dir2), info = "Cartesian directory should be created")
    expect_true(file.exists(file.path(results_dir2, "Distancia_NormEN.csv")))

    # Clean up
    unlink(temp_dir, recursive = TRUE)
})

test_that("calculate_distances handles input correctly", {
    # Test with invalid inputs
    expect_error(calculate_distances(NULL))
    expect_error(calculate_distances(matrix(1:10, ncol = 1))) # Missing columns
})

test_that("preprocess_image handles missing files", {
    expect_error(preprocess_image("non_existent_file.jpg"))
})
