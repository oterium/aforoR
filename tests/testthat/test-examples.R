library(testthat)
library(aforoR)

test_that("examples from documentation work as expected", {
    skip_if_not_installed("EBImage")
    skip_if_not_installed("Momocs")

    # 1. Preprocess Image Example
    image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")

    # Fallback if the package is not installed but we are running tests locally
    if (image_path == "") {
        image_path <- file.path("../../inst/extdata", "otolith.jpg")
    }

    expect_true(file.exists(image_path))

    processed_img <- preprocess_image(image_path)
    expect_s4_class(processed_img, "Image")

    # 2. Extract Contour Example
    contour <- extract_contour(processed_img)
    expect_true(is.matrix(contour))
    expect_gt(nrow(contour), 200)

    # 3. Calculate Distances Example
    dists <- calculate_distances(contour, n_points = 512)
    expect_type(dists, "list")
    expect_equal(length(dists$polar$normalized), 512)
    expect_equal(length(dists$perimeter$normalized), 512)

    # 4. Calculate Morphometrics Example
    metrics <- calculate_morphometrics(contour, pixels_per_mm = 100)
    expect_equal(metrics$Units, "mm")
    expect_true(all(c("Area", "Roundness", "Circularity") %in% names(metrics)))

    # 5. Wavelet Analysis Example (using distances from step 3)
    wavelets <- calculate_wavelets_analysis(dists, n_scales = 9)
    expect_type(wavelets, "list")
    expect_equal(nrow(wavelets$polar), 9)
    expect_equal(ncol(wavelets$polar), 512)

    # 6. Full process_images Example (Quick Start)
    temp_dir <- file.path(tempdir(), "test_examples_full")
    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir)
    file.copy(image_path, file.path(temp_dir, "otolith.jpg"))

    # process_images uses cat and progress bars, so we capture the output
    suppressWarnings(
        capture.output(
            process_images(
                folder = temp_dir,
                detect_scale = TRUE,
                wavelets = TRUE,
                ef = TRUE,
                save = TRUE,
                testing = FALSE
            )
        )
    )

    # Verify output folders
    expect_true(dir.exists(file.path(temp_dir, "Polar")))
    expect_true(dir.exists(file.path(temp_dir, "Cartesian")))
    expect_true(file.exists(file.path(temp_dir, "Polar", "MorphometricsEN.csv")))

    unlink(temp_dir, recursive = TRUE)
})
