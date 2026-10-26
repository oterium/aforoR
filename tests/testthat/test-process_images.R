test_that("process_images works with example image", {
    # Skip if dependencies are not installed
    skip_if_not_installed("EBImage")

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
    # We suppress warnings because Contour extraction might produce some warnings about contours
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

test_that("preprocess_image handles missing files and validates blur params", {
    expect_error(preprocess_image("non_existent_file.jpg"))
    skip_if_not_installed("EBImage")
    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    # Invalid blur parameter
    expect_error(preprocess_image(example_image, blur_size = "invalid"))
    expect_error(preprocess_image(example_image, blur_sigma = "invalid"))

    # Auto blur mode
    bin_auto <- preprocess_image(example_image, blur_size = "auto", blur_sigma = "auto")
    expect_s4_class(bin_auto, "Image")

    # Explicit numeric with even blur_size (should adjust to odd automatically)
    bin_num <- preprocess_image(example_image, blur_size = 20, blur_sigma = 5)
    expect_s4_class(bin_num, "Image")
})


test_that("save_visualization renders lightweight JPEG files over real image and blank canvas", {
    skip_if_not_installed("EBImage")

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    img <- EBImage::readImage(example_image)
    bin <- preprocess_image(img)
    cnt <- extract_contour(bin)
    dists <- calculate_distances(cnt, n_points = 512)
    wavs <- calculate_wavelets_analysis(dists, n_scales = 5)

    temp_out <- file.path(tempdir(), "test_viz_dir")
    if (dir.exists(temp_out)) unlink(temp_out, recursive = TRUE)
    dir.create(temp_out)

    # Test with Image object
    save_visualization(image = img, distances = dists, wavelets = wavs,
                       image_name = "test_img.jpg", output_dir = temp_out,
                       plot_width = 800, plot_height = 600, plot_res = 100)

    zonal_file <- file.path(temp_out, "zonal_test_img.jpg")
    wavelet_file <- file.path(temp_out, "wavelet_test_img.jpg")

    expect_true(file.exists(zonal_file))
    expect_true(file.exists(wavelet_file))
    # File size should be lightweight (< 300 KB, not 2.5 MB)
    expect_lt(file.info(zonal_file)$size, 300 * 1024)

    # Test perimeter visualization
    save_visualization_perimeter(image = img, distances = dists, wavelets = wavs,
                                 image_name = "test_img.jpg", output_dir = temp_out,
                                 plot_width = 800, plot_height = 600, plot_res = 100)
    expect_true(file.exists(file.path(temp_out, "zonal_test_img.jpg")))

    # Test with image = NULL (Procrustes mode)
    unlink(temp_out, recursive = TRUE)
    dir.create(temp_out)
    expect_silent(save_visualization(image = NULL, distances = dists, wavelets = wavs,
                                     image_name = "test_null.jpg", output_dir = temp_out,
                                     plot_width = 600, plot_height = 400, plot_res = 72))
    expect_true(file.exists(file.path(temp_out, "zonal_test_null.jpg")))

    unlink(temp_out, recursive = TRUE)
})

test_that("preprocess_image downsample_threshold and extract_contour robustness", {
    skip_if_not_installed("EBImage")

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    # Fast downsample threshold mode
    bin_fast <- preprocess_image(example_image, downsample_threshold = TRUE)
    expect_s4_class(bin_fast, "Image")

    # extract_contour handles degenerate contours without throwing out-of-bounds error
    cnt <- extract_contour(bin_fast)
    expect_true(is.matrix(cnt))
    expect_gte(nrow(cnt), 200)
})

test_that("process_images works with testing=TRUE, subfolder=TRUE, and parallel cores", {
    skip_if_not_installed("EBImage")

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") example_image <- file.path("../../inst/extdata", "otolith.jpg")
    skip_if_not(file.exists(example_image), "Example image not found")

    temp_root <- file.path(tempdir(), "test_subfolders_root")
    if (dir.exists(temp_root)) unlink(temp_root, recursive = TRUE)
    sub1 <- file.path(temp_root, "stockA")
    dir.create(sub1, recursive = TRUE)

    file.copy(example_image, file.path(sub1, "specimen1.jpg"))

    # Test subfolder = TRUE with cores = 2 and testing = TRUE
    suppressMessages(
        process_images(folder = temp_root, subfolder = TRUE, cores = 2,
                       wavelets = TRUE, testing = TRUE, save = TRUE,
                       plot_width = 800, plot_height = 600, plot_res = 100)
    )

    polar_dir <- file.path(sub1, "Polar")
    expect_true(dir.exists(polar_dir))
    expect_true(file.exists(file.path(polar_dir, "Wavelet_1EN.csv")))
    expect_true(file.exists(file.path(polar_dir, "zonal_specimen1.jpg")))
    expect_true(file.exists(file.path(polar_dir, "wavelet_specimen1.jpg")))

    unlink(temp_root, recursive = TRUE)
})

