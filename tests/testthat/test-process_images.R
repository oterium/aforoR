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

test_that("preprocess_image strictly clamps FFT filter2 overshoot to [0, 1]", {
    skip_if_not_installed("EBImage")

    # Synthetic sharp step-edge that produces FFT ringing
    m <- matrix(0, nrow = 120, ncol = 120)
    m[30:90, 30:90] <- 1.0

    # Test raw filter2 produces ringing
    brush <- EBImage::makeBrush(size = 15, shape = "gaussian", sigma = 3)
    raw_filtered <- EBImage::filter2(m, brush)
    # FFT ringing produces overshoot > 1 or undershoot < 0
    expect_true(max(raw_filtered) > 1.0 || min(raw_filtered) < 0.0)

    # preprocess_image clamps post-filter2
    bin_img <- preprocess_image(EBImage::Image(m))
    expect_s4_class(bin_img, "Image")
    expect_true(is.logical(bin_img) || all(as.vector(bin_img) %in% c(0, 1, TRUE, FALSE)))
})

test_that("calculate_adaptive_threshold guards / 255 rescaling and clamps inputs", {
    skip_if_not_installed("EBImage")

    # 1. Microscopic float overshoot (1.0001) must NOT trigger / 255
    mat_overshoot <- matrix(runif(10000, min = 0.1, max = 0.9), nrow = 100, ncol = 100)
    mat_overshoot[1, 1] <- 1.0001
    th1 <- aforoR:::calculate_adaptive_threshold(mat_overshoot)
    expect_true(th1 > 0.05, info = "Float overshoot must not collapse threshold to < 0.05")

    # 2. Intermediate value near 1.0 (1.04) must NOT trigger / 255
    mat_inter <- matrix(runif(10000, min = 0.1, max = 0.8), nrow = 100, ncol = 100)
    mat_inter[1, 1] <- 1.04
    th2 <- aforoR:::calculate_adaptive_threshold(mat_inter)
    expect_true(th2 > 0.05, info = "Value 1.04 must not trigger / 255")

    # 3. 8-bit integer image (e.g. max = 220) MUST trigger / 255
    mat_int <- matrix(round(runif(10000, min = 20, max = 220)), nrow = 100, ncol = 100)
    th_int <- aforoR:::calculate_adaptive_threshold(mat_int)
    expect_true(th_int >= 0.05 && th_int <= 0.95)

    # 4. Otsu input clamping prevents crash on out-of-range breaks
    mat_out <- matrix(seq(-0.1, 1.1, length.out = 10000), nrow = 100, ncol = 100)
    expect_silent(th_out <- aforoR:::calculate_adaptive_threshold(mat_out))
    expect_true(th_out >= 0 && th_out <= 1)
})

test_that("calculate_adaptive_threshold applies safeguard on textured grey backgrounds", {
    skip_if_not_installed("EBImage")

    # Simulate textured grey background (mean ~ 0.18, stdev ~ 0.05) and specimen (mean ~ 0.8)
    set.seed(42)
    bg <- matrix(rnorm(10000, mean = 0.18, sd = 0.04), nrow = 100, ncol = 100)
    bg[bg < 0] <- 0
    # Specimen core
    bg[30:70, 30:70] <- rnorm(1681, mean = 0.85, sd = 0.05)
    bg[bg > 1] <- 1

    th <- aforoR:::calculate_adaptive_threshold(bg)
    th_otsu <- EBImage::otsu(bg)

    # Adaptive threshold must not fall into background noise (< 0.65 * Otsu)
    expect_gte(th, 0.65 * th_otsu)
    expect_gt(th, 0.22) # comfortably above the 0.18 background
})

test_that("read_and_process_image and preprocess_image strip alpha from 2-channel TIFFs", {
    skip_if_not_installed("EBImage")

    # Create synthetic 2-channel Image (Ch1: grayscale otolith, Ch2: alpha 1.0)
    arr <- array(0.02, dim = c(120, 100, 2))
    # Otolith in channel 1
    arr[30:90, 20:80, 1] <- 0.85
    # Alpha mask in channel 2
    arr[, , 2] <- 1.0
    img2ch <- EBImage::Image(arr, colormode = "Color")

    # 1. read_and_process_image strips channel 2 and sets Grayscale
    stripped <- read_and_process_image(img2ch)
    expect_equal(length(dim(stripped)), 2)
    expect_equal(EBImage::colorMode(stripped), 0L)
    # Background in corner remains 0.02, not shifted to 0.34
    expect_equal(as.numeric(stripped[1, 1]), 0.02, tolerance = 1e-5)

    # 2. preprocess_image handles 2-channel Image cleanly
    bin <- preprocess_image(img2ch)
    expect_s4_class(bin, "Image")
    expect_equal(length(dim(bin)), 2)
    # Corner background is FALSE, otolith center is TRUE
    expect_false(as.logical(bin[1, 1]))
    expect_true(as.logical(bin[50, 50]))
})

test_that("extract_contour rejects outer frame border (area >= 85%)", {
    skip_if_not_installed("EBImage")

    total_px <- 600 * 800 # 480,000
    m <- matrix(0, nrow = 600, ncol = 800)

    # 1. Outer border frame: area ~95% of total image
    m[5:595, 5:795] <- 1
    m[25:575, 25:775] <- 0

    # Internal biological specimen: area ~30,000 px
    xx <- col(m) - 400
    yy <- row(m) - 300
    m[(xx^2 / 120^2 + yy^2 / 60^2) <= 1] <- 1

    # Extract contour: must discard the 95% outer frame and pick the internal specimen
    cnt <- extract_contour(m, min_points = 150, min_area = 5000)
    expect_false(is.null(cnt))
    specimen_area <- aforoR:::.coo_area(cnt)
    expect_lt(specimen_area, 0.85 * total_px)
    expect_gt(specimen_area, 5000)

    # 2. Image containing ONLY an outer frame border returns NULL with warning
    border_only <- matrix(0, nrow = 400, ncol = 400)
    border_only[2:398, 2:398] <- 1
    border_only[20:380, 20:380] <- 0

    expect_warning(res_border <- extract_contour(border_only, min_points = 100, min_area = 1000))
    expect_null(res_border)
})

test_that("extract_contour supports small images with min_area = 1000", {
    skip_if_not_installed("EBImage")

    # Small 480x480 image with specimen area ~3,500 px
    m <- matrix(0, nrow = 480, ncol = 480)
    xx <- col(m) - 240
    yy <- row(m) - 240
    rx <- 55
    ry <- 20 # pi * 55 * 20 ~ 3456 px
    m[(xx^2 / rx^2 + yy^2 / ry^2) <= 1] <- 1

    # Fails default min_area = 5000
    expect_warning(res_fail <- extract_contour(m, min_points = 200, min_area = 5000))
    expect_null(res_fail)

    # Succeeds with min_area = 1000
    res_ok <- extract_contour(m, min_points = 200, min_area = 1000)
    expect_false(is.null(res_ok))
    area_ok <- aforoR:::.coo_area(res_ok)
    expect_gte(area_ok, 1000)
    expect_lt(area_ok, 5000)
    expect_true(all(res_ok[, 1] >= 1 & res_ok[, 1] <= 480))
    expect_true(all(res_ok[, 2] >= 1 & res_ok[, 2] <= 480))
})

test_that("preprocess_image defensively validates threshold_factor", {
    skip_if_not_installed("EBImage")

    img <- EBImage::Image(matrix(0.5, 50, 50))
    expect_error(preprocess_image(img, threshold_factor = -1), "threshold_factor must be a positive numeric value")
    expect_error(preprocess_image(img, threshold_factor = 0), "threshold_factor must be a positive numeric value")
    expect_error(preprocess_image(img, threshold_factor = "invalid"), "threshold_factor must be a positive numeric value")
    expect_error(preprocess_image(img, threshold_factor = NA), "threshold_factor must be a positive numeric value")
    expect_error(preprocess_image(img, threshold_factor = NULL), "threshold_factor must be a positive numeric value")
    expect_error(preprocess_image(img, threshold_factor = c(1, 2)), "threshold_factor must be a positive numeric value")

    # Positive numeric threshold_factor works
    expect_silent(res <- preprocess_image(img, threshold_factor = 1.2))
    expect_true(inherits(res, "Image"))
})



