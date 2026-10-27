test_that("scale detection works on example image", {
    skip_if_not_installed("EBImage")

    temp_dir <- file.path(tempdir(), "test_scale")
    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir)

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") {
        example_image <- file.path("../../inst/extdata", "otolith.jpg")
    }
    file.copy(example_image, file.path(temp_dir, "otolith.jpg"))

    suppressWarnings(process_images(temp_dir, wavelets = FALSE, ef = FALSE, save = TRUE, detect_scale = TRUE))

    results_file <- file.path(temp_dir, "Polar", "MorphometricsEN.csv")
    expect_true(file.exists(results_file))

    data <- read.table(results_file, header = TRUE, sep = ";", dec = ".")
    expect_equal(data$Units[1], "mm")
    expect_lt(data$Area[1], 1000)

    # Check diagnostic visualization output
    zonal_file <- file.path(temp_dir, "Polar", "zonal_otolith.jpg")
    expect_true(file.exists(zonal_file))
})


test_that("find_scale_bar_length accurately detects horizontal scale bar within <= 2 px (tolerance <= 0.02)", {
    skip_if_not_installed("EBImage")

    # Create synthetic 1200x2000 binary image with otolith and 200x10 rectangular bar
    m <- matrix(0, nrow = 1200, ncol = 2000)
    xx <- col(m) - 1000
    yy <- row(m) - 600
    m[(xx^2 / 400^2 + yy^2 / 200^2) <= 1] <- 1

    # Horizontal scale bar: x from 1500 to 1699 (200 px), y from 1100 to 1109 (10 px)
    m[1100:1109, 1500:1699] <- 1

    px_per_mm <- aforoR:::find_scale_bar_length(m, scale_mm = 2)
    expect_false(is.null(px_per_mm))

    scale_info <- attr(px_per_mm, "scale_info")
    expect_false(is.null(scale_info))

    # True pixel span is 200 px, eliminating +15.5% overestimation (~231 px)
    expect_equal(scale_info$bar_length_px, 200, tolerance = 0.01)
    expect_equal(as.numeric(px_per_mm), 100, tolerance = 0.01)
    expect_equal(scale_info$px_per_mm, 100, tolerance = 0.01)
    expect_equal(scale_info$pixels_per_mm, 100, tolerance = 0.01)
    expect_equal(scale_info$scale_mm, 2)

    # Verify endpoints coinciding with tips
    expect_false(is.null(scale_info$endpoints))
    expect_equal(scale_info$endpoints$x1, 1500)
    expect_equal(scale_info$endpoints$x2, 1699)
    expect_true(scale_info$endpoints$y1 >= 1100 && scale_info$endpoints$y1 <= 1109)
    expect_true(scale_info$endpoints$y2 >= 1100 && scale_info$endpoints$y2 <= 1109)
})


test_that("find_scale_bar_length accurately detects vertical scale bar with correct endpoints", {
    skip_if_not_installed("EBImage")

    # Create synthetic image with specimen and vertical scale bar (e.g. ShrarifJemaa)
    m <- matrix(0, nrow = 1200, ncol = 1600)
    xx <- col(m) - 800
    yy <- row(m) - 600
    m[(xx^2 / 350^2 + yy^2 / 250^2) <= 1] <- 1

    # Vertical scale bar: x from 100 to 114 (15 px), y from 400 to 671 (272 px)
    m[400:671, 100:114] <- 1

    px_per_mm <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
    expect_false(is.null(px_per_mm))

    scale_info <- attr(px_per_mm, "scale_info")
    expect_false(is.null(scale_info))

    # True vertical span is 272 px
    expect_equal(scale_info$bar_length_px, 272, tolerance = 0.01)
    expect_equal(as.numeric(px_per_mm), 272, tolerance = 0.01)

    # Endpoints top and bottom tips
    expect_false(is.null(scale_info$endpoints))
    expect_equal(scale_info$endpoints$y1, 400)
    expect_equal(scale_info$endpoints$y2, 671)
    expect_true(scale_info$endpoints$x1 >= 100 && scale_info$endpoints$x1 <= 114)
    expect_true(scale_info$endpoints$x2 >= 100 && scale_info$endpoints$x2 <= 114)
})


test_that("find_scale_bar_length detects chromatic red scale bar on dark background without discarding isolated bar", {
    skip_if_not_installed("EBImage")

    # Create synthetic 3-channel RGB image (Canaries_JosePajuelo pattern)
    # Background: black
    # Specimen: white otolith (R=0.9, G=0.9, B=0.9)
    # Scale bar: pure red (R=1.0, G=0.0, B=0.0)
    w <- 1200
    h <- 800
    arr <- array(0, dim = c(w, h, 3))

    # Otolith ellipse in center
    xx <- row(arr[, , 1]) - 600
    yy <- col(arr[, , 1]) - 400
    otolith_mask <- (xx^2 / 300^2 + yy^2 / 180^2) <= 1
    arr[, , 1][otolith_mask] <- 0.95
    arr[, , 2][otolith_mask] <- 0.95
    arr[, , 3][otolith_mask] <- 0.95

    # Red scale bar: x from 800 to 999 (200 px), y from 700 to 707 (8 px)
    # Pure saturated red
    arr[800:999, 700:707, 1] <- 1.0
    arr[800:999, 700:707, 2] <- 0.05
    arr[800:999, 700:707, 3] <- 0.05

    img_rgb <- EBImage::Image(arr, colormode = "Color")

    # Detect scale bar on color image
    px_per_mm <- aforoR:::find_scale_bar_length(img_rgb, scale_mm = 1)
    expect_false(is.null(px_per_mm))

    scale_info <- attr(px_per_mm, "scale_info")
    expect_false(is.null(scale_info))
    expect_equal(scale_info$bar_length_px, 200, tolerance = 0.01)
    expect_equal(as.numeric(px_per_mm), 200, tolerance = 0.01)

    expect_false(is.null(scale_info$endpoints))
    expect_equal(scale_info$endpoints$x1, 800)
    expect_equal(scale_info$endpoints$x2, 999)
})


test_that("find_scale_bar_length accurately detects vertical chromatic red scale bar without coordinate transposition", {
    skip_if_not_installed("EBImage")

    w <- 1200
    h <- 800
    arr <- array(0, dim = c(w, h, 3))

    xx <- row(arr[, , 1]) - 600
    yy <- col(arr[, , 1]) - 400
    arr[, , 1][(xx^2 / 300^2 + yy^2 / 180^2) <= 1] <- 0.9
    arr[, , 2][(xx^2 / 300^2 + yy^2 / 180^2) <= 1] <- 0.9
    arr[, , 3][(xx^2 / 300^2 + yy^2 / 180^2) <= 1] <- 0.9

    # Vertical red scale bar at X in [100, 107] (8 px), Y in [200, 399] (200 px)
    arr[100:107, 200:399, 1] <- 1.0
    arr[100:107, 200:399, 2] <- 0.05
    arr[100:107, 200:399, 3] <- 0.05

    img_rgb <- EBImage::Image(arr, colormode = "Color")
    px_per_mm <- aforoR:::find_scale_bar_length(img_rgb, scale_mm = 1)
    expect_false(is.null(px_per_mm))

    scale_info <- attr(px_per_mm, "scale_info")
    expect_false(is.null(scale_info))
    expect_equal(scale_info$bar_length_px, 200, tolerance = 0.01)

    # Coordinates must NOT be transposed
    expect_equal(scale_info$min_x, 100)
    expect_equal(scale_info$max_x, 107)
    expect_equal(scale_info$min_y, 200)
    expect_equal(scale_info$max_y, 399)

    expect_false(is.null(scale_info$endpoints))
    expect_true(scale_info$endpoints$x1 >= 100 && scale_info$endpoints$x1 <= 107)
    expect_true(scale_info$endpoints$x2 >= 100 && scale_info$endpoints$x2 <= 107)
    expect_equal(scale_info$endpoints$y1, 200)
    expect_equal(scale_info$endpoints$y2, 399)
})



test_that("find_scale_bar_length supports plaque-style scale bars (caballa/PHOTOS)", {
    skip_if_not_installed("EBImage")

    # Create synthetic image with otolith and plaque-style rectangular bar (aspect ratio ~3.5)
    m <- matrix(0, nrow = 800, ncol = 1000)
    xx <- col(m) - 500
    yy <- row(m) - 400
    m[(xx^2 / 200^2 + yy^2 / 150^2) <= 1] <- 1

    # Plaque-style bar: x from 20 to 143 (124 px), y from 720 to 752 (33 px)
    m[720:752, 20:143] <- 1

    px_per_mm <- aforoR:::find_scale_bar_length(m, scale_mm = 1)
    expect_false(is.null(px_per_mm))

    scale_info <- attr(px_per_mm, "scale_info")
    expect_false(is.null(scale_info))

    # True physical coordinate span is 124 px
    expect_equal(scale_info$bar_length_px, 124, tolerance = 0.01)
    expect_equal(scale_info$endpoints$x1, 20)
    expect_equal(scale_info$endpoints$x2, 143)
})


test_that("save_visualization and save_visualization_perimeter render endpoint markers without error", {
    skip_if_not_installed("EBImage")

    temp_dir <- file.path(tempdir(), "test_viz_endpoints")
    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir)

    # Mock polar and perimeter distance structures
    theta <- seq(0, 2 * pi, length.out = 128)
    coords <- cbind(200 + 100 * cos(theta), 200 + 80 * sin(theta))
    radii <- sqrt((coords[, 1] - 200)^2 + (coords[, 2] - 200)^2)

    polar_dist <- list(
        polar = list(
            radii = radii,
            coord = coords,
            centroid = c(200, 200),
            rostrum = c(coords[1, 1], coords[1, 2])
        ),
        reordered_coords = list(x = coords[, 1], y = coords[, 2])
    )

    perimeter_dist <- list(
        perimeter = list(
            dist = radii,
            coords = coords,
            centroid = c(200, 200),
            rostrum = c(coords[1, 1], coords[1, 2])
        )
    )

    scale_info <- list(
        min_x = 50, max_x = 250, min_y = 350, max_y = 360,
        x1 = 50, y1 = 355, x2 = 250, y2 = 355,
        endpoints = list(x1 = 50, y1 = 355, x2 = 250, y2 = 355),
        bar_length_px = 201, scale_mm = 1, px_per_mm = 201,
        pixels_per_mm = 201
    )

    # Test save_visualization
    expect_silent(
        save_visualization(
            image = NULL, distances = polar_dist, wavelets = NULL,
            image_name = "sample_test", output_dir = temp_dir,
            scale_info = scale_info
        )
    )
    expect_true(file.exists(file.path(temp_dir, "zonal_sample_test.jpg")))

    # Test save_visualization_perimeter
    expect_silent(
        save_visualization_perimeter(
            image = NULL, distances = perimeter_dist, wavelets = NULL,
            image_name = "sample_perim", output_dir = temp_dir,
            scale_info = scale_info
        )
    )
    expect_true(file.exists(file.path(temp_dir, "zonal_sample_perim.jpg")))
})
