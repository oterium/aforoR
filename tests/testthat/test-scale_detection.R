test_that("scale detection works on example image", {
    skip_if_not_installed("EBImage")

    # Load internal helper if possible, or source it.
    # Since find_scale_bar_length is internal (noRd), we might need to source process_images.R again or test via process_images.

    # Setup temp dir
    temp_dir <- file.path(tempdir(), "test_scale")
    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir)

    example_image <- system.file("extdata", "otolith.jpg", package = "aforoR")
    if (example_image == "") {
        example_image <- file.path("../../inst/extdata", "otolith.jpg")
    }
    file.copy(example_image, file.path(temp_dir, "otolith.jpg"))

    # Run process_images with detection
    # We expect the morphometrics to be in mm.
    # Based on manual inspection, the bar is ~56 pixels.
    # So Area in pixels (approx 98000) -> Area in mm^2 = 98000 / (56^2) = 98000 / 3136 ~ 31 mm^2.

    suppressWarnings(process_images(temp_dir, wavelets = FALSE, ef = FALSE, save = TRUE, detect_scale = TRUE))

    results_file <- file.path(temp_dir, "Procesamiento", "MorphometricsEN.csv")
    expect_true(file.exists(results_file))

    data <- read.table(results_file, header = TRUE, sep = ";", dec = ".")

    # Check units
    expect_equal(data$Units[1], "mm")

    # Check Area value roughly
    # If detection failed, pixels_per_mm would be NULL -> Units = px -> Area ~ 98000
    # If detection worked, Area should be much smaller.
    expect_lt(data$Area[1], 1000) # 31 is much less than 1000
})
