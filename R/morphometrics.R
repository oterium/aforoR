#' Calculate Basic Morphometric Measurements
#'
#' Computes standard otolith morphometric measurements (Area, Perimeter, Length, Width)
#' along with advanced shape descriptors (Feret_Max, Feret_Min, PCA_Angle).
#' The function supports scale conversion from pixels to millimeters if a `pixels_per_mm`
#' value is provided.
#'
#' @param contour A matrix or data.frame with 'X' and 'Y' columns representing the contour.
#' @param pixels_per_mm A numeric value specifying the number of pixels per millimeter.
#'                      If NULL (default), measurements are returned in pixels.
#' @return A named list containing:
#'   \itemize{
#'     \item \code{Area}: Surface area of the otolith.
#'     \item \code{Perimeter}: Length of the otolith boundary.
#'     \item \code{Length}: Maximum dimension (major axis) of the PCA-aligned bounding box.
#'     \item \code{Width}: Minimum dimension (minor axis) of the PCA-aligned bounding box.
#'     \item \code{Feret_Max}: Maximum Feret diameter (longest distance between any two boundary points).
#'     \item \code{Feret_Min}: Minimum Feret diameter (minimum distance between parallel tangency lines).
#'     \item \code{PCA_Angle}: Orientation angle of the major axis of inertia in degrees (normalized to [0, 180)).
#'     \item \code{Units}: Measurement units ("px" or "mm").
#'   }
#' @export
#' @examples
#' # Example using the built-in Aphanopus dataset
#' data(Aphanopus)
#' # Note: For internal calculation, you would typically pass the raw coordinates.
#' # Here we illustrate the function with sample coordinates from the package:
#' image_path <- system.file("extdata", "otolith.jpg", package = "aforoR")
#' binary_img <- preprocess_image(image_path)
#' contour <- extract_contour(binary_img)
#' metrics <- calculate_morphometrics(contour, pixels_per_mm = 100)
#' print(metrics)
calculate_morphometrics <- function(contour, pixels_per_mm = NULL) {
    # Input validation
    if (is.null(contour)) {
        stop("contour cannot be NULL")
    }

    if (is.matrix(contour)) {
        contour <- as.data.frame(contour)
        if (ncol(contour) == 2 && (is.null(colnames(contour)) || !all(c("X", "Y") %in% colnames(contour)))) {
            colnames(contour) <- c("X", "Y")
        }
    }

    if (!is.data.frame(contour) || !all(c("X", "Y") %in% names(contour))) {
        stop("contour must be a data.frame or matrix with X and Y columns")
    }

    if (!is.null(pixels_per_mm) && (!is.numeric(pixels_per_mm) || pixels_per_mm <= 0)) {
        stop("pixels_per_mm must be a positive numeric value")
    }

    # Convert contour for Momocs
    coo <- as.matrix(contour[, c("X", "Y")])

    # Calculate basic measurements in pixels
    # Area
    area_px <- Momocs::coo_area(coo)

    # Perimeter
    perimeter_px <- Momocs::coo_perim(coo)

    # Length (Major Axis) and Width (Minor Axis)
    # Momocs::coo_length gives the maximum dimension (Length)
    # Momocs::coo_width gives the minimum dimension (Width)
    length_px <- Momocs::coo_length(coo)
    width_px <- Momocs::coo_width(coo)

    # Feret Max (Caliper Length / Maximum Feret Diameter)
    feret_max_px <- Momocs::coo_calliper(coo)

    # Feret Min (Minimum Feret Diameter) via angular sweep (360 angles)
    angles <- seq(0, pi, length.out = 360)
    widths <- sapply(angles, function(a) {
        rot_mat <- matrix(c(cos(a), -sin(a), sin(a), cos(a)), 2, 2)
        rot_coo <- coo %*% rot_mat
        diff(range(rot_coo[, 2]))
    })
    feret_min_px <- min(widths)

    # PCA Alignment / Rotation Angle (in degrees, scale-invariant, normalized to [0, 180))
    V <- stats::var(coo)
    s <- svd(V)
    pca_angle <- (atan2(s$u[2, 1], s$u[1, 1]) * 180 / pi) %% 180

    # Apply scale if provided
    if (!is.null(pixels_per_mm)) {
        area <- area_px / (pixels_per_mm^2)
        perimeter <- perimeter_px / pixels_per_mm
        length_val <- length_px / pixels_per_mm
        width_val <- width_px / pixels_per_mm
        feret_max_val <- feret_max_px / pixels_per_mm
        feret_min_val <- feret_min_px / pixels_per_mm
        units <- "mm"
    } else {
        area <- area_px
        perimeter <- perimeter_px
        length_val <- length_px
        width_val <- width_px
        feret_max_val <- feret_max_px
        feret_min_val <- feret_min_px
        units <- "px"
    }

    # Return results
    result <- list(
        # Measurements
        Area = area,
        Perimeter = perimeter,
        Length = length_val,
        Width = width_val,
        Feret_Max = feret_max_val,
        Feret_Min = feret_min_val,
        PCA_Angle = pca_angle,
        Units = units
    )

    return(result)
}
