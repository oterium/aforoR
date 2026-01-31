#' Calculate Morphometric Indices and Measurements
#'
#' Computes standard otolith morphometric measurements (Area, Perimeter, Length, Width)
#' and shape indices (Roundness, Circularity, etc.) based on Tuset et al. (2003) formulas.
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
#'     \item \code{Length}: Maximum dimension (major axis).
#'     \item \code{Width}: Minimum dimension (minor axis).
#'     \item \code{Units}: Measurement units ("px" or "mm").
#'     \item \code{Roundness}: (4 * Area) / (pi * Length^2).
#'     \item \code{FormFactor}: (4 * pi * Area) / Perimeter^2.
#'     \item \code{Circularity}: Perimeter^2 / Area.
#'     \item \code{Rectangularity}: Area / (Length * Width).
#'     \item \code{Ellipticity}: (Length - Width) / (Length + Width).
#'     \item \code{AspectRatio}: Length / Width.
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

    # Apply scale if provided
    if (!is.null(pixels_per_mm)) {
        area <- area_px / (pixels_per_mm^2)
        perimeter <- perimeter_px / pixels_per_mm
        length_val <- length_px / pixels_per_mm
        width_val <- width_px / pixels_per_mm
        units <- "mm"
    } else {
        area <- area_px
        perimeter <- perimeter_px
        length_val <- length_px
        width_val <- width_px
        units <- "px"
    }

    # Calculate Shape Indices (Tuset et al.)
    # Formulas are unit-independent (ratios), so we can use either px or mm values

    # Roundness (Redondez): 4*Area / (pi * Length^2)
    # Note: Some definitions use Length (Major Axis) squared.
    roundness <- (4 * area) / (pi * length_val^2)

    # Form Factor: 4 * pi * Area / Perimeter^2
    form_factor <- (4 * pi * area) / (perimeter^2)

    # Circularity (Circularidad): Perimeter^2 / Area
    circularity <- (perimeter^2) / area

    # Rectangularity (Rectangularidad): Area / (Length * Width)
    rectangularity <- area / (length_val * width_val)

    # Ellipticity (Ellipticidad): (Length - Width) / (Length + Width)
    ellipticity <- (length_val - width_val) / (length_val + width_val)

    # Aspect Ratio (Relación de Aspecto): Length / Width
    aspect_ratio <- length_val / width_val

    # Return results
    result <- list(
        # Measurements
        Area = area,
        Perimeter = perimeter,
        Length = length_val,
        Width = width_val,
        Units = units,

        # Indices
        Roundness = roundness,
        FormFactor = form_factor,
        Circularity = circularity,
        Rectangularity = rectangularity,
        Ellipticity = ellipticity,
        AspectRatio = aspect_ratio
    )

    return(result)
}
