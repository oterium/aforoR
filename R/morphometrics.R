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
        if (ncol(contour) >= 2) {
            contour <- as.data.frame(contour[, 1:2, drop = FALSE])
            colnames(contour) <- c("X", "Y")
        } else {
            stop("contour must have at least 2 columns")
        }
    } else if (is.data.frame(contour)) {
        if (all(c("X", "Y") %in% colnames(contour))) {
            contour <- contour[, c("X", "Y")]
        } else if (all(c("x", "y") %in% colnames(contour))) {
            contour <- contour[, c("x", "y")]
            colnames(contour) <- c("X", "Y")
        } else if (ncol(contour) >= 2) {
            contour <- contour[, 1:2]
            colnames(contour) <- c("X", "Y")
        } else {
            stop("contour must have at least 2 columns")
        }
    } else {
        stop("contour must be a data.frame or matrix with X and Y columns")
    }

    if (nrow(contour) < 3) {
        stop("contour must contain at least 3 points")
    }

    if (!is.numeric(contour$X) || !is.numeric(contour$Y)) {
        stop("contour coordinates must be numeric")
    }

    if (any(is.na(contour$X)) || any(is.na(contour$Y))) {
        stop("contour coordinates cannot contain NA values")
    }

    if (any(is.infinite(contour$X)) || any(is.infinite(contour$Y))) {
        stop("contour coordinates cannot contain infinite values")
    }

    if (all(contour$X == contour$X[1]) && all(contour$Y == contour$Y[1])) {
        stop("contour is degenerate with identical coordinates")
    }

    if (!is.null(pixels_per_mm) && (!is.numeric(pixels_per_mm) || length(pixels_per_mm) != 1 || is.na(pixels_per_mm) || is.infinite(pixels_per_mm) || pixels_per_mm <= 0)) {
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


#' Calculate Relative Otolith Length (ROL)
#'
#' Computes the Relative Otolith Length (ROL) index relating otolith length to
#' fish standard or total length: \eqn{ROL = (OL \times 100) / FL}.
#'
#' @param otolith_length Numeric vector of otolith length measurements.
#' @param fish_length Numeric vector of fish length measurements (e.g. standard length or total length).
#' @return A numeric vector containing the ROL values (percentage).
#' @export
#' @references
#' Tuset VM, Olivar MP, Manjabacas A, Lombarte A. (2026). Two-dimensional otolith size trends in myctophid fishes.
#' Journal of Fish Biology, 109(1), 323-336. \doi{10.1111/jfb.70412}.
#' @examples
#' calculate_rol(otolith_length = 5.2, fish_length = 55)
calculate_rol <- function(otolith_length, fish_length) {
    if (!is.numeric(otolith_length) || !is.numeric(fish_length)) {
        stop("Both otolith_length and fish_length must be numeric vectors")
    }
    if (length(otolith_length) == 0 || length(fish_length) == 0) {
        stop("otolith_length and fish_length cannot be empty")
    }
    if (length(otolith_length) != length(fish_length) &&
        length(otolith_length) != 1 && length(fish_length) != 1) {
        stop("Lengths of otolith_length and fish_length must match or one must be of length 1")
    }
    if (any(is.infinite(otolith_length)) || any(is.infinite(fish_length))) {
        stop("Inputs cannot contain infinite values")
    }
    if (any(otolith_length < 0, na.rm = TRUE)) {
        stop("otolith_length cannot be negative")
    }
    if (any(fish_length <= 0, na.rm = TRUE)) {
        stop("fish_length must be strictly positive")
    }
    (otolith_length * 100) / fish_length
}

#' Calculate Relative Otolith Surface / Size (ROS)
#'
#' Computes the Relative Otolith Surface/Size (ROS) index relating otolith surface
#' area to fish body volume proxy: \eqn{ROS = (OA \times 1000) / (FL^3)}.
#'
#' @param otolith_area Numeric vector of otolith area measurements.
#' @param fish_length Numeric vector of fish length measurements (e.g. standard length or total length).
#' @return A numeric vector containing the ROS values.
#' @export
#' @references
#' Tuset VM, Olivar MP, Manjabacas A, Lombarte A. (2026). Two-dimensional otolith size trends in myctophid fishes.
#' Journal of Fish Biology, 109(1), 323-336. \doi{10.1111/jfb.70412}.
#' @examples
#' calculate_ros(otolith_area = 14.5, fish_length = 55)
calculate_ros <- function(otolith_area, fish_length) {
    if (!is.numeric(otolith_area) || !is.numeric(fish_length)) {
        stop("Both otolith_area and fish_length must be numeric vectors")
    }
    if (length(otolith_area) == 0 || length(fish_length) == 0) {
        stop("otolith_area and fish_length cannot be empty")
    }
    if (length(otolith_area) != length(fish_length) &&
        length(otolith_area) != 1 && length(fish_length) != 1) {
        stop("Lengths of otolith_area and fish_length must match or one must be of length 1")
    }
    if (any(is.infinite(otolith_area)) || any(is.infinite(fish_length))) {
        stop("Inputs cannot contain infinite values")
    }
    if (any(otolith_area < 0, na.rm = TRUE)) {
        stop("otolith_area cannot be negative")
    }
    if (any(fish_length <= 0, na.rm = TRUE)) {
        stop("fish_length must be strictly positive")
    }
    (otolith_area * 1000) / (fish_length^3)
}

#' Calculate Relative Otolith Indices (ROL and ROS)
#'
#' Computes both Relative Otolith Length (ROL) and Relative Otolith Surface (ROS)
#' indices in a structured data frame.
#'
#' @param otolith_length Numeric vector of otolith length.
#' @param otolith_area Numeric vector of otolith area.
#' @param fish_length Numeric vector of fish length.
#' @return A data.frame containing columns \code{ROL} and \code{ROS}.
#' @export
#' @examples
#' calculate_relative_indices(otolith_length = 5.2, otolith_area = 14.5, fish_length = 55)
calculate_relative_indices <- function(otolith_length, otolith_area, fish_length) {
    if (length(otolith_length) != length(otolith_area) &&
        length(otolith_length) != 1 && length(otolith_area) != 1) {
        stop("Lengths of otolith_length and otolith_area must match or one must be of length 1")
    }
    rol <- calculate_rol(otolith_length, fish_length)
    ros <- calculate_ros(otolith_area, fish_length)
    data.frame(ROL = rol, ROS = ros)
}

