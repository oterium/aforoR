#' Aphanopus morphometric dataset
#'
#' A dataset containing morphometric measurements and biological data for two species
#' of Aphanopus (A. carbo and A. intermedius).
#'
#' @format A data frame with 14 variables:
#' \describe{
#'   \item{Species}{Species name (Acar or Aint)}
#'   \item{ID}{Unique identifier for the otolith (from image filename)}
#'   \item{SL}{Standard Length of the fish (mm)}
#'   \item{OA}{Otolith area (pixels)}
#'   \item{OP}{Otolith perimeter (pixels)}
#'   \item{OL}{Otolith maximum length (pixels)}
#'   \item{OW}{Otolith maximum width (pixels)}
#'   \item{Roundness}{Shape index: (4 * Area) / (pi * Length^2)}
#'   \item{FormFactor}{Shape index: (4 * pi * Area) / Perimeter^2}
#'   \item{Circularity}{Shape index: Perimeter^2 / Area}
#'   \item{Rectangularity}{Shape index: Area / (Length * Width)}
#'   \item{Ellipticity}{Shape index: (Length - Width) / (Length + Width)}
#'   \item{AspectRatio}{Shape index: Length / Width}
#'   \item{Units}{Measurement units (px)}
#' }
"Aphanopus"

#' Aphanopus wavelet coefficients (Scale 5)
#'
#' Wavelet coefficients at scale 5 for Aphanopus carbo and A. intermedius otoliths.
#'
#' @format A data frame with 514 variables: `Species`, `ID`, and wavelet coefficients from `W1` to `W512`.
"Aphanopus_W5"

#' Aphanopus elliptic Fourier descriptors
#'
#' Elliptic Fourier descriptors for Aphanopus carbo and A. intermedius otoliths.
#'
#' @format A data frame with 129 variables: `Species`, `ID`, and Fourier coefficients `A2` through `A32`, `B1`..`B32`, `C1`..`C32`, and `D1`..`D32`.
"Aphanopus_EF"
