#' Example Metabolomic Longitudinal Dataset
#'
#' A list containing:
#' \describe{
#'   \item{n}{Integer, number of subjects.}
#'   \item{p}{Integer, number of features.}
#'   \item{dat}{List of length \code{p}, each a list of length \code{n} with elements \code{time} and \code{val}.}
#'   \item{sub_label}{Character vector of subject IDs.}
#'   \item{feature_label}{Character vector of feature names.}
#'   \item{time_grid}{Numeric vector of normalized time points (0–1).}
#'   \item{time_min}{Numeric, original minimum time.}
#'   \item{time_max}{Numeric, original maximum time.}
#' }
#' @source Generated via \code{wrap_dat()}
"dat_X"

#' Example Lipid Longitudinal Dataset
#'
#' Same structure as \code{dat_X}, but for lipid measurements.
#' @source Generated via \code{wrap_dat()}
"dat_Y"