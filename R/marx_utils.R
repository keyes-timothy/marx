# marx_utils.R
# This script contains functions for various marx helpers and utilities.


#' Calculate the magnitudes of the disease-specific and lineage-specific
#' components in a marx projection.
#'
#' @param marx_projection The output data.frame of marx_project()
#'
#' @return A data.frame containing the healthy- and disease-specific component
#' magnitudes in marx_projection.
#'
#' @export
#'
#' @importFrom dplyr ends_with
#' @importFrom dplyr select
#' @importFrom dplyr tibble
#'
marx_find_magnitudes <- function(marx_projection) {

  # find lineage magnitudes
  lineage_magnitudes <-
    marx_projection |>
    dplyr::select(dplyr::ends_with("_lineage")) |>
    as.matrix() |>
    (\(x) x^2)() |>
    rowSums() |>
    sqrt()

  # find disease magnitudes
  disease_magnitudes <-
    marx_projection |>
    dplyr::select(dplyr::ends_with("_disease")) |>
    as.matrix() |>
    (\(x) x^2)() |>
    rowSums() |>
    sqrt()

  # find total magnitude
  total_magnitudes <- sqrt(lineage_magnitudes^2 + disease_magnitudes^2)

  result <-
    dplyr::tibble(
      lineage_magnitude = lineage_magnitudes,
      disease_magnitude = disease_magnitudes,
      total_magnitude = total_magnitudes
    )

  return(result)
}
