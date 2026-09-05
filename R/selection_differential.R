# ============================================================================
# selection_differential
#
# Purpose: Calculate selection differential (S) for a single trait
#
# Definition: S = Cov(z, w)
#   where:
#     z = standardized trait value
#     w = relative fitness
#
# IMPORTANT NOTE:
#   When traits are standardized to mean 0 and SD 1, and fitness is relativized
#   to mean 1, then S = mean(z x w)
#
#   For multi-year/site studies, S should be calculated WITHIN each group
#   (e.g., within each year) to ensure individuals are compared to their
#   relevant context.
#
# Parameters:
#   data               : data frame with fitness and trait measurements
#   fitness_col        : name of the fitness column
#   trait_col          : name of the trait column
#   standardized       : if TRUE, assume trait is already standardized
#                        if FALSE, standardize the trait
#   use_relative       : if TRUE, use relative fitness (w / mean(w))
#   group              : optional grouping variable (e.g., "year", "site")
#                        When specified, S is calculated separately within
#                        each group and averaged (or returned as vector)
#   return_grouped     : if TRUE and group is specified, return a data frame
#                        with S per group; if FALSE, return overall mean S
#
# Returns:
#   - If group = NULL: single numeric value (overall S)
#   - If group is specified and return_grouped = TRUE: data frame with S per group
#   - If group is specified and return_grouped = FALSE: mean S across groups
# ============================================================================

#' Calculate selection differential (S)
#'
#' Computes the selection differential S = Cov(z, w) for a single trait.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_col A string specifying the name of the trait column.
#' @param standardized Logical indicating whether the trait is already standardized. If \code{FALSE}, it will be standardized.
#' @param use_relative Logical indicating whether to use relative fitness.
#' @param group Optional grouping variable string. Calculates S within each group if provided.
#' @param return_grouped Logical indicating whether to return a data frame with S per group (\code{TRUE}) or overall mean S (\code{FALSE}).
#'
#' @return A single numeric value representing overall S, or a data frame with S per group if \code{return_grouped = TRUE}.
#' @export
#'
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", "total_length")
#' selection_differential(prep, "survival", "total_length")
#' selection_differential(prep, "survival", "total_length", group = "sex", return_grouped = TRUE)
selection_differential <- function(data,
                                   fitness_col,
                                   trait_col,
                                   standardized = TRUE,
                                   use_relative = TRUE,
                                   group = NULL,
                                   return_grouped = FALSE) {
  # Input validation
  if (!trait_col %in% names(data)) {
    stop("Trait column '", trait_col, "' not found in data")
  }
  if (!fitness_col %in% names(data)) {
    stop("Fitness column '", fitness_col, "' not found in data")
  }

  # Helper function to calculate S for a subset of data
  calc_S <- function(subset) {
    z <- subset[[trait_col]]
    w <- subset[[fitness_col]]

    # Remove NAs
    keep <- stats::complete.cases(z, w)
    z <- z[keep]
    w <- w[keep]

    if (length(z) < 2) {
      return(NA_real_)
    }

    # Standardize trait if needed
    if (!standardized) {
      # Check for zero variance
      if (sd(z) == 0) {
        warning("Trait has zero variance in group; returning NA")
        return(NA_real_)
      }
      z <- as.numeric(scale(z))
    }

    # Calculate relative fitness if needed
    if (use_relative) {
      mu <- mean(w, na.rm = TRUE)
      if (!is.finite(mu) || mu == 0) {
        warning("Mean fitness is zero or non-finite; cannot compute relative fitness")
        return(NA_real_)
      }
      w <- w / mu
    }

    # S = Cov(z, w), the population covariance (the selection-induced change in
    # mean trait). Centring both terms makes this exact whether or not the
    # trait was pre-centred, and reduces to mean(z * w) when the trait is centred.
    mean((z - mean(z)) * (w - mean(w)))
  }

  if (!is.null(group)) {
    if (!group %in% names(data)) {
      stop("Group column '", group, "' not found in data")
    }

    # Drop rows missing trait, fitness, or the group label up front so each
    # group's n counts only the observations calc_S uses; otherwise a
    # group with more missing data (or an NA-group phantom level) would skew the
    # pooled mean below.
    usable <- stats::complete.cases(data[[trait_col]], data[[fitness_col]], data[[group]])

    # Calculate S for each group
    S_by_group <- data[usable, , drop = FALSE] %>%
      dplyr::group_by(.data[[group]]) %>%
      dplyr::summarise(
        S = calc_S(dplyr::pick(dplyr::everything())),
        n = dplyr::n(),
        .groups = "drop"
      )

    if (return_grouped) {
      return(S_by_group)
    } else {
      # Return weighted mean (by sample size), ignoring groups with no estimate.
      # If no group yields an estimate, return NA (not NaN) to match the
      # ungrouped degenerate path.
      ok <- !is.na(S_by_group$S)
      if (!any(ok)) {
        return(NA_real_)
      }
      weighted_mean <- sum(S_by_group$S[ok] * S_by_group$n[ok]) /
        sum(S_by_group$n[ok])
      return(weighted_mean)
    }
  } else {
    # No grouping: calculate overall S
    return(calc_S(data))
  }
}
