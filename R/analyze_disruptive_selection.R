# ============================================================================
# analyze_disruptive_selection
#
# Purpose: Detect disruptive or stabilizing selection on a single trait.
#
# Model: w = alpha + betaz + gammaz^2 + epsilon
#   beta  = linear (directional) selection gradient
#   gamma = quadratic selection gradient (gamma > 0: disruptive; gamma < 0: stabilizing)
#
# This is the single-trait case of the standard analysis, so it defers to
# selection_coefficients(). That keeps beta on the linear-only fit and the
# doubled gamma on the full quadratic fit, so a trait analysed here and via
# selection_coefficients() always agrees. See analyze_linear_selection() and
# analyze_nonlinear_selection() for the underlying models.
# ============================================================================

#' Analyze disruptive/stabilizing selection
#'
#' Detects disruptive or stabilizing selection on a single trait by estimating its
#' linear (beta) and quadratic (gamma) selection gradients. A positive gamma indicates
#' disruptive selection; a negative gamma indicates stabilizing selection.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_col A string specifying the name of the single trait column.
#' @param fitness_type A string indicating the fitness type: \code{"auto"} (detect from the data),
#'   \code{"binary"}, or \code{"continuous"}. Default is \code{"auto"}.
#' @param standardize Logical indicating whether to standardize the trait to mean 0 and SD 1. Default is \code{TRUE}.
#' @param group Optional string specifying a grouping variable; standardisation and
#'   relative fitness are then computed within each group.
#' @param return_grouped Logical; if \code{TRUE} and \code{group} is given, the
#'   gradients are estimated separately for each group and returned with a
#'   \code{Group} column. Default is \code{FALSE}.
#'
#' @return A data frame with one row per gradient (\code{Term}, \code{Type},
#'   \code{Beta_Coefficient}, \code{Standard_Error}, \code{P_Value}, \code{Variance}).
#' @export
analyze_disruptive_selection <- function(
  data,
  fitness_col,
  trait_col,
  fitness_type = c("auto", "binary", "continuous"),
  standardize = TRUE,
  group = NULL,
  return_grouped = FALSE
) {
  fitness_type <- match.arg(fitness_type)

  if (length(trait_col) != 1) {
    stop("analyze_disruptive_selection() takes a single trait; use selection_coefficients() for several.")
  }

  selection_coefficients(
    data           = data,
    fitness_col    = fitness_col,
    trait_cols     = trait_col,
    fitness_type   = fitness_type,
    standardize    = standardize,
    group          = group,
    return_grouped = return_grouped
  )
}
