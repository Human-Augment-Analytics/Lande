# ============================================================================
# detect_family
#
# Purpose: Automatically determine the appropriate fitness type and GLM family
#
# Detection logic:
#   - Binary: values are 0/1 (survival, presence/absence)
#   - Count: non-negative integers (e.g., offspring number)
#   - Continuous: any other numeric values, including proportions in [0,1]
#     and relative fitness (mean ~ 1)
#
# Selection gradients always come from OLS on relative fitness. Binary fitness
# additionally uses a logistic GLM for p-values; every other type uses OLS.
# The returned `family` is a suggestion for a user's own diagnostics, not a
# family the package fits automatically.
#
# Parameters:
#   y : fitness vector
#
# Returns:
#   list with:
#     type : "binary", "count", or "continuous" (proportions report "continuous")
#     family : suggested GLM family object
#     note : additional information about detection
# ============================================================================

#' Automatically detect fitness type and GLM family
#'
#' Heuristically determines whether a fitness vector represents binary, count, proportion, or continuous data.
#'
#' @param y A numeric vector representing fitness values.
#'
#' @return A list containing \code{type} (string), \code{family} (GLM family object), and \code{note} (string).
#' @export
#'
#' @examples
#' detect_family(c(0, 1, 0, 0, 1, 1))
detect_family <- function(y) {
  # Clean data
  y_clean <- y[!is.na(y)]

  if (length(y_clean) == 0) {
    stop("No non-NA values in fitness vector")
  }

  # Get unique values and summary statistics
  unique_vals <- sort(unique(y_clean))
  n_unique <- length(unique_vals)
  min_val <- min(y_clean, na.rm = TRUE)
  max_val <- max(y_clean, na.rm = TRUE)
  mean_val <- mean(y_clean, na.rm = TRUE)

  # Detect BINARY (0/1)
  is_binary <- n_unique <= 2 && all(unique_vals %in% c(0, 1))

  if (is_binary) {
    if (length(y_clean) < 10) {
      warning("Binary fitness detected, but sample size < 10 (may be unstable)")
    }

    if (all(y_clean == 0) || all(y_clean == 1)) {
      warning("Complete separation detected: all fitness values are the same")
    }

    return(list(
      type = "binary",
      family = stats::binomial("logit"),
      note = "Binary fitness (0/1) detected. Use logistic GLM for p-values."
    ))
  }

  # Detect COUNT DATA (non-negative integers)
  # Check if values are integers (within tolerance) and non-negative
  is_integer <- all(abs(y_clean - round(y_clean)) < 1e-8)
  is_non_negative <- all(y_clean >= 0)

  if (is_integer && is_non_negative && !is_binary) {
    if (length(y_clean) < 20) {
      warning("Count fitness detected, but sample size < 20 (may be unstable)")
    }

    return(list(
      type = "count",
      family = stats::poisson("log"),
      note = paste(
        "Count fitness detected (non-negative integers). Gradients and",
        "p-values use OLS on relative fitness; Poisson family suggested",
        "if you want to check p-values with a count GLM."
      )
    ))
  }

  # Detect PROPORTION DATA (values in [0,1] with >2 unique values)
  # This could be survival rates, proportions, etc.
  is_proportion <- min_val >= 0 && max_val <= 1 && !is_binary

  if (is_proportion) {
    # For proportion data, binomial with weights would be ideal,
    # but we don't have denominator info. Treat as continuous for now.
    # Alternatively, could use quasi-binomial.
    return(list(
      type = "continuous",
      family = stats::gaussian(),
      note = paste(
        "Proportion fitness detected (values in [0,1]). Using Gaussian family.",
        "Consider binomial with weights if the denominator is available."
      )
    ))
  }

  # Detect RELATIVE FITNESS (mean ~ 1)
  # Relative fitness is often used in selection analysis
  is_relative <- abs(mean_val - 1) < 0.1 && !is_binary

  if (is_relative) {
    return(list(
      type = "continuous",
      family = stats::gaussian(),
      note = "Relative fitness detected (mean ~ 1). Using Gaussian family."
    ))
  }

  # DEFAULT: CONTINUOUS
  return(list(
    type = "continuous",
    family = stats::gaussian(),
    note = "Continuous fitness detected. Using Gaussian family."
  ))
}
