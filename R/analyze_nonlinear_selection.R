# ============================================================================
# analyze_nonlinear_selection
#
# Purpose: Estimate quadratic (γ) and correlational (γ_ij) selection gradients
#
# IMPORTANT NOTE:
#   - ALWAYS use OLS to estimate selection gradients (β, γ, γ_ij)
#   - For binary fitness: OLS gives the correct gradient estimates, but p-values
#     come from logistic GLM (with Wald tests) because OLS residuals violate
#     normality assumptions.
#   - For continuous fitness: OLS provides both gradient estimates AND valid
#     p-values (via t-tests and F-tests).
#
# Model:
#   w = α + β₁z₁ + β₂z₂ + ½γ₁₁z₁² + ½γ₂₂z₂² + γ₁₂z₁z₂ + ε
#
# Where:
#   - β = linear selection gradients
#   - γ_ii = quadratic selection gradients (stabilizing/disruptive)
#   - γ_ij = correlational selection gradients (interactions)
#
# Workflow:
#   1. Build formula with linear, quadratic, and interaction terms
#   2. Fit OLS to get all gradient estimates (β, γ, γ_ij)
#   3. For binary fitness: fit logistic GLM for valid p-values
#   4. Type III ANOVA for significance tests
#   5. VIF check for multicollinearity
#
# Returns:
#   For continuous: list with model_ols, summary_ols, anova, vif, fitness_type
#   For binary: list with model_ols, model_glm, summary_ols, summary_glm,
#               anova (from GLM), vif, fitness_type
# ============================================================================

#' Analyze nonlinear selection gradients (gamma)
#'
#' Estimates quadratic and correlational selection gradients using OLS.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_cols A character vector of trait column names.
#' @param fitness_type A string indicating the fitness type: \code{"binary"} or \code{"continuous"}.
#'
#' @return A list containing the fitted nonlinear models, summaries, ANOVA tables, and VIFs.
#' @export
analyze_nonlinear_selection <- function(data, fitness_col, trait_cols, fitness_type) {
  if (length(trait_cols) < 2) {
    stop("Nonlinear selection requires at least 2 traits")
  }

  if (nrow(data) < 20) {
    warning("Small sample size (n < 20) — nonlinear estimates may be unreliable")
  }

  # Quadratic terms: I(trait1^2), I(trait2^2), ...
  quad <- paste0("I(", trait_cols, "^2)")

  # Interaction terms: trait1:trait2, trait1:trait3, ...
  inter <- combn(trait_cols, 2, FUN = function(x) paste(x, collapse = ":"), simplify = TRUE)

  rhs <- paste(c(trait_cols, quad, inter), collapse = " + ")

  resp <- if ("relative_fitness" %in% names(data)) {
    "relative_fitness"
  } else {
    fitness_col
  }


  # OLS
  fml_ols <- as.formula(paste(resp, "~", rhs))

  # Remove rows with missing data
  fit_data <- data[complete.cases(data[, c(resp, trait_cols)]), ]

  # Check if sample size is sufficient for the number of parameters
  n_params <- length(trait_cols) + length(quad) + length(inter) + 1 # +1 for intercept
  if (nrow(fit_data) < n_params * 2) {
    warning(
      "Sample size (", nrow(fit_data), ") may be too small for ",
      n_params, " parameters — results may be unreliable"
    )
  }

  # Fit OLS model
  fit_ols <- lm(fml_ols, data = fit_data)
  sm_ols <- summary(fit_ols)

  # Variance Inflation Factor (VIF) check
  vif_vals <- NULL
  if (requireNamespace("car", quietly = TRUE)) {
    vif_vals <- tryCatch(
      car::vif(fit_ols),
      error = function(e) {
        warning("VIF calculation failed: ", e$message)
        NULL
      }
    )
    if (!is.null(vif_vals) && any(vif_vals > 5)) {
      warning("High multicollinearity detected (VIF > 5) — standard errors may be inflated")
    }
  }


  if (fitness_type == "binary") {
    # LOGISTIC GLM (for valid p-values)
    fml_glm <- as.formula(paste(fitness_col, "~", rhs))
    fit_glm <- tryCatch(
      glm(fml_glm, data = fit_data, family = binomial),
      error = function(e) {
        stop("Nonlinear GLM fitting failed: ", e$message)
      }
    )
    sm_glm <- summary(fit_glm)

    # Convergence: If GLM doesn't converge, results are unreliable
    if (!fit_glm$converged) {
      warning("Nonlinear GLM did not converge — results may be unreliable")
    }

    if (any(abs(coef(fit_glm)) > 10)) {
      warning("Possible complete separation detected — large coefficients (>10)")
    }

    # Type III ANOVA
    anova_bin <- NULL
    if (requireNamespace("car", quietly = TRUE)) {
      anova_bin <- tryCatch(
        car::Anova(fit_glm, type = "III", test.statistic = "Wald"),
        error = function(e) {
          warning("Type III ANOVA for nonlinear model failed: ", e$message)
          NULL
        }
      )
    } else {
      warning("Package 'car' not installed — skipping Type III ANOVA")
    }

    return(list(
      mmodel = list(ols = fit_ols, glm = fit_glm),
      summary = list(ols = sm_ols, glm = sm_glm),
      anova = anova_bin,
      vif = vif_vals,
      fitness_type = "binary"
    ))
  } else {
    # Type III ANOVA
    anova_cont <- NULL
    if (requireNamespace("car", quietly = TRUE)) {
      anova_cont <- tryCatch(
        car::Anova(fit_ols, type = "III"),
        error = function(e) {
          warning("Type III ANOVA for nonlinear model failed: ", e$message)
          NULL
        }
      )
    } else {
      warning("Package 'car' not installed — skipping Type III ANOVA")
    }

    return(list(
      model_ols = fit_ols, # lm object (coefficients = gradients)
      summary_ols = sm_ols, # summary.lm (coefficients, SE, p-values)
      anova = anova_cont, # Type III ANOVA table
      vif = vif_vals, # variance inflation factors
      fitness_type = "continuous"
    ))
  }
}
