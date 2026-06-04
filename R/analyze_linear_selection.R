# ============================================================================
# analyze_linear_selection
#
# Purpose: Estimate linear selection gradients (β) using Lande & Arnold (1983)
#
# IMPORTANT NOTE:
#   - ALWAYS use OLS (ordinary least squares) to estimate selection gradients
#   - For binary fitness (0/1 survival): OLS gives the correct β estimates,
#     but the p-values from OLS are not trustworthy because residuals violate
#     normality assumptions. Therefore, we use logistic regression (GLM)
#     specifically to obtain valid p-values.
#   - For continuous fitness (e.g., growth rate, offspring number): OLS
#     provides both valid coefficient estimates AND valid p-values (via t-tests
#     and F-tests), because the residuals can reasonably approximate normality.
#
# Workflow:
#   Continuous fitness:
#     1. Fit OLS model → coefficients = selection gradients (β)
#     2. p-values come directly from OLS summary
#     3. Type III ANOVA (car::Anova) for significance of each term
#     4. Check VIF for multicollinearity
#
#   Binary fitness:
#     1. Fit OLS model → coefficients = selection gradients (β)
#     2. Fit logistic GLM model → p-values from Wald tests (Type III)
#     3. Check convergence and separation
#     4. Type III ANOVA on GLM for overall term significance
#     5. Check VIF using OLS (GLM VIF not directly comparable)
#
# Returns:
#   For continuous: list with model (lm), summary, anova, vif, fitness_type
#   For binary: list with model (list of ols and glm), summary (list of ols and glm),
#               anova (from GLM), vif, fitness_type
# ============================================================================

#' Analyze linear selection gradients (beta)
#'
#' Estimates linear selection gradients using OLS, and logistic regression for valid p-values when fitness is binary.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_cols A character vector of trait column names.
#' @param fitness_type A string indicating the fitness type: \code{"binary"} or \code{"continuous"}.
#'
#' @return A list containing the fitted models, summaries, ANOVA tables, and VIFs.
#' @export
analyze_linear_selection <- function(data, fitness_col, trait_cols, fitness_type) {
  # Check sample size
  if (nrow(data) < 10) {
    warning("Small sample size (n < 10) — results may be unreliable")
  }

  # Build fitness ~ trait1 + trait2 + trait3
  rhs <- paste(trait_cols, collapse = " + ")

  if (fitness_type == "continuous") {
    # Use relative fitness if available (standardized to mean = 1)
    resp <- if ("relative_fitness" %in% names(data)) {
      "relative_fitness"
    } else {
      fitness_col
    }

    # OLS
    fml <- as.formula(paste(resp, "~", rhs))

    # Remove rows with missing data
    fit_data <- data[complete.cases(data[, c(resp, trait_cols)]), ]
    fit_ols <- lm(fml, data = fit_data)
    sm_ols <- summary(fit_ols)

    # Type III ANOVA, it gives the partial regression coefficients (β) and their significance.
    anova_cont <- NULL
    if (requireNamespace("car", quietly = TRUE)) {
      anova_cont <- tryCatch(
        car::Anova(fit_ols, type = "III"),
        error = function(e) {
          warning("Type III ANOVA failed: ", e$message)
          NULL
        }
      )
    } else {
      warning("Package 'car' not installed — skipping Type III ANOVA")
    }

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

    return(list(
      model = fit_ols, # lm object (contains coefficients, etc.)
      summary = sm_ols, # summary.lm object (p-values)
      anova = anova_cont, # Type III ANOVA table
      vif = vif_vals, # variance inflation factors
      fitness_type = "continuous"
    ))
  } else {
    # For binary fitness, we use the original 0/1 values
    fml <- as.formula(paste(fitness_col, "~", rhs))

    # Remove rows with missing data
    fit_data <- data[complete.cases(data[, c(fitness_col, trait_cols)]), ]

    # OLS
    fit_ols <- lm(fml, data = fit_data)
    sm_ols <- summary(fit_ols)

    # LOGISTIC GLM (for valid p-values)
    fit_glm <- glm(fml, data = fit_data, family = binomial)
    sm_glm <- summary(fit_glm)

    # Convergence: If the GLM algorithm fails to converge, results are unreliable
    if (!fit_glm$converged) {
      warning("GLM did not converge — results may be unreliable")
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
          warning("Type III ANOVA failed: ", e$message)
          NULL
        }
      )
    } else {
      warning("Package 'car' not installed — skipping Type III ANOVA")
    }

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

    return(list(
      model = list(
        ols = fit_ols, # lm object (coefficients = β)
        glm = fit_glm # glm object (for p-values and ANOVA)
      ),
      summary = list(
        ols = sm_ols, # summary.lm (coefficients, SE)
        glm = sm_glm # summary.glm (p-values from Wald tests)
      ),
      anova = anova_bin, # Type III ANOVA from GLM
      vif = vif_vals,
      fitness_type = "binary"
    ))
  }
}
