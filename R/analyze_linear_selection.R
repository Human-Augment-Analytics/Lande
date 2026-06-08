# ============================================================================
# analyze_linear_selection
#
# Purpose: Estimate linear selection gradients (beta) using Lande & Arnold (1983)
#
# IMPORTANT NOTE:
#   - Gradients are estimated by OLS on RELATIVE fitness (w = W / mean W).
#     The caller (selection_coefficients) supplies relative fitness as
#     `fitness_col`; a direct caller may pass a raw fitness column, which is
#     relativised internally for binary data.
#   - For binary fitness (0/1 survival): OLS on relative fitness gives the
#     correct beta, but its p-values are not trustworthy because the residuals
#     violate normality. A logistic GLM on the raw 0/1 outcome supplies valid
#     p-values (via Wald tests).
#   - For continuous, count, or proportion fitness: OLS provides both the
#     gradients and valid p-values (t- and F-tests).
#
# Returns:
#   Non-binary: list with model (lm), summary, anova, vif, fitness_type
#   Binary:     list with model (ols + glm), summary (ols + glm),
#               anova (from GLM), vif, fitness_type
# ============================================================================

#' Analyze linear selection gradients (beta)
#'
#' Estimates linear selection gradients by OLS on relative fitness, using logistic regression for valid p-values when fitness is binary.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the response column for the OLS gradient model (relative fitness).
#' @param trait_cols A character vector of trait column names.
#' @param fitness_type A string indicating the fitness type: \code{"binary"}, \code{"continuous"}, \code{"count"}, or \code{"proportion"}.
#' @param binary_response_col Optional string naming the raw 0/1 column used for the logistic GLM when \code{fitness_type = "binary"}. If \code{NULL}, \code{fitness_col} is treated as the raw outcome and relativised internally.
#'
#' @return A list containing the fitted models, summaries, ANOVA tables, and VIFs.
#' @export
analyze_linear_selection <- function(data, fitness_col, trait_cols, fitness_type,
                                     binary_response_col = NULL) {
  # Check sample size
  if (nrow(data) < 10) {
    warning("Small sample size (n < 10) - results may be unreliable")
  }

  # Build fitness ~ trait1 + trait2 + trait3
  rhs <- paste(trait_cols, collapse = " + ")

  if (fitness_type == "binary") {
    # Gradients from OLS on relative fitness; p-values from a logistic GLM on
    # the raw 0/1 outcome.
    glm_col <- if (!is.null(binary_response_col)) binary_response_col else fitness_col
    fit_data <- data[complete.cases(data[, c(fitness_col, glm_col, trait_cols)]), ]

    if (is.null(binary_response_col)) {
      fit_data$.rel_fitness <- fit_data[[fitness_col]] / mean(fit_data[[fitness_col]])
      ols_resp <- ".rel_fitness"
    } else {
      ols_resp <- fitness_col
    }

    fit_ols <- lm(as.formula(paste(ols_resp, "~", rhs)), data = fit_data)
    sm_ols <- summary(fit_ols)

    fit_glm <- glm(as.formula(paste(glm_col, "~", rhs)), data = fit_data, family = binomial)
    sm_glm <- summary(fit_glm)

    # Convergence: If the GLM algorithm fails to converge, results are unreliable
    if (!fit_glm$converged) {
      warning("GLM did not converge - results may be unreliable")
    }

    if (any(abs(coef(fit_glm)) > 10)) {
      warning("Possible complete separation detected - large coefficients (>10)")
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
      warning("Package 'car' not installed - skipping Type III ANOVA")
    }

    # Variance Inflation Factor (VIF) check, computed on the OLS fit
    vif_vals <- .compute_vif(fit_ols)

    return(list(
      model = list(
        ols = fit_ols, # lm object on relative fitness (coefficients = beta)
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
  } else {
    # Continuous, count, or proportion fitness: OLS supplies gradients and
    # p-values. The caller passes relative fitness as `fitness_col`.
    fit_data <- data[complete.cases(data[, c(fitness_col, trait_cols)]), ]
    fit_ols <- lm(as.formula(paste(fitness_col, "~", rhs)), data = fit_data)
    sm_ols <- summary(fit_ols)

    # Type III ANOVA gives the partial regression coefficients (beta) and their significance.
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
      warning("Package 'car' not installed - skipping Type III ANOVA")
    }

    vif_vals <- .compute_vif(fit_ols)

    return(list(
      model = fit_ols, # lm object (contains coefficients, etc.)
      summary = sm_ols, # summary.lm object (p-values)
      anova = anova_cont, # Type III ANOVA table
      vif = vif_vals, # variance inflation factors
      fitness_type = fitness_type
    ))
  }
}
