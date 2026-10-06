# ============================================================================
# analyze_linear_selection
#
# Purpose: Estimate linear selection gradients (beta) using Lande & Arnold (1983)
#
# IMPORTANT NOTE:
#   - Gradients are estimated by OLS on relative fitness (w = W / mean W).
#     The caller (selection_coefficients) supplies relative fitness as
#     `fitness_col`; a direct caller may pass a raw fitness column, which is
#     relativised internally for binary data.
#   - For binary fitness (0/1 survival): OLS on relative fitness gives the
#     correct beta, but its p-values are not trustworthy because the residuals
#     violate normality. A logistic GLM on the raw 0/1 outcome supplies valid
#     p-values (via Wald tests).
#   - For count fitness (offspring numbers): the same split, with a Poisson
#     GLM on the raw counts for the p-values, or a negative binomial one when
#     the counts are overdispersed.
#   - For continuous or proportion fitness: OLS provides both the gradients
#     and valid p-values (t- and F-tests).
#
# Returns:
#   Continuous:   list with model (lm), summary, anova, vif, fitness_type
#   Binary/count: list with model (ols + glm), summary (ols + glm),
#                 anova (from GLM), vif, fitness_type, glm_family
# ============================================================================

#' Analyze linear selection gradients (beta)
#'
#' Estimates linear selection gradients by OLS on relative fitness. Binary
#' fitness gets its p-values from a logistic GLM and count fitness from a
#' Poisson GLM (negative binomial if overdispersed), since the OLS tests are not
#' valid for either.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the response column for the OLS gradient model (relative fitness).
#' @param trait_cols A character vector of trait column names.
#' @param fitness_type A string indicating the fitness type: \code{"binary"}, \code{"continuous"}, \code{"count"}, or \code{"proportion"}.
#' @param binary_response_col Optional string naming the raw fitness column (0/1 for binary, counts for count fitness) used for the GLM that supplies p-values. If \code{NULL}, \code{fitness_col} is treated as the raw outcome and relativised internally.
#' @param group Optional grouping column. With two or more groups the GLM that
#'   supplies the p-values gets an intercept for each, to match the
#'   standardising within groups; the least-squares fit has no group term.
#' @param se_type Standard errors of the least-squares gradients: \code{"ols"}
#'   (the default) or \code{"hc3"}, heteroscedasticity-consistent; see
#'   \code{selection_coefficients()}.
#'
#' @return A list containing the fitted models, summaries, ANOVA tables, and VIFs.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' fit <- analyze_linear_selection(prep, "survival", c("total_length", "weight"), "binary")
#' extract_linear_coefficients(c("total_length", "weight"), fit)
#' @export
analyze_linear_selection <- function(data, fitness_col, trait_cols, fitness_type,
                                     binary_response_col = NULL, group = NULL,
                                     se_type = c("ols", "hc3")) {
  se_type <- .se_type_arg(se_type)
  # Check sample size
  if (nrow(data) < 10) {
    warning("Fewer than 10 individuals, so the gradients are poorly estimated")
  }

  # Build fitness ~ trait1 + trait2 + trait3
  rhs <- paste(trait_cols, collapse = " + ")

  if (fitness_type %in% c("binary", "count")) {
    # Gradients from OLS on relative fitness; p-values from a GLM on the raw
    # outcome (logistic for 0/1, Poisson or negative binomial for counts).
    glm_col <- if (!is.null(binary_response_col)) binary_response_col else fitness_col
    fit_data <- data[complete.cases(data[, c(fitness_col, glm_col, trait_cols)]), ]

    if (is.null(binary_response_col)) {
      raw <- fit_data[[fitness_col]]
      if (!.is_raw_fitness(raw, fitness_type)) {
        stop(
          "For ", fitness_type, " fitness pass the raw ",
          if (fitness_type == "binary") "0/1" else "count",
          " column as `fitness_col` (it is relativised internally), ",
          "or name it in `binary_response_col`."
        )
      }
      fit_data$.rel_fitness <- raw / mean(raw)
      ols_resp <- ".rel_fitness"
    } else {
      ols_resp <- fitness_col
    }

    fit_ols <- lm(as.formula(paste(ols_resp, "~", rhs)), data = fit_data)
    sm_ols <- summary(fit_ols)
    if (se_type == "hc3") sm_ols <- .hc3_summary(fit_ols, sm_ols)

    fit_data <- .add_glm_group(fit_data, group)
    glm_rhs <- if (".group" %in% names(fit_data)) paste(".group +", rhs) else rhs
    fit_glm <- .fit_pvalue_glm(as.formula(paste(glm_col, "~", glm_rhs)), fit_data, fitness_type)
    sm_glm <- summary(fit_glm)

    # Convergence: If the GLM algorithm fails to converge, results are unreliable
    if (isFALSE(fit_glm$converged)) {
      warning("GLM did not converge; its p-values may be off")
    }

    # separation is about the traits: a group in which every bird survived has
    # a large intercept of its own, which is not the traits separating the outcome
    slopes <- coef(fit_glm)[!grepl("^\\(Intercept\\)$|^\\.group", names(coef(fit_glm)))]
    if (fitness_type == "binary" && any(abs(slopes) > 10, na.rm = TRUE)) {
      warning("Possible complete separation, with a logistic slope above 10 in absolute value")
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
      warning("No Type III ANOVA without the car package")
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
      fitness_type = fitness_type,
      glm_family = attr(fit_glm, "family_label")
    ))
  } else {
    # Continuous or proportion fitness: OLS supplies gradients and p-values.
    # the column given is the one fitted; if the data also has relative_fitness,
    # say so
    if (is.null(binary_response_col) && fitness_col != "relative_fitness" &&
        !grepl("_relative$", fitness_col) && "relative_fitness" %in% names(data)) {
      message("Fitting '", fitness_col, "' as given; the data also has 'relative_fitness'")
    }
    fit_data <- data[complete.cases(data[, c(fitness_col, trait_cols)]), ]
    fit_ols <- lm(as.formula(paste(fitness_col, "~", rhs)), data = fit_data)
    sm_ols <- summary(fit_ols)
    if (se_type == "hc3") sm_ols <- .hc3_summary(fit_ols, sm_ols)

    # Type III ANOVA gives the partial regression coefficients (beta) and their significance.
    anova_cont <- NULL
    if (requireNamespace("car", quietly = TRUE)) {
      anova_cont <- tryCatch(
        if (se_type == "hc3") car::Anova(fit_ols, type = "III", white.adjust = "hc3") else car::Anova(fit_ols, type = "III"),
        error = function(e) {
          warning("Type III ANOVA failed: ", e$message)
          NULL
        }
      )
    } else {
      warning("No Type III ANOVA without the car package")
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
