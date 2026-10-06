# ============================================================================
# analyze_nonlinear_selection
#
# Purpose: Estimate quadratic (gamma) and correlational (gamma_ij) selection gradients
#
# IMPORTANT NOTE:
#   - Gradients are estimated by OLS on relative fitness (w = W / mean W).
#   - For binary fitness: OLS on relative fitness gives the gradients; p-values
#     come from a logistic GLM on the raw 0/1 outcome (Wald tests).
#   - For count fitness: p-values from a Poisson GLM on the raw counts, or a
#     negative binomial one when they are overdispersed.
#   - For continuous or proportion fitness: OLS supplies both.
#
# Model:
#   w = alpha + beta1z1 + beta2z2 + 1/2gamma11z1^2 + 1/2gamma22z2^2 + gamma12z1z2 + epsilon
#
# Where:
#   - beta = linear selection gradients
#   - gamma_ii = quadratic selection gradients (stabilising or disruptive)
#   - gamma_ij = correlational selection gradients (interactions)
#
# Returns:
#   Non-binary: list with model, summary, anova, vif, fitness_type
#   Binary:     list with model (ols + glm), summary (ols + glm),
#               anova (from GLM), vif, fitness_type
# ============================================================================

#' Analyze nonlinear selection gradients (gamma)
#'
#' Estimates quadratic and correlational selection gradients by OLS on relative fitness.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the response column for the OLS gradient model (relative fitness).
#' @param trait_cols A character vector of trait column names.
#' @param fitness_type A string indicating the fitness type: \code{"binary"}, \code{"continuous"}, \code{"count"}, or \code{"proportion"}.
#' @param binary_response_col Optional string naming the raw fitness column (0/1 for binary, counts for count fitness) used for the GLM that supplies p-values. If \code{NULL}, \code{fitness_col} is treated as the raw outcome and relativised internally.
#' @param group Optional grouping column, used as in
#'   \code{analyze_linear_selection()}: a separate intercept for each group in
#'   the GLM that supplies the p-values.
#' @param se_type Standard errors of the least-squares gradients: \code{"ols"}
#'   (the default) or \code{"hc3"}, heteroscedasticity-consistent; see
#'   \code{selection_coefficients()}.
#'
#' @return A list containing the fitted nonlinear models, summaries, ANOVA tables, and VIFs.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' fit <- analyze_nonlinear_selection(prep, "survival", c("total_length", "weight"), "binary")
#' extract_quadratic_coefficients(c("total_length", "weight"), fit)
#' extract_interaction_coefficients(c("total_length", "weight"), fit)
#' @export
analyze_nonlinear_selection <- function(data, fitness_col, trait_cols, fitness_type,
                                        binary_response_col = NULL, group = NULL,
                                        se_type = c("ols", "hc3")) {
  se_type <- .se_type_arg(se_type)
  if (length(trait_cols) < 1) {
    stop("Nonlinear selection requires at least one trait")
  }

  if (nrow(data) < 20) {
    warning("Fewer than 20 individuals, so the nonlinear gradients are poorly estimated")
  }

  # Quadratic terms: I(trait1^2), I(trait2^2), ...
  quad <- paste0("I(", trait_cols, "^2)")

  # Interaction terms (correlational selection) only exist with two or more traits.
  inter <- if (length(trait_cols) >= 2) {
    combn(trait_cols, 2, FUN = function(x) paste(x, collapse = ":"), simplify = TRUE)
  } else {
    character(0)
  }

  rhs <- paste(c(trait_cols, quad, inter), collapse = " + ")
  n_params <- length(trait_cols) + length(quad) + length(inter) + 1 # +1 for intercept

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

    .warn_small_sample(nrow(fit_data), n_params)

    fit_ols <- lm(as.formula(paste(ols_resp, "~", rhs)), data = fit_data)
    sm_ols <- summary(fit_ols)
    if (se_type == "hc3") sm_ols <- .hc3_summary(fit_ols, sm_ols)
    vif_vals <- .compute_vif(fit_ols)

    fit_data <- .add_glm_group(fit_data, group)
    glm_rhs <- if (".group" %in% names(fit_data)) paste(".group +", rhs) else rhs
    fit_glm <- tryCatch(
      .fit_pvalue_glm(as.formula(paste(glm_col, "~", glm_rhs)), fit_data, fitness_type),
      error = function(e) stop("Nonlinear GLM fitting failed: ", e$message)
    )
    sm_glm <- summary(fit_glm)

    if (isFALSE(fit_glm$converged)) {
      warning("Quadratic GLM did not converge; its p-values may be off")
    }
    # only the trait terms count here; a group where everyone survived gets a
    # huge intercept without any separation by the traits
    slopes <- coef(fit_glm)[!grepl("^\\(Intercept\\)$|^\\.group", names(coef(fit_glm)))]
    if (fitness_type == "binary" && any(abs(slopes) > 10, na.rm = TRUE)) {
      warning("Possible complete separation, with a logistic slope above 10 in absolute value")
    }

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
      warning("No Type III ANOVA without the car package")
    }

    return(list(
      model = list(ols = fit_ols, glm = fit_glm),
      summary = list(ols = sm_ols, glm = sm_glm),
      anova = anova_bin,
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
    .warn_small_sample(nrow(fit_data), n_params)

    fit_ols <- lm(as.formula(paste(fitness_col, "~", rhs)), data = fit_data)
    sm_ols <- summary(fit_ols)
    if (se_type == "hc3") sm_ols <- .hc3_summary(fit_ols, sm_ols)
    vif_vals <- .compute_vif(fit_ols)

    anova_cont <- NULL
    if (requireNamespace("car", quietly = TRUE)) {
      anova_cont <- tryCatch(
        if (se_type == "hc3") car::Anova(fit_ols, type = "III", white.adjust = "hc3") else car::Anova(fit_ols, type = "III"),
        error = function(e) {
          warning("Type III ANOVA for nonlinear model failed: ", e$message)
          NULL
        }
      )
    } else {
      warning("No Type III ANOVA without the car package")
    }

    return(list(
      model = fit_ols, # lm object (coefficients = gradients)
      summary = sm_ols, # summary.lm (coefficients, SE, p-values)
      anova = anova_cont, # Type III ANOVA table
      vif = vif_vals, # variance inflation factors
      fitness_type = fitness_type
    ))
  }
}
