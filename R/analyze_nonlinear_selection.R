# ============================================================================
# analyze_nonlinear_selection
#
# Purpose: Estimate quadratic (gamma) and correlational (gamma_ij) selection gradients
#
# IMPORTANT NOTE:
#   - Gradients are estimated by OLS on RELATIVE fitness (w = W / mean W).
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
#   - gamma_ii = quadratic selection gradients (stabilizing/disruptive)
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
#'
#' @return A list containing the fitted nonlinear models, summaries, ANOVA tables, and VIFs.
#' @export
analyze_nonlinear_selection <- function(data, fitness_col, trait_cols, fitness_type,
                                        binary_response_col = NULL) {
  if (length(trait_cols) < 1) {
    stop("Nonlinear selection requires at least one trait")
  }

  if (nrow(data) < 20) {
    warning("Small sample size (n < 20) - nonlinear estimates may be unreliable")
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
    vif_vals <- .compute_vif(fit_ols)

    fit_glm <- tryCatch(
      .fit_pvalue_glm(as.formula(paste(glm_col, "~", rhs)), fit_data, fitness_type),
      error = function(e) stop("Nonlinear GLM fitting failed: ", e$message)
    )
    sm_glm <- summary(fit_glm)

    if (isFALSE(fit_glm$converged)) {
      warning("Nonlinear GLM did not converge - results may be unreliable")
    }
    if (fitness_type == "binary" && any(abs(coef(fit_glm)) > 10, na.rm = TRUE)) {
      warning("Possible complete separation detected - large coefficients (>10)")
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
      warning("Package 'car' not installed - skipping Type III ANOVA")
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
    # Continuous or proportion fitness: OLS supplies gradients and
    # p-values. A direct caller working from prepare_selection_data() output
    # gets its relative_fitness column picked up automatically.
    if (is.null(binary_response_col) && fitness_col != "relative_fitness" &&
        "relative_fitness" %in% names(data)) {
      message("Using 'relative_fitness' as the response for gradient estimation")
      fitness_col <- "relative_fitness"
    }
    fit_data <- data[complete.cases(data[, c(fitness_col, trait_cols)]), ]
    .warn_small_sample(nrow(fit_data), n_params)

    fit_ols <- lm(as.formula(paste(fitness_col, "~", rhs)), data = fit_data)
    sm_ols <- summary(fit_ols)
    vif_vals <- .compute_vif(fit_ols)

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
      warning("Package 'car' not installed - skipping Type III ANOVA")
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
