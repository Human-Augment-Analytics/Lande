# ============================================================================
# analyze_disruptive_selection
#
# Purpose: Detect disruptive or stabilizing selection on a single trait
#
# Model: w = α + βz + γz² + ε
#   where:
#     β = linear selection gradient (directional selection)
#     γ = quadratic selection gradient (γ > 0: disruptive; γ < 0: stabilizing)
#
# IMPORTANT NOTES:
#   For ALL fitness types, selection gradients (β and γ) come from OLS.
#   P-values come from:
#     - Continuous fitness: OLS (t-tests are valid)
#     - Binary fitness: Logistic GLM (Wald tests)
#
#   Quadratic coefficients are multiplied by 2 to obtain γ following
#   Lande & Arnold (1983): γ_ii = 2 × β_quad
#
#   Standardization is done WITHIN each group (e.g., year) to ensure
#   individuals are compared relative to their relevant context.
#
# Parameters:
#   data         : data frame with fitness and trait measurements
#   fitness_col  : name of the fitness column
#   trait_col    : name of the trait column (single trait)
#   fitness_type : "binary" or "continuous"
#   standardize  : if TRUE, standardize trait to mean 0, SD 1
#   group        : optional grouping variable (e.g., "year", "site")
#                  When specified, standardization and relative fitness are
#                  calculated separately within each group.
#
# Returns:
#   Data frame with:
#     Term              : trait name and "trait²"
#     Type              : "Linear" or "Quadratic"
#     Beta_Coefficient  : selection gradient (β or γ)
#     Standard_Error    : standard error of estimate
#     P_Value           : statistical significance (from appropriate source)
#     Variance          : squared standard error
# ============================================================================

#' Analyze disruptive/stabilizing selection
#'
#' Detects disruptive or stabilizing selection on a single trait by estimating linear and quadratic gradients.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_col A string specifying the name of the single trait column.
#' @param fitness_type A string indicating the fitness type: \code{"binary"} or \code{"continuous"}.
#' @param standardize Logical indicating whether to standardize the trait to mean 0 and SD 1. Default is \code{TRUE}.
#' @param group Optional string specifying a grouping variable.
#'
#' @return A data frame containing selection coefficients and statistics.
#' @export
analyze_disruptive_selection <- function(
  data,
  fitness_col,
  trait_col,
  fitness_type = c("binary", "continuous"),
  standardize = TRUE,
  group = NULL
) {
  # Input validation
  fitness_type <- match.arg(fitness_type)

  df <- prepare_selection_data(
    data = data,
    fitness_col = fitness_col,
    trait_cols = trait_col,
    standardize = standardize,
    group = group,
    add_relative = (fitness_type == "continuous"),
    na_action = "warn"
  )

  fml <- as.formula(
    paste(fitness_col, "~", trait_col, "+ I(", trait_col, "^2)")
  )

  # OLS
  fit_ols <- lm(fml, data = df)
  coef_ols <- summary(fit_ols)$coefficients

  quad_term <- paste0("I(", trait_col, "^2)")

  # Linear term (β) from OLS
  beta_linear <- coef_ols[trait_col, "Estimate"]
  se_linear <- coef_ols[trait_col, "Std. Error"]

  # Quadratic term (γ = 2 × β_quad) from OLS
  if (quad_term %in% rownames(coef_ols)) {
    gamma_quad <- 2 * coef_ols[quad_term, "Estimate"]
    se_quad <- 2 * coef_ols[quad_term, "Std. Error"]
  } else {
    gamma_quad <- NA_real_
    se_quad <- NA_real_
  }


  # GET P-VALUES BASED ON FITNESS TYPE
  if (fitness_type == "continuous") {
    # For continuous fitness:
    # OLS p-values (t-tests) are valid because residuals approximate normality
    p_linear <- coef_ols[trait_col, "Pr(>|t|)"]

    if (quad_term %in% rownames(coef_ols)) {
      p_quad <- coef_ols[quad_term, "Pr(>|t|)"]
    } else {
      p_quad <- NA_real_
    }
  } else {
    # For binary fitness (0/1 survival):
    # OLS p-values are NOT valid (residuals violate normality)
    # Use logistic GLM for valid p-values
    fit_glm <- glm(fml, data = df, family = binomial)
    coef_glm <- summary(fit_glm)$coefficients

    p_linear <- coef_glm[trait_col, "Pr(>|z|)"]

    if (quad_term %in% rownames(coef_glm)) {
      p_quad <- coef_glm[quad_term, "Pr(>|z|)"]
    } else {
      p_quad <- NA_real_
    }
  }

  results <- data.frame(
    Term = c(trait_col, paste0(trait_col, "²")),
    Type = c("Linear", "Quadratic"),
    Beta_Coefficient = c(beta_linear, gamma_quad),
    Standard_Error = c(se_linear, se_quad),
    P_Value = c(p_linear, p_quad),
    stringsAsFactors = FALSE
  )

  results$Variance <- results$Standard_Error^2

  return(results)
}
