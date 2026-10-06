# ============================================================================
# selection_coefficients
#
# Purpose: Main wrapper function for selection analysis following Lande & Arnold (1983)
#
# Workflow:
#   1. Data preparation (standardization, relative fitness, optional grouping)
#   2. Automatic fitness type detection
#   3. Linear selection analysis (beta)
#   4. Nonlinear selection analysis (gamma, gamma_ij)
#   5. Extract and combine all coefficients
#
# IMPORTANT NOTE:
#   - OLS is always used to estimate selection gradients (beta, gamma, gamma_ij)
#   - For binary fitness: OLS gives gradients, logistic GLM gives p-values
#   - For continuous fitness: OLS gives both gradients and valid p-values
#   - Standardization and relative fitness can be done within groups (e.g., year)
#
# Parameters:
#   data                 : data frame with fitness and trait measurements
#   fitness_col          : name of the fitness column
#   trait_cols           : vector of trait column names
#   fitness_type         : "auto", "binary", or "continuous"
#   standardize          : if TRUE, standardize traits to mean 0, SD 1
#   group                : optional grouping variable (e.g., "year", "site")
#                          When specified, standardization and relative fitness
#                          are calculated separately within each group.
#   use_relative_for_fit : if TRUE (default), gradients are estimated on relative
#                          fitness W / mean(W) for every fitness type; FALSE
#                          gives coefficients on the absolute fitness scale
#
# Returns:
#   Data frame with columns:
#     Term               : coefficient name (e.g., "size", "size^2", "sizexcolor")
#     Type               : "Linear", "Quadratic", or "Correlational"
#     Beta_Coefficient   : estimated selection gradient (beta or gamma)
#     Standard_Error     : standard error of estimate
#     P_Value            : statistical significance
#     Variance           : square of standard error
# ============================================================================
#' Calculate selection coefficients
#'
#' Main wrapper function for selection analysis following Lande & Arnold (1983).
#' Extracts linear (beta), quadratic (gamma), and correlational selection gradients.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_cols A character vector of trait column names.
#' @param fitness_type A string indicating the fitness type: \code{"auto"}, \code{"binary"}, \code{"count"}, or \code{"continuous"}. Binary and count fitness take their p-values from a GLM (logistic, or Poisson and negative binomial) on the raw values; the gradients always come from OLS on relative fitness.
#' @param standardize Logical indicating whether to standardize traits to mean 0 and SD 1. Default is \code{TRUE}.
#' @param group Optional string specifying a grouping variable (e.g., "year", "site").
#'   Traits and relative fitness are standardised within each group, and for
#'   binary and count fitness the GLM that supplies the p-values gets a
#'   separate intercept for each group.
#' @param use_relative_for_fit Logical; if \code{TRUE} (default) the gradients are estimated on relative fitness \eqn{W / \bar{W}} for every fitness type, as in Lande & Arnold (1983). Set \code{FALSE} only to reproduce coefficients on the absolute fitness scale.
#' @param return_grouped Logical indicating whether to return results grouped if a \code{group} is specified.
#' @param se_type Standard errors of the gradients: \code{"ols"}, the usual
#'   least-squares ones (the default), or \code{"hc3"}, leave-one-out standard
#'   errors that allow for residual spread changing with the traits.
#'   Mitchell-Olds and Shaw (1987) suggested the jackknife for selection
#'   gradients when the residuals are not normal; HC3 (MacKinnon and White
#'   1985) gets much the same in closed form, from how far the estimates move
#'   when each individual is left out. In the package's simulation it recovered
#'   most of the coverage beta loses when the model leaves out curvature, from
#'   88 to 93\% with normal traits, 67 to 86\% with log-normal ones and 77 to
#'   91\% with heavy-tailed symmetric ones. It did nothing for the coverage lost
#'   to estimating a skewed trait's SD, and with survival and counts its
#'   intervals covered a little less, down to 91\% (see
#'   \code{check_selection_assumptions()}). For continuous fitness the p-values
#'   follow the chosen errors; for survival and counts they come from the GLM
#'   with either \code{se_type}.
#'
#' @return A data frame containing selection coefficients (Term, Type,
#'   Beta_Coefficient, Standard_Error, P_Value, Variance), with the standard
#'   errors used in the attribute \code{"se_type"}.
#' @references MacKinnon, J. G. and White, H. (1985) Some
#'   heteroskedasticity-consistent covariance matrix estimators with improved
#'   finite sample properties. Journal of Econometrics 29, 305-325.
#'   Mitchell-Olds, T. and Shaw, R. G. (1987) Regression analysis of natural
#'   selection: statistical inference and biological interpretation. Evolution
#'   41, 1149-1161.
#' @export
#'
#' @examples
#' selection_coefficients(bumpus, "survival", c("total_length", "weight"), fitness_type = "binary")
#'
#' # one set of gradients per sex, each sex standardised on its own
#' selection_coefficients(bumpus, "survival", c("total_length", "weight"),
#'                        group = "sex", return_grouped = TRUE)
selection_coefficients <- function(data,
                                   fitness_col,
                                   trait_cols,
                                   fitness_type = c("auto", "binary", "count", "continuous"),
                                   standardize = TRUE,
                                   group = NULL,
                                   use_relative_for_fit = TRUE,
                                   return_grouped = FALSE,
                                   se_type = c("ols", "hc3")) {
  fitness_type <- match.arg(fitness_type)
  se_type <- .se_type_arg(se_type)

  # ======================================================
  # CASE 1: return by group
  # ======================================================
  if (return_grouped && !is.null(group)) {
    if (!group %in% names(data)) {
      stop("Group column '", group, "' not found in data")
    }
    groups <- unique(data[[group]])
    if (anyNA(groups)) {
      warning(
        sum(is.na(data[[group]])), " row(s) have a missing '", group,
        "' label and are excluded from the per-group results"
      )
      groups <- groups[!is.na(groups)]
    }
    results_list <- list()

    for (g in groups) {
      data_g <- data[!is.na(data[[group]]) & data[[group]] == g, , drop = FALSE]

      res <- selection_coefficients(
        data = data_g,
        fitness_col = fitness_col,
        trait_cols = trait_cols,
        fitness_type = fitness_type,
        standardize = standardize,
        group = NULL,
        use_relative_for_fit = use_relative_for_fit,
        return_grouped = FALSE,
        se_type = se_type
      )

      res$Group <- g
      results_list[[as.character(g)]] <- res
    }

    all_results <- do.call(rbind, results_list)
    attr(all_results, "grouped") <- TRUE
    attr(all_results, "groups") <- groups
    attr(all_results, "se_type") <- se_type
    return(all_results)
  }

  # ======================================================
  # CASE 2: group all together
  # ======================================================

  # Determine relative fitness column name
  rel_col <- paste0(fitness_col, "_relative")

  df <- prepare_selection_data(
    data          = data,
    fitness_col   = fitness_col,
    trait_cols    = trait_cols,
    standardize   = standardize,
    group         = group,
    add_relative  = TRUE,
    na_action     = "warn",
    name_relative = rel_col
  )

  # Detect fitness type if auto
  det <- detect_family(df[[fitness_col]])
  if (fitness_type == "auto") {
    fitness_type <- det$type
  }

  # Selection gradients come from OLS on relative fitness (or on absolute
  # fitness if use_relative_for_fit = FALSE). Binary and count fitness
  # also use the raw column for the GLM that supplies p-values.
  ols_response_col <- if (use_relative_for_fit) {
    if (!rel_col %in% names(df)) {
      stop(
        "Relative fitness column '", rel_col, "' not found; ",
        "prepare_selection_data(add_relative = TRUE) adds it."
      )
    }
    rel_col
  } else {
    fitness_col
  }

  binary_response_col <- if (fitness_type %in% c("binary", "count")) fitness_col else NULL

  # Run analyses
  linear_result <- analyze_linear_selection(
    data                = df,
    fitness_col         = ols_response_col,
    trait_cols          = trait_cols,
    fitness_type        = fitness_type,
    binary_response_col = binary_response_col,
    group               = group,
    se_type             = se_type
  )

  nonlinear_result <- analyze_nonlinear_selection(
    data                = df,
    fitness_col         = ols_response_col,
    trait_cols          = trait_cols,
    fitness_type        = fitness_type,
    binary_response_col = binary_response_col,
    group               = group,
    se_type             = se_type
  )

  # Extract coefficients
  linear_coefs <- extract_linear_coefficients(trait_cols, linear_result)
  quadratic_coefs <- extract_quadratic_coefficients(trait_cols, nonlinear_result)
  interaction_coefs <- extract_interaction_coefficients(trait_cols, nonlinear_result)

  # Combine all coefficients
  all_coefs <- rbind(linear_coefs, quadratic_coefs, interaction_coefs)

  # Compute variance
  all_coefs$Variance <- all_coefs$Standard_Error^2

  # Add attributes
  attr(all_coefs, "fitness_type_detected") <- det$type
  attr(all_coefs, "fitness_type_used") <- fitness_type
  attr(all_coefs, "model_family_used") <- linear_result$glm_family %||% "gaussian"
  # the quadratic model can switch to a negative binomial on its own
  attr(all_coefs, "model_family_quadratic") <- nonlinear_result$glm_family %||% "gaussian"
  attr(all_coefs, "model_fitness_col") <- ols_response_col
  attr(all_coefs, "relative_available") <- rel_col %in% names(df)
  attr(all_coefs, "group_used") <- group
  attr(all_coefs, "se_type") <- se_type

  return(all_coefs)
}
