# ============================================================================
# selection_coefficients
#
# Purpose: Main wrapper function for selection analysis following Lande & Arnold (1983)
#
# Workflow:
#   1. Data preparation (standardization, relative fitness, optional grouping)
#   2. Automatic fitness type detection
#   3. Linear selection analysis (β)
#   4. Nonlinear selection analysis (γ, γ_ij)
#   5. Extract and combine all coefficients
#
# IMPORTANT NOTE:
#   - OLS is ALWAYS used to estimate selection gradients (β, γ, γ_ij)
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
#   use_relative_for_fit : if TRUE, use relative fitness for continuous data
#
# Returns:
#   Data frame with columns:
#     Term               : coefficient name (e.g., "size", "size²", "size×color")
#     Type               : "Linear", "Quadratic", or "Correlational"
#     Beta_Coefficient   : estimated selection gradient (β or γ)
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
#' @param fitness_type A string indicating the fitness type: \code{"auto"}, \code{"binary"}, or \code{"continuous"}.
#' @param standardize Logical indicating whether to standardize traits to mean 0 and SD 1. Default is \code{TRUE}.
#' @param group Optional string specifying a grouping variable (e.g., "year", "site").
#' @param use_relative_for_fit Logical indicating whether to use relative fitness for continuous data. Default is \code{TRUE}.
#' @param return_grouped Logical indicating whether to return results grouped if a \code{group} is specified.
#'
#' @return A data frame containing selection coefficients (Term, Type, Beta_Coefficient, Standard_Error, P_Value, Variance).
#' @export
#'
#' @examples
#' \dontrun{
#' coefs <- selection_coefficients(my_data, "fitness", c("trait1", "trait2"))
#' }
selection_coefficients <- function(data,
                                   fitness_col,
                                   trait_cols,
                                   fitness_type = c("auto", "binary", "continuous"),
                                   standardize = TRUE,
                                   group = NULL,
                                   use_relative_for_fit = TRUE,
                                   return_grouped = FALSE) {
  fitness_type <- match.arg(fitness_type)

  # ======================================================
  # CASE 1: return by group
  # ======================================================
  if (return_grouped && !is.null(group)) {
    groups <- unique(data[[group]])
    results_list <- list()

    for (g in groups) {
      data_g <- data[data[[group]] == g, ]

      res <- selection_coefficients(
        data = data_g,
        fitness_col = fitness_col,
        trait_cols = trait_cols,
        fitness_type = fitness_type,
        standardize = standardize,
        group = NULL,
        use_relative_for_fit = use_relative_for_fit,
        return_grouped = FALSE
      )

      res$Group <- g
      results_list[[as.character(g)]] <- res
    }

    all_results <- do.call(rbind, results_list)
    attr(all_results, "grouped") <- TRUE
    attr(all_results, "groups") <- groups
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

  # Determine which fitness column to model
  model_fitness_col <- if (fitness_type == "binary") {
    if (use_relative_for_fit) {
      message("Binary fitness detected: modeling on absolute 0/1")
    }
    fitness_col
  } else {
    if (use_relative_for_fit) {
      if (!rel_col %in% names(df)) {
        stop(
          "Relative fitness column '", rel_col, "' not found. ",
          "Ensure prepare_selection_data(add_relative=TRUE) creates it."
        )
      }
      rel_col
    } else {
      fitness_col
    }
  }

  # Run analyses
  linear_result <- analyze_linear_selection(
    data         = df,
    fitness_col  = model_fitness_col,
    trait_cols   = trait_cols,
    fitness_type = fitness_type
  )

  nonlinear_result <- analyze_nonlinear_selection(
    data         = df,
    fitness_col  = model_fitness_col,
    trait_cols   = trait_cols,
    fitness_type = fitness_type
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
  attr(all_coefs, "model_family_used") <- if (fitness_type == "binary") "binomial(logit)" else "gaussian"
  attr(all_coefs, "model_fitness_col") <- model_fitness_col
  attr(all_coefs, "relative_available") <- rel_col %in% names(df)
  attr(all_coefs, "group_used") <- group

  return(all_coefs)
}
