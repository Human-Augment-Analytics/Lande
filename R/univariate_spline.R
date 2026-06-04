# ======================================================
# univariate_spline.R
# Estimate univariate correlated fitness function
#
# Important Concept Explanation:
# This function calculates a UNIVARIATE CORRELATED FITNESS FUNCTION
# Definition: Individual fitness ~ Individual phenotype (single trait)
# Formula: w ~ f(z)
#
# This is a special case of correlated fitness surface (1D instead of 2D)
# Adaptive landscape would require: Mean fitness ~ Population mean phenotype
#
# IMPORTANT NOTES
#   - Traits should be standardized BEFORE calling this function
#   - Use prepare_selection_data() with standardize = TRUE and optional group
#   - The GAM smooth term estimates the shape of the fitness function
#   - DO NOT use scale() within this function (double standardization)
# ======================================================


#' Estimate univariate correlated fitness function
#'
#' This function calculates a univariate correlated fitness function using a GAM smooth term.
#' The formula is \code{w ~ f(z)}. Traits should be standardized BEFORE calling this function.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_col A string specifying the name of the trait column (must be numeric and standardized).
#' @param fitness_type A string indicating the fitness type: either \code{"binary"} or \code{"continuous"}.
#' @param group Optional string specifying a grouping variable. If provided, group fixed effects are included.
#' @param k Integer specifying the basis dimension for the GAM smooth term. Default is 10.
#'
#' @return A list of class \code{"univariate_fitness"} containing the fitted GAM model, a prediction grid, and metadata.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- univariate_spline(data = my_data, fitness_col = "survival", trait_col = "size")
#' }
univariate_spline <- function(data,
                              fitness_col,
                              trait_col,
                              fitness_type = c("binary", "continuous"),
                              group = NULL,
                              k = 10) {
  fitness_type <- match.arg(fitness_type)

  # Input validation
  if (length(trait_col) != 1L) {
    stop("`trait_col` must be a single column name.")
  }
  if (!trait_col %in% names(data)) {
    stop("`trait_col` not found in `data`.")
  }
  if (!is.numeric(data[[trait_col]])) {
    stop("`trait_col` must be numeric (standardize upstream if needed).")
  }

  # Check if group column exists
  if (!is.null(group) && !group %in% names(data)) {
    stop("Group column '", group, "' not found in data")
  }

  # CHECK FOR DOUBLE STANDARDIZATION
  # Warn users not to use scale() within this function
  message("IMPORTANT: Traits should already be standardized (mean = 0, SD = 1).")
  message("           Do NOT apply scale() again within this function.")

  # Check if trait appears to be standardized
  z_mean <- mean(data[[trait_col]], na.rm = TRUE)
  z_sd <- sd(data[[trait_col]], na.rm = TRUE)

  if (abs(z_mean) > 0.1 || abs(z_sd - 1) > 0.1) {
    warning(
      "Trait '", trait_col, "' does not appear standardized ",
      "(mean = ", round(z_mean, 3), ", SD = ", round(z_sd, 3), "). ",
      "Consider using prepare_selection_data() first."
    )
  } else {
    message("Trait appears standardized (mean ≈ 0, SD ≈ 1)")
  }

  if (fitness_type == "continuous") {
    # For continuous fitness, use relative fitness if available
    if ("relative_fitness" %in% names(data)) {
      y <- data[["relative_fitness"]]
      fit_note <- "Using relative_fitness column"
    } else {
      # Compute relative fitness on the fly (warning: not group-specific)
      if (!is.null(group)) {
        warning(
          "Group specified but no relative_fitness column found. ",
          "Consider using prepare_selection_data() first."
        )
      }
      y <- data[[fitness_col]] / mean(data[[fitness_col]], na.rm = TRUE)
      fit_note <- "Computed relative fitness on the fly (pooled)"
    }
    fam <- stats::gaussian()
    family_name <- "gaussian"
  } else {
    # Binary fitness - use original 0/1 values
    y <- data[[fitness_col]]
    fam <- stats::binomial("logit")
    family_name <- "binomial(logit)"

    # Check if really binary
    unique_vals <- unique(y[!is.na(y)])
    if (!all(unique_vals %in% c(0, 1))) {
      warning(
        "fitness_type = 'binary' but values are not all 0/1. ",
        "Proceeding but results may be unreliable."
      )
    }
    fit_note <- "Using original binary fitness"
  }

  df <- data
  df[[".y"]] <- y

  # Adjust k based on unique values (avoid GAM errors)
  n_unique <- length(unique(df[[trait_col]][complete.cases(df[[trait_col]])]))
  n_obs <- sum(complete.cases(df[[trait_col]], y))

  # k should not exceed n_unique - 1
  k_adj <- min(k, max(3, n_unique - 1))
  if (k_adj < k) {
    warning(
      "Reducing k from ", k, " to ", k_adj,
      " (only ", n_unique, " unique values)"
    )
    k <- k_adj
  }

  # Check if sample size is sufficient
  if (n_obs < k * 2) {
    warning(
      "Sample size (", n_obs, ") may be insufficient for k = ", k,
      ". Consider reducing k."
    )
  }

  # Create formula with smooth term
  if (!is.null(group)) {
    fml <- stats::as.formula(paste0(".y ~ ", group, " + s(", trait_col, ", k = ", k, ")"))
    cat("Including group fixed effect: '", group, "'\n")
  } else {
    fml <- stats::as.formula(paste0(".y ~ s(", trait_col, ", k = ", k, ")"))
  }

  # Fit GAM
  fit <- tryCatch(
    {
      mgcv::gam(fml,
        data = df,
        family = fam,
        method = "REML",
        na.action = stats::na.omit
      )
    },
    error = function(e) {
      stop(
        "GAM fitting failed: ", e$message,
        "\nTry reducing k (currently k = ", k, ")"
      )
    }
  )

  # Check convergence
  if (!is.null(fit$converged) && !fit$converged) {
    warning("GAM algorithm did not fully converge")
  }

  # Create prediction grid across observed trait range
  rng <- range(df[[trait_col]], na.rm = TRUE)
  grid <- data.frame(seq(rng[1], rng[2], length.out = 200))
  names(grid) <- trait_col

  # If group was specified, predictions need a reference group
  if (!is.null(group)) {
    group_levels <- unique(df[[group]])
    # Use the median group (or first) as reference
    # For factor groups, use the most common level
    if (is.factor(df[[group]])) {
      ref_group <- names(sort(table(df[[group]]), decreasing = TRUE))[1]
    } else {
      ref_group <- group_levels[which.min(abs(group_levels - median(group_levels)))]
    }
    grid[[group]] <- ref_group
    cat("Predictions use group = '", ref_group, "' as reference\n")
  }

  # Predict on link scale, then transform to response scale
  pr <- stats::predict(fit, newdata = grid, se.fit = TRUE, type = "link")
  linkinv <- fit$family$linkinv

  grid$fit <- linkinv(pr$fit)
  grid$lwr <- linkinv(pr$fit - 1.96 * pr$se.fit)
  grid$upr <- linkinv(pr$fit + 1.96 * pr$se.fit)

  # Ensure confidence bounds stay within [0,1] for binary fitness
  if (fitness_type == "binary") {
    grid$lwr <- pmax(grid$lwr, 0)
    grid$upr <- pmin(grid$upr, 1)
  }

  result <- list(
    model = fit,
    grid = grid,
    trait = trait_col,
    fitness_type = fitness_type,
    family = family_name,
    k = k,
    n_obs = n_obs,
    fit_note = fit_note,
    group_used = group,
    trait_mean = z_mean,
    trait_sd = z_sd,
    surface_type = "correlated_fitness_univariate",
    note = "Univariate correlated fitness function (individual fitness)"
  )

  class(result) <- "univariate_fitness"

  return(result)
}
