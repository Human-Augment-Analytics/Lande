# ======================================================
# adaptive_landscape.R
# Calculate Adaptive Landscape (Sewall Wright's concept)
#
# IMPORTANT CONCEPT:
# Adaptive Landscape = Mean fitness ~ Population mean phenotype
# Formula: W̄ ~ z̄₁ + z̄₂
#
# This is DIFFERENT from correlated fitness surface:
# - Correlated fitness: w ~ z (individual fitness)
# - Adaptive landscape: W̄ ~ z̄ (population mean fitness)
#
# KEY PRINCIPLE:
# - Traits MUST already be standardized (mean = 0, SD = 1)
# - Use prepare_selection_data() with standardize = TRUE before calling
# - DO NOT standardize again within this function
#
# Requires:
# - A fitted fitness model (from GAM or TPS) that predicts individual fitness
# - Individual-level data to estimate within-population variance
# ======================================================

#' Calculate Adaptive Landscape
#'
#' Computes the adaptive landscape (mean fitness as a function of population mean phenotype)
#' using a fitted individual-level fitness model.
#'
#' @param data A data frame containing the original trait and fitness data.
#' @param fitness_model A fitted model object (e.g., GAM or Tps) predicting individual fitness.
#' @param trait_cols A character vector of length 2 specifying the trait column names.
#' @param group_col Optional character string specifying a grouping variable.
#' @param population_variance Optional covariance matrix for the traits. Estimated from data if \code{NULL}.
#' @param simulation_n Integer specifying the number of individuals to simulate per grid point. Default is 1000.
#' @param grid_n Integer specifying the resolution of the population mean grid. Default is 50.
#' @param custom_range Optional list specifying custom ranges for the traits.
#' @return An object of class \code{"adaptive_landscape"}.
#' @export
adaptive_landscape <- function(
  data,
  fitness_model,
  trait_cols,
  group_col = NULL,
  population_variance = NULL,
  simulation_n = 1000,
  grid_n = 50,
  custom_range = NULL
) {
    # Input validation
    stopifnot(length(trait_cols) == 2L)
    stopifnot(inherits(fitness_model, "gam") || inherits(fitness_model, "Tps"))

    # DOUBLE STANDARDIZATION WARNING
    message("IMPORTANT: Traits should already be standardized (mean = 0, SD = 1).")
    message("           Use prepare_selection_data() before calling this function.")
    message("           Do NOT standardize again within this function.")

    # Check if traits appear standardized
    for (t in trait_cols) {
        if (t %in% names(data)) {
            z_mean <- mean(data[[t]], na.rm = TRUE)
            z_sd <- sd(data[[t]], na.rm = TRUE)
            if (abs(z_mean) > 0.1 || abs(z_sd - 1) > 0.1) {
                warning(
                    "Trait '", t, "' does not appear standardized ",
                    "(mean = ", round(z_mean, 3), ", SD = ", round(z_sd, 3), "). ",
                    "Consider using prepare_selection_data() first."
                )
            }
        }
    }


    x1 <- data[[trait_cols[1]]]
    x2 <- data[[trait_cols[2]]]

    # Determine range for population mean grid
    if (!is.null(custom_range)) {
        x1_range <- custom_range[[trait_cols[1]]]
        x2_range <- custom_range[[trait_cols[2]]]
    } else {
        # Expand range slightly to include possible evolutionary space
        x1_range <- range(x1, na.rm = TRUE)
        x2_range <- range(x2, na.rm = TRUE)
        x1_range <- x1_range + c(-0.2, 0.2) * diff(x1_range)
        x2_range <- x2_range + c(-0.2, 0.2) * diff(x2_range)
    }

    cat("Population mean grid ranges:\n")
    cat("  ", trait_cols[1], ":", round(x1_range, 2), "\n")
    cat("  ", trait_cols[2], ":", round(x2_range, 2), "\n")

    # Create grid of population mean phenotypes
    g1 <- seq(x1_range[1], x1_range[2], length.out = grid_n)
    g2 <- seq(x2_range[1], x2_range[2], length.out = grid_n)

    population_grid <- expand.grid(g1, g2)
    names(population_grid) <- trait_cols

    # Estimate within-population variance if not provided
    if (is.null(population_variance)) {
        # Use all data to estimate phenotypic variance
        trait_data <- data[, trait_cols, drop = FALSE]
        trait_data <- trait_data[complete.cases(trait_data), ]

        if (nrow(trait_data) > 1) {
            population_variance <- var(trait_data, na.rm = TRUE)
        } else {
            stop("Not enough data to estimate population variance")
        }

        # If group_col provided, use average within-group variance
        if (!is.null(group_col) && group_col %in% names(data)) {
            group_variances <- by(
                trait_data,
                data[[group_col]][complete.cases(trait_data)],
                function(x) if (nrow(x) > 1) var(x) else NULL
            )
            group_variances <- group_variances[!sapply(group_variances, is.null)]
            if (length(group_variances) > 0) {
                population_variance <- Reduce(`+`, group_variances) / length(group_variances)
                cat("Using average within-group variance\n")
            }
        }

        cat("Estimated within-population variance-covariance:\n")
        print(round(population_variance, 4))
    }

    # Ensure variance matrix is positive definite
    if (inherits(try(chol(population_variance), silent = TRUE), "try-error")) {
        cat("Variance matrix not positive definite. Applying correction...\n")
        if (requireNamespace("Matrix", quietly = TRUE)) {
            population_variance <- as.matrix(Matrix::nearPD(population_variance)$mat)
        } else {
            # Simple correction: add small constant to diagonal
            population_variance <- population_variance + diag(1e-6, nrow(population_variance))
        }
    }

    mean_fitness <- numeric(nrow(population_grid))

    # Progress indicator
    cat("Calculating mean fitness for", nrow(population_grid), "grid points...\n")

    for (i in seq_len(nrow(population_grid))) {
        pop_mean <- population_grid[i, ]

        # Simulate individuals around this population mean
        simulated <- MASS::mvrnorm(
            n = simulation_n,
            mu = as.numeric(pop_mean),
            Sigma = population_variance
        )

        colnames(simulated) <- trait_cols
        sim_df <- as.data.frame(simulated)

        # Predict individual fitness
        if (inherits(fitness_model, "gam")) {
            ind_fitness <- predict(fitness_model, newdata = sim_df, type = "response")
        } else if (inherits(fitness_model, "Tps")) {
            ind_fitness <- predict(fitness_model, as.matrix(sim_df))
        }

        # Mean population fitness
        mean_fitness[i] <- mean(ind_fitness, na.rm = TRUE)
    }


    # Add mean fitness to grid
    population_grid$.mean_fit <- mean_fitness

    # Find optimum (maximum mean fitness)
    optimum <- population_grid[which.max(mean_fitness), ]
    cat("\nOptimal population mean phenotype:\n")
    print(optimum[, trait_cols, drop = FALSE])
    cat("Mean fitness at optimum:", round(optimum$.mean_fit, 4), "\n")

    # Calculate actual population means if group_col provided
    actual_means <- NULL
    if (!is.null(group_col) && group_col %in% names(data)) {
        actual_means <- aggregate(
            data[, trait_cols, drop = FALSE],
            by = list(group = data[[group_col]]),
            FUN = mean,
            na.rm = TRUE
        )
        names(actual_means)[1] <- group_col
        cat("\nActual population means:\n")
        print(actual_means)
    }

    result <- list(
        grid = population_grid,
        trait_cols = trait_cols,
        population_variance = population_variance,
        simulation_n = simulation_n,
        optimum = optimum,
        actual_population_means = actual_means,
        fitness_model_class = class(fitness_model)[1],
        data_summary = list(
            n_individuals = nrow(data),
            trait_ranges = list(
                x1 = range(x1, na.rm = TRUE),
                x2 = range(x2, na.rm = TRUE)
            ),
            traits_standardized = all(
                sapply(trait_cols, function(t) {
                    abs(mean(data[[t]], na.rm = TRUE)) < 0.1 &&
                        abs(sd(data[[t]], na.rm = TRUE) - 1) < 0.1
                })
            )
        ),
        surface_type = "adaptive_landscape",
        note = "Adaptive landscape: mean fitness as function of population mean phenotype"
    )

    class(result) <- "adaptive_landscape"
    return(result)
}


#' Print method for adaptive landscape
#'
#' @param x An object of class \code{"adaptive_landscape"}.
#' @param ... Additional arguments passed to \code{print}.
#' @export
print.adaptive_landscape <- function(x, ...) {
    cat("\nAdaptive Landscape Object\n")
    cat("========================\n")
    cat("Traits:", paste(x$trait_cols, collapse = " and "), "\n")
    cat("Grid size:", nrow(x$grid), "population means\n")
    cat("Simulations per point:", x$simulation_n, "individuals\n")
    cat("Fitness model:", x$fitness_model_class, "\n")

    # Check if traits were standardized
    if (!is.null(x$data_summary$traits_standardized)) {
        if (x$data_summary$traits_standardized) {
            cat("Traits: standardized (mean ≈ 0, SD ≈ 1)\n")
        } else {
            cat("WARNING: Traits may not be standardized!\n")
        }
    }

    cat("\nOptimal population mean phenotype:\n")
    print(x$optimum[, x$trait_cols, drop = FALSE])
    cat("Mean fitness at optimum:", round(x$optimum$.mean_fit, 4), "\n")

    if (!is.null(x$actual_population_means)) {
        cat("\nActual population means:\n")
        print(x$actual_population_means)
    }
    cat("\n", x$note, "\n")
}
