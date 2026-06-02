# ======================================================
# test_bumpus_sparrows.R
# Integration test using Bumpus sparrow dataset
# ======================================================

cat("\n========================================\n")
cat("BUMPUS SPARROW FUNCTION TEST\n")
cat("========================================\n")

cat("Working directory:", getwd(), "\n")

# ------------------------------------------------------
# 1 Initialize environment
# ------------------------------------------------------

if (file.exists("R/scripts/0.0_initialize.R")) {
    source("R/scripts/0.0_initialize.R")
}

# ======================================================
# 2 Load scripts
# ======================================================

cat("\nLoading script files...\n")

script_files <- list.files(
    "R/scripts",
    pattern = "\\.R$",
    full.names = TRUE
)

for (f in script_files) {
    if (basename(f) != "0.0_initialize.R") {
        source(f)
        cat("Loaded script:", basename(f), "\n")
    }
}

# ======================================================
# 3 Load functions
# ======================================================

cat("\nLoading function files...\n")

fn_files <- list.files(
    "R/functions",
    pattern = "\\.R$",
    full.names = TRUE
)

for (f in fn_files) {
    source(f)
    cat("Loaded:", basename(f), "\n")
}

# ======================================================
# 4 Load plotting
# ======================================================

cat("\nLoading plotting functions...\n")

plot_files <- list.files(
    "R/plotting",
    pattern = "\\.R$",
    full.names = TRUE
)

for (f in plot_files) {
    source(f)
    cat("Loaded plot:", basename(f), "\n")
}

# ======================================================
# 5 Output directories (defined AFTER scripts)
# ======================================================

library(here)

rm(list = intersect(ls(), c("output_dir", "figure_dir", "table_dir", "model_dir")))

output_dir <- here("R", "results", "bumpus_sparrows_results")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
model_dir <- file.path(output_dir, "models")

dirs <- c(output_dir, figure_dir, table_dir, model_dir)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

cat("\nResults saved to:\n")
cat("  Figures:", figure_dir, "\n")
cat("  Tables: ", table_dir, "\n")
cat("  Models: ", model_dir, "\n")


# ------------------------------------------------------
# 4 Load Bumpus data
# ------------------------------------------------------

cat("\nLoading Bumpus dataset...\n")

data_path <- here("R", "data", "Bumpus_data.csv")

if (!file.exists(data_path)) {
    stop("Cannot find data/Bumpus_data.csv")
}

bumpus <- read.csv(data_path)

bumpus$Survived <- as.numeric(bumpus$surv == "alive")
bumpus$TL <- bumpus$totlen
bumpus$HL <- bumpus$humer
bumpus$WT <- bumpus$wgt
bumpus$KL <- bumpus$stern

# Filter for adults to match original study (if AG column exists)
if ("AG" %in% names(bumpus)) {
    cat("Filtering for adults only (AG == 'adult')...\n")
    bumpus <- bumpus[bumpus$AG == "adult", ]
}

cat("Rows:", nrow(bumpus), "\n")
cat("Survival rate:", mean(bumpus$Survived), "\n")

# ------------------------------------------------------
# 5 Define traits
# ------------------------------------------------------

TRAITS <- c("TL", "HL", "WT", "KL")
FITNESS <- "Survived"

# ------------------------------------------------------
# 6 selection_coefficients
# ------------------------------------------------------

cat("\nTesting selection_coefficients...\n")

selection_results <- list()

for (k in seq_along(TRAITS)) {
    traits <- TRAITS[1:k]
    name <- paste(traits, collapse = "_")

    if (k == 1) {
        cat("\n---", name, "(skipped: single trait requires at least 2 traits) ---\n")
        next
    }

    tryCatch(
        {
            res <- selection_coefficients(
                data = bumpus,
                fitness_col = FITNESS,
                trait_cols = traits,
                fitness_type = "binary",
                standardize = TRUE
            )

            selection_results[[name]] <- res

            write.csv(
                res,
                file.path(table_dir, paste0("selection_", name, ".csv")),
                row.names = FALSE
            )

            sig <- res[res$P_Value < 0.05, ]

            if (nrow(sig) > 0) {
                cat("Significant terms:", paste(sig$Term, collapse = ", "), "\n")
            } else {
                cat("No significant terms detected\n")
            }
        },
        error = function(e) {
            cat("FAILED:", name, "-", e$message, "\n")
        }
    )
}

# ------------------------------------------------------
# 7 selection_differential
# ------------------------------------------------------

cat("\nTesting selection_differential...\n")

sel_diffs <- list()

for (trait in TRAITS) {
    tryCatch(
        {
            d <- selection_differential(
                data = bumpus,
                fitness_col = FITNESS,
                trait_col = trait,
                standardized = FALSE
            )

            sel_diffs[[trait]] <- d
            cat(trait, "=", round(d, 4), "\n")
        },
        error = function(e) {
            cat("FAILED:", trait, "-", e$message, "\n")
        }
    )
}

sel_df <- data.frame(
    trait = names(sel_diffs),
    differential = unlist(sel_diffs)
)

write.csv(
    sel_df,
    file.path(table_dir, "selection_differentials.csv"),
    row.names = FALSE
)

# ------------------------------------------------------
# 8 analyze_disruptive_selection
# ------------------------------------------------------

cat("\nTesting disruptive selection...\n")

disruptive_results <- list()

for (trait in TRAITS) {
    tryCatch(
        {
            res <- analyze_disruptive_selection(
                data = bumpus,
                fitness_col = FITNESS,
                trait_col = trait,
                fitness_type = "binary",
                standardize = TRUE
            )

            disruptive_results[[trait]] <- res

            write.csv(
                res,
                file.path(table_dir, paste0("disruptive_", trait, ".csv")),
                row.names = FALSE
            )
        },
        error = function(e) {
            cat("FAILED:", trait, "-", e$message, "\n")
        }
    )
}


if (length(disruptive_results) > 0) {
    summary_df <- data.frame(
        trait = names(disruptive_results),
        linear = sapply(disruptive_results, function(x) x$Beta_Coefficient[1]),
        quadratic = sapply(disruptive_results, function(x) x$Beta_Coefficient[2]),
        p_linear = sapply(disruptive_results, function(x) x$P_Value[1]),
        p_quadratic = sapply(disruptive_results, function(x) x$P_Value[2])
    )

    write.csv(
        summary_df,
        file.path(table_dir, "disruptive_selection_summary.csv"),
        row.names = FALSE
    )
}


# ------------------------------------------------------
# 9 univariate_spline + plotting
# ------------------------------------------------------
cat("\nTesting univariate spline + plotting...\n")

prepared_bumpus <- prepare_selection_data(
    data = bumpus,
    fitness_col = FITNESS,
    trait_cols = TRAITS,
    standardize = TRUE,
    group = NULL,
    add_relative = TRUE,
    na_action = "warn"
)

for (trait in TRAITS) {
    tryCatch(
        {
            spline <- univariate_spline(
                data = prepared_bumpus,
                fitness_col = FITNESS,
                trait_col = trait,
                fitness_type = "binary"
            )

            saveRDS(
                spline,
                file.path(model_dir, paste0("spline_", trait, ".rds"))
            )

            if (exists("plot_univariate_fitness")) {
                p <- plot_univariate_fitness(
                    uni = spline,
                    trait_col = trait,
                    title = paste("Survival vs", trait)
                )

                ggplot2::ggsave(
                    file.path(figure_dir, paste0("univariate_", trait, ".png")),
                    p,
                    width = 7,
                    height = 5,
                    dpi = 300
                )

                cat("Plot saved for:", trait, "\n")
            } else {
                cat("plot_univariate_fitness function not found\n")
            }
        },
        error = function(e) {
            cat("FAILED:", trait, "-", e$message, "\n")
        }
    )
}

# ------------------------------------------------------
# 9 correlated_fitness_surface + plots
# ------------------------------------------------------
cat("\nTesting correlated fitness surface + plotting...\n")

trait_pairs <- combn(TRAITS, 2, simplify = FALSE)

cfs_results <- list()
cfs_plots <- list()

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    tryCatch(
        {
            cfs <- correlated_fitness_surface(
                data = prepared_bumpus,
                fitness_col = FITNESS,
                trait_cols = pair,
                grid_n = 60,
                method = "auto",
                scale_traits = FALSE,
                group = NULL,
                k = 30
            )

            cfs_results[[pair_name]] <- cfs

            saveRDS(cfs, file.path(model_dir, paste0("cfs_", pair_name, ".rds")))

            write.csv(cfs$grid,
                file.path(table_dir, paste0("cfs_", pair_name, "_grid.csv")),
                row.names = FALSE
            )

            if (exists("plot_correlated_fitness")) {
                plot <- plot_correlated_fitness(
                    tps = cfs,
                    trait_cols = pair
                )

                cfs_plots[[pair_name]] <- plot

                ggplot2::ggsave(
                    file.path(figure_dir, paste0("cfs_", pair_name, ".png")),
                    plot,
                    width = 8,
                    height = 6,
                    dpi = 300
                )
            }

            if (exists("plot_correlated_fitness_enhanced")) {
                plot_enhanced <- plot_correlated_fitness_enhanced(
                    tps = cfs,
                    trait_cols = pair,
                    original_data = prepared_bumpus,
                    fitness_col = FITNESS
                )

                ggplot2::ggsave(
                    file.path(figure_dir, paste0("cfs_", pair_name, "_enhanced.png")),
                    plot_enhanced,
                    width = 10,
                    height = 8,
                    dpi = 300
                )
            }
        },
        error = function(e) {
            cat("FAILED:", trait, "-", e$message, "\n")
        }
    )
}

if (length(cfs_results) > 0) {
    summary_df <- data.frame(
        pair = names(cfs_results),
        method = sapply(cfs_results, function(x) x$method),
        data_type = sapply(cfs_results, function(x) x$data_type),
        grid_points = sapply(cfs_results, function(x) nrow(x$grid))
    )

    write.csv(
        summary_df,
        file.path(table_dir, "correlated_fitness_surfaces_summary.csv"),
        row.names = FALSE
    )
}

# ------------------------------------------------------
# 11 adaptive landscape
# ------------------------------------------------------
if (exists("adaptive_landscape")) {
    tryCatch(
        {
            trait_combinations <- list(
                c("TL", "HL"),
                c("TL", "WT"),
                c("TL", "KL"),
                c("HL", "WT"),
                c("HL", "KL"),
                c("WT", "KL")
            )

            for (traits in trait_combinations) {
                pair_name <- paste(traits, collapse = "_")

                formula_str <- paste("Survived ~ s(", traits[1], ", ", traits[2], ")")
                fitness_model <- mgcv::gam(
                    as.formula(formula_str),
                    family = binomial,
                    data = prepared_bumpus
                )

                landscape <- adaptive_landscape(
                    data = prepared_bumpus,
                    fitness_model = fitness_model,
                    trait_cols = traits,
                    simulation_n = 500,
                    grid_n = 60
                )

                saveRDS(landscape, file.path(model_dir, paste0("adaptive_landscape_", pair_name, ".rds")))

                if (exists("plot_adaptive_landscape")) {
                    points_data <- prepared_bumpus[, c("Survived", traits)]

                    p2d <- plot_adaptive_landscape(
                        landscape = landscape,
                        trait_cols = traits,
                        original_data = points_data,
                        bins = 12
                    )

                    ggplot2::ggsave(
                        file.path(figure_dir, paste0("adaptive_landscape_2d_", pair_name, ".png")),
                        p2d,
                        width = 8,
                        height = 6,
                        dpi = 300
                    )
                }

                # 3D plot
                if (exists("plot_adaptive_landscape_3d")) {
                    png(
                        file.path(figure_dir, paste0("adaptive_landscape_3d_", pair_name, ".png")),
                        width = 8,
                        height = 6,
                        units = "in",
                        res = 300
                    )

                    plot_adaptive_landscape_3d(
                        landscape = landscape,
                        trait_cols = traits,
                        theta = -30,
                        phi = 30,
                        grid_n = 200
                    )

                    dev.off()
                }

                # Save summary for this pair
                summary_df <- data.frame(
                    Trait1 = traits[1],
                    Trait2 = traits[2],
                    Optimal_Trait1 = round(landscape$optimum[[traits[1]]], 3),
                    Optimal_Trait2 = round(landscape$optimum[[traits[2]]], 3),
                    Optimal_Fitness = round(landscape$optimum$.mean_fit, 3)
                )

                # Append to master summary
                if (!exists("master_summary")) {
                    master_summary <- summary_df
                } else {
                    master_summary <- rbind(master_summary, summary_df)
                }
            }

            write.csv(master_summary, file.path(table_dir, "adaptive_landscape_summary_all.csv"), row.names = FALSE)

            # Create 3D plots from multiple angles for each pair
            if (exists("plot_adaptive_landscape_3d")) {
                angles <- list(
                    c(-30, 30), # default
                    c(0, 30), # front view
                    c(-60, 20), # side view
                    c(30, 40), # back-right view
                    c(-90, 30) # left side view
                )

                for (traits in trait_combinations) {
                    pair_name <- paste(traits, collapse = "_")
                    landscape <- readRDS(file.path(model_dir, paste0("adaptive_landscape_", pair_name, ".rds")))

                    for (i in seq_along(angles)) {
                        png(
                            file.path(figure_dir, paste0("adaptive_landscape_3d_", pair_name, "_angle", i, ".png")),
                            width = 8,
                            height = 6,
                            units = "in",
                            res = 300
                        )

                        plot_adaptive_landscape_3d(
                            landscape = landscape,
                            trait_cols = traits,
                            theta = angles[[i]][1],
                            phi = angles[[i]][2],
                            grid_n = 200
                        )

                        dev.off()
                    }
                }
            }
        },
        error = function(e) {}
    )
}

# ------------------------------------------------------
# 12 compare_fitness_surfaces
# ------------------------------------------------------

cat("\nTesting fitness surface comparison...\n")

if (length(cfs_results) > 0) {
    tryCatch(
        {
            trait_pairs <- names(cfs_results)
            comparison_results <- list()

            for (pair_name in trait_pairs) {
                trait_cols <- cfs_results[[pair_name]]$trait_cols

                # Load corresponding adaptive landscape
                landscape_file <- file.path(model_dir, paste0("adaptive_landscape_", pair_name, ".rds"))

                if (file.exists(landscape_file)) {
                    landscape <- readRDS(landscape_file)

                    comparison_data <- compare_fitness_surfaces_data(
                        correlated_surface = cfs_results[[pair_name]],
                        adaptive_landscape = landscape,
                        trait_cols = trait_cols
                    )

                    plots <- plot_fitness_surfaces_comparison(
                        comparison_data = comparison_data,
                        bins = 10,
                        title = paste("Comparison:", paste(trait_cols, collapse = " vs "))
                    )

                    # Save side-by-side plot
                    if (!is.null(plots$side_by_side)) {
                        ggplot2::ggsave(
                            file.path(figure_dir, paste0("comparison_", pair_name, "_side_by_side.png")),
                            plots$side_by_side,
                            width = 12,
                            height = 5,
                            dpi = 300
                        )
                    }

                    # Save overlay plot
                    if (!is.null(plots$overlay)) {
                        ggplot2::ggsave(
                            file.path(figure_dir, paste0("comparison_", pair_name, "_overlay.png")),
                            plots$overlay,
                            width = 8,
                            height = 6,
                            dpi = 300
                        )
                    }

                    # Calculate correlation between surfaces
                    cor_df <- comparison_data$combined_data
                    individual_fit <- cor_df$fitness[cor_df$type == "Correlated Fitness (Individual)"]
                    population_fit <- cor_df$fitness[cor_df$type == "Adaptive Landscape (Population)"]

                    cor_val <- cor(individual_fit, population_fit, use = "complete.obs")

                    comparison_results[[pair_name]] <- list(
                        correlation = cor_val,
                        comparison_data = comparison_data
                    )

                    # Save correlation
                    write.csv(
                        data.frame(
                            trait1 = trait_cols[1],
                            trait2 = trait_cols[2],
                            correlation = cor_val,
                            n_points = length(individual_fit)
                        ),
                        file.path(table_dir, paste0("comparison_correlation_", pair_name, ".csv")),
                        row.names = FALSE
                    )
                }
            }

            # Save master summary
            if (length(comparison_results) > 0) {
                summary_df <- do.call(rbind, lapply(names(comparison_results), function(name) {
                    data.frame(
                        Trait_Pair = name,
                        Correlation = round(comparison_results[[name]]$correlation, 4)
                    )
                }))

                write.csv(summary_df, file.path(table_dir, "comparison_summary.csv"), row.names = FALSE)
            }
        },
        error = function(e) {
            cat("FAILED: comparison -", e$message, "\n")
        }
    )
} else {
    cat("Skipping comparison (required functions or data not found)\n")
}
