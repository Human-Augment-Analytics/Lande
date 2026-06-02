# ======================================================
# test_birds.R
# Integration test using bird dataset (Darwin's finches)
# ======================================================

cat("\n========================================\n")
cat("MEDIUM GROUND FINCH DATASET FUNCTION TEST\n")
cat("========================================\n")

cat("Working directory:", getwd(), "\n")

# ------------------------------------------------------
# 1 Initialize environment
# ------------------------------------------------------

library(dplyr)

# Set a base directory resilient to being run from root or tests/
base_dir <- if (basename(getwd()) == "tests") ".." else "R"

if (file.exists(file.path(base_dir, "scripts", "0.0_initialize.R"))) {
    source(file.path(base_dir, "scripts", "0.0_initialize.R"))
}

# ======================================================
# 2 Load scripts
# ======================================================

cat("\nLoading script files...\n")

script_files <- list.files(file.path(base_dir, "scripts"),
    pattern = "\\.R$",
    full.names = TRUE)

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

fn_files <- list.files(file.path(base_dir, "functions"),
    pattern = "\\.R$",
    full.names = TRUE)

for (f in fn_files) {
    source(f)
    cat("Loaded:", basename(f), "\n")
}

# ======================================================
# 4 Load plotting
# ======================================================

cat("\nLoading plotting functions...\n")

plot_files <- list.files(file.path(base_dir, "plotting"),
    pattern = "\\.R$",
    full.names = TRUE)

for (f in plot_files) {
    source(f)
    cat("Loaded plot:", basename(f), "\n")
}

# ======================================================
# 5 Output directories
# ======================================================

rm(list = intersect(ls(), c("output_dir", "figure_dir", "table_dir", "model_dir")))

output_dir <- file.path(base_dir, "results", "medium_grouped_finch_results")
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
# 6 Load bird data
# ------------------------------------------------------

cat("\nLoading bird dataset...\n")

data_path <- file.path(base_dir, "data", "bird.data.RData")

# Load the data into a temporary environment and assign the first object to bird_data
# This prevents 'object not found' errors if the saved object is named differently
temp_env <- new.env()
load(data_path, envir = temp_env)
loaded_objs <- ls(temp_env)
bird_data <- temp_env[[loaded_objs[1]]]

# Remove duplicate X.y.* columns if present
bird_data <- bird_data[, !grepl("^X\\.y\\.", names(bird_data))]

cat("Rows:", nrow(bird_data), "\n")
cat("Columns:", ncol(bird_data), "\n")

# ------------------------------------------------------
# 7 Define traits and fitness
# ------------------------------------------------------

# Build Beak_PC1 from raw beak traits
beak_cols <- c("MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth")
if (all(beak_cols %in% names(bird_data))) {
    cat("\nBuilding Beak_PC1 via PCA...\n")

    beak_complete <- bird_data[, beak_cols]
    ok <- complete.cases(beak_complete)

    if (sum(ok) > 10) {
        pca_fit <- prcomp(beak_complete[ok, ], scale. = TRUE)
        pc1_scores <- rep(NA_real_, nrow(bird_data))
        pc1_scores[ok] <- pca_fit$x[, 1]

        # Ensure Beak_PC1 is positively correlated with depth
        cor_depth <- suppressWarnings(cor(pc1_scores[ok], beak_complete[ok, "MedianBeakDepth"]))
        if (!is.na(cor_depth) && cor_depth < 0) {
            pc1_scores <- -pc1_scores
        }

        bird_data$Beak_PC1 <- pc1_scores
        cat("  Added Beak_PC1\n")
    }
}

# Define traits for analysis
TRAITS <- c("Beak_PC1", "PC.body1")
TRAITS <- TRAITS[TRAITS %in% names(bird_data)]
FITNESS <- "Survived"

cat("\nTraits for analysis:", paste(TRAITS, collapse = ", "), "\n")

# ------------------------------------------------------
# 8 Data preparation function
# ------------------------------------------------------

prepare_survival_data <- function(year, last_year = 2011) {
    current_col <- paste0("y.", year)
    later_cols <- paste0("y.", (year + 1):last_year)

    needed_cols <- c(current_col, later_cols)
    missing_cols <- needed_cols[!needed_cols %in% names(bird_data)]
    if (length(missing_cols) > 0) {
        return(NULL)
    }

    survival_data <- bird_data %>%
        dplyr::filter(.data[[current_col]] == 1) %>%
        dplyr::mutate(
            Survived = as.numeric(rowSums(dplyr::across(dplyr::all_of(later_cols)), na.rm = TRUE) > 0),
            Year = year
        )

    available_traits <- TRAITS[TRAITS %in% names(survival_data)]
    id_col <- if ("BANDFINAL" %in% names(survival_data)) "BANDFINAL" else NULL

    survival_data <- survival_data %>%
        dplyr::select(dplyr::all_of(c("Year", "Survived", available_traits, id_col)))

    return(survival_data)
}

# ------------------------------------------------------
# 9 Check data availability across years
# ------------------------------------------------------

cat("\nChecking data availability across years...\n")

test_years <- 2003:2011
available_years <- c()

for (year in test_years) {
    test_data <- prepare_survival_data(year)
    if (!is.null(test_data) && nrow(test_data) >= 20) {
        available_years <- c(available_years, year)
        cat(
            "  Year", year, ": n =", nrow(test_data),
            ", survival =", round(mean(test_data$Survived), 3), "\n"
        )
    }
}

if (length(available_years) == 0) {
    stop("No years with sufficient data found.")
}

cat(
    "\nAnalyzing", length(available_years), "years:",
    paste(available_years, collapse = ", "), "\n"
)

# ------------------------------------------------------
# 10 Multi-year analysis
# ------------------------------------------------------

yearly_results <- list()
beak_trait <- TRAITS[1]

for (year in available_years) {
    cat("\n", strrep("=", 50), "\n")
    cat("YEAR:", year, "\n")
    cat(strrep("=", 50), "\n")

    year_data_raw <- prepare_survival_data(year)

    if (is.null(year_data_raw) || nrow(year_data_raw) < 20) {
        cat("  Skipped: insufficient sample size\n")
        next
    }

    year_data <- prepare_selection_data(
        data = year_data_raw,
        fitness_col = "Survived",
        trait_cols = TRAITS,
        standardize = TRUE,
        add_relative = FALSE,
        na_action = "warn"
    )

    cat("  Sample size:", nrow(year_data), "\n")
    cat("  Survival rate:", round(mean(year_data$Survived), 3), "\n")

    year_result <- list(
        year = year,
        n = nrow(year_data),
        survival_rate = mean(year_data$Survived)
    )

    # 10.1 selection_coefficients
    cat("\n  --- Selection Coefficients ---\n")
    tryCatch(
        {
            res <- selection_coefficients(
                data = year_data,
                fitness_col = "Survived",
                trait_cols = TRAITS,
                fitness_type = "binary",
                standardize = FALSE
            )
            year_result$selection <- res
            write.csv(res, file.path(table_dir, paste0("selection_", year, ".csv")), row.names = FALSE)

            quad_row <- res[res$Term == paste0(beak_trait, "²") & res$Type == "Quadratic", ]
            if (nrow(quad_row) > 0) {
                cat("    Gamma (", beak_trait, "^2):", round(quad_row$Beta_Coefficient, 4),
                    "(p =", round(quad_row$P_Value, 4), ")\n",
                    sep = ""
                )
            }
        },
        error = function(e) {
            cat("    FAILED:", e$message, "\n")
        }
    )

    # 10.2 analyze_disruptive_selection
    cat("\n  --- Disruptive Selection ---\n")
    tryCatch(
        {
            disrupt <- analyze_disruptive_selection(
                data = year_data,
                fitness_col = "Survived",
                trait_col = beak_trait,
                fitness_type = "binary",
                standardize = FALSE
            )
            year_result$disruptive <- disrupt
            write.csv(disrupt, file.path(table_dir, paste0("disruptive_", year, ".csv")), row.names = FALSE)

            gamma <- disrupt$Beta_Coefficient[2]
            p_val <- disrupt$P_Value[2]
            cat("Linear (Beta):", round(disrupt$Beta_Coefficient[1], 4), "\n")
            cat("Quadratic (Gamma):", round(gamma, 4), "\n")
            cat("P-value:", round(p_val, 4), "\n")

            if (gamma > 0 && p_val < 0.05) {
                cat("Significant DISRUPTIVE selection!\n")
            } else if (gamma < 0 && p_val < 0.05) {
                cat("Significant STABILIZING selection!\n")
            } else {
                cat("No significant nonlinear selection\n")
            }
        },
        error = function(e) {
            cat("    FAILED:", e$message, "\n")
        }
    )

    # 10.3 univariate_spline
    cat("\n  --- Univariate Spline ---\n")
    tryCatch(
        {
            spline <- univariate_spline(
                data = year_data,
                fitness_col = "Survived",
                trait_col = beak_trait,
                fitness_type = "binary"
            )
            year_result$spline <- spline
            saveRDS(spline, file.path(model_dir, paste0("spline_", year, ".rds")))

            if (exists("plot_univariate_fitness")) {
                p <- plot_univariate_fitness(
                    uni = spline,
                    trait_col = beak_trait,
                    title = paste("Fitness Function - Year", year)
                )
                ggsave(file.path(figure_dir, paste0("univariate_", year, ".png")), p, width = 7, height = 5, dpi = 300)
            }
        },
        error = function(e) {
            cat("FAILED:", e$message, "\n")
        }
    )

    # 10.4 correlated_fitness_surface
    if (length(TRAITS) >= 2) {
        cat("\n  --- Correlated Fitness Surface ---\n")
        tryCatch(
            {
                cfs <- correlated_fitness_surface(
                    data = year_data,
                    fitness_col = "Survived",
                    trait_cols = TRAITS[1:2],
                    grid_n = 50,
                    method = "gam"
                )
                year_result$cfs <- cfs
                saveRDS(cfs, file.path(model_dir, paste0("cfs_", year, ".rds")))

                if (exists("plot_correlated_fitness")) {
                    p <- plot_correlated_fitness(cfs, TRAITS[1:2]) + labs(title = paste("Fitness Surface - Year", year))
                    ggsave(file.path(figure_dir, paste0("cfs_", year, ".png")), p, width = 8, height = 6, dpi = 300)
                }
            },
            error = function(e) {
                cat("FAILED:", e$message, "\n")
            }
        )
    }

    # 10.5 adaptive landscape (for significant years)
    if (length(TRAITS) >= 2 && exists("adaptive_landscape")) {
        cat("\n  --- Adaptive Landscape ---\n")
        tryCatch(
            {
                # Fit GAM for adaptive landscape
                gam_model <- mgcv::gam(
                    Survived ~ s(Beak_PC1, PC.body1),
                    family = binomial,
                    data = year_data
                )

                landscape <- adaptive_landscape(
                    data = year_data,
                    fitness_model = gam_model,
                    trait_cols = c("Beak_PC1", "PC.body1"),
                    simulation_n = 500,
                    grid_n = 50
                )
                year_result$landscape <- landscape
                saveRDS(landscape, file.path(model_dir, paste0("landscape_", year, ".rds")))

                # ============================================
                # 2D Adaptive Landscape Plot
                # ============================================
                if (exists("plot_adaptive_landscape")) {
                    p2d <- plot_adaptive_landscape(
                        landscape = landscape,
                        trait_cols = c("Beak_PC1", "PC.body1"),
                        original_data = year_data,
                        bins = 12
                    )
                    ggsave(file.path(figure_dir, paste0("adaptive_landscape_", year, ".png")), p2d, width = 8, height = 6, dpi = 300)
                }

                # ============================================
                # 3D Adaptive Landscape Plot
                # ============================================
                if (exists("plot_adaptive_landscape_3d")) {
                    png(file.path(figure_dir, paste0("adaptive_landscape_3d_", year, ".png")),
                        width = 8, height = 6, units = "in", res = 300
                    )
                    plot_adaptive_landscape_3d(
                        landscape = landscape,
                        trait_cols = c("Beak_PC1", "PC.body1"),
                        theta = -30,
                        phi = 30,
                        grid_n = 200
                    )
                    dev.off()
                    cat("    3D plot saved\n")


                    angles <- list(
                        c(-30, 30),
                        c(0, 30),
                        c(-60, 20),
                        c(30, 40),
                        c(-90, 30)
                    )

                    for (i in seq_along(angles)) {
                        png(file.path(figure_dir, paste0("adaptive_landscape_3d_", year, "_angle", i, ".png")),
                            width = 8, height = 6, units = "in", res = 300
                        )
                        plot_adaptive_landscape_3d(
                            landscape = landscape,
                            trait_cols = c("Beak_PC1", "PC.body1"),
                            theta = angles[[i]][1],
                            phi = angles[[i]][2],
                            grid_n = 200
                        )
                        dev.off()
                    }
                    cat("    Multiple angles saved\n")
                }
            },
            error = function(e) {
                cat("FAILED:", e$message, "\n")
            }
        )
    }

    yearly_results[[as.character(year)]] <- year_result
}

# ------------------------------------------------------
# 11 Multi-year summary
# ------------------------------------------------------

cat("\n\n", strrep("=", 60), "\n")
cat("MULTI-YEAR SUMMARY\n")
cat(strrep("=", 60), "\n")

summary_df <- data.frame()

for (yr in names(yearly_results)) {
    res <- yearly_results[[yr]]
    if (!is.null(res$disruptive) && nrow(res$disruptive) >= 2) {
        gamma <- res$disruptive$Beta_Coefficient[2]
        p_val <- res$disruptive$P_Value[2]
        se <- res$disruptive$Standard_Error[2]

        summary_df <- rbind(summary_df, data.frame(
            Year = as.numeric(yr),
            N = res$n,
            Survival_Rate = round(res$survival_rate, 3),
            Gamma = round(gamma, 4),
            SE = round(se, 4),
            P_Value = round(p_val, 4),
            Significant = p_val < 0.05
        ))
        cat(
            "Year", yr, ": Gamma =", round(gamma, 4), "(p =", round(p_val, 4),
            ifelse(p_val < 0.05,
                ifelse(gamma > 0, "Disruptive", "Stabilizing"),
                "Not significant"
            ), "\n"
        )
    }
}

if (nrow(summary_df) > 0) {
    cat("\n=== Summary Table ===\n")
    print(summary_df)

    write.csv(summary_df, file.path(table_dir, "multi_year_summary.csv"), row.names = FALSE)
    cat("\nSummary saved to:", file.path(table_dir, "multi_year_summary.csv"), "\n")

    # Temporal plot
    if (nrow(summary_df) > 1) {
        p_temporal <- ggplot2::ggplot(summary_df, ggplot2::aes(x = Year, y = Gamma)) +
            ggplot2::geom_point(ggplot2::aes(size = N, color = Significant), alpha = 0.7) +
            ggplot2::geom_line(alpha = 0.5, color = "gray50") +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = Gamma - 1.96 * SE, ymax = Gamma + 1.96 * SE), width = 0.5, alpha = 0.5) +
            ggplot2::labs(
                title = "Temporal Variation in Disruptive Selection",
                subtitle = paste("Quadratic selection on", beak_trait),
                y = expression(paste("Quadratic Selection Coefficient (", gamma, ")")),
                x = "Year",
                color = "Significant (p < 0.05)",
                size = "Sample Size"
            ) +
            ggplot2::scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            ggplot2::theme_minimal(base_size = 12)

        print(p_temporal)
        ggplot2::ggsave(file.path(figure_dir, "temporal_disruptive_selection.png"), p_temporal, width = 10, height = 6, dpi = 300)
    }
} else {
    cat("\nNo disruptive selection results available.\n")
}

# ------------------------------------------------------
# 12 Final summary
# ------------------------------------------------------

cat("\n\n", strrep("=", 60), "\n")
cat("TESTING COMPLETE\n")
cat(strrep("=", 60), "\n")

cat("\nResults saved to:\n")
cat("  Tables:", table_dir, "\n")
cat("  Figures:", figure_dir, "\n")
cat("  Models:", model_dir, "\n")

cat("\nGenerated files:\n")
cat("\nTables:\n")
if (dir.exists(table_dir)) print(list.files(table_dir))
cat("\nFigures:\n")
if (dir.exists(figure_dir)) print(list.files(figure_dir))
cat("\nModels:\n")
if (dir.exists(model_dir)) print(list.files(model_dir))

cat("\n\n========================================\n")




cor(bird_data$Beak_PC1, bird_data$PC.body1, use = "complete.obs")


# ------------------------------------------------------
# 13 compare_fitness_surfaces (2009 only)
# ------------------------------------------------------

# Re-prepare the 2009 data
year_data_raw_2009 <- prepare_survival_data(2009)
year_2009 <- prepare_selection_data(
    data = year_data_raw_2009,
    fitness_col = "Survived",
    trait_cols = TRAITS,
    standardize = TRUE,
    add_relative = FALSE,
    na_action = "warn"
)

# Extract the correlated fitness surface for 2009 from the loop results
cfs_2009 <- yearly_results[["2009"]]$cfs

gam_model_2009 <- mgcv::gam(
    Survived ~ s(Beak_PC1, PC.body1),
    family = binomial,
    data = year_2009
)

landscape_2009 <- adaptive_landscape(
    data = year_2009,
    fitness_model = gam_model_2009,
    trait_cols = c("Beak_PC1", "PC.body1"),
    simulation_n = 500,
    grid_n = 50
)

if (exists("compare_fitness_surfaces_data") &&
    exists("plot_fitness_surfaces_comparison") &&
    exists("adaptive_landscape")) {
    cat("\n=== Testing Surface Comparison (2009) ===\n")

    trait_cols <- c("Beak_PC1", "PC.body1")

    comparison_data <- compare_fitness_surfaces_data(
        correlated_surface = cfs_2009,
        adaptive_landscape = landscape_2009,
        trait_cols = trait_cols
    )

    plots <- plot_fitness_surfaces_comparison(
        comparison_data = comparison_data,
        bins = 10,
        title = "2009 Comparison: Beak Size vs Body Size"
    )

    if (!is.null(plots$side_by_side)) {
        ggsave(
            file.path(figure_dir, "comparison_2009_side_by_side.png"),
            plots$side_by_side,
            width = 12,
            height = 5,
            dpi = 300
        )
        cat("  Side-by-side plot saved\n")
    }

    if (!is.null(plots$overlay)) {
        ggsave(
            file.path(figure_dir, "comparison_2009_overlay.png"),
            plots$overlay,
            width = 8,
            height = 6,
            dpi = 300
        )
        cat("  Overlay plot saved\n")
    }

    cor_df <- comparison_data$combined_data
    individual_fit <- cor_df$fitness[cor_df$type == "Correlated Fitness (Individual)"]
    population_fit <- cor_df$fitness[cor_df$type == "Adaptive Landscape (Population)"]

    cor_val <- cor(individual_fit, population_fit, use = "complete.obs")
    cat("\n  Correlation between surfaces:", round(cor_val, 4), "\n")

    write.csv(
        data.frame(
            Year = 2009,
            trait1 = trait_cols[1],
            trait2 = trait_cols[2],
            correlation = cor_val,
            n_points = length(individual_fit)
        ),
        file.path(table_dir, "comparison_correlation_2009.csv"),
        row.names = FALSE
    )

    cat("\n=== Comparison Summary (2009) ===\n")
    cat("Trait pair:", paste(trait_cols, collapse = " × "), "\n")
    cat("Correlation:", round(cor_val, 4), "\n")
} else {
    cat("\nRequired functions not found\n")
}
