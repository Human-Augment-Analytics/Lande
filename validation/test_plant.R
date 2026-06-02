# ======================================================
# test_plants.R
# Integration test using Aster plant dataset
# ======================================================

cat("\n========================================\n")
cat("ASTER PLANT DATASET FUNCTION TEST\n")
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
# 5 Output directories
# ======================================================

library(here)

output_dir <- here("R", "results", "plant_results")
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
# 6 Load Aster data
# ------------------------------------------------------

cat("\nLoading Aster dataset...\n")

data_dir <- here("R", "data")

data1 <- read.delim(file.path(data_dir, "Aster_analyses_2011_Cohort.txt"), sep = "\t")
data2 <- read.delim(file.path(data_dir, "Aster_analyses_2012_Cohort_full.txt"), sep = "\t")
data3 <- read.delim(file.path(data_dir, "Aster_analyses_2011_Cohort_full.txt"), sep = "\t")
data4 <- read.delim(file.path(data_dir, "Aster_analyses_2012_Cohort.txt"), sep = "\t")

cat("2011 Cohort:", dim(data1)[1], "rows,", dim(data1)[2], "columns\n")
cat("2012 Full:", dim(data2)[1], "rows,", dim(data2)[2], "columns\n")
cat("2011 Full:", dim(data3)[1], "rows,", dim(data3)[2], "columns\n")
cat("2012 Cohort:", dim(data4)[1], "rows,", dim(data4)[2], "columns\n")

data_2011_full <- data3
data_2011_fitness <- data1
data_2012_full <- data2
data_2012_fitness <- data4

# ------------------------------------------------------
# 7 Prepare 2011 data
# ------------------------------------------------------

cat("\nPreparing 2011 data...\n")

traits_2011 <- data_2011_full %>%
    dplyr::select(geno, elevation, season, SLA, deltaC13, FDsnow, height) %>%
    dplyr::distinct(geno, season, .keep_all = TRUE)

fitness_2011 <- data_2011_fitness %>%
    dplyr::select(geno, season, varb, resp) %>%
    dplyr::filter(!is.na(resp))

fitness_wide_2011 <- fitness_2011 %>%
    tidyr::pivot_wider(
        names_from = c(season, varb),
        names_sep = "_",
        names_prefix = "Year",
        values_from = resp,
        values_fn = mean
    )

analysis_2011 <- traits_2011 %>%
    dplyr::left_join(fitness_wide_2011, by = "geno") %>%
    dplyr::filter(!is.na(SLA), !is.na(deltaC13), !is.na(FDsnow), !is.na(height))

cat("2011 data: n =", nrow(analysis_2011), "\n")

# ------------------------------------------------------
# 8 Prepare 2012 data
# ------------------------------------------------------

cat("\nPreparing 2012 data...\n")

traits_2012 <- data_2012_full %>%
    dplyr::select(geno, elevation, season, SLA, deltaC, Fdsnow, height) %>%
    dplyr::rename(deltaC13 = deltaC, FDsnow = Fdsnow) %>%
    dplyr::distinct(geno, season, .keep_all = TRUE)

fitness_2012 <- data_2012_fitness %>%
    dplyr::select(geno, season, varb, resp) %>%
    dplyr::filter(!is.na(resp))

fitness_wide_2012 <- fitness_2012 %>%
    tidyr::pivot_wider(
        names_from = c(season, varb),
        names_sep = "_",
        names_prefix = "Year",
        values_from = resp,
        values_fn = mean
    )

analysis_2012 <- traits_2012 %>%
    dplyr::left_join(fitness_wide_2012, by = "geno") %>%
    dplyr::filter(!is.na(SLA), !is.na(deltaC13), !is.na(FDsnow), !is.na(height))

cat("2012 data: n =", nrow(analysis_2012), "\n")

# ------------------------------------------------------
# 9 Prepare selection data (standardize)
# ------------------------------------------------------

cat("\nPreparing selection data...\n")

fitness_2011_col <- grep("fecund", names(analysis_2011), value = TRUE)[1]
fitness_2012_col <- grep("fecund", names(analysis_2012), value = TRUE)[1]

traits <- c("SLA", "deltaC13", "FDsnow", "height")

prepared_2011 <- prepare_selection_data(
    data = analysis_2011,
    fitness_col = fitness_2011_col,
    trait_cols = traits,
    standardize = TRUE,
    add_relative = TRUE,
    na_action = "drop"
)

prepared_2012 <- prepare_selection_data(
    data = analysis_2012,
    fitness_col = fitness_2012_col,
    trait_cols = traits,
    standardize = TRUE,
    add_relative = TRUE,
    na_action = "drop"
)

cat("\nPrepared data:\n")
cat("  2011: n =", nrow(prepared_2011), "\n")
cat("  2012: n =", nrow(prepared_2012), "\n")

# ------------------------------------------------------
# 10 Selection Differential (Separate)
# ------------------------------------------------------

cat("\n=== Selection Differential (Separate) ===\n")

sel_diff_2011 <- list()
sel_diff_2012 <- list()

for (trait in traits) {
    sel_diff_2011[[trait]] <- tryCatch(
        selection_differential(
            data = prepared_2011,
            fitness_col = "relative_fitness",
            trait_col = trait,
            standardized = FALSE,
            use_relative = FALSE
        ),
        error = function(e) NA
    )

    sel_diff_2012[[trait]] <- tryCatch(
        selection_differential(
            data = prepared_2012,
            fitness_col = "relative_fitness",
            trait_col = trait,
            standardized = FALSE,
            use_relative = FALSE
        ),
        error = function(e) NA
    )
}

sel_diff_df <- data.frame(
    Trait = traits,
    S_2011 = unlist(sel_diff_2011),
    S_2012 = unlist(sel_diff_2012)
)

write.csv(sel_diff_df, file.path(table_dir, "selection_differentials.csv"), row.names = FALSE)


# ------------------------------------------------------
# 11 Selection Coefficients (Separate)
# ------------------------------------------------------

cat("\n=== Selection Coefficients (Separate) ===\n")

selection_2011 <- selection_coefficients(
    data = prepared_2011,
    fitness_col = "relative_fitness",
    trait_cols = traits,
    fitness_type = "continuous",
    standardize = TRUE
)

write.csv(selection_2011, file.path(table_dir, "selection_2011.csv"), row.names = FALSE)

cat("\n2011 Significant:\n")
sig_2011 <- selection_2011[selection_2011$P_Value < 0.05, ]
if (nrow(sig_2011) > 0) print(sig_2011) else cat("  None\n")


selection_2012 <- selection_coefficients(
    data = prepared_2012,
    fitness_col = "relative_fitness",
    trait_cols = traits,
    fitness_type = "continuous",
    standardize = TRUE
)

write.csv(selection_2012, file.path(table_dir, "selection_2012.csv"), row.names = FALSE)

cat("\n2012 Significant:\n")
sig_2012 <- selection_2012[selection_2012$P_Value < 0.05, ]
if (nrow(sig_2012) > 0) print(sig_2012) else cat("  None\n")

# ------------------------------------------------------
# 12 Combined Analysis (2011 + 2012)
# ------------------------------------------------------

cat("\n=== Combined Analysis (2011 + 2012) ===\n")

combined_data <- rbind(
    prepared_2011[, c("SLA", "deltaC13", "FDsnow", "height", "relative_fitness", "geno")],
    prepared_2012[, c("SLA", "deltaC13", "FDsnow", "height", "relative_fitness", "geno")]
)

cat("Combined sample size:", nrow(combined_data), "\n")
cat("  2011: n =", nrow(prepared_2011), "\n")
cat("  2012: n =", nrow(prepared_2012), "\n")

# 12.1 Selection Differential (Combined)
cat("\n--- Selection Differential (Combined) ---\n")

sel_diff_combined <- list()
for (trait in traits) {
    d <- tryCatch(
        selection_differential(
            data = combined_data,
            fitness_col = "relative_fitness",
            trait_col = trait,
            standardized = FALSE,
            use_relative = FALSE
        ),
        error = function(e) NA
    )
    sel_diff_combined[[trait]] <- d
    cat(trait, ": S =", round(d, 4), "\n")
}

sel_diff_combined_df <- data.frame(
    Trait = traits,
    S_Combined = unlist(sel_diff_combined)
)
write.csv(sel_diff_combined_df, file.path(table_dir, "selection_differentials_combined.csv"), row.names = FALSE)

# 12.2 Selection Coefficients (Combined)
cat("\n--- Selection Coefficients (Combined) ---\n")

selection_combined <- selection_coefficients(
    data = combined_data,
    fitness_col = "relative_fitness",
    trait_cols = traits,
    fitness_type = "continuous",
    standardize = TRUE
)

write.csv(selection_combined, file.path(table_dir, "selection_combined.csv"), row.names = FALSE)

cat("\nSelection Coefficients (Combined):\n")
print(selection_combined)


sig_combined <- selection_combined[selection_combined$P_Value < 0.05, ]
if (nrow(sig_combined) > 0) {
    cat("\n=== Significant Selection (Combined) ===\n")
    print(sig_combined)
} else {
    cat("\nNo significant selection detected (combined n =", nrow(combined_data), ")\n")
}

# ------------------------------------------------------
# 13 Comparison Summary
# ------------------------------------------------------

cat("\n=== Comparison: Separate vs Combined ===\n")

cat("\n2011 (n =", nrow(prepared_2011), "):\n")
if (nrow(sig_2011) > 0) print(sig_2011) else cat("  No significant selection\n")

cat("\n2012 (n =", nrow(prepared_2012), "):\n")
if (nrow(sig_2012) > 0) print(sig_2012) else cat("  No significant selection\n")

cat("\nCombined (n =", nrow(combined_data), "):\n")
if (nrow(sig_combined) > 0) print(sig_combined) else cat("  No significant selection\n")

# ------------------------------------------------------
# 14 Univariate spline
# ------------------------------------------------------

cat("\n=== Univariate Fitness Functions ===\n")

for (trait in traits) {
    cat("  ", trait, "...")

    spline <- tryCatch(
        {
            univariate_spline(
                data = combined_data,
                fitness_col = "relative_fitness",
                trait_col = trait,
                fitness_type = "continuous",
                k = 6
            )
        },
        error = function(e) NULL
    )

    if (!is.null(spline) && exists("plot_univariate_fitness")) {
        p <- plot_univariate_fitness(
            uni = spline,
            trait_col = trait,
            title = paste("Fitness Function:", trait)
        )
        ggsave(file.path(figure_dir, paste0("univariate_", trait, ".png")), p, width = 7, height = 5, dpi = 300)
    } else {
        cat("\n")
    }
}

# ------------------------------------------------------
# 15 Correlated fitness
# ------------------------------------------------------

cat("\n=== Correlated Fitness Surfaces ===\n")

trait_pairs <- combn(traits, 2, simplify = FALSE)

cfs_results <- list()

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")
    cat("\n---", pair_name, "---\n")

    tryCatch(
        {
            cfs <- correlated_fitness_surface(
                data = combined_data,
                fitness_col = "relative_fitness",
                trait_cols = pair,
                grid_n = 30,
                method = "tps"
            )

            cfs_results[[pair_name]] <- cfs

            saveRDS(cfs, file.path(model_dir, paste0("cfs_", pair_name, ".rds")))

            if (exists("plot_correlated_fitness")) {
                p <- plot_correlated_fitness(cfs, pair) +
                    labs(title = paste("Fitness Surface:", pair_name))
                ggsave(file.path(figure_dir, paste0("cfs_", pair_name, ".png")),
                    p,
                    width = 8, height = 6, dpi = 300
                )
                cat("  Basic plot saved\n")
            }

            if (exists("plot_correlated_fitness_enhanced")) {
                p_enhanced <- plot_correlated_fitness_enhanced(
                    cfs, pair,
                    original_data = combined_data,
                    fitness_col = "relative_fitness"
                )
                ggsave(file.path(figure_dir, paste0("cfs_", pair_name, "_enhanced.png")),
                    p_enhanced,
                    width = 8, height = 6, dpi = 300
                )
                cat("  Enhanced plot saved\n")
            }
        },
        error = function(e) {
            cat("FAILED:", e$message, "\n")
        }
    )
}

# ------------------------------------------------------
# 15 disruptive
# ------------------------------------------------------

# ------------------------------------------------------
# 17 Disruptive Selection (Univariate)
# ------------------------------------------------------

cat("\n=== Disruptive Selection ===\n")

disruptive_results <- list()

for (trait in traits) {
    cat("\n---", trait, "---\n")

    res <- tryCatch(
        {
            analyze_disruptive_selection(
                data = combined_data,
                fitness_col = "relative_fitness",
                trait_col = trait,
                fitness_type = "continuous",
                standardize = TRUE
            )
        },
        error = function(e) NULL
    )

    if (!is.null(res) && nrow(res) >= 2) {
        disruptive_results[[trait]] <- res

        gamma <- res$Beta_Coefficient[2]
        p_val <- res$P_Value[2]

        cat("  Linear (β):", round(res$Beta_Coefficient[1], 4), "\n")
        cat("  Quadratic (γ):", round(gamma, 4), "\n")
        cat("  P-value:", round(p_val, 4), "\n")

        if (p_val < 0.05) {
            if (gamma < 0) {
                cat("STABILIZING selection\n")
            } else {
                cat("DISRUPTIVE selection\n")
            }
        } else {
            cat("No significant nonlinear selection\n")
        }
    } else {
        cat("Failed\n")
    }
}

if (length(disruptive_results) > 0) {
    disruptive_summary <- data.frame(
        Trait = names(disruptive_results),
        Beta = sapply(disruptive_results, function(x) x$Beta_Coefficient[1]),
        Gamma = sapply(disruptive_results, function(x) x$Beta_Coefficient[2]),
        P_Linear = sapply(disruptive_results, function(x) x$P_Value[1]),
        P_Quadratic = sapply(disruptive_results, function(x) x$P_Value[2])
    )

    write.csv(disruptive_summary, file.path(table_dir, "disruptive_selection.csv"), row.names = FALSE)

    cat("\n=== Disruptive Selection Summary ===\n")
    print(disruptive_summary)
}


# ------------------------------------------------------
# 18 Adaptive Landscape
# ------------------------------------------------------

if (exists("adaptive_landscape") && length(cfs_results) > 0) {
    cat("\n=== Adaptive Landscape ===\n")

    for (pair_name in names(cfs_results)) {
        cat("\n---", pair_name, "---\n")

        trait_pair <- strsplit(pair_name, "_")[[1]]

        complete_cases <- sum(complete.cases(combined_data[, c("relative_fitness", trait_pair)]))
        if (complete_cases < 20) {
            cat("  SKIP: Insufficient complete cases (", complete_cases, ")\n", sep = "")
            next
        }

        tryCatch(
            {
                gam_model <- mgcv::gam(
                    as.formula(paste("relative_fitness", "~ s(", trait_pair[1], ",", trait_pair[2], ")")),
                    family = gaussian,
                    data = combined_data
                )

                landscape <- adaptive_landscape(
                    data = combined_data,
                    fitness_model = gam_model,
                    trait_cols = trait_pair,
                    simulation_n = 500,
                    grid_n = 50
                )

                saveRDS(landscape, file.path(model_dir, paste0("adaptive_landscape_", pair_name, ".rds")))

                # 2D plot
                if (exists("plot_adaptive_landscape")) {
                    p2d <- plot_adaptive_landscape(
                        landscape = landscape,
                        trait_cols = trait_pair,
                        original_data = combined_data,
                        bins = 12
                    )
                    ggsave(file.path(figure_dir, paste0("adaptive_landscape_", pair_name, "_2d.png")),
                        p2d,
                        width = 8, height = 6, dpi = 300
                    )
                    cat("  2D plot saved\n")
                }

                # 3D plot
                if (exists("plot_adaptive_landscape_3d")) {
                    png(file.path(figure_dir, paste0("adaptive_landscape_", pair_name, "_3d.png")),
                        width = 8, height = 6, units = "in", res = 300
                    )
                    plot_adaptive_landscape_3d(
                        landscape, trait_pair,
                        theta = -30, phi = 30, grid_n = 200
                    )
                    dev.off()
                    cat("  3D plot saved\n")
                }
            },
            error = function(e) {
                cat("FAILED:", e$message, "\n")
            }
        )
    }
}
