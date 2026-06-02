# ======================================================
# test_fish.R
# Integration test using pupfish datasets (Crescent Pond + Little Lake)
# ======================================================

rm(list = ls())
gc()

cat("\n========================================\n")
cat("PUPFISH DATASET FUNCTION TEST\n")
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

output_dir <- here("R", "results", "pupfish_results")
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
# 6 Load pupfish datasets
# ------------------------------------------------------
cat("\n=== Loading Pupfish Datasets ===\n")

data_dirs <- c(here("R", "data"), here("R", "test_data"))

# ======================================================
# 6.1 Crescent Pond data (main dataset)
# ======================================================
crescent_files <- list.files(data_dirs,
    pattern = "Crescent.*Pond.*\\.csv$",
    full.names = TRUE, recursive = TRUE
)

crescent_data <- read.csv(crescent_files[1])
cat("\n--- Crescent Pond ---\n")
cat("File:", basename(crescent_files[1]), "\n")
cat("Rows:", nrow(crescent_data), "\n")
cat("Columns:", ncol(crescent_data), "\n")

# ======================================================
# 6.2 Little Lake data (comparison dataset)
# ======================================================
little_files <- list.files(data_dirs,
    pattern = "Little.*Lake.*\\.csv$",
    full.names = TRUE, recursive = TRUE
)

little_data <- read.csv(little_files[1])
cat("\n--- Little Lake ---\n")
cat("File:", basename(little_files[1]), "\n")
cat("Rows:", nrow(little_data), "\n")
cat("Columns:", ncol(little_data), "\n")

# ======================================================
# 6.3 Combine datasets (for lake comparison)
# ======================================================

crescent_data$lake <- "Crescent Pond"
little_data$lake <- "Little Lake"

all_pupfish_data <- rbind(crescent_data, little_data)

cat("\n=== Combined Dataset ===\n")
cat("Total rows:", nrow(all_pupfish_data), "\n")
cat("Crescent Pond:", sum(all_pupfish_data$lake == "Crescent Pond"), "\n")
cat("Little Lake:", sum(all_pupfish_data$lake == "Little Lake"), "\n")
cat("Lakes:", paste(unique(all_pupfish_data$lake), collapse = ", "), "\n")

main_data <- all_pupfish_data

# checking
cat("\nColumn names in combined data:\n")
print(names(main_data))

# ------------------------------------------------------
# 7 Define traits and fitness
# ------------------------------------------------------

FITNESS_BINARY <- "survival"
FITNESS_CONTINUOUS <- "ln.growth"
TRAITS <- c("jaw", "eye", "body", "nasal", "mouth", "SL")
GROUP <- "lake"

pupfish_data <- main_data %>%
    dplyr::select(dplyr::all_of(c(FITNESS_BINARY, FITNESS_CONTINUOUS, TRAITS, GROUP))) %>%
    dplyr::filter(complete.cases(.))

cat("\nComplete observations:", nrow(pupfish_data), "\n")
cat("Traits:", paste(TRAITS, collapse = ", "), "\n")
cat("Lakes:", paste(unique(pupfish_data[[GROUP]]), collapse = ", "), "\n")

# ------------------------------------------------------
# 8 Data preparation
# ------------------------------------------------------

prepared_binary <- prepare_selection_data(
    data = pupfish_data,
    fitness_col = FITNESS_BINARY,
    trait_cols = TRAITS,
    standardize = TRUE,
    group = NULL,
    add_relative = FALSE,
    na_action = "drop"
)

prepared_binary_group <- prepare_selection_data(
    data = pupfish_data,
    fitness_col = FITNESS_BINARY,
    trait_cols = TRAITS,
    standardize = TRUE,
    group = GROUP,
    add_relative = FALSE,
    na_action = "drop"
)

prepared_continuous <- prepare_selection_data(
    data = pupfish_data,
    fitness_col = FITNESS_CONTINUOUS,
    trait_cols = TRAITS,
    standardize = TRUE,
    group = NULL,
    add_relative = TRUE,
    na_action = "drop"
)

prepared_continuous_group <- prepare_selection_data(
    data = pupfish_data,
    fitness_col = FITNESS_CONTINUOUS,
    trait_cols = TRAITS,
    standardize = TRUE,
    group = GROUP,
    add_relative = TRUE,
    na_action = "drop"
)

cat("\n--- No grouping ---\n")
cat("Binary data (overall): n =", nrow(prepared_binary), "\n")
cat("Continuous data (overall): n =", nrow(prepared_continuous), "\n")

cat("\n--- With grouping by lake ---\n")
cat("Binary data (grouped): n =", nrow(prepared_binary_group), "\n")
cat("Continuous data (grouped): n =", nrow(prepared_continuous_group), "\n")

cat("\n--- Sample size by group ---\n")
cat("Binary grouped data:\n")
print(table(prepared_binary_group[[GROUP]]))

cat("\nContinuous grouped data:\n")
print(table(prepared_continuous_group[[GROUP]]))

# ======================================================
# 9 Selection Differential
# ======================================================

sel_diff_all <- list()


for (trait in TRAITS) {
    cat("\n--- Trait:", trait, "---\n")

    # ======================================================
    # Binary Fitness
    # ======================================================

    # Without group
    sel_diff_bin <- selection_differential(
        data = prepared_binary,
        fitness_col = FITNESS_BINARY,
        trait_col = trait,
        standardized = TRUE,
        use_relative = FALSE
    )
    cat("  Binary (no group):", round(sel_diff_bin, 4), "\n")

    # With group
    sel_diff_bin_grouped <- selection_differential(
        data = prepared_binary_group,
        fitness_col = FITNESS_BINARY,
        trait_col = trait,
        group = GROUP,
        return_grouped = TRUE,
        standardized = TRUE,
        use_relative = FALSE
    )

    cat("  Binary (by group):\n")
    print(sel_diff_bin_grouped)

    # With group (return weighted mean)
    sel_diff_bin_weighted <- selection_differential(
        data = prepared_binary_group,
        fitness_col = FITNESS_BINARY,
        trait_col = trait,
        group = GROUP,
        return_grouped = FALSE,
        standardized = TRUE,
        use_relative = FALSE
    )
    cat("  Binary (weighted mean):", round(sel_diff_bin_weighted, 4), "\n")

    # ======================================================
    # Continuous Fitness
    # ======================================================

    # Without group
    sel_diff_cont <- selection_differential(
        data = prepared_continuous,
        fitness_col = FITNESS_CONTINUOUS,
        trait_col = trait,
        standardized = TRUE,
        use_relative = TRUE
    )
    cat("  Continuous (no group):", round(sel_diff_cont, 4), "\n")

    # With group
    sel_diff_cont_grouped <- selection_differential(
        data = prepared_continuous_group,
        fitness_col = FITNESS_CONTINUOUS,
        trait_col = trait,
        group = GROUP,
        return_grouped = TRUE,
        standardized = TRUE,
        use_relative = TRUE
    )

    cat("  Continuous (by group):\n")
    print(sel_diff_cont_grouped)

    sel_diff_cont_weighted <- selection_differential(
        data = prepared_continuous_group,
        fitness_col = FITNESS_CONTINUOUS,
        trait_col = trait,
        group = GROUP,
        return_grouped = FALSE,
        standardized = TRUE,
        use_relative = TRUE
    )
    cat("  Continuous (weighted mean):", round(sel_diff_cont_weighted, 4), "\n")

    sel_diff_all[[trait]] <- list(
        binary = list(
            overall = sel_diff_bin,
            grouped = sel_diff_bin_grouped,
            weighted = sel_diff_bin_weighted
        ),
        continuous = list(
            overall = sel_diff_cont,
            grouped = sel_diff_cont_grouped,
            weighted = sel_diff_cont_weighted
        )
    )


    # Save grouped results (binary)
    write.csv(sel_diff_bin_grouped,
        file.path(table_dir, paste0("selection_differentials_binary_", trait, "_grouped.csv")),
        row.names = FALSE
    )

    # Save grouped results (continuous)
    write.csv(sel_diff_cont_grouped,
        file.path(table_dir, paste0("selection_differentials_continuous_", trait, "_grouped.csv")),
        row.names = FALSE
    )

    # Save overall results
    overall_bin_df <- data.frame(
        trait = trait,
        selection_differential = sel_diff_bin,
        fitness = "binary",
        method = "overall"
    )
    write.csv(overall_bin_df,
        file.path(table_dir, paste0("selection_differentials_binary_", trait, "_overall.csv")),
        row.names = FALSE
    )

    overall_cont_df <- data.frame(
        trait = trait,
        selection_differential = sel_diff_cont,
        fitness = "continuous",
        method = "overall"
    )
    write.csv(overall_cont_df,
        file.path(table_dir, paste0("selection_differentials_continuous_", trait, "_overall.csv")),
        row.names = FALSE
    )

    # Save weighted mean results
    weighted_bin_df <- data.frame(
        trait = trait,
        selection_differential = sel_diff_bin_weighted,
        fitness = "binary",
        method = "weighted_mean"
    )
    write.csv(weighted_bin_df,
        file.path(table_dir, paste0("selection_differentials_binary_", trait, "_weighted.csv")),
        row.names = FALSE
    )

    weighted_cont_df <- data.frame(
        trait = trait,
        selection_differential = sel_diff_cont_weighted,
        fitness = "continuous",
        method = "weighted_mean"
    )
    write.csv(weighted_cont_df,
        file.path(table_dir, paste0("selection_differentials_continuous_", trait, "_weighted.csv")),
        row.names = FALSE
    )
}


cat("\n=== Selection Differential Summary ===\n")

summary_all <- data.frame()

for (trait in TRAITS) {
    res <- sel_diff_all[[trait]]

    # Binary overall
    summary_all <- rbind(summary_all, data.frame(
        Trait = trait,
        Fitness = "Binary",
        Analysis = "Overall",
        S = res$binary$overall,
        stringsAsFactors = FALSE
    ))

    # Binary weighted
    summary_all <- rbind(summary_all, data.frame(
        Trait = trait,
        Fitness = "Binary",
        Analysis = "Weighted Mean",
        S = res$binary$weighted,
        stringsAsFactors = FALSE
    ))

    # Binary by group
    if (is.data.frame(res$binary$grouped)) {
        for (i in 1:nrow(res$binary$grouped)) {
            summary_all <- rbind(summary_all, data.frame(
                Trait = trait,
                Fitness = "Binary",
                Analysis = res$binary$grouped[[GROUP]][i],
                S = res$binary$grouped$S[i],
                stringsAsFactors = FALSE
            ))
        }
    }

    # Continuous overall
    summary_all <- rbind(summary_all, data.frame(
        Trait = trait,
        Fitness = "Continuous",
        Analysis = "Overall",
        S = res$continuous$overall,
        stringsAsFactors = FALSE
    ))

    # Continuous weighted
    summary_all <- rbind(summary_all, data.frame(
        Trait = trait,
        Fitness = "Continuous",
        Analysis = "Weighted Mean",
        S = res$continuous$weighted,
        stringsAsFactors = FALSE
    ))

    # Continuous by group
    if (is.data.frame(res$continuous$grouped)) {
        for (i in 1:nrow(res$continuous$grouped)) {
            summary_all <- rbind(summary_all, data.frame(
                Trait = trait,
                Fitness = "Continuous",
                Analysis = res$continuous$grouped[[GROUP]][i],
                S = res$continuous$grouped$S[i],
                stringsAsFactors = FALSE
            ))
        }
    }
}


write.csv(summary_all, file.path(table_dir, "selection_differentials_all.csv"), row.names = FALSE)


# ======================================================
# 10 Selection Coefficients
# ======================================================

trait_pairs <- combn(TRAITS, 2, simplify = FALSE)

selection_coef_all <- list()

# ======================================================
# 10.1 Binary - Without group
# ======================================================
cat("\n--- Binary (Overall) ---\n")

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    result <- selection_coefficients(
        data = prepared_binary,
        fitness_col = FITNESS_BINARY,
        trait_cols = pair,
        fitness_type = "binary",
        standardize = TRUE
    )

    write.csv(result, file.path(table_dir, paste0("selection_binary_", pair_name, ".csv")), row.names = FALSE)

    selection_coef_all[[paste0("binary_", pair_name)]] <- result
}

# ======================================================
# 10.2 Continuous - Without group
# ======================================================
cat("\n--- Continuous (Overall) ---\n")

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    result <- selection_coefficients(
        data = prepared_continuous,
        fitness_col = FITNESS_CONTINUOUS,
        trait_cols = pair,
        fitness_type = "continuous",
        standardize = TRUE
    )

    write.csv(result, file.path(table_dir, paste0("selection_continuous_", pair_name, ".csv")), row.names = FALSE)

    selection_coef_all[[paste0("continuous_", pair_name)]] <- result
}

# ======================================================
# 10.3 Binary - With group
# ======================================================

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    result <- selection_coefficients(
        data = prepared_binary_group,
        fitness_col = FITNESS_BINARY,
        trait_cols = pair,
        fitness_type = "binary",
        standardize = TRUE,
        group = GROUP,
        return_grouped = TRUE
    )

    write.csv(result, file.path(table_dir, paste0("selection_binary_", pair_name, "_group.csv")), row.names = FALSE)

    selection_coef_all[[paste0("binary_", pair_name, "_group")]] <- result
}

# ======================================================
# 10.4 Continuous - With group
# ======================================================

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    result <- selection_coefficients(
        data = prepared_continuous_group,
        fitness_col = FITNESS_CONTINUOUS,
        trait_cols = pair,
        fitness_type = "continuous",
        standardize = TRUE,
        group = GROUP,
        return_grouped = TRUE
    )

    write.csv(result, file.path(table_dir, paste0("selection_continuous_", pair_name, "_group.csv")), row.names = FALSE)

    selection_coef_all[[paste0("continuous_", pair_name, "_group")]] <- result
}

# ======================================================
# 10.5 Selection Coefficients Summary Table
# ======================================================

summary_coef <- data.frame()

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    # Binary overall
    binary_key <- paste0("binary_", pair_name)
    if (!is.null(selection_coef_all[[binary_key]])) {
        df <- selection_coef_all[[binary_key]]

        beta_rows <- df[df$Type == "Linear", ]
        for (i in 1:nrow(beta_rows)) {
            summary_coef <- rbind(summary_coef, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Term = beta_rows$Term[i],
                Fitness = "Binary",
                Analysis = "Overall",
                Type = "Linear",
                Coefficient = round(beta_rows$Beta_Coefficient[i], 4),
                P_Value = beta_rows$P_Value[i],
                stringsAsFactors = FALSE
            ))
        }

        gamma_rows <- df[df$Type == "Quadratic", ]
        for (i in 1:nrow(gamma_rows)) {
            summary_coef <- rbind(summary_coef, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Term = gamma_rows$Term[i],
                Fitness = "Binary",
                Analysis = "Overall",
                Type = "Quadratic",
                Coefficient = round(gamma_rows$Beta_Coefficient[i], 4),
                P_Value = gamma_rows$P_Value[i],
                stringsAsFactors = FALSE
            ))
        }

        inter_rows <- df[df$Type == "Correlational", ]
        for (i in 1:nrow(inter_rows)) {
            summary_coef <- rbind(summary_coef, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Term = inter_rows$Term[i],
                Fitness = "Binary",
                Analysis = "Overall",
                Type = "Correlational",
                Coefficient = round(inter_rows$Beta_Coefficient[i], 4),
                P_Value = inter_rows$P_Value[i],
                stringsAsFactors = FALSE
            ))
        }
    }

    # Continuous overall
    continuous_key <- paste0("continuous_", pair_name)
    if (!is.null(selection_coef_all[[continuous_key]])) {
        df <- selection_coef_all[[continuous_key]]

        beta_rows <- df[df$Type == "Linear", ]
        for (i in 1:nrow(beta_rows)) {
            summary_coef <- rbind(summary_coef, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Term = beta_rows$Term[i],
                Fitness = "Continuous",
                Analysis = "Overall",
                Type = "Linear",
                Coefficient = round(beta_rows$Beta_Coefficient[i], 4),
                P_Value = beta_rows$P_Value[i],
                stringsAsFactors = FALSE
            ))
        }

        gamma_rows <- df[df$Type == "Quadratic", ]
        for (i in 1:nrow(gamma_rows)) {
            summary_coef <- rbind(summary_coef, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Term = gamma_rows$Term[i],
                Fitness = "Continuous",
                Analysis = "Overall",
                Type = "Quadratic",
                Coefficient = round(gamma_rows$Beta_Coefficient[i], 4),
                P_Value = gamma_rows$P_Value[i],
                stringsAsFactors = FALSE
            ))
        }

        inter_rows <- df[df$Type == "Correlational", ]
        for (i in 1:nrow(inter_rows)) {
            summary_coef <- rbind(summary_coef, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Term = inter_rows$Term[i],
                Fitness = "Continuous",
                Analysis = "Overall",
                Type = "Correlational",
                Coefficient = round(inter_rows$Beta_Coefficient[i], 4),
                P_Value = inter_rows$P_Value[i],
                stringsAsFactors = FALSE
            ))
        }
    }

    # Binary with group
    binary_group_key <- paste0("binary_", pair_name, "_group")
    if (!is.null(selection_coef_all[[binary_group_key]])) {
        df <- selection_coef_all[[binary_group_key]]
        linear_rows <- df[df$Type == "Linear", ]

        for (i in 1:nrow(linear_rows)) {
            group_name <- linear_rows[[GROUP]][i]
            summary_coef <- rbind(summary_coef, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Term = linear_rows$Term[i],
                Fitness = "Binary",
                Analysis = paste0("Group: ", group_name),
                Type = "Linear",
                Coefficient = round(linear_rows$Beta_Coefficient[i], 4),
                P_Value = linear_rows$P_Value[i],
                stringsAsFactors = FALSE
            ))
        }
    }

    # Continuous with group
    continuous_group_key <- paste0("continuous_", pair_name, "_group")
    if (!is.null(selection_coef_all[[continuous_group_key]])) {
        df <- selection_coef_all[[continuous_group_key]]
        linear_rows <- df[df$Type == "Linear", ]

        for (i in 1:nrow(linear_rows)) {
            group_name <- linear_rows[[GROUP]][i]
            summary_coef <- rbind(summary_coef, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Term = linear_rows$Term[i],
                Fitness = "Continuous",
                Analysis = paste0("Group: ", group_name),
                Type = "Linear",
                Coefficient = round(linear_rows$Beta_Coefficient[i], 4),
                P_Value = linear_rows$P_Value[i],
                stringsAsFactors = FALSE
            ))
        }
    }
}

write.csv(summary_coef, file.path(table_dir, "selection_coefficients_summary_all.csv"), row.names = FALSE)

# ======================================================
# 11 analyze_disruptive_selection
# ======================================================

disruptive_binary <- list()

for (trait in TRAITS) {
    res <- analyze_disruptive_selection(
        data = prepared_binary,
        fitness_col = FITNESS_BINARY,
        trait_col = trait,
        fitness_type = "binary",
        standardize = TRUE
    )

    disruptive_binary[[trait]] <- res

    write.csv(res, file.path(table_dir, paste0("disruptive_binary_", trait, ".csv")), row.names = FALSE)

    linear_row <- res[res$Type == "Linear", ]
    quadratic_row <- res[res$Type == "Quadratic", ]
}

# ======================================================
# Binary - With group
# ======================================================
disruptive_binary_group <- list()

for (trait in TRAITS) {
    res <- analyze_disruptive_selection(
        data = prepared_binary_group,
        fitness_col = FITNESS_BINARY,
        trait_col = trait,
        fitness_type = "binary",
        standardize = TRUE,
        group = GROUP,
        return_grouped = TRUE
    )

    disruptive_binary_group[[trait]] <- res

    write.csv(res, file.path(table_dir, paste0("disruptive_binary_", trait, "_group.csv")), row.names = FALSE)
}

# ======================================================
# Continuous - Without group
# ======================================================
disruptive_continuous <- list()

for (trait in TRAITS) {
    res <- analyze_disruptive_selection(
        data = prepared_continuous,
        fitness_col = FITNESS_CONTINUOUS,
        trait_col = trait,
        fitness_type = "continuous",
        standardize = TRUE
    )

    disruptive_continuous[[trait]] <- res

    write.csv(
        res,
        file.path(table_dir, paste0("disruptive_continuous_", trait, ".csv")),
        row.names = FALSE
    )

    linear_row <- res[res$Type == "Linear", ]
    quadratic_row <- res[res$Type == "Quadratic", ]
}

# ======================================================
# Continuous - With group
# ======================================================

disruptive_continuous_group <- list()

for (trait in TRAITS) {
    res <- analyze_disruptive_selection(
        data = prepared_continuous_group,
        fitness_col = FITNESS_CONTINUOUS,
        trait_col = trait,
        fitness_type = "continuous",
        standardize = TRUE,
        group = GROUP,
        return_grouped = TRUE
    )

    disruptive_continuous_group[[trait]] <- res

    write.csv(
        res,
        file.path(table_dir, paste0("disruptive_continuous_", trait, "_group.csv")),
        row.names = FALSE
    )
}

cat("\n=== Disruptive Selection Summary Table ===\n")

summary_df <- data.frame()

for (trait in TRAITS) {
    # Binary overall
    if (!is.null(disruptive_binary[[trait]])) {
        df <- disruptive_binary[[trait]]
        linear_row <- df[df$Type == "Linear", ]
        quadratic_row <- df[df$Type == "Quadratic", ]

        summary_df <- rbind(summary_df, data.frame(
            Trait = trait,
            Fitness = "Binary",
            Analysis = "Overall",
            Beta = round(linear_row$Beta_Coefficient, 4),
            P_Beta = linear_row$P_Value,
            Gamma = round(quadratic_row$Beta_Coefficient, 4),
            P_Gamma = quadratic_row$P_Value,
            Selection = ifelse(quadratic_row$P_Value < 0.05,
                ifelse(quadratic_row$Beta_Coefficient > 0, "Disruptive", "Stabilizing"),
                "None"
            ),
            stringsAsFactors = FALSE
        ))
    }

    # Continuous overall
    if (!is.null(disruptive_continuous[[trait]])) {
        df <- disruptive_continuous[[trait]]
        linear_row <- df[df$Type == "Linear", ]
        quadratic_row <- df[df$Type == "Quadratic", ]

        summary_df <- rbind(summary_df, data.frame(
            Trait = trait,
            Fitness = "Continuous",
            Analysis = "Overall",
            Beta = round(linear_row$Beta_Coefficient, 4),
            P_Beta = linear_row$P_Value,
            Gamma = round(quadratic_row$Beta_Coefficient, 4),
            P_Gamma = quadratic_row$P_Value,
            Selection = ifelse(quadratic_row$P_Value < 0.05,
                ifelse(quadratic_row$Beta_Coefficient > 0, "Disruptive", "Stabilizing"),
                "None"
            ),
            stringsAsFactors = FALSE
        ))
    }
}

write.csv(summary_df, file.path(table_dir, "disruptive_selection_summary.csv"), row.names = FALSE)


# ======================================================
# 12 Univariate Spline + Plot
# ======================================================

cat("\n=== Univariate Spline Analysis ===\n")

# Binary - Without group
cat("\n--- Binary (Overall) ---\n")
for (trait in TRAITS) {
    spline <- univariate_spline(
        data = prepared_binary,
        fitness_col = FITNESS_BINARY,
        trait_col = trait,
        fitness_type = "binary",
        group = NULL
    )
    saveRDS(spline, file.path(model_dir, paste0("spline_binary_", trait, ".rds")))

    if (exists("plot_univariate_fitness")) {
        p <- plot_univariate_fitness(uni = spline, trait_col = trait)
        ggplot2::ggsave(file.path(figure_dir, paste0("spline_binary_", trait, ".png")), p, width = 7, height = 5, dpi = 300)
    }
}

# Binary - With group
cat("\n--- Binary (By Lake) ---\n")
for (trait in TRAITS) {
    spline_list <- univariate_spline(
        data = prepared_binary_group,
        fitness_col = FITNESS_BINARY,
        trait_col = trait,
        fitness_type = "binary",
        group = GROUP
    )

    for (lake in names(spline_list)) {
        saveRDS(spline_list[[lake]], file.path(model_dir, paste0("spline_binary_", trait, "_", lake, ".rds")))

        if (exists("plot_univariate_fitness")) {
            p <- plot_univariate_fitness(uni = spline_list[[lake]], trait_col = trait)
            ggplot2::ggsave(file.path(figure_dir, paste0("spline_binary_", trait, "_", lake, ".png")), p, width = 7, height = 5, dpi = 300)
        }
    }
}

# Continuous - Without group
cat("\n--- Continuous (Overall) ---\n")
for (trait in TRAITS) {
    spline <- univariate_spline(
        data = prepared_continuous,
        fitness_col = FITNESS_CONTINUOUS,
        trait_col = trait,
        fitness_type = "continuous",
        group = NULL
    )
    saveRDS(spline, file.path(model_dir, paste0("spline_continuous_", trait, ".rds")))

    if (exists("plot_univariate_fitness")) {
        p <- plot_univariate_fitness(uni = spline, trait_col = trait)
        ggplot2::ggsave(file.path(figure_dir, paste0("spline_continuous_", trait, ".png")), p, width = 7, height = 5, dpi = 300)
    }
}

# Continuous - With group
cat("\n--- Continuous (By Lake) ---\n")
for (trait in TRAITS) {
    spline_list <- univariate_spline(
        data = prepared_continuous_group,
        fitness_col = FITNESS_CONTINUOUS,
        trait_col = trait,
        fitness_type = "continuous",
        group = GROUP
    )

    for (lake in names(spline_list)) {
        saveRDS(spline_list[[lake]], file.path(model_dir, paste0("spline_continuous_", trait, "_", lake, ".rds")))

        if (exists("plot_univariate_fitness")) {
            p <- plot_univariate_fitness(uni = spline_list[[lake]], trait_col = trait)
            ggplot2::ggsave(file.path(figure_dir, paste0("spline_continuous_", trait, "_", lake, ".png")), p, width = 7, height = 5, dpi = 300)
        }
    }
}

# ======================================================
# 13 Correlated Fitness Surfaces
# ======================================================

cat("\n=== Correlated Fitness Surfaces ===\n")

trait_pairs <- combn(TRAITS, 2, simplify = FALSE)

# Binary - Without group
cat("\n--- Binary (Overall) ---\n")
for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    cfs <- correlated_fitness_surface(
        data = prepared_binary,
        fitness_col = FITNESS_BINARY,
        trait_cols = pair,
        grid_n = 60,
        method = "auto",
        scale_traits = FALSE,
        group = NULL
    )

    saveRDS(cfs, file.path(model_dir, paste0("cfs_binary_", pair_name, ".rds")))

    if (exists("plot_correlated_fitness")) {
        plot <- plot_correlated_fitness(tps = cfs, trait_cols = pair)
        ggplot2::ggsave(file.path(figure_dir, paste0("cfs_binary_", pair_name, ".png")), plot, width = 8, height = 6, dpi = 300)
    }
}

# Binary - With group
cat("\n--- Binary (By Lake) ---\n")
for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    cfs_list <- correlated_fitness_surface(
        data = prepared_binary_group,
        fitness_col = FITNESS_BINARY,
        trait_cols = pair,
        grid_n = 60,
        method = "auto",
        scale_traits = FALSE,
        group = GROUP,
        k = 20
    )

    for (lake in names(cfs_list)) {
        saveRDS(cfs_list[[lake]], file.path(model_dir, paste0("cfs_binary_", pair_name, "_", lake, ".rds")))

        if (exists("plot_correlated_fitness")) {
            plot <- plot_correlated_fitness(tps = cfs_list[[lake]], trait_cols = pair)
            ggplot2::ggsave(file.path(figure_dir, paste0("cfs_binary_", pair_name, "_", lake, ".png")), plot, width = 8, height = 6, dpi = 300)
        }
    }
}

# Continuous - Without group
cat("\n--- Continuous (Overall) ---\n")
for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    cfs <- correlated_fitness_surface(
        data = prepared_continuous,
        fitness_col = FITNESS_CONTINUOUS,
        trait_cols = pair,
        grid_n = 60,
        method = "auto",
        scale_traits = FALSE,
        group = NULL
    )

    saveRDS(cfs, file.path(model_dir, paste0("cfs_continuous_", pair_name, ".rds")))

    if (exists("plot_correlated_fitness")) {
        plot <- plot_correlated_fitness(tps = cfs, trait_cols = pair)
        ggplot2::ggsave(file.path(figure_dir, paste0("cfs_continuous_", pair_name, ".png")), plot, width = 8, height = 6, dpi = 300)
    }
}

# Continuous - With group
cat("\n--- Continuous (By Lake) ---\n")
for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    cfs_list <- correlated_fitness_surface(
        data = prepared_continuous_group,
        fitness_col = FITNESS_CONTINUOUS,
        trait_cols = pair,
        grid_n = 60,
        method = "auto",
        scale_traits = FALSE,
        group = GROUP,
        k = 20
    )

    for (lake in names(cfs_list)) {
        saveRDS(cfs_list[[lake]], file.path(model_dir, paste0("cfs_continuous_", pair_name, "_", lake, ".rds")))

        if (exists("plot_correlated_fitness")) {
            plot <- plot_correlated_fitness(tps = cfs_list[[lake]], trait_cols = pair)
            ggplot2::ggsave(file.path(figure_dir, paste0("cfs_continuous_", pair_name, "_", lake, ".png")), plot, width = 8, height = 6, dpi = 300)
        }
    }
}

# ======================================================
# 11 Adaptive Landscape
# ======================================================

cat("\n=== Adaptive Landscape Analysis ===\n")

sites <- unique(prepared_binary_group[[GROUP]])
site_info <- data.frame(
    name = sites,
    label = gsub(" ", "_", sites), # "Crescent Pond" -> "Crescent_Pond"
    type = ifelse(grepl("Pond", sites), "pond", ifelse(grepl("Lake", sites), "lake", "site")),
    stringsAsFactors = FALSE
)

trait_pairs <- combn(TRAITS, 2, simplify = FALSE)

# ======================================================
# 11.1 Binary - Without group
# ======================================================

cat("\n--- Binary (Overall) ---\n")

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    # Fit GAM model for this trait pair
    formula_str <- paste(FITNESS_BINARY, "~ s(", pair[1], ", ", pair[2], ")")
    fitness_model <- mgcv::gam(
        as.formula(formula_str),
        family = binomial,
        data = prepared_binary
    )

    # Create adaptive landscape
    landscape <- adaptive_landscape(
        data = prepared_binary,
        fitness_model = fitness_model,
        trait_cols = pair,
        simulation_n = 500,
        grid_n = 60
    )

    saveRDS(landscape, file.path(model_dir, paste0("adaptive_landscape_binary_", pair_name, "_overall.rds")))

    # 2D plot
    if (exists("plot_adaptive_landscape")) {
        points_data <- prepared_binary[, c(FITNESS_BINARY, pair)]

        p2d <- plot_adaptive_landscape(
            landscape = landscape,
            trait_cols = pair,
            original_data = points_data,
            bins = 12
        )

        ggplot2::ggsave(
            file.path(figure_dir, paste0("adaptive_landscape_binary_", pair_name, "_overall_2d.png")),
            p2d,
            width = 8,
            height = 6,
            dpi = 300
        )
    }

    # 3D plot
    if (exists("plot_adaptive_landscape_3d")) {
        png(
            file.path(figure_dir, paste0("adaptive_landscape_binary_", pair_name, "_overall_3d.png")),
            width = 8,
            height = 6,
            units = "in",
            res = 300
        )

        plot_adaptive_landscape_3d(
            landscape = landscape,
            trait_cols = pair,
            theta = -30,
            phi = 30,
            grid_n = 200
        )

        dev.off()
    }
}

# ======================================================
# 11.2 Binary - With group
# ======================================================

cat("\n--- Binary (By Site) ---\n")

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    for (i in 1:nrow(site_info)) {
        site <- site_info$name[i]
        site_label <- site_info$label[i]
        site_type <- site_info$type[i]

        cat("  Processing:", site, "(", site_type, ")\n")

        site_data <- prepared_binary_group[prepared_binary_group[[GROUP]] == site, ]

        # Fit GAM for this site
        formula_str <- paste(FITNESS_BINARY, "~ s(", pair[1], ", ", pair[2], ")")
        fitness_model <- mgcv::gam(
            as.formula(formula_str),
            family = binomial,
            data = site_data
        )

        # Create adaptive landscape
        landscape <- adaptive_landscape(
            data = site_data,
            fitness_model = fitness_model,
            trait_cols = pair,
            simulation_n = 500,
            grid_n = 60
        )

        saveRDS(landscape, file.path(model_dir, paste0("adaptive_landscape_binary_", pair_name, "_", site_label, ".rds")))

        # 2D plot
        if (exists("plot_adaptive_landscape")) {
            points_data <- site_data[, c(FITNESS_BINARY, pair)]

            p2d <- plot_adaptive_landscape(
                landscape = landscape,
                trait_cols = pair,
                original_data = points_data,
                bins = 12
            )

            ggplot2::ggsave(
                file.path(figure_dir, paste0("adaptive_landscape_binary_", pair_name, "_", site_label, "_2d.png")),
                p2d,
                width = 8,
                height = 6,
                dpi = 300
            )
        }

        # 3D plot
        if (exists("plot_adaptive_landscape_3d")) {
            png(
                file.path(figure_dir, paste0("adaptive_landscape_binary_", pair_name, "_", site_label, "_3d.png")),
                width = 8,
                height = 6,
                units = "in",
                res = 300
            )

            plot_adaptive_landscape_3d(
                landscape = landscape,
                trait_cols = pair,
                theta = -30,
                phi = 30,
                grid_n = 200
            )

            dev.off()
        }
    }
}

# ======================================================
# 11.3 Continuous - Without group
# ======================================================

cat("\n--- Continuous (Overall) ---\n")

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    # Fit GAM model for this trait pair
    formula_str <- paste(FITNESS_CONTINUOUS, "~ s(", pair[1], ", ", pair[2], ")")
    fitness_model <- mgcv::gam(
        as.formula(formula_str),
        family = gaussian,
        data = prepared_continuous
    )

    # Create adaptive landscape
    landscape <- adaptive_landscape(
        data = prepared_continuous,
        fitness_model = fitness_model,
        trait_cols = pair,
        simulation_n = 500,
        grid_n = 60
    )

    saveRDS(landscape, file.path(model_dir, paste0("adaptive_landscape_continuous_", pair_name, "_overall.rds")))

    # 2D plot
    if (exists("plot_adaptive_landscape")) {
        points_data <- prepared_continuous[, c(FITNESS_CONTINUOUS, pair)]

        p2d <- plot_adaptive_landscape(
            landscape = landscape,
            trait_cols = pair,
            original_data = points_data,
            bins = 12
        )

        ggplot2::ggsave(
            file.path(figure_dir, paste0("adaptive_landscape_continuous_", pair_name, "_overall_2d.png")),
            p2d,
            width = 8,
            height = 6,
            dpi = 300
        )
    }

    # 3D plot
    if (exists("plot_adaptive_landscape_3d")) {
        png(
            file.path(figure_dir, paste0("adaptive_landscape_continuous_", pair_name, "_overall_3d.png")),
            width = 8,
            height = 6,
            units = "in",
            res = 300
        )

        plot_adaptive_landscape_3d(
            landscape = landscape,
            trait_cols = pair,
            theta = -30,
            phi = 30,
            grid_n = 200
        )

        dev.off()
    }
}

# ======================================================
# 11.4 Continuous - With group
# ======================================================

cat("\n--- Continuous (By Site) ---\n")

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    for (i in 1:nrow(site_info)) {
        site <- site_info$name[i]
        site_label <- site_info$label[i]
        site_type <- site_info$type[i]

        cat("  Processing:", site, "(", site_type, ")\n")

        site_data <- prepared_continuous_group[prepared_continuous_group[[GROUP]] == site, ]

        # Fit GAM for this site
        formula_str <- paste(FITNESS_CONTINUOUS, "~ s(", pair[1], ", ", pair[2], ")")
        fitness_model <- mgcv::gam(
            as.formula(formula_str),
            family = gaussian,
            data = site_data
        )

        # Create adaptive landscape
        landscape <- adaptive_landscape(
            data = site_data,
            fitness_model = fitness_model,
            trait_cols = pair,
            simulation_n = 500,
            grid_n = 60
        )

        saveRDS(landscape, file.path(model_dir, paste0("adaptive_landscape_continuous_", pair_name, "_", site_label, ".rds")))

        # 2D plot
        if (exists("plot_adaptive_landscape")) {
            points_data <- site_data[, c(FITNESS_CONTINUOUS, pair)]

            p2d <- plot_adaptive_landscape(
                landscape = landscape,
                trait_cols = pair,
                original_data = points_data,
                bins = 12
            )

            ggplot2::ggsave(
                file.path(figure_dir, paste0("adaptive_landscape_continuous_", pair_name, "_", site_label, "_2d.png")),
                p2d,
                width = 8,
                height = 6,
                dpi = 300
            )
        }

        # 3D plot
        if (exists("plot_adaptive_landscape_3d")) {
            png(
                file.path(figure_dir, paste0("adaptive_landscape_continuous_", pair_name, "_", site_label, "_3d.png")),
                width = 8,
                height = 6,
                units = "in",
                res = 300
            )

            plot_adaptive_landscape_3d(
                landscape = landscape,
                trait_cols = pair,
                theta = -30,
                phi = 30,
                grid_n = 200
            )

            dev.off()
        }
    }
}

# ======================================================
# 11.5 Summary Table
# ======================================================

cat("\n=== Adaptive Landscape Summary ===\n")

summary_landscape <- data.frame()

# Collect binary overall optima
for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")
    file_path <- file.path(model_dir, paste0("adaptive_landscape_binary_", pair_name, "_overall.rds"))

    if (file.exists(file_path)) {
        landscape <- readRDS(file_path)
        summary_landscape <- rbind(summary_landscape, data.frame(
            Trait1 = pair[1],
            Trait2 = pair[2],
            Fitness = "Binary",
            Site = "Overall",
            Site_Type = "overall",
            Optimum_Trait1 = round(landscape$optimum[[pair[1]]], 4),
            Optimum_Trait2 = round(landscape$optimum[[pair[2]]], 4),
            Max_Fitness = round(landscape$optimum$.mean_fit, 4),
            stringsAsFactors = FALSE
        ))
    }
}

# Collect continuous overall optima
for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")
    file_path <- file.path(model_dir, paste0("adaptive_landscape_continuous_", pair_name, "_overall.rds"))

    if (file.exists(file_path)) {
        landscape <- readRDS(file_path)
        summary_landscape <- rbind(summary_landscape, data.frame(
            Trait1 = pair[1],
            Trait2 = pair[2],
            Fitness = "Continuous",
            Site = "Overall",
            Site_Type = "overall",
            Optimum_Trait1 = round(landscape$optimum[[pair[1]]], 4),
            Optimum_Trait2 = round(landscape$optimum[[pair[2]]], 4),
            Max_Fitness = round(landscape$optimum$.mean_fit, 4),
            stringsAsFactors = FALSE
        ))
    }
}

# Collect binary by site optima
for (i in 1:nrow(site_info)) {
    site <- site_info$name[i]
    site_label <- site_info$label[i]
    site_type <- site_info$type[i]

    for (pair in trait_pairs) {
        pair_name <- paste(pair, collapse = "_")
        file_path <- file.path(model_dir, paste0("adaptive_landscape_binary_", pair_name, "_", site_label, ".rds"))

        if (file.exists(file_path)) {
            landscape <- readRDS(file_path)
            summary_landscape <- rbind(summary_landscape, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Fitness = "Binary",
                Site = site,
                Site_Type = site_type,
                Optimum_Trait1 = round(landscape$optimum[[pair[1]]], 4),
                Optimum_Trait2 = round(landscape$optimum[[pair[2]]], 4),
                Max_Fitness = round(landscape$optimum$.mean_fit, 4),
                stringsAsFactors = FALSE
            ))
        }
    }
}

# Collect continuous by site optima
for (i in 1:nrow(site_info)) {
    site <- site_info$name[i]
    site_label <- site_info$label[i]
    site_type <- site_info$type[i]

    for (pair in trait_pairs) {
        pair_name <- paste(pair, collapse = "_")
        file_path <- file.path(model_dir, paste0("adaptive_landscape_continuous_", pair_name, "_", site_label, ".rds"))

        if (file.exists(file_path)) {
            landscape <- readRDS(file_path)
            summary_landscape <- rbind(summary_landscape, data.frame(
                Trait1 = pair[1],
                Trait2 = pair[2],
                Fitness = "Continuous",
                Site = site,
                Site_Type = site_type,
                Optimum_Trait1 = round(landscape$optimum[[pair[1]]], 4),
                Optimum_Trait2 = round(landscape$optimum[[pair[2]]], 4),
                Max_Fitness = round(landscape$optimum$.mean_fit, 4),
                stringsAsFactors = FALSE
            ))
        }
    }
}

write.csv(summary_landscape, file.path(table_dir, "adaptive_landscape_summary.csv"), row.names = FALSE)

# ======================================================
# 12 Compare Fitness Surfaces
# ======================================================

cat("\n=== Compare Fitness Surfaces ===\n")

trait_pairs <- combn(TRAITS, 2, simplify = FALSE)

sites <- c("Crescent Pond", "Little Lake")
site_labels <- c("Crescent_Pond", "Little_Lake")

comparison_results <- list()

# ======================================================
# 12.1 Compare Binary Fitness Surfaces
# ======================================================

cat("\n--- Binary Fitness Surfaces ---\n")

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    cfs_overall_file <- file.path(model_dir, paste0("cfs_binary_", pair_name, ".rds"))
    if (!file.exists(cfs_overall_file)) {
        cat("  Skipping", pair_name, "- missing overall CFS\n")
        next
    }
    cfs_overall <- readRDS(cfs_overall_file)

    al_overall_file <- file.path(model_dir, paste0("adaptive_landscape_binary_", pair_name, "_overall.rds"))
    if (!file.exists(al_overall_file)) {
        cat("  Skipping", pair_name, "- missing overall AL\n")
        next
    }
    al_overall <- readRDS(al_overall_file)

    cat("  Comparing:", pair_name, "(Overall)\n")

    comparison_data_overall <- compare_fitness_surfaces_data(
        correlated_surface = cfs_overall,
        adaptive_landscape = al_overall,
        trait_cols = pair
    )

    plots_overall <- plot_fitness_surfaces_comparison(
        comparison_data = comparison_data_overall,
        bins = 10,
        title = paste("Fitness Surfaces Comparison -", pair_name, "(Overall)")
    )

    if (!is.null(plots_overall$side_by_side)) {
        ggplot2::ggsave(
            file.path(figure_dir, paste0("comparison_binary_", pair_name, "_overall_side_by_side.png")),
            plots_overall$side_by_side,
            width = 12,
            height = 5,
            dpi = 300
        )
    }

    if (!is.null(plots_overall$overlay)) {
        ggplot2::ggsave(
            file.path(figure_dir, paste0("comparison_binary_", pair_name, "_overall_overlay.png")),
            plots_overall$overlay,
            width = 8,
            height = 6,
            dpi = 300
        )
    }

    cor_df_overall <- comparison_data_overall$combined_data
    individual_fit_overall <- cor_df_overall$fitness[cor_df_overall$type == "Correlated Fitness (Individual)"]
    population_fit_overall <- cor_df_overall$fitness[cor_df_overall$type == "Adaptive Landscape (Population)"]
    cor_val_overall <- cor(individual_fit_overall, population_fit_overall, use = "complete.obs")

    write.csv(
        data.frame(
            trait1 = pair[1],
            trait2 = pair[2],
            fitness = "binary",
            comparison_type = "overall",
            site = "Overall",
            correlation = cor_val_overall,
            n_points = length(individual_fit_overall)
        ),
        file.path(table_dir, paste0("comparison_binary_", pair_name, "_overall_correlation.csv")),
        row.names = FALSE
    )

    comparison_results[[paste0("binary_", pair_name, "_overall")]] <- list(
        correlation = cor_val_overall,
        site = "Overall"
    )

    for (i in seq_along(sites)) {
        site <- sites[i]
        site_label <- site_labels[i]
        site_type <- ifelse(grepl("Pond", site), "pond", "lake")

        cfs_site_file <- file.path(model_dir, paste0("cfs_binary_", pair_name, "_", site, ".rds"))
        if (!file.exists(cfs_site_file)) {
            cat("    Skipping", site, "- missing CFS\n")
            next
        }
        cfs_site <- readRDS(cfs_site_file)

        al_site_file <- file.path(model_dir, paste0("adaptive_landscape_binary_", pair_name, "_", site_label, ".rds"))
        if (!file.exists(al_site_file)) {
            cat("    Skipping", site, "- missing AL\n")
            next
        }
        al_site <- readRDS(al_site_file)

        cat("    Comparing:", pair_name, "-", site, "(", site_type, ")\n")

        comparison_data_site <- compare_fitness_surfaces_data(
            correlated_surface = cfs_site,
            adaptive_landscape = al_site,
            trait_cols = pair
        )

        plots_site <- plot_fitness_surfaces_comparison(
            comparison_data = comparison_data_site,
            title = paste("Fitness Surfaces Comparison -", pair_name, "-", site)
        )

        plots_site$side_by_side <- plots_site$side_by_side +
            ggplot2::scale_colour_manual(values = rainbow(20)) +
            ggplot2::scale_fill_manual(values = rainbow(20))

        plots_site$overlay <- plots_site$overlay +
            ggplot2::scale_colour_manual(values = rainbow(20)) +
            ggplot2::scale_fill_manual(values = rainbow(20))


        if (!is.null(plots_site$side_by_side)) {
            ggplot2::ggsave(
                file.path(figure_dir, paste0("comparison_binary_", pair_name, "_", site_label, "_side_by_side.png")),
                plots_site$side_by_side,
                width = 12,
                height = 5,
                dpi = 300
            )
        }

        if (!is.null(plots_site$overlay)) {
            ggplot2::ggsave(
                file.path(figure_dir, paste0("comparison_binary_", pair_name, "_", site_label, "_overlay.png")),
                plots_site$overlay,
                width = 8,
                height = 6,
                dpi = 300
            )
        }

        cor_df_site <- comparison_data_site$combined_data
        individual_fit_site <- cor_df_site$fitness[cor_df_site$type == "Correlated Fitness (Individual)"]
        population_fit_site <- cor_df_site$fitness[cor_df_site$type == "Adaptive Landscape (Population)"]
        cor_val_site <- cor(individual_fit_site, population_fit_site, use = "complete.obs")

        write.csv(
            data.frame(
                trait1 = pair[1],
                trait2 = pair[2],
                fitness = "binary",
                comparison_type = "site",
                site = site,
                site_type = site_type,
                correlation = cor_val_site,
                n_points = length(individual_fit_site)
            ),
            file.path(table_dir, paste0("comparison_binary_", pair_name, "_", site_label, "_correlation.csv")),
            row.names = FALSE
        )

        comparison_results[[paste0("binary_", pair_name, "_", site_label)]] <- list(
            correlation = cor_val_site,
            site = site,
            site_type = site_type
        )
    }
}

# ======================================================
# 12.2 Compare Continuous Fitness Surfaces
# ======================================================

cat("\n--- Continuous Fitness Surfaces ---\n")

options(ggplot2.discrete.colour = function(...) scale_colour_manual(values = rainbow(20)))

for (pair in trait_pairs) {
    pair_name <- paste(pair, collapse = "_")

    cfs_overall_file <- file.path(model_dir, paste0("cfs_continuous_", pair_name, ".rds"))
    if (!file.exists(cfs_overall_file)) {
        cat("  Skipping", pair_name, "- missing overall CFS\n")
        next
    }
    cfs_overall <- readRDS(cfs_overall_file)

    al_overall_file <- file.path(model_dir, paste0("adaptive_landscape_continuous_", pair_name, "_overall.rds"))
    if (!file.exists(al_overall_file)) {
        cat("  Skipping", pair_name, "- missing overall AL\n")
        next
    }
    al_overall <- readRDS(al_overall_file)

    cat("  Comparing:", pair_name, "(Overall)\n")

    comparison_data_overall <- compare_fitness_surfaces_data(
        correlated_surface = cfs_overall,
        adaptive_landscape = al_overall,
        trait_cols = pair
    )

    plots_overall <- plot_fitness_surfaces_comparison(
        comparison_data = comparison_data_overall,
        title = paste("Fitness Surfaces Comparison -", pair_name, "(Overall - Continuous)")
    )

    if (!is.null(plots_overall$side_by_side)) {
        ggplot2::ggsave(
            file.path(figure_dir, paste0("comparison_continuous_", pair_name, "_overall_side_by_side.png")),
            plots_overall$side_by_side,
            width = 12,
            height = 5,
            dpi = 300
        )
    }

    if (!is.null(plots_overall$overlay)) {
        ggplot2::ggsave(
            file.path(figure_dir, paste0("comparison_continuous_", pair_name, "_overall_overlay.png")),
            plots_overall$overlay,
            width = 8,
            height = 6,
            dpi = 300
        )
    }

    cor_df_overall <- comparison_data_overall$combined_data
    individual_fit_overall <- cor_df_overall$fitness[cor_df_overall$type == "Correlated Fitness (Individual)"]
    population_fit_overall <- cor_df_overall$fitness[cor_df_overall$type == "Adaptive Landscape (Population)"]
    cor_val_overall <- cor(individual_fit_overall, population_fit_overall, use = "complete.obs")

    write.csv(
        data.frame(
            trait1 = pair[1],
            trait2 = pair[2],
            fitness = "continuous",
            comparison_type = "overall",
            site = "Overall",
            correlation = cor_val_overall,
            n_points = length(individual_fit_overall)
        ),
        file.path(table_dir, paste0("comparison_continuous_", pair_name, "_overall_correlation.csv")),
        row.names = FALSE
    )

    for (i in seq_along(sites)) {
        site <- sites[i]
        site_label <- site_labels[i]
        site_type <- ifelse(grepl("Pond", site), "pond", "lake")

        cfs_site_file <- file.path(model_dir, paste0("cfs_continuous_", pair_name, "_", site, ".rds"))
        if (!file.exists(cfs_site_file)) {
            cat("    Skipping", site, "- missing CFS\n")
            next
        }
        cfs_site <- readRDS(cfs_site_file)

        al_site_file <- file.path(model_dir, paste0("adaptive_landscape_continuous_", pair_name, "_", site_label, ".rds"))
        if (!file.exists(al_site_file)) {
            cat("    Skipping", site, "- missing AL\n")
            next
        }
        al_site <- readRDS(al_site_file)

        cat("    Comparing:", pair_name, "-", site, "(", site_type, ")\n")

        comparison_data_site <- compare_fitness_surfaces_data(
            correlated_surface = cfs_site,
            adaptive_landscape = al_site,
            trait_cols = pair
        )

        plots_site <- plot_fitness_surfaces_comparison(
            comparison_data = comparison_data_site,
            title = paste("Fitness Surfaces Comparison -", pair_name, "-", site, "(Continuous)")
        )
        plots_site$side_by_side <- plots_site$side_by_side +
            ggplot2::scale_colour_manual(values = rainbow(20)) +
            ggplot2::scale_fill_manual(values = rainbow(20))

        plots_site$overlay <- plots_site$overlay +
            ggplot2::scale_colour_manual(values = rainbow(20)) +
            ggplot2::scale_fill_manual(values = rainbow(20))

        if (!is.null(plots_site$side_by_side)) {
            ggplot2::ggsave(
                file.path(figure_dir, paste0("comparison_continuous_", pair_name, "_", site_label, "_side_by_side.png")),
                plots_site$side_by_side,
                width = 12,
                height = 5,
                dpi = 300
            )
        }

        if (!is.null(plots_site$overlay)) {
            ggplot2::ggsave(
                file.path(figure_dir, paste0("comparison_continuous_", pair_name, "_", site_label, "_overlay.png")),
                plots_site$overlay,
                width = 8,
                height = 6,
                dpi = 300
            )
        }

        cor_df_site <- comparison_data_site$combined_data
        individual_fit_site <- cor_df_site$fitness[cor_df_site$type == "Correlated Fitness (Individual)"]
        population_fit_site <- cor_df_site$fitness[cor_df_site$type == "Adaptive Landscape (Population)"]
        cor_val_site <- cor(individual_fit_site, population_fit_site, use = "complete.obs")

        write.csv(
            data.frame(
                trait1 = pair[1],
                trait2 = pair[2],
                fitness = "continuous",
                comparison_type = "site",
                site = site,
                site_type = site_type,
                correlation = cor_val_site,
                n_points = length(individual_fit_site)
            ),
            file.path(table_dir, paste0("comparison_continuous_", pair_name, "_", site_label, "_correlation.csv")),
            row.names = FALSE
        )

        comparison_results[[paste0("continuous_", pair_name, "_", site_label)]] <- list(
            correlation = cor_val_site,
            site = site,
            site_type = site_type
        )
    }
}

# ======================================================
# 12.3 Summary Table
# ======================================================

cat("\n=== Comparison Summary ===\n")

comparison_summary <- data.frame()

for (name in names(comparison_results)) {
    res <- comparison_results[[name]]

    comparison_summary <- rbind(comparison_summary, data.frame(
        Comparison = name,
        Correlation = round(res$correlation, 4),
        Site = res$site,
        stringsAsFactors = FALSE
    ))
}

write.csv(comparison_summary, file.path(table_dir, "fitness_surfaces_comparison_summary.csv"), row.names = FALSE)
