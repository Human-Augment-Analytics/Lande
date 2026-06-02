# =============================================================================
# PART 1: SETUP AND INITIALIZATION
# =============================================================================

cat("1. Working directory and package initialization\n")

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "mgcv", "fields",
  "purrr", "patchwork", "viridis", "scales", "knitr", "here"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
  cat("Loaded package:", pkg, "\n")
}

cat("Project Root:", here(), "\n")

cat("\n2. Loading selection analysis functions\n")

function_files <- c(
  "extract_results.R",
  "selection_coefficients.R", 
  "detect_family.R", 
  "selection_differential.R",
  "univariate_spline.R", 
  "univariate_surface.R", 
  "correlational_tps.R",
  "correlation_surface.R", 
  "bootstrap_selection.R"
)

scripts_files <- c(
  "1_prepare_selection_data.R",
  "2_linear_selection_analysis.R",
  "3_nonlinear_selection_analysis.R",
  "4_disruptive_selection_analysis.R"
)

for (f in function_files) {
  file_path <- here("R","functions", f)
  if (file.exists(file_path)) {
    source(file_path)
    cat("Sourced:", f, "\n")
  } else {
    cat("File not found:", file_path, "\n")
  }
}

for (s in scripts_files) {
  file_path <- here("R","scripts", s)
  if (file.exists(file_path)) {
    source(file_path)
    cat("Sourced:", s, "\n")
  } else {
    cat("File not found:", file_path, "\n")
  }
}

# =============================================================================
# PART 2: DATA LOADING AND EXPLORATION
# =============================================================================

cat("\n3. Data loading and exploration\n")

data_files <- list(
  data1 = here("R","data", "Aster_analyses_2011_Cohort.txt"),
  data2 = here("R","data", "Aster_analyses_2012_Cohort_full.txt"),
  data3 = here("R","data", "Aster_analyses_2011_Cohort_full.txt"),
  data4 = here("R","data", "Aster_analyses_2012_Cohort.txt")
)

data1 <- read.delim(data_files$data1, sep = "\t")
data2 <- read.delim(data_files$data2, sep = "\t")
data3 <- read.delim(data_files$data3, sep = "\t")
data4 <- read.delim(data_files$data4, sep = "\t")

cat("2011 Cohort:", dim(data1)[1], "rows,", dim(data1)[2], "columns\n")
cat("2012 Full:", dim(data2)[1], "rows,", dim(data2)[2], "columns\n")
cat("2011 Full:", dim(data3)[1], "rows,", dim(data3)[2], "columns\n")
cat("2012 Cohort:", dim(data4)[1], "rows,", dim(data4)[2], "columns\n")

# =============================================================================
# PART 3: DATA PREPARATION
# =============================================================================

cat("\n4. Data preparation\n")

c13_col <- grep("c13|deltaC13", names(data3), ignore.case = TRUE, value = TRUE)[1]

traits_2011 <- data3 %>%
  select(geno, elevation, season, SLA, !!c13_col, FDsnow, height) %>%
  rename(deltaC13 = !!c13_col) %>%
  distinct(geno, season, .keep_all = TRUE)

fitness_2011 <- data1 %>%
  select(geno, season, varb, resp) %>%
  filter(!is.na(resp))

fitness_wide <- fitness_2011 %>%
  pivot_wider(
    names_from = c(season, varb),
    names_sep = "_",
    names_prefix = "Year",
    values_from = resp
  )

analysis_data <- traits_2011 %>%
  left_join(fitness_wide, by = "geno")

analysis_data_clean <- analysis_data %>%
  filter(
    !is.na(SLA),
    !is.na(deltaC13),
    !is.na(FDsnow),
    !is.na(height)
  )

names(analysis_data_clean) <- make.names(names(analysis_data_clean))

trait_cols <- c("SLA", "deltaC13", "FDsnow", "height")

cat("Families retained:", nrow(analysis_data_clean), "\n")
cat("Traits:", paste(trait_cols, collapse = ", "), "\n")

fecund_cols <- grep("_fecund", names(analysis_data_clean), value = TRUE)
surv_cols <- grep("_surv|survived", names(analysis_data_clean), value = TRUE)
flower_cols <- grep("_flr", names(analysis_data_clean), value = TRUE)

main_fecund_col <- grep("^Year2013_.*fecund",
  names(analysis_data_clean),
  value = TRUE
)[1]

if (is.na(main_fecund_col)) main_fecund_col <- fecund_cols[1]

cat("Primary fitness variable:", main_fecund_col, "\n")

# =============================================================================
# PART 4: SINGLE-YEAR SELECTION ANALYSIS
# =============================================================================

cat("\n5. Single-year selection analysis\n")

main_data <- analysis_data_clean
main_fitness <- main_fecund_col

main_data_prepared <- prepare_selection_data(
  data = main_data,
  fitness_col = main_fitness,
  trait_cols = trait_cols,
  standardize = TRUE,
  add_relative = TRUE,
  na_action = "warn"
)

cat("Sample size:", nrow(main_data_prepared), "\n")

multi_trait_result <- selection_coefficients(
  data = main_data,
  fitness_col = main_fitness,
  trait_cols = trait_cols,
  fitness_type = "continuous",
  standardize = TRUE
)

cat("Multivariate selection gradients:\n")
print(multi_trait_result)

# =============================================================================
# PART 5: TEMPORAL ANALYSIS
# =============================================================================

cat("\n6. Temporal analysis\n")

yearly_results <- list()
fecund_cols <- grep("^Year201[2-4]_.*fecund",
  names(analysis_data_clean),
  value = TRUE
)

years_to_analyze <- as.integer(sub("Year", "", sub("_.*", "", fecund_cols)))

for (i in seq_along(fecund_cols)) {
  year <- years_to_analyze[i]
  fitness_col <- fecund_cols[i]

  year_data <- analysis_data_clean %>%
    filter(!is.na(.data[[fitness_col]]))

  if (nrow(year_data) < 10) next

  res <- analyze_disruptive_selection(
    data = year_data,
    fitness_col = fitness_col,
    trait_col = "FDsnow",
    fitness_type = "continuous",
    standardize = TRUE
  )

  lin <- res %>% filter(Type == "Linear")
  quad <- res %>% filter(Type == "Quadratic")

  if (nrow(lin) == 0 | nrow(quad) == 0) next

  yearly_results[[as.character(year)]] <- data.frame(
    Year = year,
    Sample_Size = nrow(year_data),
    Linear = lin$Beta_Coefficient,
    Quadratic = quad$Beta_Coefficient
  )
}

if (length(yearly_results) > 0) {
  yearly_summary <- bind_rows(yearly_results)
  print(yearly_summary)
}

# =============================================================================
# PART 6: SAVE RESULTS
# =============================================================================

cat("\n7. Saving results\n")

output_dir <- here("R","results","plant_selection_results")
if (!dir.exists(output_dir)) dir.create(output_dir)

# 1. Save the main multivariate results to CSV
write.csv(multi_trait_result, 
          file = file.path(output_dir, "multivariate_selection_results.csv"), 
          row.names = FALSE)

# 2. Save the temporal analysis to CSV if it exists
if (exists("yearly_summary")) {
  write.csv(yearly_summary, 
            file = file.path(output_dir, "temporal_selection_results.csv"), 
            row.names = FALSE)
}

# 3. Create a simple text summary report
summary_conn <- file(file.path(output_dir, "summary_report.txt"))
writeLines(c(
  paste("Analysis Date:", Sys.time()),
  paste("Families analyzed:", nrow(main_data_prepared)),
  paste("Traits analyzed:", paste(trait_cols, collapse = ", ")),
  paste("Significant gradients:", sum(multi_trait_result$P_Value < 0.05, na.rm = TRUE))
), summary_conn)
close(summary_conn)

# 4. Save the R-specific objects for future loading
analysis_results <- list(
  date = Sys.time(),
  sample_size = nrow(main_data_prepared),
  traits = trait_cols,
  main_fitness = main_fitness,
  multi_trait = multi_trait_result,
  temporal = if (exists("yearly_summary")) yearly_summary else NULL
)

saveRDS(
  analysis_results, 
  file.path(output_dir, "analysis_results.rds")
)
save.image(file.path(output_dir, "workspace.RData"))

cat("Results saved to directory:", output_dir, "\n")
cat("  - multivariate_selection_results.csv\n")
if (exists("yearly_summary")) cat("  - temporal_selection_results.csv\n")
cat("  - summary_report.txt\n")
cat("  - analysis_results.rds (for R)\n")

# =============================================================================
# PART 7: FINAL SUMMARY
# =============================================================================

cat("\n8. Summary\n")
cat("Families analyzed:", nrow(main_data_prepared), "\n")
cat("Traits analyzed:", length(trait_cols), "\n")

if (!is.null(multi_trait_result)) {
  sig_count <- sum(multi_trait_result$P_Value < 0.05, na.rm = TRUE)
  cat("Significant gradients detected:", sig_count, "\n")
}
