# ======================================================
# test_all_functions.R
# Full integration test for evolutionary selection package
# ======================================================

cat("\n========================================\n")
cat("RUNNING FULL FUNCTION TEST SUITE\n")
cat("========================================\n")

required_packages <- c(
  "here"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
  cat("Loaded package:", pkg, "\n")
}

set.seed(123)

cat("Project Root:", here(), "\n")

# ------------------------------------------------------
# 1 Initialize environment
# ------------------------------------------------------

if (file.exists(here("R","scripts","/0.0_initialize.R"))) {
  source(here("R","scripts","/0.0_initialize.R"))
}

# ======================================================
# 2 Output directories
# ======================================================

# ======================================================
# Load scripts
# ======================================================

cat("\nLoading script files...\n")

script_files <- list.files(
  here("R","scripts"),
  pattern = "\\.R$",
  full.names = TRUE
)

for (f in script_files) {
  source(f)
  cat("Loaded script:", basename(f), "\n")
}

# ======================================================
# Load functions
# ======================================================

cat("\nLoading function files...\n")

fn_files <- list.files(
  here("R","functions"),
  pattern = "\\.R$",
  full.names = TRUE
)

for (f in fn_files) {
  source(f)
  cat("Loaded:", basename(f), "\n")
}

# ======================================================
# Load plotting
# ======================================================

cat("\nLoading plotting functions...\n")

plot_files <- list.files(
  here("R","plotting"),
  pattern = "\\.R$",
  full.names = TRUE
)

for (f in plot_files) {
  source(f)
  cat("Loaded plot:", basename(f), "\n")
}

# ======================================================
# Define Output Directories
# ======================================================

output_dir <- here("R", "results", "test_extended_tests_results")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
model_dir <- file.path(output_dir, "models")

dirs <- c(output_dir, figure_dir, table_dir, model_dir)

for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("Results saved to:", normalizePath(output_dir), "\n")

# ======================================================
# 5 Create synthetic test data
# ======================================================

cat("\nCreating synthetic test data...\n")

test_datasets <- list(
  normal_binary = data.frame(
    fitness = rbinom(50, 1, 0.6),
    size = rnorm(50, 5, 1),
    speed = rnorm(50, 10, 2),
    color = rnorm(50, 3, 0.5)
  ),
  normal_continuous = data.frame(
    fitness = rpois(50, 3) + 1,
    size = rnorm(50, 5, 1),
    speed = rnorm(50, 10, 2),
    color = rnorm(50, 3, 0.5)
  ),
  small_sample = data.frame(
    fitness = c(0, 1, 0, 1, 0),
    size = c(4.1, 5.2, 3.9, 5.5, 4.8),
    speed = c(9.5, 10.2, 8.9, 11.1, 10.5)
  ),
  all_zeros = data.frame(
    fitness = rep(0, 20),
    size = rnorm(20, 5, 1),
    speed = rnorm(20, 10, 2)
  ),
  all_ones = data.frame(
    fitness = rep(1, 20),
    size = rnorm(20, 5, 1),
    speed = rnorm(20, 10, 2)
  ),
  high_correlation = data.frame(
    fitness = rpois(30, 2) + 1,
    size = rnorm(30, 5, 1),
    speed = rnorm(30, 5, 0.1),
    color = rnorm(30, 3, 0.5)
  )
)

# ======================================================
# Prepare constants
# ======================================================

df_binary <- test_datasets$normal_binary
df_continuous <- test_datasets$normal_continuous

FITNESS <- "fitness"
TRAITS <- c("size", "speed", "color")

trait_pairs <- list(
  c("size", "speed"),
  c("size", "color"),
  c("speed", "color")
)

cfs_plots <- list()
univariate_plots <- list()

# ======================================================
# 6 selection_coefficients
# ======================================================

cat("\nTesting selection_coefficients...\n")

test_binary <- selection_coefficients(
  data = df_binary,
  fitness_col = FITNESS,
  trait_cols = c("size", "speed"),
  fitness_type = "binary",
  standardize = TRUE
)

test_continuous <- selection_coefficients(
  data = df_continuous,
  fitness_col = FITNESS,
  trait_cols = c("size", "speed"),
  fitness_type = "continuous",
  standardize = TRUE
)

write.csv(test_binary, file.path(table_dir, "selection_binary.csv"), row.names = FALSE)
write.csv(test_continuous, file.path(table_dir, "selection_continuous.csv"), row.names = FALSE)

cat("Selection coefficients Done\n")

# ======================================================
# 7 selection_differential
# ======================================================

cat("\nTesting selection_differential...\n")

sel_diff <- selection_differential(
  data = df_continuous,
  fitness_col = FITNESS,
  trait_col = "size"
)

cat("Selection differential:", round(sel_diff, 4), "\n")

# ======================================================
# 8 Univariate spline
# ======================================================

cat("\nTesting univariate surfaces...\n")

for (trait in TRAITS) {
  spline <- univariate_spline(
    data = df_continuous,
    fitness_col = FITNESS,
    trait_col = trait,
    fitness_type = "continuous",
    k = 8
  )

  if (exists("plot_univariate_fitness")) {
    p <- plot_univariate_fitness(
      uni = spline,
      trait_col = trait,
      title = paste("Univariate fitness:", trait)
    )

    univariate_plots[[trait]] <- p

    ggsave(
      file.path(figure_dir, paste0("univariate_", trait, ".png")),
      p,
      width = 8, height = 6, dpi = 300
    )
  }
}

cat("Univariate surfaces Done\n")

# ======================================================
# 9 Correlated fitness surfaces
# ======================================================

cat("\nTesting correlated fitness surfaces...\n")

for (pair in trait_pairs) {
  cfs <- correlated_fitness_surface(
    data = df_continuous,
    fitness_col = FITNESS,
    trait_cols = pair,
    grid_n = 30
  )

  plot <- plot_correlated_fitness(
    tps = cfs,
    trait_cols = pair
  )

  name <- paste(pair, collapse = "_")

  cfs_plots[[name]] <- plot

  ggsave(
    file.path(figure_dir, paste0("cfs_", name, ".png")),
    plot,
    width = 8, height = 6, dpi = 300
  )
}

cat("Correlated fitness surfaces Done\n")

# ======================================================
# 10 Bootstrap
# ======================================================

if (exists("bootstrap_selection")) {
  cat("\nTesting bootstrap_selection...\n")

  boot <- bootstrap_selection(
    data = df_continuous,
    fitness_col = FITNESS,
    trait_cols = c("size", "speed"),
    fitness_type = "continuous",
    B = 50,
    seed = 42
  )

  write.csv(boot$ci, file.path(table_dir, "bootstrap_ci.csv"), row.names = FALSE)

  cat("Bootstrap Done\n")
}

# ======================================================
# 11 Adaptive landscape
# ======================================================

if (exists("adaptive_landscape")) {
  cat("\nTesting adaptive landscape...\n")

  fitness_model <- mgcv::gam(
    fitness ~ s(size, speed),
    data = df_continuous
  )

  landscape <- adaptive_landscape(
    data = df_continuous,
    fitness_model = fitness_model,
    trait_cols = c("size", "speed"),
    simulation_n = 400,
    grid_n = 30
  )

  saveRDS(
    landscape,
    file.path(model_dir, "adaptive_landscape.rds")
  )

  if (exists("plot_adaptive_landscape")) {
    p <- plot_adaptive_landscape(
      landscape = landscape,
      trait_cols = c("size", "speed"),
      original_data = df_continuous
    )

    ggsave(
      file.path(figure_dir, "adaptive_landscape.png"),
      p,
      width = 8, height = 6, dpi = 300
    )
  }

  cat("Adaptive landscape Done\n")
}

# ======================================================
# 12 Validation
# ======================================================

validate <- function(obj, type) {
  if (is.null(obj)) {
    return(FALSE)
  }

  if (type == "data.frame") {
    return(is.data.frame(obj))
  }
  if (type == "list") {
    return(is.list(obj))
  }
  if (type == "plot") {
    return(inherits(obj, "gg") || inherits(obj, "ggplot"))
  }

  TRUE
}

tests <- c(
  validate(test_binary, "data.frame"),
  validate(test_continuous, "data.frame"),
  if (length(univariate_plots) > 0) validate(univariate_plots[[1]], "plot") else FALSE,
  if (length(cfs_plots) > 0) validate(cfs_plots[[1]], "plot") else FALSE,
  if (exists("boot")) validate(boot, "list") else FALSE,
  if (exists("landscape")) validate(landscape, "list") else FALSE
)

cat("\nValidation results:", sum(tests), "/", length(tests), "passed\n")

# ======================================================
# 13 Final summary
# ======================================================

cat("\n========================================\n")
cat("FUNCTION TESTING COMPLETE\n")
cat("========================================\n")

cat("\nFunctions loaded:", length(fn_files), "\n")
cat("Plots created:", length(univariate_plots) + length(cfs_plots), "\n")

cat("\nFiles generated:\n")
cat("Figures:", length(list.files(figure_dir)), "\n")
cat("Tables:", length(list.files(table_dir)), "\n")
cat("Models:", length(list.files(model_dir)), "\n")

cat("\nResults saved to:\n")
cat("Figures:", figure_dir, "\n")
cat("Tables:", table_dir, "\n")
cat("Models:", model_dir, "\n")

cat("\n========================================\n")
