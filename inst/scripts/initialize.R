# ======================================================
# 0.0_initialize.R
# Project initialization script
# ======================================================

cat("\n========================================\n")
cat("Initializing R Environment\n")
cat("========================================\n")


set.seed(42)


required_packages <- c(
    "mgcv", # GAM models
    "fields", # Thin plate splines (Tps)
    "MASS", # multivariate normal simulation
    "Matrix", # matrix operations
    "ggplot2", # plotting
    "dplyr", # data manipulation
    "tidyr", # data tidying
    "ggrepel", # for text labels in plots
    "viridis",
    "patchwork",
    "akima",
    "car"
)

#' Helper function to load or install a package
#' @param pkg String representing the package name.
#' @noRd
load_or_install <- function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
        library(pkg, character.only = TRUE)
    }
}

invisible(lapply(required_packages, load_or_install))
cat("Packages loaded\n")


cat("\nR version:\n")
print(R.version.string)
