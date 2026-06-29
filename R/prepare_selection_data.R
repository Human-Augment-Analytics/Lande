# ============================================================================
# prepare_selection_data
#
# Purpose:
#   1. Check fitness and traits columns
#   2. Handle missing values (including group column)
#   3. Standardize trait variables (optionally within groups)
#   4. Compute relative fitness (optionally within groups)
#
# IMPORTANT NOTE:
#   For studies spanning multiple years, sites, or populations, traits should be
#   standardized within each group (e.g., year, site) separately, not across all
#   individuals pooled together. This ensures individuals are compared relative
#   to their relevant context (same year, same environment). The same applies
#   to relative fitness - it should be calculated within the same group.
#
#   This approach is equivalent to including group fixed effects in the
#   regression models. However, it does NOT allow testing for group-by-trait
#   interactions (i.e., whether selection varies across years). For that,
#   researchers should fit separate models per group or include interaction
#   terms in the regression.
#
# Parameters:
#   data           : data frame containing fitness and trait measurements
#   fitness_col    : name of the fitness column (character)
#   trait_cols     : vector of trait column names
#   standardize    : logical; if TRUE, standardize traits to mean 0, SD 1
#   group          : optional grouping variable (e.g., "year", "site", "population")
#                    When specified, standardization and relative fitness are
#                    calculated separately within each group level.
#   add_relative   : logical; if TRUE, add a relative fitness column
#   na_action      : how to handle NAs ("warn", "drop", "none")
#   name_relative  : name for the relative fitness column
#
# Returns:
#   Modified data frame with:
#     - Standardized traits (if standardize = TRUE)
#     - Relative fitness column (if add_relative = TRUE)
#     - Potentially dropped NA rows (if na_action = "drop")
# ============================================================================

#' Prepare data for selection analysis
#'
#' Cleans, standardizes, and calculates relative fitness for trait and fitness data.
#' Standardizations and relative fitness calculations can be performed within groups.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_cols A character vector of trait column names.
#' @param standardize Logical indicating whether to standardize traits to mean 0 and SD 1. Default is \code{TRUE}.
#' @param group Optional string specifying a grouping variable.
#' @param add_relative Logical indicating whether to add a relative fitness column. Default is \code{TRUE}.
#' @param na_action A string specifying how to handle missing values: \code{"warn"}, \code{"drop"}, or \code{"none"}.
#' @param name_relative A string specifying the name for the relative fitness column. Default is \code{"relative_fitness"}.
#'
#' @return A modified data frame ready for selection analysis.
#' @export
prepare_selection_data <- function(data,
                                   fitness_col,
                                   trait_cols,
                                   standardize = TRUE,
                                   group = NULL,
                                   add_relative = TRUE,
                                   na_action = c("warn", "drop", "none"),
                                   name_relative = "relative_fitness") {
  # Input validation
  na_action <- match.arg(na_action)

  df <- data.frame(data, check.names = FALSE)

  if (!is.character(fitness_col) || length(fitness_col) != 1L) {
    stop("`fitness_col` must be a single column name (character)")
  }

  if (!all(c(fitness_col, trait_cols) %in% names(df))) {
    stop("`fitness_col` or some `trait_cols` not found in data")
  }

  # Coerce fitness to numeric if it's logical or integer
  if (is.logical(df[[fitness_col]]) || is.integer(df[[fitness_col]])) {
    df[[fitness_col]] <- as.numeric(df[[fitness_col]])
  }

  # Check if group column exists
  if (!is.null(group)) {
    if (!group %in% names(df)) {
      stop("Group column '", group, "' not found in data")
    }
    message("Standardizing and computing relative fitness within groups: '", group, "'")
  }

  # Handle missing values
  cols_check <- c(fitness_col, trait_cols)
  if (!is.null(group)) {
    cols_check <- c(cols_check, group)
  }

  na_rows <- !stats::complete.cases(df[, cols_check, drop = FALSE])

  if (any(na_rows)) {
    n_bad <- sum(na_rows)
    msg <- sprintf(
      "Detected %d row(s) with NA in: %s",
      n_bad, paste(cols_check, collapse = ", ")
    )
    if (na_action == "warn") {
      warning(msg)
    } else if (na_action == "drop") {
      df <- df[!na_rows, , drop = FALSE]
      message("Dropped ", n_bad, " row(s) with missing values")
    }
    # if na_action == "none": do nothing, but later operations may fail
  }


  # Both standardisation and relative fitness are computed from the rows that
  # are analysed (complete fitness, traits, and group), within each
  # group when `group` is set. Individuals dropped from the gradient models for
  # a missing value would otherwise shift the trait mean/SD and the mean fitness
  # and bias every gradient (Lande & Arnold standardise within the sample under
  # selection). `stat_by_group()` returns each row's group statistic over the
  # analysed rows.
  cc <- stats::complete.cases(df[, c(fitness_col, trait_cols, group), drop = FALSE])
  grp_key <- if (!is.null(group)) df[[group]] else rep(1L, nrow(df))
  stat_by_group <- function(x, FUN) {
    stats::ave(ifelse(cc, x, NA_real_), grp_key, FUN = function(v) FUN(v, na.rm = TRUE))
  }

  # Standardize traits, z = (x - mean(x)) / sd(x). A trait with no variance in
  # a group (constant, or a single observation) cannot be scaled there: it is
  # centred only, so those rows stay in the analysis with z = 0 rather than
  # becoming NaN, and the other groups are still scaled normally. A trait with
  # no variance in any group is left unscaled and named in a warning.
  if (standardize) {
    for (t in trait_cols) {
      mu <- stat_by_group(df[[t]], mean)
      sdev <- stat_by_group(df[[t]], stats::sd)
      scalable <- is.finite(sdev) & sdev > 0
      if (!any(scalable[cc])) {
        warning("Zero-variance trait left unstandardized: ", t)
        next
      }
      if (any(!scalable[cc])) {
        warning(
          "Trait '", t, "' has no variance in group(s) ",
          paste(unique(grp_key[cc & !scalable]), collapse = ", "),
          "; centred only there"
        )
      }
      df[[t]] <- (df[[t]] - mu) / ifelse(scalable, sdev, 1)
    }
  }

  # Add relative fitness, w_i = W_i / mean(W), with mean(W) over the analysed rows
  if (add_relative) {
    mean_fit <- stat_by_group(df[[fitness_col]], mean)
    if (!any(is.finite(mean_fit) & mean_fit != 0)) {
      warning("Mean fitness is zero or non-finite; cannot compute relative fitness. Skipping.")
    } else {
      df[[name_relative]] <- df[[fitness_col]] / mean_fit
    }
  }

  return(df)
}
