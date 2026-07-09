# ======================================================
# correlated_fitness_surface.R
# Important Concept Explanation:
# This function calculates the Correlated Fitness Surface.
# Definition: Individual fitness ~ Individual phenotype
# Formula: w ~ z1 + z2 + z1^2 + z2^2 + z1xz2
#
# IMPORTANT NOTES: Traits MUST be standardized BEFORE calling this function
# Use prepare_selection_data() with standardize = TRUE and optional group
# DO NOT use scale_traits = TRUE (double standardization)
# ======================================================

#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @noRd
# internal utility: the convex hull of the observed trait pairs as a closed
# polygon (first vertex repeated at the end). Predictions outside it are
# extrapolation.
.data_hull <- function(x, y) {
  pts <- cbind(x, y)
  pts <- pts[stats::complete.cases(pts), , drop = FALSE]
  h <- grDevices::chull(pts)
  pts[c(h, h[1]), , drop = FALSE]
}

#' @noRd
# internal utility: blank the grid outside the data when asked to. The full
# predictions are kept in .fit_all so the plots can draw the surface up to
# the hull edge and cover the outside with the hull polygon.
.mask_grid <- function(grid, hull, trait_cols, mask) {
  grid$.fit_all <- grid$.fit
  if (!mask) {
    grid$.inside <- TRUE
    return(grid)
  }
  grid$.inside <- mgcv::in.out(hull, cbind(grid[[trait_cols[1]]], grid[[trait_cols[2]]]))
  grid$.fit[!grid$.inside] <- NA_real_
  cat("Masked", sum(!grid$.inside), "of", nrow(grid), "grid points outside the data\n")
  grid
}

#' Calculate the Correlated Fitness Surface
#'
#' Fits a model for individual fitness based on multiple individual phenotypes (w ~ z1 + z2 + interactions).
#' Traits MUST be standardized BEFORE calling this function.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_cols A character vector of exactly length 2 specifying the trait column names.
#' @param grid_n Integer specifying the resolution of the prediction grid. Default is 60.
#' @param method A string specifying the modeling method: \code{"auto"}, \code{"gam"}, or \code{"tps"}.
#' @param scale_traits Deprecated. Logical. Set to \code{FALSE} to avoid double standardization.
#' @param group Optional string specifying a grouping variable.
#' @param k Basis dimension for the GAM smooth. The default \code{NULL} sets it
#'   from the data as \code{min(30, max(10, floor(sqrt(n1 * n2))))}, where
#'   \code{n1} and \code{n2} are the numbers of distinct values of each trait;
#'   it is always capped at one less than the number of observations. Pass an
#'   integer to override.
#' @param mask Logical; if \code{TRUE} (the default) grid points outside the
#'   convex hull of the observed trait pairs get \code{NA} fitness, so the
#'   surface is only drawn where there are data. The grid column
#'   \code{.inside} records which points were kept, \code{.fit_all} holds the
#'   unmasked predictions, and the result's \code{hull} is the polygon the
#'   plot functions use to cover the outside.
#' @param bs Basis for the GAM smooth: \code{"tp"} (thin plate, the default),
#'   \code{"cr"} or \code{"ps"}.
#' @param smoothing How the GAM's smoothing parameter is chosen: \code{"REML"}
#'   (the default), \code{"GCV.Cp"} or \code{"ML"}.
#'
#' @details By default the GAM uses a thin-plate smooth with the smoothing
#'   parameter chosen by REML; \code{bs} and \code{smoothing} are there to match
#'   another study's smoother and do not apply to \code{method = "tps"}. With
#'   \code{"cr"} or \code{"ps"}, which are one-dimensional bases, the two
#'   traits enter as a tensor product smooth.
#'   Predicting a fitted surface over the full rectangle of the grid
#'   extrapolates into corners that no individual occupies; \code{mask = TRUE}
#'   leaves those blank rather than showing a fitted value there.
#'
#' @return A list containing the fitted model, grid predictions, and metadata.
#' @export
#'
#' @examples
#' \dontrun{
#' surf <- correlated_fitness_surface(my_data, "fitness", c("trait1", "trait2"))
#' }
correlated_fitness_surface <- function(
  data,
  fitness_col,
  trait_cols,
  grid_n = 60,
  method = "auto",
  scale_traits = FALSE,
  group = NULL,
  k = NULL,
  mask = TRUE,
  bs = c("tp", "cr", "ps"),
  smoothing = c("REML", "GCV.Cp", "ML")
) {
  stopifnot(length(trait_cols) == 2L)
  bs <- match.arg(bs)
  smoothing <- match.arg(smoothing)
  need <- c(fitness_col, trait_cols)

  # Input validation
  if (!all(need %in% names(data))) {
    stop("Missing columns: ", paste(setdiff(need, names(data)), collapse = ", "))
  }

  if (!is.null(group) && !group %in% names(data)) {
    stop("Group column '", group, "' not found in data")
  }

  # DOUBLE STANDARDIZATION WARNING
  if (scale_traits) {
    warning(
      "scale_traits = TRUE is deprecated. ",
      "Traits should be standardized using prepare_selection_data() ",
      "before calling this function. Setting scale_traits = FALSE."
    )
    scale_traits <- FALSE
  }

  message("IMPORTANT: Traits should already be standardized (mean = 0, SD = 1).")
  message("           Do NOT apply scale() again within this function.")

  # Check if traits appear standardized
  for (t in trait_cols) {
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

  # Fitness and traits must be numeric (or logical); reject other types before
  # coercion, since as.numeric() on a character/factor column silently yields NA.
  for (col in need) {
    if (!is.numeric(data[[col]]) && !is.logical(data[[col]])) {
      stop("Column '", col, "' must be numeric; got ", class(data[[col]])[1])
    }
  }

  y <- as.numeric(data[[fitness_col]])
  x1 <- as.numeric(data[[trait_cols[1]]])
  x2 <- as.numeric(data[[trait_cols[2]]])

  # Remove incomplete cases
  keep <- stats::complete.cases(y, x1, x2)
  y <- y[keep]
  x1 <- x1[keep]
  x2 <- x2[keep]

  if (!is.null(group)) {
    grp <- data[[group]][keep]
  } else {
    grp <- NULL
  }

  if (length(y) < 10) stop("Too few complete cases: ", length(y), " (<10)")

  # Detect binary fitness
  uniq_y <- unique(y)
  is_binary <- length(uniq_y) == 2 && all(sort(uniq_y) == c(0, 1))

  if (method == "auto") {
    method <- if (is_binary) "gam" else "tps"
  }

  if (!method %in% c("gam", "tps")) stop("method must be 'auto' | 'gam' | 'tps'")

  # Check trait variation
  if (length(unique(x1)) < 3 || length(unique(x2)) < 3) {
    stop(
      "Too few unique trait values: ",
      trait_cols[1], " has ", length(unique(x1)),
      " unique values; ", trait_cols[2], " has ", length(unique(x2)), " unique values."
    )
  }

  n1 <- length(unique(x1))
  n2 <- length(unique(x2))
  if (is.null(k)) {
    k <- min(30, max(10, floor(sqrt(n1 * n2))))
  }

  cat("Data type:", ifelse(is_binary, "binary", "continuous"), "\n")
  cat("Selected method:", method, "\n")
  cat("Data points:", length(y), "\n")
  cat("Basis dimension k:", k, "\n")
  if (!is.null(group)) {
    cat("Grouping variable:", group, "\n")
    cat("Number of groups:", length(unique(grp)), "\n")
  }

  x1s <- x1
  x2s <- x2

  # For prediction grid (in original units, not standardized)
  # But since traits are standardized, original = standardized
  g1 <- seq(min(x1, na.rm = TRUE), max(x1, na.rm = TRUE), length.out = grid_n)
  g2 <- seq(min(x2, na.rm = TRUE), max(x2, na.rm = TRUE), length.out = grid_n)

  grid <- expand.grid(g1, g2, KEEP.OUT.ATTRS = FALSE)
  names(grid) <- trait_cols

  # Grid is already in standardized units
  grid_scaled <- grid

  hull <- .data_hull(x1, x2)
  hull_df <- if (mask) stats::setNames(as.data.frame(hull), trait_cols) else NULL

  if (method == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("mgcv package required. Please install.packages('mgcv')")
    }

    fam <- if (is_binary) stats::binomial("logit") else stats::gaussian()

    # Prepare data frame
    df_fit <- data.frame(
      .y = as.numeric(y),
      trait1 = as.numeric(x1s),
      trait2 = as.numeric(x2s)
    )
    names(df_fit)[2:3] <- trait_cols

    if (!is.null(group)) {
      df_fit[[group]] <- grp
    }

    df_fit <- df_fit[complete.cases(df_fit), ]

    cat("GAM fitting with", nrow(df_fit), "observations\n")

    # Build formulas
    k_adj <- min(k, nrow(df_fit) - 1)

    # cr and ps are one-dimensional bases, so with those the two-trait smooth is
    # a tensor product; tp handles both traits in one isotropic smooth
    joint <- if (bs == "tp") {
      paste0("s(", trait_cols[1], ", ", trait_cols[2], ", bs = 'tp', k = ", k_adj, ")")
    } else {
      paste0("te(", trait_cols[1], ", ", trait_cols[2], ", bs = '", bs, "', k = ", max(3, floor(sqrt(k_adj))), ")")
    }
    k1 <- paste0("min(floor(", k_adj, "/2), nrow(df_fit) - 1)")
    separate <- paste0("s(", trait_cols[1], ", bs = '", bs, "', k = ", k1, ") + s(", trait_cols[2], ", bs = '", bs, "', k = ", k1, ")")
    lhs <- if (!is.null(group)) paste0(".y ~ ", group, " + ") else ".y ~ "
    fml <- as.formula(paste0(lhs, joint))
    fml_alt1 <- as.formula(paste0(lhs, separate))
    fml_alt2 <- as.formula(paste0(lhs, trait_cols[1], " + ", trait_cols[2]))

    try_formulas <- list(
      main = fml,
      alt1 = fml_alt1,
      alt2 = fml_alt2
    )

    fit <- NULL
    formula_used <- NULL

    for (form_name in names(try_formulas)) {
      if (is.null(fit)) {
        tryCatch(
          {
            cat("  Trying formula:", form_name, "\n")
            fit <- mgcv::gam(try_formulas[[form_name]],
              data = df_fit,
              family = fam,
              method = smoothing
            )
            formula_used <- form_name
            cat("Success with formula:", form_name, "\n")
            break
          },
          error = function(e) {
            cat("Failed with formula", form_name, ":", e$message, "\n")
          }
        )
      }
    }

    if (is.null(fit)) {
      stop("All GAM formula attempts failed")
    }

    # Predict on grid
    newdat <- grid_scaled[, trait_cols, drop = FALSE]
    names(newdat) <- trait_cols

    if (!is.null(group)) {
      ref_group <- .reference_group(df_fit[[group]])
      newdat[[group]] <- ref_group
      cat("Predictions use group = '", ref_group, "' as reference\n")
    }

    .fit <- tryCatch(
      {
        as.numeric(stats::predict(fit, newdata = newdat, type = "response"))
      },
      error = function(e) {
        cat("Prediction failed, using mean:", e$message, "\n")
        rep(mean(y, na.rm = TRUE), nrow(newdat))
      }
    )

    grid$.fit <- .fit

    if (anyNA(grid$.fit)) {
      warning("NA predictions detected, using mean imputation")
      grid$.fit[is.na(grid$.fit)] <- mean(grid$.fit, na.rm = TRUE)
    }

    cat("Success Predictions range:", round(range(grid$.fit), 4), "\n")
    grid <- .mask_grid(grid, hull, trait_cols, mask)

    result <- list(
      model = fit,
      grid = grid,
      method = "gam",
      formula_used = formula_used,
      k = k_adj,
      basis = bs,
      smoothing = smoothing,
      mask = mask,
      hull = hull_df,
      data_type = ifelse(is_binary, "binary", "continuous"),
      trait_cols = trait_cols,
      fitness_col = fitness_col,
      group_used = group,
      surface_type = "correlated_fitness",
      note = "Correlated fitness surface (individual fitness)"
    )
    class(result) <- "correlated_fitness"
    return(result)
  }

  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("For continuous fitness with tps method, install fields: install.packages('fields')")
  }

  if (is_binary) {
    warning("Binary data detected but method='tps' chosen. Using Tps on binary outcomes (not ideal).")
  }

  Xs <- cbind(as.numeric(x1s), as.numeric(x2s))

  tps_model <- tryCatch(
    fields::Tps(Xs, as.numeric(y)),
    error = function(e) {
      warning("Tps failed, retrying with m=2: ", e$message)
      fields::Tps(Xs, as.numeric(y), m = 2)
    }
  )

  grid_scaled_mat <- cbind(
    as.numeric(grid_scaled[[trait_cols[1]]]),
    as.numeric(grid_scaled[[trait_cols[2]]])
  )

  .fit <- as.numeric(stats::predict(tps_model, grid_scaled_mat))
  grid$.fit <- .fit

  if (anyNA(grid$.fit)) {
    warning("NA predictions, using mean imputation")
    grid$.fit[is.na(grid$.fit)] <- mean(grid$.fit, na.rm = TRUE)
  }
  grid <- .mask_grid(grid, hull, trait_cols, mask)

  result <- list(
    model = tps_model,
    grid = grid,
    method = "tps",
    mask = mask,
    hull = hull_df,
    data_type = ifelse(is_binary, "binary", "continuous"),
    trait_cols = trait_cols,
    fitness_col = fitness_col,
    group_used = group,
    surface_type = "correlated_fitness",
    note = "Correlated fitness surface (individual fitness)"
  )
  class(result) <- "correlated_fitness"
  return(result)
}
