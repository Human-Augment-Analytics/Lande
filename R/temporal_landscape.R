# ======================================================
# temporal_landscape.R
# Fitness functions, surfaces and adaptive landscapes fitted separately for
# each level of a time column (usually year), so that changes in selection
# over time can be seen. Each period is fitted on its own; the traits are
# not restandardised per period, so the periods share one trait axis.
# ======================================================

#' @noRd
# internal utility: number of interior peaks of a fitted curve
.count_peaks <- function(y) {
  y <- y[!is.na(y)]
  if (length(y) < 3) return(NA_integer_)
  s <- sign(diff(y))
  s <- s[s != 0]
  sum(diff(s) == -2)
}

#' Fitness functions and adaptive landscapes over time
#'
#' Fits the fitness function (one trait) or fitness surface (two traits)
#' separately for each level of a time column, usually year, and optionally
#' the adaptive landscape for each, so that how selection changes over time
#' can be seen. Beausoleil et al. (2019) did this for Darwin's finches, one
#' year at a time.
#'
#' @param data A data frame with fitness, trait and time columns. Standardise
#'   the traits once, over all periods, with \code{prepare_selection_data()}
#'   before calling; the periods must share one trait axis, so nothing is
#'   restandardised per period.
#' @param fitness_col Name of the fitness column.
#' @param trait_cols One or two trait column names.
#' @param time_col Name of the column giving the period of each row.
#' @param fitness_type Passed to \code{univariate_spline()} for one trait; the
#'   surface detects the type itself.
#' @param min_n Periods with fewer complete rows than this are skipped and
#'   listed in the result's \code{skipped}. Default is 20.
#' @param landscape Logical; also compute the adaptive landscape for each
#'   period. Default is \code{TRUE}.
#' @param k,bs,smoothing Basis size, basis and smoothing criterion for the
#'   period fits. \code{NULL} uses the defaults of \code{univariate_spline()}
#'   or \code{correlated_fitness_surface()}.
#' @param bootstrap,n_boot For one trait, a bootstrap band around each
#'   period's fitness function.
#' @param grid_n Grid resolution for the surfaces and landscapes.
#' @param simulation_n Simulated individuals per grid point in the landscapes.
#' @param mask,too_far Passed to \code{correlated_fitness_surface()} for two
#'   traits.
#'
#' @return An object of class \code{"temporal_landscape"}: \code{fits} and
#'   \code{landscapes}, one per period; \code{summary}, one row per period
#'   with n, mean fitness, trait means, the smooth's degrees of freedom, the
#'   position and height of the highest fitted fitness and whether it sits at
#'   the edge of the data (\code{optimum_edge}), the number of interior
#'   peaks (one trait) and the landscape optimum; \code{grid}, the fitted
#'   values of every period stacked with a \code{time} column; for one trait
#'   \code{heat}, each period's fitness function on one common trait grid,
#'   \code{NA} outside that period's data; \code{landscape_grid}, the
#'   landscapes stacked; \code{data}, the rows used; and \code{skipped}.
#' @export
#'
#' @examples
#' finch <- read.csv(system.file("extdata", "finch_yearly.csv", package = "RforEvolution"))
#' prep <- prepare_selection_data(finch, "survived", "beak_pc1")
#' years <- temporal_landscape(prep, "survived", "beak_pc1", "year", landscape = FALSE)
#' years$summary
#' plot_temporal_landscape(years, type = "heatmap")
temporal_landscape <- function(
  data,
  fitness_col,
  trait_cols,
  time_col,
  fitness_type = c("auto", "binary", "count", "continuous"),
  min_n = 20,
  landscape = TRUE,
  k = NULL,
  bs = NULL,
  smoothing = NULL,
  bootstrap = FALSE,
  n_boot = 200,
  grid_n = 60,
  simulation_n = 300,
  mask = TRUE,
  too_far = NULL
) {
  fitness_type <- match.arg(fitness_type)
  stopifnot(length(trait_cols) %in% c(1L, 2L))
  need <- c(fitness_col, trait_cols, time_col)
  absent <- setdiff(need, names(data))
  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))
  one <- length(trait_cols) == 1L

  d <- data[stats::complete.cases(data[, need, drop = FALSE]), , drop = FALSE]
  if (!nrow(d)) stop("No complete rows")
  times <- sort(unique(d[[time_col]]))
  ranges <- stats::setNames(lapply(trait_cols, function(t) range(d[[t]])), trait_cols)
  common <- if (one) seq(ranges[[1]][1], ranges[[1]][2], length.out = 200) else NULL

  # the period fits print their own notes and warn that the traits are not
  # standardised, which is expected here: they were standardised once, over
  # all periods, so that the periods share one axis
  quiet <- function(expr) {
    out <- NULL
    utils::capture.output(out <- withCallingHandlers(
      suppressMessages(expr),
      warning = function(w) {
        if (grepl("standardized|standardised", conditionMessage(w))) invokeRestart("muffleWarning")
      }
    ))
    out
  }

  fits <- list()
  lands <- list()
  rows <- list()
  grids <- list()
  heats <- list()
  land_grids <- list()
  hulls <- list()
  skipped <- character()

  for (tm in times) {
    key <- as.character(tm)
    sub <- d[d[[time_col]] == tm, , drop = FALSE]
    if (nrow(sub) < min_n) {
      skipped <- c(skipped, key)
      next
    }
    fit <- if (one) {
      quiet(univariate_spline(sub, fitness_col, trait_cols, fitness_type = fitness_type,
                              k = k %||% 10, bs = bs %||% "cr", smoothing = smoothing %||% "GCV.Cp",
                              bootstrap = bootstrap, n_boot = n_boot))
    } else {
      quiet(correlated_fitness_surface(sub, fitness_col, trait_cols, method = "gam", grid_n = grid_n,
                                       k = k, mask = mask, too_far = too_far,
                                       bs = bs %||% "tp", smoothing = smoothing %||% "REML"))
    }
    fits[[key]] <- fit

    if (one) {
      g <- fit$grid
      g$time <- tm
      grids[[key]] <- g
      inside <- common >= min(sub[[trait_cols]]) & common <= max(sub[[trait_cols]])
      h <- stats::setNames(data.frame(common), trait_cols)
      h$fit <- NA_real_
      h$fit[inside] <- as.numeric(stats::predict(fit$model, newdata = h[inside, , drop = FALSE], type = "response"))
      h$time <- tm
      heats[[key]] <- h
      i <- which.max(g$fit)
      opt <- c(g[[trait_cols]][i], NA_real_)
      opt_fit <- g$fit[i]
      peaks <- .count_peaks(g$fit)
      at_edge <- i == 1L || i == nrow(g)
    } else {
      g <- fit$grid
      g$time <- tm
      grids[[key]] <- g
      hulls[[key]] <- fit$hull
      i <- which.max(g$.fit)
      opt <- c(g[[trait_cols[1]]][i], g[[trait_cols[2]]][i])
      opt_fit <- g$.fit[i]
      peaks <- NA_integer_
      at_edge <- .cell_states(g, trait_cols)$on_edge[i]
    }

    land_opt <- c(NA_real_, NA_real_)
    if (landscape) {
      land <- quiet(adaptive_landscape(sub, fit$model, trait_cols, simulation_n = simulation_n,
                                       grid_n = if (one) 60 else 30, custom_range = ranges))
      lands[[key]] <- land
      lg <- land$grid
      lg$time <- tm
      land_grids[[key]] <- lg
      land_opt[seq_along(trait_cols)] <- as.numeric(land$optimum[1, trait_cols])
    }

    edf <- if (inherits(fit$model, "gam")) sum(fit$model$edf) - 1 else NA_real_
    row <- data.frame(time = tm, n = nrow(sub), mean_fitness = mean(sub[[fitness_col]]),
                      stringsAsFactors = FALSE)
    for (t in trait_cols) row[[paste0("mean_", t)]] <- mean(sub[[t]])
    row$edf <- edf
    row[[paste0("optimum_", trait_cols[1])]] <- opt[1]
    if (!one) row[[paste0("optimum_", trait_cols[2])]] <- opt[2]
    row$optimum_fit <- opt_fit
    row$optimum_edge <- at_edge
    if (one) row$peaks <- peaks
    if (landscape) {
      row[[paste0("landscape_optimum_", trait_cols[1])]] <- land_opt[1]
      if (!one) row[[paste0("landscape_optimum_", trait_cols[2])]] <- land_opt[2]
    }
    rows[[key]] <- row
    message(sprintf("%s: n = %d, mean fitness %.3f, edf %.1f%s%s", key, nrow(sub), mean(sub[[fitness_col]]), edf,
                    if (one) sprintf(", %d interior peak%s", peaks, if (peaks == 1) "" else "s") else "",
                    if (at_edge) ", highest fitness at the edge of the data" else ""))
  }

  if (!length(fits)) stop("No period has at least ", min_n, " complete rows")
  if (length(skipped)) message("Skipped (fewer than ", min_n, " rows): ", paste(skipped, collapse = ", "))

  summary_df <- do.call(rbind, rows)
  rownames(summary_df) <- NULL
  result <- list(
    fits = fits,
    landscapes = if (landscape) lands else NULL,
    summary = summary_df,
    grid = do.call(rbind, grids),
    heat = if (one) do.call(rbind, heats) else NULL,
    landscape_grid = if (landscape) do.call(rbind, land_grids) else NULL,
    hulls = if (!one) hulls else NULL,
    data = d[d[[time_col]] %in% as.numeric(names(fits)) | as.character(d[[time_col]]) %in% names(fits), need, drop = FALSE],
    trait_cols = trait_cols,
    fitness_col = fitness_col,
    time_col = time_col,
    times = times[as.character(times) %in% names(fits)],
    skipped = skipped,
    too_far = too_far
  )
  class(result) <- "temporal_landscape"
  result
}

#' @export
print.temporal_landscape <- function(x, ...) {
  cat("Fitness", if (length(x$trait_cols) == 1) "function" else "surface", "by", x$time_col,
      "for", paste(x$trait_cols, collapse = " and "), "\n")
  print(x$summary, row.names = FALSE, digits = 3)
  edge <- x$summary$optimum_edge
  if (!is.null(edge) && any(edge)) {
    cat("Highest fitness at the edge of the data in", sum(edge), "of", length(edge), "periods\n")
  }
  if (length(x$skipped)) cat("Skipped:", paste(x$skipped, collapse = ", "), "\n")
  invisible(x)
}
