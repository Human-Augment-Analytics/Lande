# ============================================================================
# bootstrap_selection
#
# Bootstrap standard errors and confidence intervals for selection gradients.
# Individuals are resampled with replacement and the whole estimation procedure
# (standardisation, relativisation, OLS) is repeated on each resample, so the
# intervals reflect the full pipeline rather than a single fitted model.
# ============================================================================

#' Bootstrap selection gradients
#'
#' Resamples individuals with replacement and re-estimates the selection
#' gradients on each resample to obtain bootstrap standard errors and percentile
#' confidence intervals for every coefficient (beta, gamma, gamma_ij).
#'
#' @inheritParams selection_coefficients
#' @param n_boot Integer number of bootstrap resamples. Default is 1000.
#' @param conf Confidence level for the percentile interval. Default is 0.95.
#'
#' @details Each resample goes through the whole procedure again: with
#'   \code{standardize = TRUE} the traits are restandardised, and fitness is made
#'   relative to the resample's own mean, before both models are refitted, so
#'   the intervals carry the uncertainty in the standardisation as well as in
#'   the fit. Resamples whose fit fails are dropped, and \code{N_Boot} counts
#'   the ones that were used.
#'
#'   When \code{group} is given, individuals are resampled within each
#'   group so every resample keeps the original group sizes. A resample in
#'   which a trait has no variance within a group (so it could only be centred)
#'   is discarded rather than fitted on a degenerate value. The point estimate
#'   is the ordinary \code{selection_coefficients()} fit on the full data and
#'   its warnings are reported as usual.
#'
#' @return A data frame with one row per coefficient and columns \code{Term},
#'   \code{Type}, \code{Estimate} (point estimate on the full data),
#'   \code{Boot_SE}, \code{CI_lower}, \code{CI_upper}, \code{P_Value}, and
#'   \code{N_Boot} (usable resamples for that coefficient). \code{P_Value} is
#'   the parametric p-value from the point-estimate fit (OLS t-test, or the
#'   logistic-GLM Wald test for binary fitness), not a bootstrap p-value.
#' @export
#'
#' @examples
#' set.seed(1)
#' bootstrap_selection(bumpus, "survival", c("total_length", "weight"), n_boot = 50)
bootstrap_selection <- function(data,
                                fitness_col,
                                trait_cols,
                                fitness_type = c("auto", "binary", "count", "continuous"),
                                standardize = TRUE,
                                group = NULL,
                                use_relative_for_fit = TRUE,
                                n_boot = 1000,
                                conf = 0.95) {
  fitness_type <- match.arg(fitness_type)
  if (!is.numeric(conf) || length(conf) != 1L || is.na(conf) || conf <= 0 || conf >= 1) {
    stop("`conf` must be a single number between 0 and 1")
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) ||
      n_boot < 2 || n_boot != round(n_boot)) {
    stop("`n_boot` must be a single whole number of at least 2")
  }
  n_boot <- as.integer(n_boot)
  if (!is.null(group) && !group %in% names(data)) {
    stop("Group column '", group, "' not found in data")
  }

  fit_once <- function(d) {
    selection_coefficients(
      d, fitness_col, trait_cols,
      fitness_type = fitness_type,
      standardize = standardize,
      group = group,
      use_relative_for_fit = use_relative_for_fit
    )
  }

  # Point estimate on the full data; its warnings are the user's to see.
  point <- suppressMessages(fit_once(data))
  keys <- paste(point$Term, point$Type)

  # A resample is unusable if it errors, or if a trait lost its variance in
  # some group (prepare_selection_data then centres only, which would put a
  # fabricated z = 0 into the fit).
  fit_resample <- function(d) {
    degenerate <- FALSE
    res <- tryCatch(
      withCallingHandlers(
        suppressMessages(fit_once(d)),
        warning = function(w) {
          if (grepl("no variance|Zero-variance", conditionMessage(w))) degenerate <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )
    if (degenerate) NULL else res
  }

  # Resample within groups when there are groups, so each resample keeps the
  # original group sizes.
  n <- nrow(data)
  strata <- if (is.null(group)) rep(1L, n) else data[[group]]
  strata_rows <- split(seq_len(n), strata, drop = TRUE)
  resample_rows <- function() {
    unlist(lapply(strata_rows, function(r) r[sample.int(length(r), length(r), replace = TRUE)]),
           use.names = FALSE)
  }

  boot <- matrix(NA_real_, nrow = length(keys), ncol = n_boot,
                 dimnames = list(keys, NULL))

  for (b in seq_len(n_boot)) {
    res <- fit_resample(data[resample_rows(), , drop = FALSE])
    if (is.null(res)) next
    boot[, b] <- res$Beta_Coefficient[match(keys, paste(res$Term, res$Type))]
  }

  # Usable resamples per coefficient. A term that drops out of some fits (e.g. a
  # collinear interaction) has fewer usable resamples than the rest, so count
  # each row on its own instead of reading only the first.
  n_used <- rowSums(!is.na(boot))
  n_min <- min(n_used)
  if (n_min < 2) {
    stop("Bootstrap failed: fewer than 2 usable resamples for at least one term")
  }
  if (n_min < n_boot) {
    warning(n_boot - n_min, " of ", n_boot,
            " resamples were unusable for at least one term and were dropped")
  }

  a <- (1 - conf) / 2
  out <- data.frame(
    Term = point$Term,
    Type = point$Type,
    Estimate = point$Beta_Coefficient,
    Boot_SE = apply(boot, 1, stats::sd, na.rm = TRUE),
    CI_lower = apply(boot, 1, stats::quantile, probs = a, na.rm = TRUE),
    CI_upper = apply(boot, 1, stats::quantile, probs = 1 - a, na.rm = TRUE),
    P_Value = point$P_Value,
    N_Boot = as.integer(n_used),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out) <- NULL
  attr(out, "n_boot") <- n_min
  attr(out, "conf") <- conf
  out
}
