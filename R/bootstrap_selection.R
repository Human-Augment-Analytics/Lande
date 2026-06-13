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
#' @return A data frame with one row per coefficient and columns \code{Term},
#'   \code{Type}, \code{Estimate} (point estimate on the full data),
#'   \code{Boot_SE}, \code{CI_lower}, \code{CI_upper}, and \code{P_Value}.
#' @export
#'
#' @examples
#' \dontrun{
#' bootstrap_selection(my_data, "fitness", c("trait1", "trait2"), n_boot = 500)
#' }
bootstrap_selection <- function(data,
                                fitness_col,
                                trait_cols,
                                fitness_type = c("auto", "binary", "continuous"),
                                standardize = TRUE,
                                group = NULL,
                                use_relative_for_fit = TRUE,
                                n_boot = 1000,
                                conf = 0.95) {
  fitness_type <- match.arg(fitness_type)
  if (conf <= 0 || conf >= 1) {
    stop("`conf` must be between 0 and 1")
  }

  run <- function(d) {
    suppressWarnings(suppressMessages(
      selection_coefficients(
        d, fitness_col, trait_cols,
        fitness_type = fitness_type,
        standardize = standardize,
        group = group,
        use_relative_for_fit = use_relative_for_fit
      )
    ))
  }

  # Point estimate on the full data
  point <- run(data)
  keys <- paste(point$Term, point$Type)

  n <- nrow(data)
  boot <- matrix(NA_real_, nrow = length(keys), ncol = n_boot,
                 dimnames = list(keys, NULL))

  for (b in seq_len(n_boot)) {
    res <- tryCatch(run(data[sample.int(n, n, replace = TRUE), , drop = FALSE]),
                    error = function(e) NULL)
    if (is.null(res)) next
    boot[, b] <- res$Beta_Coefficient[match(keys, paste(res$Term, res$Type))]
  }

  n_ok <- sum(!is.na(boot[1, ]))
  if (n_ok < 2) {
    stop("Bootstrap failed: fewer than 2 usable resamples")
  }
  if (n_ok < n_boot) {
    warning(n_boot - n_ok, " of ", n_boot, " resamples failed and were dropped")
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
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out) <- NULL
  attr(out, "n_boot") <- n_ok
  attr(out, "conf") <- conf
  out
}
