# ======================================================
# check_selection_assumptions.R
# The checks a selection analysis rests on: multivariate normality of the
# traits (Lande & Arnold's gradients equal the covariance only then),
# collinearity, rows per term, and the residuals or dispersion of the
# gradient models. Reported, not enforced.
# ======================================================

#' @noRd
# Mardia's (1970) multivariate skewness and kurtosis, with the small-sample
# factor for the skewness statistic. Rows are subsampled above `max_n`
# because the statistics need an n by n matrix.
.mardia <- function(X, max_n = 2000) {
  X <- as.matrix(X)
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) > max_n) X <- X[sample.int(nrow(X), max_n), , drop = FALSE]
  n <- nrow(X)
  p <- ncol(X)
  Xc <- scale(X, scale = FALSE)
  S <- crossprod(Xc) / n
  Sinv <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(Sinv)) return(NULL)
  D <- Xc %*% Sinv %*% t(Xc)
  b1p <- sum(D^3) / n^2
  b2p <- sum(diag(D)^2) / n
  k <- (p + 1) * (n + 1) * (n + 3) / (n * ((n + 1) * (p + 1) - 6))
  skew <- n * k * b1p / 6
  p_skew <- stats::pchisq(skew, df = p * (p + 1) * (p + 2) / 6, lower.tail = FALSE)
  kurt <- (b2p - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
  p_kurt <- 2 * stats::pnorm(-abs(kurt))
  list(n = n, skewness = b1p, p_skewness = p_skew, kurtosis = b2p, p_kurtosis = p_kurt)
}

#' @noRd
# Breusch-Pagan test of the squared residuals on the fitted values
.breusch_pagan <- function(fit) {
  r2 <- stats::residuals(fit)^2
  aux <- stats::lm(r2 ~ stats::fitted(fit))
  stat <- length(r2) * summary(aux)$r.squared
  c(statistic = stat, p = stats::pchisq(stat, df = 1, lower.tail = FALSE))
}

#' Check the assumptions behind a selection analysis
#'
#' Runs the checks a Lande and Arnold analysis rests on and returns them in
#' one table: multivariate normality of the traits, normality of each trait,
#' collinearity, rows per model term, and for the gradient models a residual
#' normality and heteroscedasticity test (continuous fitness), a separation
#' check (binary fitness) or the dispersion ratio (count fitness). Nothing is
#' enforced; the table is there to report alongside the gradients, as the
#' protocol of Palacio et al. (2019) asks.
#'
#' @details Lande and Arnold (1983) showed that the regression gradients equal
#'   the selection differentials adjusted for the phenotypic covariance only
#'   when the traits are multivariate normal, which is why the normality of the
#'   traits, not of fitness, is the assumption that matters. Mardia's (1970)
#'   skewness and kurtosis tests are used for that; each trait is also tested
#'   on its own with Shapiro and Wilk's test when there are 5000 rows or fewer.
#'   Collinearity is the largest variance inflation factor from the linear
#'   gradient model. The rows-per-term check counts the individuals, or for
#'   binary fitness the rarer outcome, per term of the quadratic model, with
#'   ten as the working minimum. With more than 2000 individuals Mardia's
#'   statistics use a random subsample of 2000. If the \code{performance}
#'   package is installed its heteroscedasticity and overdispersion tests are
#'   added beside the package's own, along with an R squared for the gradient
#'   model.
#'
#' @inheritParams selection_coefficients
#' @return A data frame of class \code{"selection_assumptions"} with one row
#'   per check: \code{check}, \code{statistic}, \code{p_value} and
#'   \code{note}.
#' @references Lande, R. and Arnold, S. J. (1983) The measurement of selection
#'   on correlated characters. Evolution 37, 1210-1226. Mardia, K. V. (1970)
#'   Measures of multivariate skewness and kurtosis with applications.
#'   Biometrika 57, 519-530. Palacio, F. X., Ordano, M. and Benitez-Vieyra, S.
#'   (2019) Measuring natural selection on multivariate phenotypic traits: a
#'   protocol for verifiable and reproducible analyses of natural selection.
#'   Israel Journal of Ecology and Evolution 65, 130-136.
#' @export
#' @examples
#' check_selection_assumptions(bumpus, "survival", c("total_length", "weight", "humerus"))
check_selection_assumptions <- function(data,
                                        fitness_col,
                                        trait_cols,
                                        fitness_type = c("auto", "binary", "count", "continuous"),
                                        standardize = TRUE,
                                        group = NULL) {
  fitness_type <- match.arg(fitness_type)
  need <- c(fitness_col, trait_cols, group)
  absent <- setdiff(need, names(data))
  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))

  prep <- suppressWarnings(suppressMessages(prepare_selection_data(
    data, fitness_col, trait_cols, standardize = standardize, group = group,
    add_relative = TRUE, na_action = "drop", name_relative = ".w"
  )))
  if (fitness_type == "auto") fitness_type <- detect_family(prep[[fitness_col]])$type
  if (!fitness_type %in% c("binary", "count")) fitness_type <- "continuous"
  # a group with zero mean fitness gets NA relative fitness; those rows are
  # outside the gradient models, so they are outside the counts here too
  prep <- prep[is.finite(prep$.w), , drop = FALSE]
  n <- nrow(prep)
  p <- length(trait_cols)
  rows <- list()
  add <- function(check, statistic = NA_real_, p_value = NA_real_, note = "") {
    rows[[length(rows) + 1]] <<- data.frame(check = check, statistic = statistic, p_value = p_value, note = note,
                                            stringsAsFactors = FALSE)
  }

  # traits
  X <- prep[, trait_cols, drop = FALSE]
  if (p >= 2) {
    m <- .mardia(X)
    if (!is.null(m)) {
      sub <- if (m$n < n) sprintf(" (random %d of %d rows)", m$n, n) else ""
      add("Multivariate normality: Mardia skewness", m$skewness, m$p_skewness,
          paste0(if (m$p_skewness < 0.05) "skewed" else "no evidence against normality", sub))
      add("Multivariate normality: Mardia kurtosis", m$kurtosis, m$p_kurtosis,
          if (m$p_kurtosis < 0.05) "tails differ from normal" else "no evidence against normality")
    }
  }
  for (t in trait_cols) {
    x <- X[[t]]
    if (n >= 3 && n <= 5000) {
      sw <- stats::shapiro.test(x)
      add(paste0("Normality of ", t, " (Shapiro-Wilk)"), unname(sw$statistic), sw$p.value,
          if (sw$p.value < 0.05) "not normal" else "")
    } else {
      z <- (x - mean(x)) / stats::sd(x)
      add(paste0("Skewness of ", t), mean(z^3), NA_real_, "Shapiro-Wilk needs 3 to 5000 rows")
    }
  }

  # collinearity
  if (p >= 2) {
    fml <- stats::as.formula(paste(".w ~", paste(trait_cols, collapse = " + ")))
    vif <- .compute_vif(stats::lm(fml, data = prep))
    if (!is.null(vif) && all(is.finite(vif))) {
      add("Collinearity: largest VIF", max(vif), NA_real_,
          if (max(vif) > 5) paste0("above 5 for ", trait_cols[which.max(vif)]) else "below 5")
    }
  }

  # rows per term of the quadratic model
  terms <- 2 * p + p * (p - 1) / 2
  events <- if (fitness_type == "binary") min(sum(prep[[fitness_col]] == 1), sum(prep[[fitness_col]] == 0)) else n
  per_term <- events / terms
  add(if (fitness_type == "binary") "Rarer outcome per quadratic term" else "Rows per quadratic term", per_term, NA_real_,
      if (per_term < 10) sprintf("%d %s for %d terms; treat gamma with care", events,
                                 if (fitness_type == "binary") "of the rarer outcome" else "rows", terms) else "")

  # the gradient models; the performance package adds its own tests and an R2 when installed
  has_perf <- requireNamespace("performance", quietly = TRUE)
  lin <- paste(trait_cols, collapse = " + ")
  if (fitness_type == "continuous") {
    fit <- stats::lm(stats::as.formula(paste(".w ~", lin)), data = prep)
    r <- stats::residuals(fit)
    if (length(r) >= 3 && length(r) <= 5000) {
      sw <- stats::shapiro.test(r)
      add("Linear model residuals: normality (Shapiro-Wilk)", unname(sw$statistic), sw$p.value,
          if (sw$p.value < 0.05) "residuals not normal; bootstrap the intervals" else "")
    }
    bp <- .breusch_pagan(fit)
    add("Linear model residuals: heteroscedasticity (Breusch-Pagan)", bp[["statistic"]], bp[["p"]],
        if (bp[["p"]] < 0.05) "variance changes with the fitted value; bootstrap the intervals" else "")
    if (has_perf) {
      ph <- tryCatch(as.numeric(performance::check_heteroscedasticity(fit)), error = function(e) NA_real_)
      add("Linear model residuals: heteroscedasticity (performance)", NA_real_, ph,
          if (isTRUE(ph < 0.05)) "performance agrees the variance is not constant" else "")
    }
  } else if (fitness_type == "binary") {
    fit <- suppressWarnings(stats::glm(stats::as.formula(paste(fitness_col, "~", lin)), data = prep, family = stats::binomial()))
    coefs <- stats::coef(fit)[-1]
    extreme <- max(abs(coefs), na.rm = TRUE)
    fitted_edge <- mean(stats::fitted(fit) < 1e-6 | stats::fitted(fit) > 1 - 1e-6)
    add("Logistic model: separation", extreme, NA_real_,
        if (extreme > 10 || fitted_edge > 0.05) "coefficients or fitted values at the edge; possible complete separation" else "no sign of separation")
  } else {
    fit <- suppressWarnings(stats::glm(stats::as.formula(paste(fitness_col, "~", lin)), data = prep, family = stats::poisson()))
    disp <- sum(stats::residuals(fit, type = "pearson")^2) / stats::df.residual(fit)
    p_disp <- NA_real_
    if (has_perf) {
      od <- tryCatch(performance::check_overdispersion(fit), error = function(e) NULL)
      if (!is.null(od)) p_disp <- as.numeric(od$p_value)
    }
    add("Poisson model: dispersion ratio", disp, p_disp,
        if (disp > 1.5) "overdispersed; the gradients use a negative binomial model" else "")
  }
  if (has_perf) {
    r2 <- tryCatch(performance::r2(fit), error = function(e) NULL)
    if (!is.null(r2) && length(r2)) {
      add(paste0("Model fit: ", names(r2)[1], " (performance)"), as.numeric(r2[[1]]), NA_real_, "")
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "fitness_type") <- fitness_type
  attr(out, "n") <- n
  class(out) <- c("selection_assumptions", "data.frame")
  out
}

#' Print the assumption checks
#'
#' @param x An object of class \code{"selection_assumptions"}.
#' @param ... Additional arguments (ignored).
#' @return The input object \code{x}, invisibly.
#' @export
print.selection_assumptions <- function(x, ...) {
  cat("Assumption checks (", attr(x, "fitness_type"), " fitness, n = ", attr(x, "n"), ")\n\n", sep = "")
  out <- data.frame(
    check = x$check,
    statistic = ifelse(is.na(x$statistic), "", formatC(x$statistic, digits = 3, format = "g")),
    p = ifelse(is.na(x$p_value), "", ifelse(x$p_value < 0.001, "< 0.001", formatC(x$p_value, digits = 3, format = "f"))),
    note = x$note,
    stringsAsFactors = FALSE
  )
  print(out, row.names = FALSE, right = FALSE)
  invisible(x)
}
