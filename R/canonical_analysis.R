# ======================================================
# canonical_analysis.R
# Canonical analysis of the gamma matrix (Phillips & Arnold 1989): rotate the
# quadratic and correlational gradients to axes with a single curvature each.
# ======================================================

#' @noRd
# internal utility: the gamma matrix from a selection_coefficients() table.
# The quadratic rows are gamma_ii already doubled, the correlational rows
# gamma_ij as they stand.
.gamma_matrix <- function(coefs, trait_cols) {
  p <- length(trait_cols)
  G <- matrix(NA_real_, p, p, dimnames = list(trait_cols, trait_cols))
  est <- intersect(c("Beta_Coefficient", "Estimate"), names(coefs))[1]
  for (i in seq_len(p)) {
    row <- coefs[coefs$Type == "Quadratic" & coefs$Term == paste0(trait_cols[i], "\u00b2"), ]
    if (nrow(row) == 1L) G[i, i] <- row[[est]]
    for (j in seq_len(p)) {
      if (j <= i) next
      both <- c(paste(trait_cols[i], "\u00d7", trait_cols[j]), paste(trait_cols[j], "\u00d7", trait_cols[i]))
      row <- coefs[coefs$Type == "Correlational" & coefs$Term %in% both, ]
      if (nrow(row) == 1L) G[i, j] <- G[j, i] <- row[[est]]
    }
  }
  G
}

#' @noRd
# internal utility: eigen decomposition with a fixed sign, the largest loading
# of each axis positive, so the axes come out the same way every time
.canonical_axes <- function(G) {
  e <- eigen(G, symmetric = TRUE)
  M <- e$vectors
  for (k in seq_len(ncol(M))) {
    if (M[which.max(abs(M[, k])), k] < 0) M[, k] <- -M[, k]
  }
  dimnames(M) <- list(rownames(G), paste0("m", seq_len(ncol(M))))
  list(lambda = e$values, M = M)
}

#' @noRd
# internal utility: the double-regression F (t squared) of each ordered
# eigenvalue, from least squares of fitness `w` on the traits `Z`. gamma comes
# from the full quadratic model (squared terms doubled), its axes from eigen(),
# and the F of axis k from refitting on the canonical scores and their squares.
.canonical_f <- function(Z, w, qx, pairs) {
  p <- ncol(Z)
  b <- qr.coef(qx, w)
  G <- diag(2 * b[(p + 2):(2 * p + 1)], p)
  for (k in seq_len(ncol(pairs))) G[pairs[1, k], pairs[2, k]] <- G[pairs[2, k], pairs[1, k]] <- b[2 * p + 1 + k]
  S <- Z %*% eigen(G, symmetric = TRUE)$vectors
  D <- cbind(1, S, S^2)
  fit <- stats::lm.fit(D, w)
  sigma2 <- sum(fit$residuals^2) / (nrow(D) - ncol(D))
  sq <- (p + 2):(2 * p + 1)
  fit$coefficients[sq]^2 / (sigma2 * diag(solve(crossprod(D)))[sq])
}

#' @noRd
# internal utility: permutation p-values for the ordered eigenvalues of gamma
# (Reynolds et al. 2010). Fitness is shuffled among individuals, within groups
# when there are groups, gamma and its axes are estimated again each time, and
# axis k's p-value is the share of shuffles whose F for the k-th eigenvalue is
# at least the observed one, counting the observed data as one of them.
.canonical_permutation <- function(Z, w, g, n_perm) {
  p <- ncol(Z)
  pairs <- utils::combn(p, 2)
  qx <- qr(cbind(1, Z, Z^2, apply(pairs, 2, function(ij) Z[, ij[1]] * Z[, ij[2]])))
  observed <- .canonical_f(Z, w, qx, pairs)
  blocks <- if (is.null(g)) list(seq_along(w)) else split(seq_along(w), g)
  hits <- numeric(p)
  for (b in seq_len(n_perm)) {
    idx <- seq_along(w)
    for (i in blocks) idx[i] <- i[sample.int(length(i))]
    hits <- hits + (.canonical_f(Z, w[idx], qx, pairs) >= observed - 1e-12)
  }
  list(p_value = (hits + 1) / (n_perm + 1), f = observed)
}

#' Canonical analysis of the gamma matrix
#'
#' Rotates the matrix of quadratic and correlational selection gradients to
#' its canonical axes (Phillips and Arnold 1989; Blows and Brooks 2003). Each
#' axis is a combination of the traits with a single curvature, its
#' eigenvalue: negative for stabilising and positive for disruptive selection
#' along that axis. Correlational selection that is spread over many
#' \eqn{\gamma_{ij}} shows up as curvature on a few axes.
#'
#' @inheritParams selection_coefficients
#' @param bootstrap Logical; add percentile intervals for the eigenvalues by
#'   resampling individuals, within groups when \code{group} is given. Default
#'   is \code{FALSE}.
#' @param n_boot Number of bootstrap resamples. Default is 200.
#' @param conf Confidence level of the bootstrap interval. Default is 0.95.
#' @param test Where the p-values come from: \code{"permutation"} (the
#'   default), Reynolds et al.'s (2010) permutation test, or
#'   \code{"double_regression"}, the double-regression tests, which are far too
#'   liberal.
#' @param n_perm Number of permutations. Default is 999.
#'
#' @details The gamma matrix \eqn{\gamma} comes from
#'   \code{selection_coefficients()}, with the quadratic gradients doubled and
#'   the correlational gradients as they stand. Its eigenvectors are the
#'   columns of \code{M} and its eigenvalues the curvatures \eqn{\lambda};
#'   directional selection along the axes is \eqn{\theta = M^T \beta}.
#'   Standard errors come from the double regression of Bisgaard and Ankenman
#'   (1996): fitness is refitted on the canonical scores and their squares,
#'   without cross-products, and twice the coefficient of a squared score is
#'   that axis's \eqn{\lambda}.
#'
#'   The double-regression tests treat the axes as known when they were
#'   estimated from the same data. With 136 individuals and fitness unrelated
#'   to five traits, the largest eigenvalue came out significant in about three
#'   samples in four, so by default
#'   the p-values come from the permutation test of Reynolds et al. (2010):
#'   fitness is shuffled among individuals, within groups when \code{group} is
#'   given, the canonical analysis and the double regression are repeated, and
#'   each axis's p-value is the share of shuffles, counting the data as one,
#'   whose F for the eigenvalue of the same rank is at least the observed one.
#'   It tests whether an axis curves more than it would with fitness unrelated
#'   to the traits; once one axis is strongly curved, the tests of the others
#'   are only rough. Every fitness type uses least squares on relative fitness,
#'   and the result depends on the random seed; call \code{set.seed()} first. With
#'   \code{test = "double_regression"} the p-values come from the double
#'   regression itself, for survival and counts from the logistic or count
#'   model on the same terms. The standard errors treat the axes as known
#'   either way.
#'
#'   Sampling error in \eqn{\gamma} also spreads its eigenvalues, so the largest
#'   curvatures are overestimated, especially with many traits and few
#'   individuals (Reynolds et al. 2010), and the estimated axes lean towards
#'   directions of phenotype with little variance (Morrissey 2014). Read the
#'   eigenvalues together with the fitness surface along the same axes, which
#'   \code{plot_canonical_axes()} draws, and with the bootstrap intervals. With
#'   the bootstrap, each resample's axes are matched to the original ones by
#'   their largest absolute cosine, since order and sign change between
#'   resamples.
#'
#' @return An object of class \code{"canonical_analysis"}: \code{gamma}, the
#'   gamma matrix; \code{beta}; \code{M}, the loadings with one column per
#'   axis; \code{axes}, a data frame with \code{lambda}, its standard error and
#'   p-value, \code{theta}, and with the bootstrap \code{CI_lower},
#'   \code{CI_upper} and \code{N_Boot}; \code{scores}, the prepared data with
#'   the canonical scores \code{m1}, \code{m2}, ... added; the fitness type,
#'   \code{n}, \code{test} and \code{n_perm}.
#' @references Bisgaard, S. and Ankenman, B. (1996) Standard errors for the
#'   eigenvalues in second-order response surface models. Technometrics 38,
#'   238-246. Blows, M. W. and Brooks, R. (2003) Measuring nonlinear
#'   selection. The American Naturalist 162, 815-820. Morrissey, M. B. (2014)
#'   In search of the best methods for multivariate selection analysis.
#'   Methods in Ecology and Evolution 5, 1095-1109. Phillips, P. C. and
#'   Arnold, S. J. (1989) Visualizing multivariate selection. Evolution 43,
#'   1209-1222. Reynolds, R. J., Childers, D. K. and Pajewski, N. M. (2010) The
#'   distribution and hypothesis testing of eigenvalues from the canonical
#'   analysis of the gamma matrix of quadratic and correlational selection
#'   gradients. Evolution 64, 1076-1085.
#' @export
#'
#' @examples
#' ca <- canonical_analysis(bumpus, "survival", c("total_length", "weight", "humerus"))
#' ca
canonical_analysis <- function(data,
                               fitness_col,
                               trait_cols,
                               fitness_type = c("auto", "binary", "count", "continuous"),
                               standardize = TRUE,
                               group = NULL,
                               bootstrap = FALSE,
                               n_boot = 200,
                               conf = 0.95,
                               test = c("permutation", "double_regression"),
                               n_perm = 999) {
  fitness_type <- match.arg(fitness_type)
  test <- match.arg(test)
  if (length(trait_cols) < 2L) stop("Canonical analysis needs at least two traits")
  if (test == "permutation" && (!is.numeric(n_perm) || length(n_perm) != 1 || !is.finite(n_perm) ||
                                n_perm < 1 || n_perm != round(n_perm))) {
    stop("`n_perm` must be a whole number of 1 or more")
  }

  gamma_of <- function(d) {
    coefs <- suppressMessages(selection_coefficients(d, fitness_col, trait_cols, fitness_type = fitness_type,
                                                     standardize = standardize, group = group))
    list(G = .gamma_matrix(coefs, trait_cols), coefs = coefs)
  }
  full <- gamma_of(data)
  G <- full$G
  if (anyNA(G)) stop("The quadratic model did not estimate every gradient, so the gamma matrix is incomplete")
  coefs <- full$coefs
  linear <- coefs[coefs$Type == "Linear", ]
  beta <- stats::setNames(linear$Beta_Coefficient[match(trait_cols, linear$Term)], trait_cols)
  used_type <- attr(coefs, "fitness_type_used")

  ax <- .canonical_axes(G)
  M <- ax$M
  axis_names <- colnames(M)

  # canonical scores on the prepared traits
  rel_col <- ".w_rel"
  prep <- suppressWarnings(suppressMessages(prepare_selection_data(
    data, fitness_col, trait_cols, standardize = standardize, group = group,
    add_relative = TRUE, na_action = "warn", name_relative = rel_col
  )))
  keep <- stats::complete.cases(prep[, c(fitness_col, trait_cols), drop = FALSE]) & is.finite(prep[[rel_col]])
  prep <- prep[keep, , drop = FALSE]
  scores <- as.matrix(prep[, trait_cols, drop = FALSE]) %*% M
  colnames(scores) <- axis_names
  prep[axis_names] <- as.data.frame(scores)

  # double regression: the scores and their squares, no cross-products
  rhs <- paste(c(axis_names, paste0("I(", axis_names, "^2)")), collapse = " + ")
  ols <- stats::lm(stats::as.formula(paste(rel_col, "~", rhs)), data = prep)
  est <- summary(ols)$coefficients
  sq <- paste0("I(", axis_names, "^2)")
  p_from <- est
  if (used_type %in% c("binary", "count")) {
    prep <- .add_glm_group(prep, group)
    glm_rhs <- if (".group" %in% names(prep)) paste(".group +", rhs) else rhs
    glm_fit <- tryCatch(
      suppressWarnings(.fit_pvalue_glm(stats::as.formula(paste(fitness_col, "~", glm_rhs)), prep, used_type)),
      error = function(e) NULL
    )
    if (!is.null(glm_fit)) p_from <- summary(glm_fit)$coefficients
  }
  p_col <- .p_col_from_summary(p_from)
  axes <- data.frame(
    axis = axis_names,
    lambda = ax$lambda,
    se = 2 * est[sq, "Std. Error"],
    p_value = as.numeric(p_from[sq, p_col]),
    theta = as.numeric(crossprod(M, beta)),
    stringsAsFactors = FALSE
  )
  if (test == "permutation") {
    g <- if (!is.null(group)) droplevels(addNA(factor(prep[[group]]), ifany = TRUE)) else NULL
    perm <- .canonical_permutation(as.matrix(prep[, trait_cols, drop = FALSE]), prep[[rel_col]], g, round(n_perm))
    axes$p_value <- perm$p_value
  }

  if (bootstrap) {
    grp <- if (!is.null(group)) addNA(factor(data[[group]]), ifany = TRUE) else NULL
    draws <- matrix(NA_real_, n_boot, length(axis_names))
    for (b in seq_len(n_boot)) {
      rows <- if (is.null(grp)) {
        sample.int(nrow(data), replace = TRUE)
      } else {
        unlist(lapply(split(seq_len(nrow(data)), grp), function(i) i[sample.int(length(i), replace = TRUE)]), use.names = FALSE)
      }
      Gb <- tryCatch(suppressWarnings(gamma_of(data[rows, , drop = FALSE])$G), error = function(e) NULL)
      if (is.null(Gb) || anyNA(Gb)) next
      eb <- .canonical_axes(Gb)
      # match this resample's axes to the original by largest |cosine|
      cosines <- abs(crossprod(M, eb$M))
      taken <- integer(0)
      for (k in seq_along(axis_names)) {
        free <- setdiff(seq_along(axis_names), taken)
        j <- free[which.max(cosines[k, free])]
        draws[b, k] <- eb$lambda[j]
        taken <- c(taken, j)
      }
    }
    alpha <- (1 - conf) / 2
    axes$CI_lower <- apply(draws, 2, stats::quantile, probs = alpha, na.rm = TRUE)
    axes$CI_upper <- apply(draws, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
    axes$N_Boot <- colSums(!is.na(draws))
  }

  result <- list(
    gamma = G,
    beta = beta,
    M = M,
    axes = axes,
    scores = prep,
    fitness_col = fitness_col,
    trait_cols = trait_cols,
    fitness_type = used_type,
    group = group,
    n = nrow(prep),
    test = test,
    n_perm = if (test == "permutation") round(n_perm) else NA_integer_
  )
  class(result) <- "canonical_analysis"
  result
}

#' @export
print.canonical_analysis <- function(x, ...) {
  cat("Canonical analysis of gamma for", paste(x$trait_cols, collapse = ", "), "on", x$n, "individuals\n\n")
  cat("Loadings (columns are the canonical axes):\n")
  print(round(x$M, 3))
  cat("\nCurvature along each axis (negative is consistent with stabilising selection, positive with disruptive):\n")
  print(x$axes, row.names = FALSE, digits = 3)
  if (identical(x$test, "double_regression")) {
    cat("\np-values from the double regression treat the axes as known although they come from",
        "\nthese data, and are much too small; use test = \"permutation\".\n")
  } else {
    cat("\np-values from", x$n_perm, "permutations of fitness (Reynolds et al. 2010), which allow for the",
        "\naxes coming from these data. The standard errors treat the axes as known.\n")
  }
  invisible(x)
}

#' Fitness along the canonical axes
#'
#' Draws the fitness surface in the space of two canonical axes, or the
#' fitness function along one, from the scores of \code{canonical_analysis()}.
#' A negative \eqn{\lambda} is stabilising selection along an axis only when
#' fitness along that axis has a peak inside the data, otherwise it is
#' only curvature.
#'
#' @param ca Output of \code{canonical_analysis()}.
#' @param which One or two axis numbers. Default is \code{1:2}.
#' @param grid_n Grid resolution of the surface. Default is 40.
#' @param ... Passed to \code{plot_correlated_fitness()} for two axes or to
#'   \code{plot_univariate_fitness()} for one.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' ca <- canonical_analysis(bumpus, "survival", c("total_length", "weight", "humerus"))
#' plot_canonical_axes(ca, which = c(1, 3), grid_n = 25)
plot_canonical_axes <- function(ca, which = 1:2, grid_n = 40, ...) {
  stopifnot(inherits(ca, "canonical_analysis"), length(which) %in% 1:2)
  axes <- colnames(ca$M)[which]
  if (anyNA(axes)) stop("`which` must pick from the ", ncol(ca$M), " canonical axes")
  quiet <- function(expr) {
    out <- NULL
    utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))
    out
  }
  if (length(axes) == 1L) {
    fit <- quiet(univariate_spline(ca$scores, ca$fitness_col, axes))
    return(plot_univariate_fitness(fit, axes, ...))
  }
  surface <- quiet(correlated_fitness_surface(ca$scores, ca$fitness_col, axes, grid_n = grid_n, method = "gam"))
  plot_correlated_fitness(surface, axes, ...)
}
