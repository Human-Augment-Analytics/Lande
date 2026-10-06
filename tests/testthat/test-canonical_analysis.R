test_that("the axes are orthonormal and rebuild gamma", {
  G <- matrix(c(-0.016, -0.028, 0.103,
                -0.028, 0.004, 0.030,
                0.103, 0.030, -0.052), 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  ax <- .canonical_axes(G)
  expect_equal(unname(crossprod(ax$M)), diag(3), tolerance = 1e-10)
  expect_equal(unname(ax$M %*% diag(ax$lambda) %*% t(ax$M)), unname(G), tolerance = 1e-10)
  expect_equal(sum(ax$lambda), sum(diag(G)))
  # the largest loading of every axis is positive
  expect_true(all(apply(ax$M, 2, function(m) m[which.max(abs(m))]) > 0))
})

test_that("canonical analysis finds a known curvature that the traits share", {
  set.seed(14)
  n <- 600
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  # stabilising selection along z1 + z2 and nothing along z1 - z2
  u <- (z1 + z2) / sqrt(2)
  df <- data.frame(z1 = z1, z2 = z2, w = 2 - 0.4 * u^2 + rnorm(n, 0, 0.2))
  ca <- quiet(canonical_analysis(df, "w", c("z1", "z2")))
  expect_s3_class(ca, "canonical_analysis")
  expect_equal(names(ca$axes), c("axis", "lambda", "se", "p_value", "theta"))
  # invariants: the trace, orthonormal loadings, gamma rebuilt from its axes
  expect_equal(sum(ca$axes$lambda), sum(diag(ca$gamma)))
  expect_equal(unname(crossprod(ca$M)), diag(2), tolerance = 1e-10)
  expect_equal(unname(ca$M %*% diag(ca$axes$lambda) %*% t(ca$M)), unname(ca$gamma), tolerance = 1e-10)
  # the strongly curved axis is the last one and loads on both traits equally
  last <- ca$axes[2, ]
  expect_lt(last$lambda, -0.2)
  # no shuffle of fitness curves as much, so the smallest p a permutation test can give
  expect_equal(ca$test, "permutation")
  expect_equal(last$p_value, 1 / 1000)
  expect_equal(abs(ca$M[, 2]), c(z1 = 1, z2 = 1) / sqrt(2), tolerance = 0.1)
  expect_lt(abs(ca$axes$lambda[1]), 0.1)
  # the double regression on the scores gives the same curvature as the eigenvalue
  fit <- lm(.w_rel ~ m1 + m2 + I(m1^2) + I(m2^2), data = ca$scores)
  expect_equal(unname(2 * coef(fit)[c("I(m1^2)", "I(m2^2)")]), ca$axes$lambda, tolerance = 1e-6)
  expect_output(print(ca), "permutations of fitness")
  # the permutation statistic is the double regression's F, t squared, for each axis
  Z <- as.matrix(ca$scores[, c("z1", "z2")])
  pairs <- utils::combn(2, 2)
  qx <- qr(cbind(1, Z, Z^2, Z[, 1] * Z[, 2]))
  expect_equal(unname(.canonical_f(Z, ca$scores$.w_rel, qx, pairs)),
               unname(summary(fit)$coefficients[c("I(m1^2)", "I(m2^2)"), "t value"]^2), tolerance = 1e-8)
  # the double-regression tests on request, and the same seed gives the same p-values
  dr <- quiet(canonical_analysis(df, "w", c("z1", "z2"), test = "double_regression"))
  expect_lt(dr$axes$p_value[2], 0.001)
  expect_output(print(dr), "much too small")
  set.seed(8)
  a <- quiet(canonical_analysis(df, "w", c("z1", "z2"), n_perm = 99))
  set.seed(8)
  b <- quiet(canonical_analysis(df, "w", c("z1", "z2"), n_perm = 99))
  expect_identical(a$axes$p_value, b$axes$p_value)
  expect_equal(a$n_perm, 99)
  expect_error(canonical_analysis(df, "w", c("z1", "z2"), n_perm = 0), "n_perm")
})

test_that("with survival unrelated to the traits the double regression finds curvature and the permutation test does not", {
  set.seed(33)
  d <- as.data.frame(matrix(rnorm(136 * 5), 136, 5, dimnames = list(NULL, paste0("t", 1:5))))
  d$alive <- rbinom(136, 1, 0.5)
  dr <- quiet(canonical_analysis(d, "alive", paste0("t", 1:5), test = "double_regression"))
  pm <- quiet(canonical_analysis(d, "alive", paste0("t", 1:5), n_perm = 199))
  expect_true(all(dr$axes$p_value[c(1, 5)] < 0.05))
  expect_true(all(pm$axes$p_value > 0.1))
})

test_that("survival works within groups, the double regression uses the logistic model, and the bootstrap adds intervals", {
  set.seed(2)
  ca <- quiet(canonical_analysis(bumpus, "survival", c("total_length", "weight", "humerus"),
                                 group = "sex", bootstrap = TRUE, n_boot = 30))
  expect_equal(ca$fitness_type, "binary")
  dr <- quiet(canonical_analysis(bumpus, "survival", c("total_length", "weight", "humerus"),
                                 group = "sex", test = "double_regression"))
  sc <- .add_glm_group(dr$scores, "sex")
  g <- glm(survival ~ .group + m1 + m2 + m3 + I(m1^2) + I(m2^2) + I(m3^2), family = binomial, data = sc)
  expect_equal(dr$axes$p_value, unname(summary(g)$coefficients[c("I(m1^2)", "I(m2^2)", "I(m3^2)"), 4]),
               tolerance = 1e-8)
  expect_equal(nrow(ca$axes), 3)
  expect_true(all(c("CI_lower", "CI_upper", "N_Boot") %in% names(ca$axes)))
  expect_true(all(ca$axes$CI_lower <= ca$axes$CI_upper))
  expect_true(all(ca$axes$N_Boot > 20))
  expect_true(all(c("m1", "m2", "m3") %in% names(ca$scores)))
  expect_true(all(diff(ca$axes$lambda) <= 0))
  expect_error(canonical_analysis(bumpus, "survival", "weight"), "at least two")

  expect_s3_class(plot_canonical_axes(ca, which = c(1, 3), grid_n = 20), "ggplot")
  expect_s3_class(plot_canonical_axes(ca, which = 3), "ggplot")
  expect_error(plot_canonical_axes(ca, which = 5), "canonical axes")
})

test_that("n_perm must be one whole number of 1 or more", {
  set.seed(4)
  d <- data.frame(z1 = rnorm(60), z2 = rnorm(60))
  d$w <- 1 + 0.1 * d$z1 + rnorm(60, 0, 0.1)
  for (bad in list(NA_real_, 0.4, 2.5, Inf, c(10, 20), "9")) {
    expect_error(canonical_analysis(d, "w", c("z1", "z2"), n_perm = bad), "whole number of 1 or more")
  }
})
