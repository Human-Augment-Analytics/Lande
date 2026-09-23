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
  expect_lt(last$p_value, 0.001)
  expect_equal(abs(ca$M[, 2]), c(z1 = 1, z2 = 1) / sqrt(2), tolerance = 0.1)
  expect_lt(abs(ca$axes$lambda[1]), 0.1)
  # the double regression on the scores gives the same curvature as the eigenvalue
  fit <- lm(.w_rel ~ m1 + m2 + I(m1^2) + I(m2^2), data = ca$scores)
  expect_equal(unname(2 * coef(fit)[c("I(m1^2)", "I(m2^2)")]), ca$axes$lambda, tolerance = 1e-6)
  expect_output(print(ca), "anticonservative")
})

test_that("survival takes its p-values from the logistic model and the bootstrap adds intervals", {
  set.seed(2)
  ca <- quiet(canonical_analysis(bumpus, "survival", c("total_length", "weight", "humerus"),
                                 group = "sex", bootstrap = TRUE, n_boot = 30))
  expect_equal(ca$fitness_type, "binary")
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
