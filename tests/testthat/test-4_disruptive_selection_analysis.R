library(testthat)

test_that("analyze_disruptive_selection computes linear and quadratic gradients", {
  set.seed(42)
  df <- data.frame(
    w = rnorm(50, 1, 0.1),
    z = rnorm(50)
  )
  
  res <- analyze_disruptive_selection(df, "w", "z", "continuous", standardize = FALSE)
  expect_true(all(c("Linear", "Quadratic") %in% res$Type))
  expect_equal(nrow(res), 2)
})

test_that("disruptive beta is the linear-only (Lande & Arnold) gradient, not the full-model coefficient", {
  # A skewed trait makes cov(z, z^2) != 0, so the linear-only and full-model
  # linear coefficients diverge. beta must come from the linear-only fit.
  set.seed(1)
  z <- as.numeric(scale(rexp(300)))
  df <- data.frame(w = 1 + 0.5 * z + 0.3 * z^2 + rnorm(300, sd = 0.5), z = z)

  res <- analyze_disruptive_selection(df, "w", "z", "continuous", standardize = FALSE)
  beta_dis <- res$Beta_Coefficient[res$Type == "Linear"]

  df$rel <- df$w / mean(df$w)
  beta_linear_only <- unname(coef(lm(rel ~ z, data = df))["z"])
  beta_full_model  <- unname(coef(lm(rel ~ z + I(z^2), data = df))["z"])

  expect_equal(beta_dis, beta_linear_only, tolerance = 1e-8)
  expect_gt(abs(beta_dis - beta_full_model), 0.1)
})

test_that("disruptive agrees with selection_coefficients on the same single trait", {
  set.seed(7)
  df <- data.frame(w = rnorm(80, 1, 0.1), z = rnorm(80))

  dis <- analyze_disruptive_selection(df, "w", "z", "continuous", standardize = FALSE)
  sc  <- selection_coefficients(df, "w", "z", fitness_type = "continuous", standardize = FALSE)

  expect_equal(dis$Beta_Coefficient, sc$Beta_Coefficient, tolerance = 1e-10)
})