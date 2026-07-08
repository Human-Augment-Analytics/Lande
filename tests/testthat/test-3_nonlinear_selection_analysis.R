library(testthat)

test_that("analyze_nonlinear_selection fits full quadratic model", {
  set.seed(42)
  df <- data.frame(
    w = rnorm(50, 1, 0.1),
    z1 = rnorm(50),
    z2 = rnorm(50)
  )
  
  res <- analyze_nonlinear_selection(df, "w", c("z1", "z2"), "continuous")
  expect_true(any(grepl("I\\(z1\\^2\\)", rownames(res$summary$coefficients))))
})

test_that("nonlinear terms for count fitness carry GLM p-values", {
  set.seed(12)
  n <- 200
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  df <- data.frame(kids = rpois(n, exp(0.5 - 0.2 * z1^2)), z1 = z1, z2 = z2)

  res <- analyze_nonlinear_selection(df, "kids", c("z1", "z2"), "count")
  expect_equal(res$fitness_type, "count")
  expect_true(!is.null(res$summary$glm))
  quad <- extract_quadratic_coefficients(c("z1", "z2"), res)
  expect_equal(quad$Beta_Coefficient[1], 2 * coef(res$model$ols)[["I(z1^2)"]])
  expect_equal(quad$P_Value[1], coef(summary(res$model$glm))["I(z1^2)", "Pr(>|z|)"])
})