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