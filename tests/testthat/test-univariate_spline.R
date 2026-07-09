library(testthat)

test_that("univariate_spline generates valid splines", {
  set.seed(42)
  df <- data.frame(
    w = rnorm(50, 1, 0.1),
    z = rnorm(50)
  )
  
  df$z <- as.numeric(scale(df$z))
  
  res <- univariate_spline(df, "w", "z", fitness_type = "continuous", k = 3)
  expect_s3_class(res, "univariate_fitness")
  expect_true(all(c("fit", "lwr", "upr") %in% names(res$grid)))
})

test_that("count fitness is fitted on the raw counts with a Poisson family", {
  set.seed(5)
  df <- data.frame(z = as.numeric(scale(rnorm(120))))
  df$kids <- rpois(120, exp(0.3 + 0.4 * df$z))

  res <- suppressMessages(univariate_spline(df, "kids", "z", k = 4))
  expect_equal(res$fitness_type, "count")
  expect_equal(res$family, "poisson(log)")
  expect_true(all(res$grid$fit >= 0))
  expect_true(all(res$grid$lwr >= 0))
})