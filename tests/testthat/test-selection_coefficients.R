library(testthat)

test_that("selection_coefficients extracts linear, quadratic, and correlational terms", {
  set.seed(42)
  df <- data.frame(
    w = rnorm(50, 5, 1),
    z1 = rnorm(50),
    z2 = rnorm(50)
  )
  
  res <- selection_coefficients(df, "w", c("z1", "z2"), fitness_type = "continuous")

  expect_true(any(res$Type == "Linear"))
  expect_true(all(c("Beta_Coefficient", "Standard_Error", "P_Value") %in% names(res)))
})

test_that("count fitness is detected and tested with a count GLM", {
  set.seed(9)
  n <- 200
  df <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  df$kids <- rpois(n, exp(0.5 + 0.3 * df$z1))

  res <- suppressMessages(suppressWarnings(selection_coefficients(df, "kids", c("z1", "z2"))))
  expect_equal(attr(res, "fitness_type_used"), "count")
  expect_equal(attr(res, "model_family_used"), "poisson(log)")
  expect_true(res$P_Value[res$Term == "z1"] < 0.05)

  # the gradient is the OLS slope on relative fitness, whatever supplies the p-value
  z <- as.data.frame(scale(df[, c("z1", "z2")]))
  z$w <- df$kids / mean(df$kids)
  expect_equal(res$Beta_Coefficient[res$Term == "z1"], unname(coef(lm(w ~ z1 + z2, data = z))["z1"]))
})