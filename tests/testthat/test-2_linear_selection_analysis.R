library(testthat)

test_that("analyze_linear_selection works with continuous and binary data", {
  set.seed(42)
  df <- data.frame(
    w_cont = rnorm(50, 1, 0.1),
    w_bin = rbinom(50, 1, 0.5),
    z1 = rnorm(50),
    z2 = rnorm(50)
  )
  
  res_cont <- analyze_linear_selection(df, "w_cont", c("z1", "z2"), "continuous")
  expect_equal(res_cont$fitness_type, "continuous")
  
  res_bin <- analyze_linear_selection(df, "w_bin", c("z1", "z2"), "binary")
  expect_equal(res_bin$fitness_type, "binary")
})

test_that("count fitness takes its p-values from a Poisson or negative binomial GLM", {
  set.seed(11)
  n <- 200
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  mu <- exp(0.4 + 0.3 * z1)
  pois <- data.frame(kids = rpois(n, mu), z1 = z1, z2 = z2)

  res <- analyze_linear_selection(pois, "kids", c("z1", "z2"), "count")
  expect_equal(res$fitness_type, "count")
  expect_equal(res$glm_family, "poisson(log)")
  expect_s3_class(res$model$ols, "lm")
  expect_s3_class(res$model$glm, "glm")
  # gradients are the OLS slopes on relative fitness
  rel <- pois$kids / mean(pois$kids)
  expect_equal(unname(coef(res$model$ols)["z1"]), unname(coef(lm(rel ~ z1 + z2, data = pois))["z1"]))
  # p-values are the GLM's
  tab <- extract_linear_coefficients(c("z1", "z2"), res)
  expect_equal(tab$P_Value[tab$Term == "z1"], coef(summary(res$model$glm))["z1", "Pr(>|z|)"])

  over <- data.frame(kids = rnbinom(n, mu = mu, size = 0.8), z1 = z1, z2 = z2)
  res_nb <- analyze_linear_selection(over, "kids", c("z1", "z2"), "count")
  expect_equal(res_nb$glm_family, "negative binomial")
  expect_s3_class(res_nb$model$glm, "negbin")

  expect_error(
    analyze_linear_selection(data.frame(kids = c(0.5, 1.2, 2, 3, 1, 0.1, 2.2, 1.1, 0.4, 3.3), z1 = rnorm(10)), "kids", "z1", "count"),
    "count column"
  )
})