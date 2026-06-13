# Tests for the reporting, bootstrap, and cubic-spline surface additions.

test_that("selection_report returns a consistent standardised table", {
  rep <- suppressWarnings(suppressMessages(
    selection_report(bumpus, "survival", c("weight", "total_length"),
                     fitness_type = "binary")
  ))
  expect_s3_class(rep, "selection_report")
  expect_true(all(c("Term", "Type", "Estimate", "Std_Error", "P_Value") %in% names(rep)))
  expect_true(any(rep$Type == "Differential"))
  expect_true(any(rep$Type == "Linear"))
  expect_true(any(rep$Type == "Quadratic"))
  expect_output(print(rep), "Selection analysis")
})

test_that("bootstrap_selection returns bootstrap SEs and percentile intervals", {
  set.seed(1)
  b <- suppressWarnings(suppressMessages(
    bootstrap_selection(bumpus, "survival", c("weight", "total_length"),
                        fitness_type = "binary", n_boot = 100)
  ))
  expect_true(all(c("Term", "Type", "Estimate", "Boot_SE", "CI_lower", "CI_upper", "P_Value")
                  %in% names(b)))
  expect_true(all(b$Boot_SE >= 0))
  expect_true(all(b$CI_lower <= b$CI_upper))
  expect_equal(attr(b, "n_boot"), 100)
})

test_that("univariate_spline uses a cubic spline with a bootstrapped ribbon", {
  set.seed(1)
  d <- data.frame(z = as.numeric(scale(rnorm(120))), surv = rbinom(120, 1, 0.5))
  u <- suppressWarnings(suppressMessages(
    univariate_spline(d, "surv", "z", fitness_type = "binary",
                      bootstrap = TRUE, n_boot = 40)
  ))
  expect_match(u$spline_type, "cubic")
  expect_match(u$ci_method, "bootstrap")
  expect_true(all(c("fit", "lwr", "upr") %in% names(u$grid)))
  expect_true(all(u$grid$lwr <= u$grid$upr))
})
