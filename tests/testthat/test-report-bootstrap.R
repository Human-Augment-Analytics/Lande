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
  # the p-values come from the logistic model, and the print says so
  expect_equal(unname(attr(rep, "p_model")), rep("binomial(logit)", 2))
  expect_output(print(rep), "p-values are from a logistic model on the same terms")
  cont <- suppressWarnings(suppressMessages(selection_report(bumpus, "weight", c("total_length", "humerus"))))
  expect_false(any(grepl("p-values", utils::capture.output(print(cont)))))

  # counts where the linear model needs a negative binomial and the quadratic one does not
  set.seed(5)
  d <- data.frame(z1 = rnorm(400), z2 = rnorm(400))
  d$kids <- rpois(400, exp(1.6 - 0.7 * d$z1^2 + 0.1 * d$z2))
  mixed <- suppressWarnings(suppressMessages(selection_report(d, "kids", c("z1", "z2"))))
  expect_equal(unname(attr(mixed, "p_model")), c("negative binomial", "poisson(log)"))
  expect_output(print(mixed), "negative binomial model for the linear gradients and a Poisson model")
})

test_that("bootstrap_selection returns bootstrap SEs and percentile intervals", {
  set.seed(1)
  b <- suppressWarnings(suppressMessages(
    bootstrap_selection(bumpus, "survival", c("weight", "total_length"),
                        fitness_type = "binary", n_boot = 100)
  ))
  expect_true(all(c("Term", "Type", "Estimate", "Boot_SE", "CI_lower", "CI_upper", "P_Value", "N_Boot")
                  %in% names(b)))
  expect_true(all(b$Boot_SE > 0))
  expect_true(all(b$CI_lower < b$CI_upper))
  # N_Boot is per-coefficient: an integer count no larger than the requested draws.
  expect_true(all(b$N_Boot <= 100))
  expect_true(all(b$N_Boot >= 2))
  expect_equal(attr(b, "n_boot"), min(b$N_Boot))

  # The point estimate is the ordinary fit, and the percentile interval covers it.
  ref <- suppressWarnings(suppressMessages(
    selection_coefficients(bumpus, "survival", c("weight", "total_length"), fitness_type = "binary")))
  expect_equal(b$Estimate, ref$Beta_Coefficient)
  expect_equal(b$P_Value, ref$P_Value)
  expect_true(all(b$CI_lower <= b$Estimate & b$Estimate <= b$CI_upper))
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

test_that("rows with no group label are resampled, as they are prepared", {
  set.seed(11)
  d <- data.frame(g = rep(c("a", "b", NA), c(60, 60, 20)), z1 = rnorm(140), z2 = rnorm(140))
  d$w <- 1 + 0.4 * d$z1 + rnorm(140, 0, 0.3)
  labelled <- d
  labelled$g[is.na(labelled$g)] <- "z"
  boot <- function(x) {
    set.seed(3)
    suppressWarnings(suppressMessages(bootstrap_selection(x, "w", c("z1", "z2"), group = "g", n_boot = 40)))
  }
  # a missing label behaves exactly like a third label that sorts last
  expect_equal(boot(d)$Boot_SE, boot(labelled)$Boot_SE)
  expect_equal(boot(d)$Estimate, boot(labelled)$Estimate)
})
