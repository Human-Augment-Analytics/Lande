library(testthat)

test_that("the assumption checks cover traits, collinearity, sample size and the model", {
  chk <- check_selection_assumptions(bumpus, "survival", c("total_length", "weight", "humerus"))
  expect_s3_class(chk, "selection_assumptions")
  expect_true(all(c("check", "statistic", "p_value", "note") %in% names(chk)))
  expect_true(any(grepl("Mardia skewness", chk$check)))
  expect_true(any(grepl("Mardia kurtosis", chk$check)))
  expect_equal(sum(grepl("Shapiro-Wilk", chk$check)), 3)
  expect_true(any(grepl("largest VIF", chk$check)))
  expect_true(any(grepl("Rarer outcome per quadratic term", chk$check)))
  expect_true(any(grepl("separation", chk$check)))
  p <- chk$p_value[!is.na(chk$p_value)]
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(attr(chk, "fitness_type"), "binary")
  expect_output(print(chk), "Assumption checks")
})

test_that("continuous fitness gets residual checks and heteroscedasticity is picked up", {
  set.seed(21)
  n <- 300
  d <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  # variance grows with the fitted value
  d$w <- 1 + 0.4 * d$z1 + rnorm(n, 0, 0.1 + 0.4 * (d$z1 - min(d$z1)))
  chk <- check_selection_assumptions(d, "w", c("z1", "z2"))
  expect_equal(attr(chk, "fitness_type"), "continuous")
  bp <- chk[grepl("Breusch-Pagan", chk$check), ]
  expect_equal(nrow(bp), 1)
  expect_lt(bp$p_value, 0.05)
  expect_true(any(grepl("residuals: normality", chk$check)))
  expect_true(any(grepl("Rows per quadratic term", chk$check)))
})

test_that("count fitness reports the dispersion ratio and normal traits pass Mardia", {
  set.seed(22)
  n <- 400
  d <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  d$kids <- MASS::rnegbin(n, mu = exp(0.5 + 0.3 * d$z1), theta = 1.2)
  chk <- check_selection_assumptions(d, "kids", c("z1", "z2"))
  expect_equal(attr(chk, "fitness_type"), "count")
  disp <- chk[grepl("dispersion ratio", chk$check), ]
  expect_gt(disp$statistic, 1.5)
  expect_match(disp$note, "overdispersed")
  m <- chk[grepl("Mardia", chk$check), ]
  expect_true(all(m$p_value > 0.01))
})

test_that("performance adds its tests and an R squared when installed", {
  skip_if_not_installed("performance")
  chk <- check_selection_assumptions(bumpus, "survival", c("total_length", "weight"))
  expect_true(any(grepl("Model fit: .*\\(performance\\)", chk$check)))
  set.seed(24)
  n <- 300
  d <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  d$kids <- MASS::rnegbin(n, mu = exp(0.5 + 0.3 * d$z1), theta = 1.2)
  cnt <- check_selection_assumptions(d, "kids", c("z1", "z2"))
  disp <- cnt[grepl("dispersion ratio", cnt$check), ]
  expect_false(is.na(disp$p_value))
  d$w <- 1 + 0.4 * d$z1 + rnorm(n, 0, 0.1 + 0.4 * (d$z1 - min(d$z1)))
  cont <- check_selection_assumptions(d, "w", c("z1", "z2"))
  expect_true(any(grepl("heteroscedasticity \\(performance\\)", cont$check)))
})

test_that("skewed traits fail Mardia and one trait skips the multivariate checks", {
  set.seed(23)
  n <- 300
  d <- data.frame(z1 = rexp(n), z2 = rexp(n))
  d$w <- rbinom(n, 1, 0.5)
  chk <- check_selection_assumptions(d, "w", c("z1", "z2"))
  expect_lt(chk$p_value[grepl("Mardia skewness", chk$check)], 0.001)
  one <- check_selection_assumptions(d, "w", "z1")
  expect_false(any(grepl("Mardia|VIF", one$check)))
  expect_error(check_selection_assumptions(d, "w", "z9"), "Missing columns")
})
