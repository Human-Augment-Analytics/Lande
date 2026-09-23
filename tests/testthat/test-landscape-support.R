test_that("the landscape counts how much of each simulated population leaves the data", {
  set.seed(12)
  n <- 200
  df <- data.frame(z1 = as.numeric(scale(rnorm(n))), z2 = as.numeric(scale(rnorm(n))))
  df$w <- 1 + 0.3 * df$z1 + rnorm(n, 0, 0.05)
  fit <- mgcv::gam(w ~ s(z1, z2, k = 10), data = df)

  # tiny variance: all inside at the centre, all outside at a corner
  tight <- suppressMessages(adaptive_landscape(df, fit, c("z1", "z2"), population_variance = diag(1e-4, 2),
                                               simulation_n = 50, grid_n = 9))
  g <- tight$grid
  expect_equal(g$.outside[which.min(g$z1^2 + g$z2^2)], 0)
  expect_equal(g$.outside[which.max(g$z1 + g$z2)], 1)

  # with the observed covariance a population at the data mean mostly stays inside
  land <- suppressMessages(adaptive_landscape(df, fit, c("z1", "z2"), simulation_n = 300, grid_n = 9))
  g <- land$grid
  expect_lt(g$.outside[which.min(g$z1^2 + g$z2^2)], 0.15)
  expect_gt(max(g$.outside), 0.9)
  expect_equal(names(land$support), c("outside", "at_optimum", "warn"))
  expect_equal(land$support$outside, mean(g$.outside))
  # fitness keeps rising, so the optimum is outside the data
  expect_gt(land$support$at_optimum, 0.25)
  expect_output(print(land), "extrapolation")
  expect_message(adaptive_landscape(df, fit, c("z1", "z2"), simulation_n = 50, grid_n = 5), "outside the data")
  expect_s3_class(plot_adaptive_landscape(land, c("z1", "z2"), show_support = TRUE), "ggplot")
})

test_that("one trait counts individuals beyond the observed range", {
  set.seed(5)
  df <- data.frame(z = as.numeric(scale(rnorm(150))))
  df$w <- 1 - 0.3 * df$z^2 + rnorm(150, 0, 0.05)
  fit <- mgcv::gam(w ~ s(z, k = 6), data = df)
  land <- suppressMessages(adaptive_landscape(df, fit, "z", simulation_n = 300, grid_n = 21))
  g <- land$grid
  expect_lt(g$.outside[which.min(abs(g$z))], 0.05)
  expect_gt(g$.outside[1], 0.5)
  # optimum in the middle of the data
  expect_lt(land$support$at_optimum, 0.25)
  expect_s3_class(plot_adaptive_landscape(land, "z", show_support = TRUE), "ggplot")
})
