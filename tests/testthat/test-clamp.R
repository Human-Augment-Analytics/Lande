test_that("clamping holds fitness inside the range of its type", {
  expect_equal(.clamp_fitness(c(-0.2, 0.5, 1.3), "binary"), c(0, 0.5, 1))
  expect_equal(.clamp_fitness(c(-2, 3), "count"), c(0, 3))
  expect_equal(.clamp_fitness(c(-2, 3), "continuous"), c(-2, 3))
})

test_that("the thin-plate surface and landscape stay within 0 and 1 for survival", {
  skip_if_not_installed("fields")
  set.seed(3)
  df <- data.frame(z1 = as.numeric(scale(rnorm(120))), z2 = as.numeric(scale(rnorm(120))))
  df$w <- rbinom(120, 1, plogis(2 * df$z1))

  held <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 12, method = "tps", mask = FALSE))
  free <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 12, method = "tps", mask = FALSE, clamp = FALSE))
  expect_true(held$clamp)
  expect_false(free$clamp)
  expect_true(all(held$grid$.fit >= 0 & held$grid$.fit <= 1))
  expect_true(any(free$grid$.fit < 0 | free$grid$.fit > 1))

  set.seed(3)
  land <- quiet(adaptive_landscape(df, held$model, c("z1", "z2"), simulation_n = 50, grid_n = 6))
  set.seed(3)
  raw <- quiet(adaptive_landscape(df, held$model, c("z1", "z2"), simulation_n = 50, grid_n = 6, clamp = FALSE))
  expect_true(all(land$grid$.mean_fit >= 0 & land$grid$.mean_fit <= 1))
  expect_gt(land$clipped, 0)
  expect_equal(raw$clipped, 0L)
  expect_gte(max(raw$grid$.mean_fit), max(land$grid$.mean_fit))
})

test_that("a GAM surface is not clamped", {
  set.seed(5)
  df <- data.frame(z1 = as.numeric(scale(rnorm(80))), z2 = as.numeric(scale(rnorm(80))))
  df$w <- rbinom(80, 1, 0.4)
  fit <- mgcv::gam(w ~ s(z1, z2, k = 10), family = binomial(), data = df)
  land <- suppressMessages(adaptive_landscape(df, fit, c("z1", "z2"), simulation_n = 20, grid_n = 4))
  expect_false(land$clamp)
  expect_equal(land$clipped, 0L)
})
