library(testthat)

test_that("adaptive_landscape calculates mean fitness on a grid", {
  set.seed(42)
  df <- data.frame(
    w = rnorm(50, 1, 0.1),
    z1 = rnorm(50),
    z2 = rnorm(50)
  )
  
  df$z1 <- as.numeric(scale(df$z1))
  df$z2 <- as.numeric(scale(df$z2))
  
  fit_model <- mgcv::gam(w ~ s(z1, z2), data = df)
  
  res <- adaptive_landscape(df, fit_model, c("z1", "z2"), simulation_n = 10, grid_n = 5)
  expect_s3_class(res, "adaptive_landscape")
  expect_true(".mean_fit" %in% names(res$grid))
  expect_equal(nrow(res$grid), 25)
})

test_that("a single trait gives a landscape curve with the fitness function alongside", {
  set.seed(7)
  df <- data.frame(z = as.numeric(scale(rnorm(80))))
  df$w <- 1 - 0.4 * df$z^2 + rnorm(80, 0, 0.1)
  fit_model <- mgcv::gam(w ~ s(z, k = 5), data = df)

  res <- adaptive_landscape(df, fit_model, "z", simulation_n = 20, grid_n = 12)
  expect_s3_class(res, "adaptive_landscape")
  expect_equal(nrow(res$grid), 12)
  expect_true(all(c("z", ".mean_fit", ".ind_fit") %in% names(res$grid)))
  expect_equal(names(res$data_summary$trait_ranges), "x1")
  # averaging over within-population variance flattens the peak
  expect_lt(max(res$grid$.mean_fit), max(res$grid$.ind_fit))
  expect_output(print(res), "Trait: z")

  p <- plot_adaptive_landscape(res, "z")
  expect_s3_class(p, "ggplot")
  expect_error(plot_adaptive_landscape_3d(res, "z"), "two traits")
})