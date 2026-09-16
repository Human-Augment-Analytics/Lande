library(testthat)

test_that("temporal_landscape fits each period and skips small ones", {
  set.seed(11)
  n <- 120
  d <- data.frame(year = rep(c(2001, 2002, 2003), each = n))
  d$z <- as.numeric(scale(rnorm(3 * n)))
  # the optimum moves to the right over the years
  opt <- c(-0.6, 0, 0.6)[match(d$year, c(2001, 2002, 2003))]
  d$w <- rbinom(3 * n, 1, plogis(1.5 - 2.5 * (d$z - opt)^2))
  d <- rbind(d, data.frame(year = 2004, z = rnorm(6), w = 1))

  tl <- suppressMessages(temporal_landscape(d, "w", "z", "year", min_n = 20, landscape = TRUE,
                                            simulation_n = 50, grid_n = 20))
  expect_s3_class(tl, "temporal_landscape")
  expect_equal(tl$summary$time, c(2001, 2002, 2003))
  expect_equal(tl$skipped, "2004")
  expect_equal(tl$summary$n, rep(n, 3))
  expect_true(all(c("mean_z", "edf", "optimum_z", "optimum_fit", "peaks", "landscape_optimum_z") %in% names(tl$summary)))
  expect_gt(tl$summary$optimum_z[3], tl$summary$optimum_z[1])
  expect_true(all(tl$summary$peaks >= 0))
  expect_true(all(c("time", "z", "fit") %in% names(tl$grid)))
  expect_equal(sort(unique(tl$heat$time)), c(2001, 2002, 2003))
  # the common grid is blank where a year has no individuals
  expect_true(anyNA(tl$heat$fit))
  expect_false(anyNA(tl$heat$fit[tl$heat$time == 2001 & tl$heat$z > -0.5 & tl$heat$z < 0.5]))
  expect_true(all(c(".mean_fit", "time") %in% names(tl$landscape_grid)))
  expect_output(print(tl), "Fitness function by year")

  for (ty in c("panels", "heatmap")) {
    p <- plot_temporal_landscape(tl, type = ty)
    expect_s3_class(p, "ggplot")
    expect_silent(b <- ggplot2::ggplot_build(p))
  }
  p <- plot_temporal_landscape(tl, show_landscape = FALSE, show_points = FALSE)
  expect_s3_class(p, "ggplot")

  expect_error(temporal_landscape(d, "w", "z", "year", min_n = 1000), "No period")
  expect_error(temporal_landscape(d, "w", "z", "season"), "Missing columns")
})

test_that("two traits give a surface per period", {
  set.seed(12)
  n <- 90
  d <- data.frame(year = rep(1:2, each = n))
  d$z1 <- as.numeric(scale(rnorm(2 * n)))
  d$z2 <- as.numeric(scale(rnorm(2 * n)))
  d$w <- 1 + 0.3 * d$z1 - 0.2 * d$z2^2 + rnorm(2 * n, 0, 0.2)

  tl <- suppressMessages(temporal_landscape(d, "w", c("z1", "z2"), "year", landscape = FALSE, grid_n = 12))
  expect_equal(nrow(tl$summary), 2)
  expect_true(all(c("optimum_z1", "optimum_z2", "mean_z1", "mean_z2") %in% names(tl$summary)))
  expect_false("peaks" %in% names(tl$summary))
  expect_true(all(c("z1", "z2", ".fit", "time") %in% names(tl$grid)))
  expect_null(tl$heat)
  expect_null(tl$landscape_grid)
  expect_length(tl$hulls, 2)

  p <- plot_temporal_landscape(tl)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_gt(length(b$data), 2)
  expect_error(plot_temporal_landscape(tl, type = "heatmap"), "one trait")
})
