test_that("a monotone surface has no interior peak", {
  set.seed(4)
  n <- 300
  df <- data.frame(z1 = as.numeric(scale(rnorm(n))), z2 = as.numeric(scale(rnorm(n))))
  df$g <- rep(c("a", "b"), each = n / 2)
  df$w <- 1 + 0.4 * df$z1 + rnorm(n, 0, 0.02)
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 25, method = "gam", group = "g"))
  expect_true(all(c("peaks", "groups") %in% names(s)))
  expect_equal(names(s$peaks), c("z1", "z2", "fit", "interior"))
  expect_gte(nrow(s$peaks), 1)
  expect_equal(sum(s$peaks$interior), 0)
  expect_true(all(c("peak_interior", "peak_edge") %in% names(s$groups)))
  expect_false(any(s$groups$peak_interior))
  expect_output(print(s), "0 interior peaks")
  # open marks, since nothing here is a peak
  p <- plot_correlated_fitness(s, c("z1", "z2"))
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_true(any(vapply(built$data, function(d) "shape" %in% names(d) && any(d$shape %in% c(5, 2)), logical(1))))
})

test_that("two bumps give two interior peaks", {
  set.seed(17)
  n <- 400
  df <- data.frame(z1 = runif(n, -2, 2), z2 = runif(n, -2, 2))
  df$w <- exp(-((df$z1 - 1)^2 + df$z2^2) / 0.5) + exp(-((df$z1 + 1)^2 + df$z2^2) / 0.5) + rnorm(n, 0, 0.02)
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 30, method = "gam", k = 30))
  expect_equal(sum(s$peaks$interior), 2)
  expect_equal(sort(round(s$peaks$z1[s$peaks$interior])), c(-1, 1))
  expect_true(s$peaks$interior[1])
  expect_output(print(s), "2 interior peaks")
  expect_s3_class(plot_correlated_fitness(s, c("z1", "z2")), "ggplot")
})

test_that("the finch community has three interior peaks and an edge maximum", {
  finch <- read.csv(system.file("extdata", "finch_community.csv", package = "Lande"))
  tr <- c("beak_length", "beak_depth")
  prep <- quiet(prepare_selection_data(finch, "recaptures", tr))
  s <- quiet(correlated_fitness_surface(prep, "recaptures", tr, grid_n = 50, too_far = 0.15,
                                        group = "species", group_effect = FALSE))
  g <- s$groups
  # fuliginosa's shallow high comes and goes with the grid, so it is not asserted
  expect_true(all(g$peak_interior[g$group %in% c("fortis large", "fortis small", "scandens")]))
  expect_false(g$peak_interior[g$group == "magnirostris"])
  expect_gte(sum(s$peaks$interior), 3)
  # magnirostris: the highest cell is not among its birds
  expect_true(any(!s$peaks$interior))
  expect_s3_class(plot_correlated_fitness_enhanced(s, tr, original_data = prep, fitness_col = "recaptures"), "ggplot")
})

test_that("temporal summaries say when the highest fitness sits at the edge", {
  set.seed(9)
  n <- 200
  d <- data.frame(year = rep(c(1, 2), each = n))
  d$z <- as.numeric(scale(rnorm(2 * n)))
  # year 1 is a slope, year 2 has a peak
  d$w <- ifelse(d$year == 1, rbinom(2 * n, 1, plogis(1.5 * d$z)), rbinom(2 * n, 1, plogis(1 - 3 * d$z^2)))
  tl <- quiet(temporal_landscape(d, "w", "z", "year", landscape = FALSE))
  expect_true("optimum_edge" %in% names(tl$summary))
  expect_true(tl$summary$optimum_edge[1])
  expect_false(tl$summary$optimum_edge[2])
  expect_output(print(tl), "edge of the data in 1 of 2")
  for (ty in c("panels", "heatmap")) expect_s3_class(plot_temporal_landscape(tl, type = ty), "ggplot")
  expect_s3_class(plot_temporal_landscape(tl, type = "heatmap", connect = TRUE), "ggplot")
})
