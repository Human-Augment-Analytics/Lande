test_that("a GAM surface carries its standard error and band", {
  set.seed(8)
  n <- 300
  df <- data.frame(z1 = as.numeric(scale(rnorm(n))), z2 = as.numeric(scale(rnorm(n))))
  df$w <- 1 + 0.3 * df$z1 - 0.2 * df$z2^2 + rnorm(n, 0, 0.3)
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 25, method = "gam"))
  g <- s$grid
  kept <- !is.na(g$.fit)
  expect_true(all(c(".se", ".fit_lo", ".fit_hi", ".se_all", ".fit_lo_all", ".fit_hi_all") %in% names(g)))
  expect_true(all(g$.se[kept] > 0))
  expect_true(all(is.na(g$.se[!kept])))
  expect_false(anyNA(g$.se_all))
  expect_true(all(g$.fit_lo[kept] <= g$.fit[kept] & g$.fit[kept] <= g$.fit_hi[kept]))
  # larger away from the middle
  d2 <- g$z1^2 + g$z2^2
  centre <- kept & d2 < 0.25
  edge <- kept & d2 > stats::quantile(d2[kept], 0.9)
  expect_gt(mean(g$.se[edge]), mean(g$.se[centre]))
  # a higher level gives a wider band
  expect_equal(s$level, 0.95)
  s99 <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 25, method = "gam", level = 0.99))
  expect_true(all((s99$grid$.fit_hi - s99$grid$.fit_lo)[kept] > (g$.fit_hi - g$.fit_lo)[kept]))
  expect_output(print(s), "Standard error")

  # survival: the band stays between 0 and 1
  df$y <- rbinom(n, 1, plogis(0.8 * df$z1))
  b <- quiet(correlated_fitness_surface(df, "y", c("z1", "z2"), grid_n = 20, method = "gam"))
  kb <- !is.na(b$grid$.fit)
  expect_true(all(b$grid$.fit_lo[kb] >= 0 & b$grid$.fit_hi[kb] <= 1))
})

test_that("the plots draw the uncertainty", {
  set.seed(3)
  n <- 250
  df <- data.frame(z1 = as.numeric(scale(rnorm(n))), z2 = as.numeric(scale(rnorm(n))))
  df$w <- 1 + 0.3 * df$z1 - 0.2 * df$z2^2 + rnorm(n, 0, 0.3)
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam"))
  panels <- function(p) length(unique(ggplot2::ggplot_build(p)$layout$layout$PANEL))

  p_se <- plot_correlated_fitness(s, c("z1", "z2"), uncertainty = "se")
  expect_s3_class(p_se, "ggplot")
  expect_equal(panels(p_se), 1)
  expect_equal(length(p_se$layers), length(plot_correlated_fitness(s, c("z1", "z2"))$layers) + 1)
  expect_equal(panels(plot_correlated_fitness(s, c("z1", "z2"), uncertainty = "band")), 3)
  p_en <- quiet(plot_correlated_fitness_enhanced(s, c("z1", "z2"), original_data = df, fitness_col = "w", uncertainty = "band"))
  expect_equal(panels(p_en), 3)
  expect_s3_class(quiet(plot_correlated_fitness_enhanced(s, c("z1", "z2"), uncertainty = "se")), "ggplot")

  # no standard errors from the thin-plate spline, so only the fit is drawn
  skip_if_not_installed("fields")
  tp <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 15, method = "tps"))
  expect_true(all(is.na(tp$grid$.se)))
  expect_warning(p_tp <- plot_correlated_fitness(tp, c("z1", "z2"), uncertainty = "band"), "standard errors")
  expect_equal(panels(p_tp), 1)
})

test_that("peak_difference compares two bumps with the dip between them", {
  set.seed(26)
  n <- 360
  df <- data.frame(z1 = runif(n, -2, 2), z2 = runif(n, -2, 2))
  bump <- function(at) exp(-((df$z1 - at)^2 + df$z2^2) / 0.5)
  df$w <- bump(-1) + bump(1) + rnorm(n, 0, 0.03)
  df$g <- ifelse(df$z1 < 0, "a", "b")
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 30, method = "gam", k = 30,
                                        group = "g", group_effect = FALSE))

  d <- peak_difference(s, c(-1, 0), c(1, 0), valley = TRUE)
  expect_equal(names(d), c("comparison", "fit_a", "fit_b", "difference", "se", "z", "p_value"))
  expect_equal(nrow(d), 3)
  # bumps of the same height, both well above the dip between them
  expect_lt(abs(d$difference[1]), 0.1)
  expect_true(all(d$difference[2:3] > 0.5))
  expect_true(all(d$z[2:3] > 5))
  # a point sits above the valley on its own route by construction, so no p-value there
  expect_true(all(is.na(d$p_value[2:3])))
  expect_false(is.na(d$p_value[1]))
  expect_lt(abs(attr(d, "valley")[["z1"]]), 0.3)

  # group names stand for the groups' highest points
  by_name <- peak_difference(s, "a", "b", valley = TRUE)
  expect_equal(nrow(by_name), 3)
  expect_match(by_name$comparison[1], "a - b")
  expect_equal(nrow(peak_difference(s, "a", "b")), 1)
  expect_error(peak_difference(s, "a", "nobody"), "No group")

  # with the group in the model the comparison is made at the reference level
  eff <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam", k = 30, group = "g"))
  expect_equal(nrow(peak_difference(eff, c(-1, 0), c(1, 0))), 1)

  # no dip on a slope, whichever grid: a point that sits above the cell it
  # starts from is not a valley
  df$lin <- 1 + 0.4 * df$z1 + rnorm(n, 0, 0.02)
  for (g in c(17, 18, 20, 21, 22)) {
    slope <- quiet(correlated_fitness_surface(df, "lin", c("z1", "z2"), grid_n = g, method = "gam"))
    expect_message(one <- peak_difference(slope, c(-1, 0), c(1, 0), valley = TRUE), "No dip")
    expect_equal(nrow(one), 1)
  }

  skip_if_not_installed("fields")
  tp <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 15, method = "tps"))
  expect_error(peak_difference(tp, c(-1, 0), c(1, 0)), "gam")
})

test_that("the pass goes round a hollow that the straight line drops into", {
  set.seed(47)
  n <- 900
  df <- data.frame(z1 = runif(n, -2, 2), z2 = runif(n, -2, 2))
  r <- sqrt(df$z1^2 + df$z2^2)
  th <- atan2(df$z2, df$z1)
  # a ring, highest at (-1.2, 0) and (1.2, 0), lower where it crosses z2 = 0
  df$w <- 1 + exp(-(r - 1.2)^2 / 0.1) * (1 + 0.5 * cos(2 * th)) + rnorm(n, 0, 0.05)
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 30, method = "gam", k = 30))
  pass <- peak_difference(s, c(-1.2, 0), c(1.2, 0), valley = TRUE)
  line <- peak_difference(s, c(-1.2, 0), c(1.2, 0), valley = TRUE, route = "line")
  expect_lt(sqrt(sum(attr(line, "valley")^2)), 0.3)
  expect_gt(abs(attr(pass, "valley")[["z2"]]), 0.9)
  expect_gt(pass$fit_b[2], line$fit_b[2] + 0.5)
  expect_true(all(pass$difference[2:3] < line$difference[2:3]))

  # the hollow in the middle is lower than any pass, so there is no dip to compare it with
  expect_message(hollow <- peak_difference(s, c(-1.2, 0), c(0, 0), valley = TRUE), "No dip")
  expect_equal(nrow(hollow), 1)
  # a point well off the grid starts from the nearest kept cell, and says so
  said <- testthat::capture_messages(peak_difference(s, c(-1.2, 0), c(5, 0), valley = TRUE))
  expect_true(any(grepl("away from the kept surface", said)))

  # two clouds of birds with nothing kept between them have no pass
  far <- data.frame(z1 = c(rnorm(150, -1.5, 0.25), rnorm(150, 1.5, 0.25)), z2 = rnorm(300, 0, 0.3))
  far$w <- 1 + exp(-far$z2^2) + rnorm(300, 0, 0.1)
  split <- quiet(correlated_fitness_surface(far, "w", c("z1", "z2"), grid_n = 30, method = "gam", too_far = 0.08))
  expect_message(apart <- peak_difference(split, c(-1.5, 0), c(1.5, 0), valley = TRUE), "not joined")
  expect_equal(nrow(apart), 1)
})
