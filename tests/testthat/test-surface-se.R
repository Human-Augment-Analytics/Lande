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
  expect_true(all(d$p_value[2:3] < 0.001))
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

  # no dip on a slope
  df$lin <- 1 + 0.4 * df$z1 + rnorm(n, 0, 0.02)
  slope <- quiet(correlated_fitness_surface(df, "lin", c("z1", "z2"), grid_n = 20, method = "gam"))
  expect_message(one <- peak_difference(slope, c(-1, 0), c(1, 0), valley = TRUE), "No dip")
  expect_equal(nrow(one), 1)

  skip_if_not_installed("fields")
  tp <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 15, method = "tps"))
  expect_error(peak_difference(tp, c(-1, 0), c(1, 0)), "gam")
})
