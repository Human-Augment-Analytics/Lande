library(testthat)

test_that("correlated_fitness_surface computes GAM grids", {
  set.seed(42)
  df <- data.frame(
    w = rnorm(50, 1, 0.1),
    z1 = rnorm(50),
    z2 = rnorm(50)
  )
  
  df$z1 <- as.numeric(scale(df$z1))
  df$z2 <- as.numeric(scale(df$z2))
  
  res <- correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 10, method = "gam")
  expect_equal(res$method, "gam")
  expect_true(".fit" %in% names(res$grid))
})

test_that("the basis dimension follows the data unless overridden", {
  set.seed(3)
  n <- 120
  df <- data.frame(z1 = as.numeric(scale(rnorm(n))), z2 = as.numeric(scale(rnorm(n))))
  df$w <- 1 + 0.3 * df$z1 - 0.2 * df$z1^2 + rnorm(n, 0, 0.2)

  auto <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 8, method = "gam"))
  expect_equal(auto$k, min(30, max(10, floor(sqrt(n * n)))))

  fixed <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 8, method = "gam", k = 12))
  expect_equal(fixed$k, 12)
})

test_that("the surface basis and smoothing method can be chosen", {
  set.seed(8)
  n <- 160
  df <- data.frame(z1 = as.numeric(scale(rnorm(n))), z2 = as.numeric(scale(rnorm(n))))
  df$w <- 1 + 0.3 * df$z1 - 0.2 * df$z2 - 0.15 * df$z1 * df$z2 + rnorm(n, 0, 0.2)

  default <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 10, method = "gam"))
  expect_equal(default$basis, "tp")
  expect_equal(default$smoothing, "REML")
  expect_equal(default$model$method, "REML")

  cr <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 10, method = "gam", bs = "cr", smoothing = "GCV.Cp"))
  expect_equal(cr$basis, "cr")
  expect_true(grepl("te\\(", deparse(formula(cr$model))[1]))
  ok <- default$grid$.inside
  expect_gt(cor(default$grid$.fit[ok], cr$grid$.fit[ok]), 0.98)

  expect_error(correlated_fitness_surface(df, "w", c("z1", "z2"), method = "gam", smoothing = "banana"), "should be one of")
})

test_that("grid points outside the data are masked", {
  set.seed(4)
  n <- 150
  # a diagonal cloud leaves the corners of the grid empty
  z1 <- rnorm(n)
  z2 <- as.numeric(scale(0.9 * z1 + rnorm(n, 0, 0.3)))
  df <- data.frame(z1 = as.numeric(scale(z1)), z2 = z2)
  df$w <- 1 + 0.2 * df$z1 + rnorm(n, 0, 0.2)

  masked <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 15, method = "gam"))
  expect_true(all(c(".inside", ".fit_all") %in% names(masked$grid)))
  expect_true(any(!masked$grid$.inside))
  expect_true(all(is.na(masked$grid$.fit[!masked$grid$.inside])))
  expect_false(anyNA(masked$grid$.fit[masked$grid$.inside]))
  expect_false(anyNA(masked$grid$.fit_all))
  # the hull is a closed polygon in trait units
  expect_equal(names(masked$hull), c("z1", "z2"))
  expect_equal(unlist(masked$hull[1, ]), unlist(masked$hull[nrow(masked$hull), ]))
  # the corners are extrapolation
  corner <- masked$grid$z1 == min(masked$grid$z1) & masked$grid$z2 == max(masked$grid$z2)
  expect_false(masked$grid$.inside[corner])

  open <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 15, method = "gam", mask = FALSE))
  expect_false(anyNA(open$grid$.fit))
  expect_true(all(open$grid$.inside))
  expect_null(open$hull)

  p <- plot_correlated_fitness(masked, c("z1", "z2"))
  expect_s3_class(p, "ggplot")
})
test_that("a distance rule blanks grid points far from any individual", {
  set.seed(5)
  n <- 120
  # two separate clusters: the hull bridges the gap between them, a distance rule does not
  z1 <- c(rnorm(n / 2, -1.5, 0.3), rnorm(n / 2, 1.5, 0.3))
  z2 <- c(rnorm(n / 2, -1.5, 0.3), rnorm(n / 2, 1.5, 0.3))
  df <- data.frame(z1 = as.numeric(scale(z1)), z2 = as.numeric(scale(z2)))
  df$w <- 1 + 0.2 * df$z1 + rnorm(n, 0, 0.2)

  hull_only <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam"))
  far <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam", too_far = 0.15))
  expect_true(".dist" %in% names(far$grid))
  expect_true(all(far$grid$.dist >= 0))
  expect_equal(far$too_far, 0.15)
  # the middle of the gap is inside the hull but far from every individual
  mid <- which.min(far$grid$z1^2 + far$grid$z2^2)
  expect_true(hull_only$grid$.inside[mid])
  expect_false(far$grid$.inside[mid])
  expect_true(is.na(far$grid$.fit[mid]))
  expect_lt(sum(far$grid$.inside), sum(hull_only$grid$.inside))
  # the points that are kept carry the same predictions
  ok <- far$grid$.inside
  expect_equal(far$grid$.fit[ok], hull_only$grid$.fit[ok])
  # a loose rule adds nothing to the hull
  loose <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam", too_far = 2))
  expect_equal(loose$grid$.inside, hull_only$grid$.inside)
  # the rule works on its own
  only_far <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam", mask = FALSE, too_far = 0.15))
  expect_null(only_far$hull)
  expect_true(any(is.na(only_far$grid$.fit)))
  expect_error(correlated_fitness_surface(df, "w", c("z1", "z2"), method = "gam", too_far = -1), "too_far")

  for (s in list(far, only_far)) {
    p <- plot_correlated_fitness(s, c("z1", "z2"))
    expect_s3_class(p, "ggplot")
    b <- ggplot2::ggplot_build(p)
    expect_gt(length(b$data), 2)
  }
})

test_that("a grouped surface reports each group's mean and local peak", {
  set.seed(6)
  n <- 160
  sp <- rep(c("a", "b"), each = n / 2)
  z1 <- c(rnorm(n / 2, -1, 0.4), rnorm(n / 2, 1, 0.4))
  z2 <- c(rnorm(n / 2, 1, 0.4), rnorm(n / 2, -1, 0.4))
  df <- data.frame(sp = sp, z1 = as.numeric(scale(z1)), z2 = as.numeric(scale(z2)))
  df$w <- 1 - 0.3 * df$z2^2 + rnorm(n, 0, 0.1)

  s <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 25, method = "gam", group = "sp"))
  g <- s$groups
  expect_equal(names(g), c("group", "n", "mean_z1", "mean_z2", "peak_z1", "peak_z2", "peak_fit"))
  expect_equal(g$group, c("a", "b"))
  expect_equal(g$n, c(n / 2, n / 2))
  expect_equal(g$mean_z1, as.numeric(tapply(df$z1, df$sp, mean)))
  expect_equal(g$mean_z2, as.numeric(tapply(df$z2, df$sp, mean)))
  # each peak sits in its own group's cloud, on the kept part of the surface
  expect_false(anyNA(g$peak_fit))
  expect_lt(g$peak_z1[1], 0)
  expect_gt(g$peak_z1[2], 0)
  expect_true(all(g$peak_fit <= max(s$grid$.fit, na.rm = TRUE)))
  # the analysed rows come back with the group
  expect_equal(nrow(s$original_data), n)
  expect_true(all(c("z1", "z2", "w", "sp") %in% names(s$original_data)))
  # no group, no table
  u <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 25, method = "gam"))
  expect_null(u$groups)

  p <- plot_correlated_fitness(s, c("z1", "z2"), show_points = TRUE)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  labelled <- vapply(b$data, function(d) "label" %in% names(d) && all(c("a", "b") %in% d$label), logical(1))
  expect_true(any(labelled))
  p0 <- plot_correlated_fitness(s, c("z1", "z2"), show_groups = FALSE)
  b0 <- ggplot2::ggplot_build(p0)
  expect_lt(length(b0$data), length(b$data))
  # the dashed line from mean to peak can be dropped on its own
  seg <- function(p) sum(vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1)))
  expect_equal(seg(p), 1)
  expect_equal(seg(plot_correlated_fitness(s, c("z1", "z2"), group_lines = FALSE)), 0)
})

test_that("the group can mark the surface without entering the model", {
  set.seed(7)
  n <- 160
  sp <- rep(c("a", "b"), each = n / 2)
  z1 <- c(rnorm(n / 2, -1, 0.4), rnorm(n / 2, 1, 0.4))
  z2 <- rnorm(n)
  df <- data.frame(sp = sp, z1 = as.numeric(scale(z1)), z2 = as.numeric(scale(z2)))
  df$w <- 1 + 0.3 * df$z1 - 0.2 * df$z2^2 + rnorm(n, 0, 0.1)

  with_effect <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam", group = "sp"))
  pooled <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam", group = "sp", group_effect = FALSE))
  expect_true(with_effect$group_effect)
  expect_false(pooled$group_effect)
  expect_true(grepl("sp", deparse(formula(with_effect$model))[1]))
  expect_false(grepl("sp", deparse(formula(pooled$model))[1]))
  # the pooled fit is the same as fitting without a group at all
  none <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam"))
  expect_equal(pooled$grid$.fit, none$grid$.fit)
  # the overlay is there either way
  expect_equal(nrow(pooled$groups), 2)
  expect_equal(pooled$groups$group, with_effect$groups$group)
  expect_null(none$group_effect)
})

test_that("count fitness gets a Poisson surface", {
  set.seed(9)
  n <- 200
  df <- data.frame(z1 = as.numeric(scale(rnorm(n))), z2 = as.numeric(scale(rnorm(n))))
  df$w <- rpois(n, exp(0.2 + 0.4 * df$z1 - 0.3 * df$z2^2))

  s <- suppressMessages(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 15))
  expect_equal(s$data_type, "count")
  expect_equal(s$method, "gam")
  expect_equal(s$model$family$family, "poisson")
  expect_true(all(s$grid$.fit[s$grid$.inside] >= 0))

  # relative fitness is not a count
  rel <- df
  rel$w <- rel$w / mean(rel$w)
  cont <- suppressMessages(correlated_fitness_surface(rel, "w", c("z1", "z2"), grid_n = 15, method = "gam"))
  expect_equal(cont$data_type, "continuous")
  expect_equal(cont$model$family$family, "gaussian")

  expect_warning(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 15, method = "tps"), "count")
})

test_that("a count fitness with two values is still fitted as a count", {
  set.seed(25)
  n <- 150
  df <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  df$kids <- ifelse(df$z1 + rnorm(n) > 0, 5L, 2L)
  s <- suppressWarnings(suppressMessages(correlated_fitness_surface(df, "kids", c("z1", "z2"), grid_n = 15)))
  expect_equal(s$data_type, "count")
  expect_equal(detect_family(df$kids)$type, "count")
})
