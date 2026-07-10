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
