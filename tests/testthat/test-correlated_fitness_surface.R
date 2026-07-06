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