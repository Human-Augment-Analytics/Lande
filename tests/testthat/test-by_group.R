test_that("by_group fits each group on its own", {
  set.seed(6)
  n <- 240
  df <- data.frame(g = rep(c("a", "b"), each = n / 2), z1 = rnorm(n), z2 = rnorm(n))
  # group a peaks in z1, group b slopes
  df$w <- ifelse(df$g == "a", 1 - 0.5 * df$z1^2, 1 + 0.5 * df$z1) + rnorm(n, 0, 0.1)
  prep <- quiet(prepare_selection_data(df, "w", c("z1", "z2"), group = "g"))

  fits <- quiet(univariate_spline(prep, "w", "z1", group = "g", by_group = TRUE))
  expect_type(fits, "list")
  expect_equal(names(fits), c("a", "b"))
  expect_equal(attr(fits, "group"), "g")
  expect_true(all(vapply(fits, inherits, logical(1), "univariate_fitness")))
  # each fit saw only its own rows, and the shapes differ
  expect_equal(nrow(fits$a$model$model), n / 2)
  peak_a <- fits$a$grid$z1[which.max(fits$a$grid$fit)]
  peak_b <- fits$b$grid$z1[which.max(fits$b$grid$fit)]
  expect_lt(abs(peak_a), 0.6)
  expect_gt(peak_b, 1)
  # the same as fitting the group's rows by hand
  by_hand <- quiet(univariate_spline(prep[prep$g == "a", ], "w", "z1"))
  expect_equal(fits$a$grid$fit, by_hand$grid$fit)

  surfaces <- quiet(correlated_fitness_surface(prep, "w", c("z1", "z2"), grid_n = 15, method = "gam",
                                               group = "g", by_group = TRUE))
  expect_equal(names(surfaces), c("a", "b"))
  expect_true(all(vapply(surfaces, function(s) "peaks" %in% names(s), logical(1))))
  expect_null(surfaces$a$groups)
  expect_equal(nrow(surfaces$a$original_data), n / 2)
  for (s in surfaces) expect_s3_class(plot_correlated_fitness(s, c("z1", "z2")), "ggplot")

  expect_error(univariate_spline(prep, "w", "z1", by_group = TRUE), "needs a `group`")
  expect_error(correlated_fitness_surface(prep, "w", c("z1", "z2"), by_group = TRUE), "needs a `group`")
})

test_that("a group that cannot be fitted is left out by name", {
  set.seed(11)
  df <- data.frame(g = c(rep("big", 100), rep("tiny", 2)), z = rnorm(102))
  df$w <- 1 + 0.3 * df$z + rnorm(102, 0, 0.1)
  said <- character()
  fits <- withCallingHandlers(
    suppressMessages(univariate_spline(df, "w", "z", group = "g", by_group = TRUE)),
    warning = function(w) {
      said <<- c(said, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("Group 'tiny' left out", said)))
  expect_equal(names(fits), "big")
})
