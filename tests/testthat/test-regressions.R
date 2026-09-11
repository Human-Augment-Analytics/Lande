# Regression tests for behaviour that broke, or nearly broke, during the
# package restructure. Each block names the failure it guards against.

sim_data <- function(n = 200, seed = 11) {
  set.seed(seed)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  data.frame(
    w = 5 + 1.5 * z1 + 0.4 * z1^2 + rnorm(n, 0, 0.5),
    z1 = z1, z2 = z2,
    grp = rep(c("A", "B"), length.out = n)
  )
}

test_that("direct analyzer calls pick up relative fitness from prepared data", {
  df <- sim_data()
  prep <- suppressWarnings(suppressMessages(prepare_selection_data(df, "w", c("z1", "z2"))))

  ref <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "w", c("z1", "z2"), fitness_type = "continuous")))
  lin <- suppressWarnings(suppressMessages(
    analyze_linear_selection(prep, "w", c("z1", "z2"), "continuous")))
  nl <- suppressWarnings(suppressMessages(
    analyze_nonlinear_selection(prep, "w", c("z1", "z2"), "continuous")))

  expect_equal(unname(coef(lin$model)["z1"]),
               ref$Beta_Coefficient[ref$Term == "z1"], tolerance = 1e-10)
  expect_equal(2 * unname(coef(nl$model)["I(z1^2)"]),
               ref$Beta_Coefficient[ref$Term == "z1²"], tolerance = 1e-10)
})

test_that("a direct binary call with a non-0/1 column is rejected clearly", {
  d <- data.frame(surv = rbinom(60, 1, 0.5), z = rnorm(60))
  prep <- suppressWarnings(suppressMessages(prepare_selection_data(d, "surv", "z")))
  expect_error(
    analyze_linear_selection(prep, "relative_fitness", "z", "binary"),
    "raw 0/1 column"
  )
})

test_that("single-trait analyses do not emit a VIF warning", {
  df <- sim_data()
  msgs <- character(0)
  withCallingHandlers(
    analyze_disruptive_selection(df, "w", "z1", fitness_type = "continuous"),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_false(any(grepl("VIF", msgs)))
})

test_that("rows with a missing group label are kept and standardised", {
  df <- sim_data(n = 40)
  df$grp[c(3, 7, 11)] <- NA
  out <- suppressWarnings(suppressMessages(
    prepare_selection_data(df, "w", c("z1", "z2"), group = "grp", na_action = "warn")))
  expect_equal(nrow(out), 40)
  expect_false(anyNA(out$z1))
  expect_false(anyNA(out$relative_fitness))
  # the labelled groups are still standardised on their own
  expect_equal(sd(out$z1[out$grp %in% "A"]), 1, tolerance = 1e-8)
})

test_that("a group with zero mean fitness is named and gets NA relative fitness", {
  set.seed(3)
  df <- data.frame(
    surv = c(rbinom(30, 1, 0.6), rep(0, 20)),
    z = rnorm(50),
    yr = rep(c(2001, 2002), c(30, 20))
  )
  expect_warning(
    out <- suppressMessages(prepare_selection_data(df, "surv", "z", group = "yr")),
    "2002"
  )
  expect_true(all(is.na(out$relative_fitness[out$yr == 2002])))
  expect_false(anyNA(out$relative_fitness[out$yr == 2001]))
})

test_that("return_grouped skips rows whose group label is missing", {
  df <- sim_data(n = 120)
  df$grp[1:2] <- NA
  expect_warning(
    res <- suppressMessages(selection_coefficients(
      df, "w", c("z1", "z2"), fitness_type = "continuous", group = "grp", return_grouped = TRUE)),
    "missing"
  )
  expect_setequal(unique(res$Group), c("A", "B"))
})

test_that("a missing group label does not break the reference-group choice", {
  df <- sim_data(n = 120)
  df$year <- rep(c(2001, 2002), length.out = 120)
  df$year[120] <- NA
  prep <- suppressWarnings(suppressMessages(
    prepare_selection_data(df, "w", c("z1", "z2"), na_action = "none")))

  u <- suppressWarnings(suppressMessages(
    univariate_spline(prep, "w", "z1", fitness_type = "continuous", group = "year")))
  expect_s3_class(u, "univariate_fitness")
  expect_false(anyNA(u$grid$fit))

  s <- suppressWarnings(suppressMessages(
    correlated_fitness_surface(prep, "w", c("z1", "z2"), method = "gam", group = "year")))
  expect_false(anyNA(s$grid$.fit[s$grid$.inside]))
})

test_that("univariate_spline detects continuous fitness by default", {
  df <- sim_data(n = 80)
  df$z1 <- as.numeric(scale(df$z1))
  u <- suppressWarnings(suppressMessages(univariate_spline(df, "w", "z1")))
  expect_equal(u$fitness_type, "continuous")
  expect_match(u$ci_method, "parametric")
})

test_that("adaptive_landscape can use a surface fitted with a group term", {
  df <- sim_data(n = 120)
  prep <- suppressWarnings(suppressMessages(prepare_selection_data(df, "w", c("z1", "z2"), group = "grp")))
  surf <- suppressWarnings(suppressMessages(
    correlated_fitness_surface(prep, "w", c("z1", "z2"), method = "gam", group = "grp")))
  land <- suppressWarnings(suppressMessages(capture.output(
    al <- adaptive_landscape(prep, surf$model, c("z1", "z2"), grid_n = 6, simulation_n = 40))))
  expect_s3_class(al, "adaptive_landscape")
  expect_false(anyNA(al$grid$.mean_fit))
})

test_that("the 3D landscape is drawn from the grid without interpolation", {
  skip_if_not_installed("fields")
  g <- expand.grid(z1 = seq(-1, 1, length.out = 7), z2 = seq(-1, 1, length.out = 7))
  g$.mean_fit <- 1 - g$z1^2 - g$z2^2
  land <- structure(list(grid = g, trait_cols = c("z1", "z2")), class = "adaptive_landscape")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot_adaptive_landscape_3d(land, c("z1", "z2")))
})

test_that("the comparison overlay renders without optimum points", {
  g <- expand.grid(z1 = 1:5, z2 = 1:5)
  df <- rbind(
    cbind(g, fitness = runif(25), type = "Correlated Fitness (Individual)"),
    cbind(g, fitness = runif(25), type = "Adaptive Landscape (Population)")
  )
  plots <- plot_fitness_surfaces_comparison(list(combined_data = df, trait_cols = c("z1", "z2")))
  expect_no_error(ggplot2::ggplot_build(plots$overlay))
  plots2 <- plot_fitness_surfaces_comparison(
    list(combined_data = df, trait_cols = c("z1", "z2"),
         optimum_individual = data.frame(z1 = 2, z2 = 3),
         optimum_population = data.frame(z1 = 3, z2 = 2)),
    show_optima = FALSE)
  expect_no_error(ggplot2::ggplot_build(plots2$overlay))
})

test_that("selection_report differentials use the rows the gradients use", {
  df <- sim_data(n = 150)
  df$z2[1:30] <- NA # these individuals never reach the gradient models
  rep <- suppressWarnings(suppressMessages(
    selection_report(df, "w", c("z1", "z2"), fitness_type = "continuous")))

  cc <- df[stats::complete.cases(df[, c("w", "z1", "z2")]), ]
  z <- as.numeric(scale(cc$z1))
  w <- cc$w / mean(cc$w)
  S_ref <- mean((z - mean(z)) * (w - mean(w)))
  expect_equal(rep$Estimate[rep$Type == "Differential" & rep$Term == "z1"], S_ref, tolerance = 1e-10)
})

test_that("selection_report labels the fitness type actually used", {
  d <- data.frame(surv = rbinom(80, 1, 0.5), z1 = rnorm(80), z2 = rnorm(80))
  rep <- suppressWarnings(suppressMessages(
    selection_report(d, "surv", c("z1", "z2"), fitness_type = "continuous")))
  expect_equal(attr(rep, "fitness_type"), "continuous")
  expect_output(print(rep), "Fitness type: continuous")
})

test_that("bootstrap_selection validates n_boot and resamples within groups", {
  df <- sim_data(n = 60)
  expect_error(bootstrap_selection(df, "w", c("z1", "z2"), n_boot = NA), "n_boot")
  expect_error(bootstrap_selection(df, "w", c("z1", "z2"), n_boot = 25.5), "n_boot")

  set.seed(2)
  b <- suppressWarnings(suppressMessages(
    bootstrap_selection(df, "w", c("z1", "z2"), fitness_type = "continuous",
                        group = "grp", n_boot = 30)))
  expect_true(all(b$N_Boot <= 30))
  ref <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "w", c("z1", "z2"), fitness_type = "continuous", group = "grp")))
  expect_equal(b$Estimate, ref$Beta_Coefficient)
})
