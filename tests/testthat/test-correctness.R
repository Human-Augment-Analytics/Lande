# Known-answer checks for the Lande-Arnold estimators. These pin the numeric
# behaviour that structural smoke tests miss (relativisation, the factor-of-2
# convention, and the fitness-type routing).

make_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  list(z1 = z1, z2 = z2)
}

test_that("continuous fits return linear, quadratic and correlational terms", {
  d <- make_data()
  w <- 1 + 0.4 * d$z1 - 0.3 * d$z2 + 0.2 * d$z1^2 + 0.15 * d$z1 * d$z2 + rnorm(400, 0, 0.2)
  df <- data.frame(w = w, z1 = d$z1, z2 = d$z2)

  res <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "w", c("z1", "z2"), fitness_type = "continuous")
  ))

  expect_true(any(res$Type == "Quadratic"))
  expect_true(any(res$Type == "Correlational"))
})

test_that("quadratic gradients are twice the OLS quadratic coefficient", {
  d <- make_data()
  w <- 1 + 0.3 * d$z1 + 0.25 * d$z1^2 + rnorm(400, 0, 0.2)
  df <- data.frame(w = w, z1 = d$z1, z2 = d$z2)

  res <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "w", c("z1", "z2"), fitness_type = "continuous")
  ))

  # Reference: OLS on standardised traits and relative fitness, gamma = 2 * b_quad
  zz1 <- as.numeric(scale(d$z1))
  zz2 <- as.numeric(scale(d$z2))
  wrel <- w / mean(w)
  fit <- lm(wrel ~ zz1 + zz2 + I(zz1^2) + I(zz2^2) + zz1:zz2)
  ref_gamma <- 2 * coef(fit)["I(zz1^2)"]

  got <- res$Beta_Coefficient[res$Term == "z1²"]
  expect_equal(unname(got), unname(ref_gamma), tolerance = 1e-6)

  # The standard error is doubled alongside the estimate (Stinchcombe 2008), so
  # a refactor that dropped the SE doubling would be caught here.
  ref_gamma_se <- 2 * summary(fit)$coefficients["I(zz1^2)", "Std. Error"]
  got_se <- res$Standard_Error[res$Term == "z1²"]
  expect_equal(unname(got_se), unname(ref_gamma_se), tolerance = 1e-6)
})

test_that("binary gradients are estimated on relative fitness, not absolute 0/1", {
  d <- make_data()
  surv <- rbinom(400, 1, plogis(0.8 * d$z1 - 0.5 * d$z2))
  df <- data.frame(surv = surv, z1 = d$z1, z2 = d$z2)

  res <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "surv", c("z1", "z2"), fitness_type = "binary")
  ))

  zz1 <- as.numeric(scale(d$z1))
  zz2 <- as.numeric(scale(d$z2))
  rel <- coef(lm((surv / mean(surv)) ~ zz1 + zz2))[c("zz1", "zz2")]
  abs <- coef(lm(surv ~ zz1 + zz2))[c("zz1", "zz2")]
  got <- res$Beta_Coefficient[match(c("z1", "z2"), res$Term)]

  expect_equal(unname(got), unname(rel), tolerance = 1e-6)
  # And they must differ from the absolute-scale gradients by the survival rate
  expect_false(isTRUE(all.equal(unname(got), unname(abs), tolerance = 1e-3)))

  # P-values must come from the logistic GLM on the raw 0/1 outcome
  # (Janzen & Stern 1998), not from the OLS fit on relative fitness.
  glm_p <- coef(summary(glm(surv ~ zz1 + zz2, family = binomial)))[c("zz1", "zz2"), "Pr(>|z|)"]
  ols_p <- coef(summary(lm((surv / mean(surv)) ~ zz1 + zz2)))[c("zz1", "zz2"), "Pr(>|t|)"]
  got_p <- res$P_Value[match(c("z1", "z2"), res$Term)]

  expect_equal(unname(got_p), unname(glm_p), tolerance = 1e-6)
  # The GLM p-values are not the OLS ones. Compare on a relative scale, since
  # both can be tiny and an absolute tolerance would call them equal.
  rel_gap <- abs(got_p - ols_p) / pmax(abs(ols_p), .Machine$double.eps)
  expect_true(any(rel_gap > 0.05))
})

test_that("count fitness is analysed without crashing", {
  d <- make_data()
  fec <- rpois(400, lambda = exp(0.3 * d$z1))
  df <- data.frame(fec = fec, z1 = d$z1, z2 = d$z2)

  res <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "fec", c("z1", "z2"), fitness_type = "auto")
  ))

  expect_s3_class(res, "data.frame")
  expect_true(any(res$Type == "Linear"))
})

test_that("selection_differential equals the covariance of trait and fitness", {
  d <- make_data()
  z <- as.numeric(scale(d$z1))
  w <- 1 + 0.2 * z + rnorm(400, 0, 0.1)
  df <- data.frame(z = z, w = w)

  S <- selection_differential(df, "w", "z", standardized = TRUE, use_relative = FALSE)
  expect_equal(S, mean((z - mean(z)) * (w - mean(w))), tolerance = 1e-9)
})

test_that("prepare_selection_data standardizes traits and centres relative fitness at 1", {
  d <- make_data(n = 100)
  df <- data.frame(w = runif(100, 1, 10), z1 = d$z1[1:100], z2 = d$z2[1:100])

  out <- suppressWarnings(suppressMessages(
    prepare_selection_data(df, "w", c("z1", "z2"), name_relative = "w_relative")
  ))

  expect_equal(mean(out$z1), 0, tolerance = 1e-8)
  expect_equal(sd(out$z1), 1, tolerance = 1e-8)
  expect_equal(mean(out$w_relative), 1, tolerance = 1e-8)
})

test_that("traits are standardized on the analyzed sample when fitness is missing", {
  d <- make_data(n = 300)
  w <- 1 + 0.3 * d$z1 - 0.2 * d$z2 + rnorm(300, 0, 0.2)
  df <- data.frame(w = w, z1 = d$z1, z2 = d$z2)
  df$w[1:60] <- NA # 60 individuals measured for traits but not fitness

  res <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "w", c("z1", "z2"), fitness_type = "continuous")))
  got <- res$Beta_Coefficient[match(c("z1", "z2"), res$Term)]

  # Reference: standardise and relativise on the complete-fitness rows only
  cc <- df[!is.na(df$w), ]
  zz1 <- as.numeric(scale(cc$z1)); zz2 <- as.numeric(scale(cc$z2))
  ref <- coef(lm((cc$w / mean(cc$w)) ~ zz1 + zz2))[c("zz1", "zz2")]
  expect_equal(unname(got), unname(ref), tolerance = 1e-8)
})

test_that("relative fitness is taken over the analyzed sample when a trait is missing", {
  d <- make_data(n = 300)
  w <- 1 + 0.3 * d$z1 - 0.2 * d$z2 + rnorm(300, 0, 0.2)
  df <- data.frame(w = w, z1 = d$z1, z2 = d$z2)
  # Drop a trait for the 60 fittest individuals, so mean(W) over everyone with
  # fitness differs from mean(W) over the rows that reach the models.
  df$z2[order(-df$w)[1:60]] <- NA

  prep <- suppressWarnings(suppressMessages(
    prepare_selection_data(df, "w", c("z1", "z2"), name_relative = "w_rel")))
  cc <- complete.cases(prep[, c("w", "z1", "z2")])
  expect_equal(mean(prep$w_rel[cc]), 1, tolerance = 1e-10)

  res <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "w", c("z1", "z2"), fitness_type = "continuous")))
  got <- res$Beta_Coefficient[match(c("z1", "z2"), res$Term)]
  d2 <- df[complete.cases(df), ]
  zz1 <- as.numeric(scale(d2$z1)); zz2 <- as.numeric(scale(d2$z2))
  ref <- coef(lm((d2$w / mean(d2$w)) ~ zz1 + zz2))[c("zz1", "zz2")]
  expect_equal(unname(got), unname(ref), tolerance = 1e-8)
})

test_that("a single-observation group does not switch off standardization elsewhere", {
  d <- make_data(n = 61)
  df <- data.frame(w = runif(61, 1, 5), z = d$z1, grp = c(rep("A", 30), rep("B", 30), "C"))

  expect_warning(
    out <- prepare_selection_data(df, "w", "z", group = "grp", name_relative = "w_rel"),
    "no variance in group"
  )
  expect_equal(sd(out$z[out$grp == "A"]), 1, tolerance = 1e-8)
  expect_equal(sd(out$z[out$grp == "B"]), 1, tolerance = 1e-8)
  expect_equal(out$z[out$grp == "C"], 0) # centred only, kept in the data
})

test_that("a zero-variance trait warns and does not drop every row", {
  d <- make_data(n = 100)
  df <- data.frame(w = runif(100, 1, 10), z1 = d$z1[1:100], flat = 5)

  expect_warning(
    out <- prepare_selection_data(df, "w", c("z1", "flat"), name_relative = "w_rel"),
    "zero-variance", ignore.case = TRUE
  )
  expect_equal(nrow(out), 100) # rows preserved, not silently dropped
  expect_equal(sd(out$z1), 1, tolerance = 1e-8) # the good trait is still standardized
  expect_true(all(out$flat == 5)) # the constant trait is left as-is
})

test_that("a trait constant within a group is caught by the grouped zero-variance guard", {
  set.seed(3)
  # `flat` varies globally (1 vs 2) but is constant within each group, so a
  # global-only check would miss it while per-group scale() would make it NaN.
  df <- data.frame(
    w = runif(60, 1, 5),
    z = rnorm(60),
    flat = rep(c(1, 2), each = 30),
    grp = rep(c("A", "B"), each = 30)
  )

  expect_warning(
    out <- prepare_selection_data(df, "w", c("z", "flat"), group = "grp", name_relative = "w_rel"),
    "zero-variance", ignore.case = TRUE
  )
  expect_equal(nrow(out), 60) # no group silently dropped
  expect_false(any(is.nan(out$flat))) # flat left unstandardized, not NaN
  expect_true(all(is.finite(out$z))) # the good trait is still usable
})

test_that("detect_family classifies the common fitness shapes", {
  expect_equal(detect_family(c(0, 1, 0, 1, 1, 0, 1, 0, 1, 0))$type, "binary")
  expect_equal(suppressWarnings(detect_family(c(0, 2, 3, 5, 1, 4, 2, 6, 3, 1))$type), "count")
  expect_equal(detect_family(rnorm(50, 5, 1))$type, "continuous")
})
