# Bias and interval coverage of the gradients when the truth is known.
# Two traits, n = 200, six scenarios: (a) continuous fitness on correlated
# normal traits, (b) survival at about 11%, (c) overdispersed counts,
# (d) log-normal traits, (e) log-normal traits with a straight surface,
# (f) symmetric traits with heavy tails (Laplace), to tell the tails from the
# skew.
# Truth: the least-squares coefficients of expected relative fitness on the
# standardised traits, averaged over ten draws of two million individuals.
# Each term is reported on its own, with the Monte Carlo SE of its bias. Two
# diagnostics for every term: traits standardised by the population mean and
# SD instead of the sample's, and HC3 standard errors. Every replicate sets its
# own seed, so the results do not depend on the number of cores.
#
# Rscript validation/simulation.R [replicates] [resamples] [cores]
suppressPackageStartupMessages(library(Lande))

args <- as.numeric(commandArgs(trailingOnly = TRUE))
reps <- if (is.na(args[1])) 1000 else args[1]
n_boot <- if (is.na(args[2])) 200 else args[2]
cores <- if (is.na(args[3])) parallel::detectCores() - 2 else args[3]
n <- 200
terms <- c("z1 Linear", "z2 Linear", "z1² Quadratic", "z2² Quadratic", "z1 × z2 Correlational")

traits_normal <- function(m) {
  z1 <- rnorm(m)
  data.frame(z1 = z1, z2 = 0.6 * z1 + sqrt(1 - 0.36) * rnorm(m))
}
traits_lognormal <- function(m) {
  z <- traits_normal(m)
  data.frame(z1 = exp(0.5 * z$z1), z2 = exp(0.5 * z$z2))
}
# z1 is Laplace, no skew and kurtosis 6; z2 mixes in a second draw, so its
# kurtosis is about 4.6 (the log-normal traits above have skewness 1.75 and
# kurtosis 8.9); same correlation as the normal traits
traits_laplace <- function(m) {
  lap <- function(k) (rexp(k) - rexp(k)) / sqrt(2)
  z1 <- lap(m)
  data.frame(z1 = z1, z2 = 0.6 * z1 + sqrt(1 - 0.36) * lap(m))
}
# population standardisation, so the surface is defined on one scale whatever the sample
standardise <- function(d, centre, scale) data.frame(z1 = (d$z1 - centre[1]) / scale[1], z2 = (d$z2 - centre[2]) / scale[2])
surface <- function(s) 0.3 * s$z1 - 0.2 * s$z2 - 0.15 * s$z1^2 + 0.1 * s$z1 * s$z2
straight <- function(s) 0.3 * s$z1 - 0.2 * s$z2

scenarios <- list(
  a = list(label = "continuous fitness, normal traits", traits = traits_normal,
           expected = function(s) 2 + surface(s),
           draw = function(mu) mu + rnorm(length(mu), 0, 0.5)),
  b = list(label = "survival, about 11% survive", traits = traits_normal,
           expected = function(s) plogis(-2 + 2 * surface(s)),
           draw = function(mu) rbinom(length(mu), 1, mu)),
  c = list(label = "overdispersed counts", traits = traits_normal,
           expected = function(s) exp(0.5 + surface(s)),
           draw = function(mu) rnbinom(length(mu), size = 1.5, mu = mu)),
  d = list(label = "continuous fitness, log-normal traits", traits = traits_lognormal,
           expected = function(s) 2 + surface(s),
           draw = function(mu) mu + rnorm(length(mu), 0, 0.5)),
  e = list(label = "continuous fitness, log-normal traits, straight surface", traits = traits_lognormal,
           expected = function(s) 2 + straight(s),
           draw = function(mu) mu + rnorm(length(mu), 0, 0.5)),
  f = list(label = "continuous fitness, heavy-tailed symmetric traits", traits = traits_laplace,
           expected = function(s) 2 + surface(s),
           draw = function(mu) mu + rnorm(length(mu), 0, 0.5))
)

# the five gradients from least squares on relative fitness, with plain or HC3
# standard errors; squared-term coefficients and their errors doubled
gradients <- function(z, w, hc3 = FALSE) {
  z$w <- w / mean(w)
  se_of <- function(fit) {
    if (!hc3) return(summary(fit)$coefficients[, 2])
    X <- model.matrix(fit)
    e <- residuals(fit) / (1 - hatvalues(fit))
    B <- solve(crossprod(X))
    setNames(sqrt(diag(B %*% crossprod(X * e) %*% B)), colnames(X))
  }
  lin <- lm(w ~ z1 + z2, data = z)
  quad <- lm(w ~ z1 + z2 + I(z1^2) + I(z2^2) + z1:z2, data = z)
  sl <- se_of(lin)
  sq <- se_of(quad)
  list(est = setNames(c(coef(lin)[c("z1", "z2")], 2 * coef(quad)[c("I(z1^2)", "I(z2^2)")], coef(quad)["z1:z2"]), terms),
       se = setNames(c(sl[c("z1", "z2")], 2 * sq[c("I(z1^2)", "I(z2^2)")], sq["z1:z2"]), terms))
}
# the HC3 errors above must be the package's own; checked on fixed data, so no
# random numbers are drawn
local({
  i <- seq_len(n)
  z <- data.frame(z1 = sin(i), z2 = cos(0.7 * i))
  w <- 3 + z$z1 - 0.5 * z$z2^2 + 0.3 * z$z1 * z$z2 + sin(1.3 * i)
  pk <- suppressWarnings(suppressMessages(
    selection_coefficients(data.frame(z, W = w), "W", c("z1", "z2"), se_type = "hc3")))
  stopifnot(isTRUE(all.equal(unname(gradients(as.data.frame(scale(z)), w, hc3 = TRUE)$se),
                             pk$Standard_Error, tolerance = 1e-8)))
})

# the truth: regressions of expected relative fitness, averaged over ten draws
# of two million individuals; truth_se is the standard error of that average
truth <- function(sc, key) {
  draws <- lapply(1:10, function(j) {
    set.seed(10000 * match(key, names(scenarios)) + j)
    raw <- sc$traits(2e6)
    centre <- colMeans(raw)
    scale <- apply(raw, 2, sd)
    list(centre = centre, scale = scale, value = gradients(standardise(raw, centre, scale), sc$expected(standardise(raw, centre, scale)))$est)
  })
  values <- sapply(draws, `[[`, "value")
  list(centre = rowMeans(sapply(draws, `[[`, "centre")), scale = rowMeans(sapply(draws, `[[`, "scale")),
       value = rowMeans(values), se = apply(values, 1, sd) / sqrt(ncol(values)))
}

cover <- function(est, se, target) abs(est - target) <= 1.96 * se

one <- function(sc, tr) {
  raw <- sc$traits(n)
  d <- raw
  d$w <- sc$draw(sc$expected(standardise(raw, tr$centre, tr$scale)))
  fit <- tryCatch(suppressWarnings(suppressMessages(selection_coefficients(d, "w", c("z1", "z2")))), error = function(e) NULL)
  boot <- tryCatch(suppressWarnings(suppressMessages(bootstrap_selection(d, "w", c("z1", "z2"), n_boot = n_boot))), error = function(e) NULL)
  if (is.null(fit) || is.null(boot)) return(NULL)
  i <- match(terms, paste(fit$Term, fit$Type))
  j <- match(terms, paste(boot$Term, boot$Type))
  est <- fit$Beta_Coefficient[i]
  known <- gradients(standardise(raw, tr$centre, tr$scale), d$w)
  robust <- gradients(as.data.frame(scale(raw)), d$w, hc3 = TRUE)
  data.frame(term = terms, truth = tr$value, estimate = est,
             covered_parametric = cover(est, fit$Standard_Error[i], tr$value),
             covered_bootstrap = boot$CI_lower[j] <= tr$value & tr$value <= boot$CI_upper[j],
             covered_known_sd = cover(known$est, known$se, tr$value),
             covered_hc3 = cover(robust$est, robust$se, tr$value))
}

out <- list()
for (key in names(scenarios)) {
  sc <- scenarios[[key]]
  tr <- truth(sc, key)
  seed0 <- 1e6 * match(key, names(scenarios))
  runs <- parallel::mclapply(seq_len(reps), function(r) {
    set.seed(seed0 + r)
    one(sc, tr)
  }, mc.cores = cores)
  failed <- vapply(runs, function(x) is.null(x) || inherits(x, "try-error") || anyNA(x$estimate), logical(1))
  all_runs <- do.call(rbind, runs[!failed])
  for (tm in terms) {
    x <- all_runs[all_runs$term == tm, ]
    err <- x$estimate - x$truth
    out[[paste(key, tm)]] <- data.frame(
      scenario = key, description = sc$label, term = tm, truth = tr$value[[tm]], truth_se = tr$se[[tm]],
      bias = mean(err), bias_mcse = sd(err) / sqrt(length(err)), rmse = sqrt(mean(err^2)),
      coverage_parametric = mean(x$covered_parametric), coverage_bootstrap = mean(x$covered_bootstrap, na.rm = TRUE),
      coverage_known_sd = mean(x$covered_known_sd), coverage_hc3 = mean(x$covered_hc3),
      failed = mean(failed), replicates = sum(!failed))
  }
  message(key, " done: ", sum(!failed), " of ", reps, " replicates")
}
res <- do.call(rbind, out)
rownames(res) <- NULL
print(res[, -2], digits = 3)
write.csv(res, file.path("validation", "simulation_results.csv"), row.names = FALSE)
