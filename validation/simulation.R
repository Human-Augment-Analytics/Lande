# Bias and interval coverage of the gradients when the truth is known.
# Two traits, n = 200, four scenarios: (a) continuous fitness on correlated
# normal traits, (b) survival at about 15%, (c) overdispersed counts,
# (d) log-normal traits. Truth: the regressions of expected relative fitness
# on the standardised traits in two million individuals.
#
# Rscript validation/simulation.R [replicates] [resamples] [cores]
suppressPackageStartupMessages(library(Lande))

args <- as.numeric(commandArgs(trailingOnly = TRUE))
reps <- if (is.na(args[1])) 500 else args[1]
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
# population standardisation, so the surface is defined on one scale whatever the sample
standardise <- function(d, centre, scale) data.frame(z1 = (d$z1 - centre[1]) / scale[1], z2 = (d$z2 - centre[2]) / scale[2])
surface <- function(s) 0.3 * s$z1 - 0.2 * s$z2 - 0.15 * s$z1^2 + 0.1 * s$z1 * s$z2

scenarios <- list(
  a = list(label = "continuous fitness, normal traits", traits = traits_normal,
           expected = function(s) 2 + surface(s),
           draw = function(mu) mu + rnorm(length(mu), 0, 0.5)),
  b = list(label = "survival, about 15% survive", traits = traits_normal,
           expected = function(s) plogis(-2 + 2 * surface(s)),
           draw = function(mu) rbinom(length(mu), 1, mu)),
  c = list(label = "overdispersed counts", traits = traits_normal,
           expected = function(s) exp(0.5 + surface(s)),
           draw = function(mu) rnbinom(length(mu), size = 1.5, mu = mu)),
  d = list(label = "continuous fitness, log-normal traits", traits = traits_lognormal,
           expected = function(s) 2 + surface(s),
           draw = function(mu) mu + rnorm(length(mu), 0, 0.5))
)

# the truth: regressions of expected relative fitness in two million individuals
truth <- function(sc) {
  raw <- sc$traits(2e6)
  centre <- colMeans(raw)
  scale <- apply(raw, 2, sd)
  s <- standardise(raw, centre, scale)
  w <- sc$expected(s)
  w <- w / mean(w)
  lin <- coef(lm(w ~ z1 + z2, data = s))
  quad <- coef(lm(w ~ z1 + z2 + I(z1^2) + I(z2^2) + z1:z2, data = s))
  list(centre = centre, scale = scale,
       value = setNames(c(lin["z1"], lin["z2"], 2 * quad["I(z1^2)"], 2 * quad["I(z2^2)"], quad["z1:z2"]), terms))
}

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
  se <- fit$Standard_Error[i]
  data.frame(term = terms, truth = tr$value, estimate = est,
             covered_parametric = abs(est - tr$value) <= 1.96 * se,
             covered_bootstrap = boot$CI_lower[j] <= tr$value & tr$value <= boot$CI_upper[j])
}

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
out <- list()
for (key in names(scenarios)) {
  sc <- scenarios[[key]]
  tr <- truth(sc)
  runs <- parallel::mclapply(seq_len(reps), function(r) one(sc, tr), mc.cores = cores, mc.set.seed = TRUE)
  failed <- vapply(runs, function(x) is.null(x) || inherits(x, "try-error") || anyNA(x$estimate), logical(1))
  all_runs <- do.call(rbind, runs[!failed])
  kind <- ifelse(grepl("Linear", all_runs$term), "beta", "gamma")
  for (k in c("beta", "gamma")) {
    x <- all_runs[kind == k, ]
    out[[paste(key, k)]] <- data.frame(
      scenario = key, description = sc$label, gradient = k,
      mean_abs_truth = mean(abs(tr$value[grepl(if (k == "beta") "Linear" else "Quadratic|Correlational", terms)])),
      bias = mean(x$estimate - x$truth), rmse = sqrt(mean((x$estimate - x$truth)^2)),
      coverage_parametric = mean(x$covered_parametric), coverage_bootstrap = mean(x$covered_bootstrap, na.rm = TRUE),
      failed = mean(failed), replicates = sum(!failed))
  }
  message(key, " done: ", sum(!failed), " of ", reps, " replicates")
}
res <- do.call(rbind, out)
rownames(res) <- NULL
print(res[, -2], digits = 3)
write.csv(res, file.path("validation", "simulation_results.csv"), row.names = FALSE)
