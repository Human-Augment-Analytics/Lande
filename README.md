# RforEvolution

Stroud lab's Lande-Arnold toolkit for measuring phenotypic selection.

The package estimates selection differentials and linear (beta), quadratic
(gamma) and correlational (gamma_ij) selection gradients, fits spline fitness
functions and two-trait fitness surfaces, and computes the adaptive landscape
of mean fitness against the population mean. Each of these can be repeated by
year, and several groups can sit on one surface.

Traits are standardised to mean 0 and SD 1 and fitness is relativised (w / mean w)
before fitting. Quadratic gradients and their standard errors are doubled
(Stinchcombe et al. 2008); correlational gradients are not. For binary fitness the
gradients come from OLS on relative fitness and the p-values from a logistic GLM;
count fitness gets its p-values from a Poisson GLM, or a negative binomial one when
the counts are overdispersed.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Human-Augment-Analytics/R-for-Evolution@sean")
```

Drop the `@sean` if and when it is merged into main.

## Gradients

```r
library(RforEvolution)

traits <- c("weight", "total_length", "humerus")

# Standardise traits and add relative fitness
prepared <- prepare_selection_data(bumpus, "survival", traits)

# Selection gradients
selection_coefficients(bumpus, "survival", traits, fitness_type = "binary")

# Differentials and gradients in one table
selection_report(bumpus, "survival", traits, fitness_type = "binary")

# Bootstrap standard errors and intervals
bootstrap_selection(bumpus, "survival", traits, fitness_type = "binary")

# Assumption checks: normality of the traits, VIF, rows per term, residuals
check_selection_assumptions(bumpus, "survival", traits)

# One set of gradients per sex, each sex standardised on its own
selection_coefficients(bumpus, "survival", traits, group = "sex", return_grouped = TRUE)
```

## Fitness functions and surfaces

```r
# cubic spline, Wald band by default, bootstrap band with bootstrap = TRUE
uni <- univariate_spline(prepared, "survival", "weight")
plot_univariate_fitness(uni, "weight")

# two-trait surface, blank outside the convex hull of the data
surf <- correlated_fitness_surface(prepared, "survival", c("weight", "total_length"))
plot_correlated_fitness(surf, c("weight", "total_length"), show_points = TRUE)
```

`too_far` also blanks cells farther than a share of the axis range from any
individual, as in Beausoleil et al. (2023). With `group` set and
`group_effect = FALSE` one surface is fitted to everyone and each group's mean
and local peak are marked, which puts several species on one surface.

## Adaptive landscapes

```r
land <- adaptive_landscape(prepared, surf$model, c("weight", "total_length"))
plot_adaptive_landscape(land, c("weight", "total_length"))

# one trait, from the spline
land1 <- adaptive_landscape(prepared, uni$model, "weight")
plot_adaptive_landscape(land1, "weight")
```

## Over time

```r
finch <- read.csv(system.file("extdata", "finch_yearly.csv", package = "RforEvolution"))
finch <- prepare_selection_data(finch, "survived", "beak_pc1")
years <- temporal_landscape(finch, "survived", "beak_pc1", "year")
years$summary
plot_temporal_landscape(years)
plot_temporal_landscape(years, type = "heatmap")
```

## App

`run_app()` opens the whole workflow in a browser: bundled or uploaded data,
gradients, fitness functions, surface, landscape and per-group results, with
the fitting options in Advanced settings. Needs `shiny`; `plotly` adds the
rotatable 3D landscape. The app lives in `inst/app`; see its README there.

## Data

`bumpus`: Bumpus's 1898 house sparrows, 136 birds, nine traits and survival.
In `inst/extdata`: the Crescent Pond and Little Lake pupfish of Martin (2016),
the yearly medium ground finch data of Beausoleil et al. (2019), and the
five-group finch community of Beausoleil et al. (2023) with recapture years as
fitness.
Sources and the dataset DOI are listed under `?RforEvolution`.

## Notes

`detect_family()` classifies fitness as binary, count or continuous and picks
the model used for the p-values; the gradients themselves always come from OLS
on relative fitness. The spline uses a cubic regression basis with the
smoothing chosen by GCV and the surface a thin-plate basis chosen by REML;
both can be changed with `bs` and `smoothing`.

See `vignette("evolutionary-selection-analysis")` for the maths and worked
examples.
