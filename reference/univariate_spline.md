# Estimate univariate correlated fitness function

This function calculates a univariate correlated fitness function using
a GAM smooth term. The formula is `w ~ f(z)`. Traits should be
standardized BEFORE calling this function.

## Usage

``` r
univariate_spline(
  data,
  fitness_col,
  trait_col,
  fitness_type = c("auto", "binary", "count", "continuous"),
  group = NULL,
  relative_col = NULL,
  k = 10,
  bs = c("cr", "tp", "ps"),
  smoothing = c("GCV.Cp", "REML", "ML"),
  bootstrap = FALSE,
  n_boot = 1000,
  by_group = FALSE
)
```

## Arguments

- data:

  A data frame containing fitness and trait measurements.

- fitness_col:

  A string specifying the name of the fitness column.

- trait_col:

  A string specifying the name of the trait column (must be numeric and
  standardized).

- fitness_type:

  A string indicating the fitness type: `"auto"` (detect from the data,
  the default), `"binary"`, `"count"`, or `"continuous"`. Binary and
  count fitness are fitted on the raw values with a binomial or Poisson
  family; continuous fitness on relative fitness with a Gaussian one.

- group:

  Optional string specifying a grouping variable. If provided, group
  fixed effects are included.

- relative_col:

  Optional string naming a pre-computed relative fitness column to use
  for continuous fitness (e.g. one produced within groups by
  `prepare_selection_data`).

- k:

  Integer specifying the basis dimension for the smooth term. Default is
  10.

- bs:

  Spline basis: `"cr"` (cubic regression spline, the default), `"tp"`
  (thin plate) or `"ps"` (P-spline).

- smoothing:

  How the smoothing parameter is chosen: `"GCV.Cp"` (generalised
  cross-validation, the default), `"REML"` or `"ML"`.

- bootstrap:

  Logical; if `TRUE` the 95% ribbon is obtained by resampling
  individuals and refitting (Schluter 1988). The default `FALSE` uses
  the parametric Wald interval, which is instant; the bootstrap refits
  the spline `n_boot` times.

- n_boot:

  Integer number of bootstrap resamples used when `bootstrap = TRUE`.
  Default is 1000.

- by_group:

  Logical; if `TRUE` fit each level of `group` on its own rows and
  return a named list of fits, one per group. The default `FALSE` fits
  one curve with the group as a fixed effect. Prepare the data with the
  same `group` first, so that each group is standardised on its own.

## Value

A list of class `"univariate_fitness"` containing the fitted GAM model,
a prediction grid, and metadata.

## Details

By default the fitness function is a penalised cubic regression spline
with the smoothing parameter chosen by generalised cross-validation,
following Schluter (1988). `bs` and `smoothing` are there to match the
smoother of another study; they do not change how much the curve is
smoothed, which is always chosen from the data. If mgcv's check suggests
the basis dimension was too small a warning says so. With
`bootstrap = TRUE` the result depends on the random seed; call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for a
reproducible ribbon.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", "total_length")
uni <- univariate_spline(prep, "survival", "total_length")
#> Fitness type detected: binary
#> IMPORTANT: Traits should already be standardized (mean = 0, SD = 1).
#>            Do NOT apply scale() again within this function.
#> Trait appears standardized (mean ~ 0, SD ~ 1)
#> Warning: k = 10 may be too small for 'total_length' (mgcv k-index 0.84); try a larger k
plot_univariate_fitness(uni, "total_length")


# a thin-plate basis with REML, to match another study's smoother
univariate_spline(prep, "survival", "total_length", bs = "tp", smoothing = "REML")$spline_type
#> Fitness type detected: binary
#> IMPORTANT: Traits should already be standardized (mean = 0, SD = 1).
#>            Do NOT apply scale() again within this function.
#> Trait appears standardized (mean ~ 0, SD ~ 1)
#> Warning: k = 10 may be too small for 'total_length' (mgcv k-index 0.77); try a larger k
#> [1] "thin-plate spline (REML)"
```
