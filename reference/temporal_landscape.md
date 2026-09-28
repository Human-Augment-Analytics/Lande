# Fitness functions and adaptive landscapes over time

Fits the fitness function (one trait) or fitness surface (two traits)
separately for each level of a time column, usually year, and optionally
the adaptive landscape for each, so that how selection changes over time
can be seen. Beausoleil et al. (2019) did this for Darwin's finches, one
year at a time.

## Usage

``` r
temporal_landscape(
  data,
  fitness_col,
  trait_cols,
  time_col,
  fitness_type = c("auto", "binary", "count", "continuous"),
  min_n = 20,
  landscape = TRUE,
  k = NULL,
  bs = NULL,
  smoothing = NULL,
  bootstrap = FALSE,
  n_boot = 200,
  grid_n = 60,
  simulation_n = 300,
  mask = TRUE,
  too_far = NULL
)
```

## Arguments

- data:

  A data frame with fitness, trait and time columns. Standardise the
  traits once, over all periods, with
  [`prepare_selection_data()`](https://human-augment-analytics.github.io/Lande/reference/prepare_selection_data.md)
  before calling; the periods must share one trait axis, so nothing is
  restandardised per period.

- fitness_col:

  Name of the fitness column.

- trait_cols:

  One or two trait column names.

- time_col:

  Name of the column giving the period of each row.

- fitness_type:

  Passed to
  [`univariate_spline()`](https://human-augment-analytics.github.io/Lande/reference/univariate_spline.md)
  for one trait; the surface detects the type itself.

- min_n:

  Periods with fewer complete rows than this are skipped and listed in
  the result's `skipped`. Default is 20.

- landscape:

  Logical; also compute the adaptive landscape for each period. Default
  is `TRUE`.

- k, bs, smoothing:

  Basis size, basis and smoothing criterion for the period fits. `NULL`
  uses the defaults of
  [`univariate_spline()`](https://human-augment-analytics.github.io/Lande/reference/univariate_spline.md)
  or
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/Lande/reference/correlated_fitness_surface.md).

- bootstrap, n_boot:

  For one trait, a bootstrap band around each period's fitness function.

- grid_n:

  Grid resolution for the surfaces and landscapes.

- simulation_n:

  Simulated individuals per grid point in the landscapes.

- mask, too_far:

  Passed to
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/Lande/reference/correlated_fitness_surface.md)
  for two traits.

## Value

An object of class `"temporal_landscape"`: `fits` and `landscapes`, one
per period; `summary`, one row per period with n, mean fitness, trait
means, the smooth's degrees of freedom, the position and height of the
highest fitted fitness and whether it sits at the edge of the data
(`optimum_edge`), the number of interior peaks (one trait) and the
landscape optimum; `grid`, the fitted values of every period stacked
with a `time` column; for one trait `heat`, each period's fitness
function on one common trait grid, `NA` outside that period's data;
`landscape_grid`, the landscapes stacked; `data`, the rows used; and
`skipped`.

## Examples

``` r
prep <- prepare_selection_data(finch_yearly, "survived", "beak_pc1")
years <- temporal_landscape(prep, "survived", "beak_pc1", "year", landscape = FALSE)
#> 2004: n = 110, mean fitness 0.273, edf 2.5, 1 interior peak, highest fitness at the edge of the data
#> 2005: n = 185, mean fitness 0.205, edf 3.9, 2 interior peaks
#> 2006: n = 233, mean fitness 0.086, edf 1.0, 0 interior peaks, highest fitness at the edge of the data
#> 2007: n = 61, mean fitness 0.262, edf 1.7, 1 interior peak
#> 2008: n = 127, mean fitness 0.205, edf 6.4, 3 interior peaks
#> 2009: n = 196, mean fitness 0.153, edf 4.0, 1 interior peak, highest fitness at the edge of the data
#> 2010: n = 175, mean fitness 0.120, edf 1.0, 0 interior peaks, highest fitness at the edge of the data
years$summary
#>   time   n mean_fitness mean_beak_pc1      edf optimum_beak_pc1 optimum_fit
#> 1 2004 110   0.27272727    0.15262095 2.514761       2.56126184   0.4175333
#> 2 2005 185   0.20540541    0.01659862 3.891967       1.51821556   0.2515484
#> 3 2006 233   0.08583691   -0.07931978 1.000098       2.70000995   0.1085371
#> 4 2007  61   0.26229508    0.12403720 1.672904       0.70740350   0.3342915
#> 5 2008 127   0.20472441   -0.12870140 6.385272      -0.02345104   0.3931988
#> 6 2009 196   0.15306122   -0.08518852 4.007243       2.84197912   0.6114977
#> 7 2010 175   0.12000000    0.13770410 1.000048      -2.10396884   0.2646446
#>   optimum_edge peaks
#> 1         TRUE     1
#> 2        FALSE     2
#> 3         TRUE     0
#> 4        FALSE     1
#> 5        FALSE     3
#> 6         TRUE     1
#> 7         TRUE     0
plot_temporal_landscape(years, type = "heatmap")
```
