# Compare the fitted fitness at two points of a surface

Takes two points on a GAM surface, given as group names or trait values,
and returns the difference in fitted fitness on the scale of the link,
with its standard error. Two peaks are separate only if each is higher
than the dip between them, which `valley = TRUE` checks.

## Usage

``` r
peak_difference(surface, from, to, valley = FALSE, n_path = 50)
```

## Arguments

- surface:

  Output of
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/Lande/reference/correlated_fitness_surface.md)
  fitted with `method = "gam"`.

- from, to:

  A group name, standing for that group's highest point in
  `surface$groups`, or a numeric vector of the two trait values.

- valley:

  Logical; also find the lowest fitted fitness on the straight line
  between the two points and compare each end with it. Default is
  `FALSE`.

- n_path:

  Number of points along that line. Default is 50.

## Value

A data frame with one row per comparison: the fitted fitness at the two
points (`fit_a`, `fit_b`), their difference on the link scale, its
standard error, `z` and the two-sided `p_value`. With `valley = TRUE`
the attribute `"valley"` holds the trait values of the lowest point on
the line.

## Details

Differences are on the scale of the link, log fitness for counts and log
odds for survival. The standard error comes from the covariance of the
model's coefficients, corrected for smoothing parameter uncertainty when
the fit was by REML or ML. The valley is picked as the lowest fitted
point on the line, so its comparisons describe the fitted surface and
are not planned tests.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#> IMPORTANT: Traits should already be standardized (mean = 0, SD = 1).
#>            Do NOT apply scale() again within this function.
#> Data type: binary; method: gam; n = 136; k = 29
#> GAM fitting with 136 observations
#>   Trying formula: main
#> Success with formula: main
#> Predictions range: 0.0472 to 0.703
#> Masked 424 of 900 grid points outside the data
peak_difference(surf, c(-1, -1), c(1, 1))
#>          comparison     fit_a     fit_b difference        se        z
#> 1 (-1, -1) - (1, 1) 0.6883016 0.3811798   1.276735 0.6329478 2.017125
#>      p_value
#> 1 0.04368246
```
