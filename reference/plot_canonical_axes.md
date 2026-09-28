# Fitness along the canonical axes

Draws the fitness surface in the space of two canonical axes, or the
fitness function along one, from the scores of
[`canonical_analysis()`](https://human-augment-analytics.github.io/Lande/reference/canonical_analysis.md).
A negative \\\lambda\\ is stabilising selection along an axis only when
fitness along that axis has a peak inside the data, otherwise it is only
curvature.

## Usage

``` r
plot_canonical_axes(ca, which = 1:2, grid_n = 40, ...)
```

## Arguments

- ca:

  Output of
  [`canonical_analysis()`](https://human-augment-analytics.github.io/Lande/reference/canonical_analysis.md).

- which:

  One or two axis numbers. Default is `1:2`.

- grid_n:

  Grid resolution of the surface. Default is 40.

- ...:

  Passed to
  [`plot_correlated_fitness()`](https://human-augment-analytics.github.io/Lande/reference/plot_correlated_fitness.md)
  for two axes or to
  [`plot_univariate_fitness()`](https://human-augment-analytics.github.io/Lande/reference/plot_univariate_fitness.md)
  for one.

## Value

A `ggplot` object.

## Examples

``` r
ca <- canonical_analysis(bumpus, "survival", c("total_length", "weight", "humerus"))
#> Warning: High multicollinearity detected (VIF > 5) - standard errors may be inflated
plot_canonical_axes(ca, which = c(1, 3), grid_n = 25)
```
