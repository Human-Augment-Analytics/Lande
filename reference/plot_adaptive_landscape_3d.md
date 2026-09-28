# Plot Adaptive Landscape (3D Perspective)

Plot Adaptive Landscape (3D Perspective)

## Usage

``` r
plot_adaptive_landscape_3d(
  landscape,
  trait_cols,
  theta = -30,
  phi = 30,
  grid_n = 200,
  color_palette = NULL,
  ...
)
```

## Arguments

- landscape:

  Output object of class `"adaptive_landscape"`.

- trait_cols:

  Character vector of length 2 specifying the trait column names.

- theta:

  Numeric azimuthal viewing angle. Default is -30.

- phi:

  Numeric colatitude viewing angle. Default is 30.

- grid_n:

  Integer specifying the resolution of the interpolation grid. Default
  is 200.

- color_palette:

  Optional vector of colors for the surface. Defaults to viridis plasma.

- ...:

  Additional arguments passed to
  [`fields::drape.plot()`](https://rdrr.io/pkg/fields/man/drape.plot.html).

## Value

A 3D plot produced by
[`fields::drape.plot()`](https://rdrr.io/pkg/fields/man/drape.plot.html).

## Examples

``` r
if (requireNamespace("fields", quietly = TRUE)) {
  prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
  surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
  land <- adaptive_landscape(prep, surf$model, c("total_length", "weight"),
                             simulation_n = 100, grid_n = 15)
  plot_adaptive_landscape_3d(land, c("total_length", "weight"))
}
#> IMPORTANT: Traits should already be standardized (mean = 0, SD = 1).
#>            Do NOT apply scale() again within this function.
#> Data type: binary; method: gam; n = 136; k = 29
#> GAM fitting with 136 observations
#>   Trying formula: main
#> Success with formula: main
#> Predictions range: 0.0472 to 0.703
#> Masked 424 of 900 grid points outside the data
#> IMPORTANT: Traits should already be standardized (mean = 0, SD = 1).
#>            Use prepare_selection_data() before calling this function.
#>            Do NOT standardize again within this function.
#> Population mean grid ranges:
#>   total_length: -2.96 to 2.94
#>   weight: -3.12 to 4.85
#> Estimated within-population variance-covariance:
#>              total_length weight
#> total_length       1.0000 0.5839
#> weight             0.5839 1.0000
#> Calculating mean fitness for 225 grid points
#> Optimal population mean phenotype:
#>   total_length    weight
#> 8   -0.0123897 -3.121579
#> Mean fitness at optimum: 0.6808
#> 75% of simulated individuals fell outside the data (96% at the optimum)
#> The optimum rests on extrapolation: more than 25% of the population simulated there lies outside the data
```
