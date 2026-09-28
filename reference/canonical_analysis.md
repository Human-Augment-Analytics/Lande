# Canonical analysis of the gamma matrix

Rotates the matrix of quadratic and correlational selection gradients to
its canonical axes (Phillips and Arnold 1989; Blows and Brooks 2003).
Each axis is a combination of the traits with a single curvature, its
eigenvalue: negative for stabilising and positive for disruptive
selection along that axis. Correlational selection that is spread over
many \\\gamma\_{ij}\\ shows up as curvature on a few axes.

## Usage

``` r
canonical_analysis(
  data,
  fitness_col,
  trait_cols,
  fitness_type = c("auto", "binary", "count", "continuous"),
  standardize = TRUE,
  group = NULL,
  bootstrap = FALSE,
  n_boot = 200,
  conf = 0.95
)
```

## Arguments

- data:

  A data frame containing fitness and trait measurements.

- fitness_col:

  A string specifying the name of the fitness column.

- trait_cols:

  A character vector of trait column names.

- fitness_type:

  A string indicating the fitness type: `"auto"`, `"binary"`, `"count"`,
  or `"continuous"`. Binary and count fitness take their p-values from a
  GLM (logistic, or Poisson and negative binomial) on the raw values;
  the gradients always come from OLS on relative fitness.

- standardize:

  Logical indicating whether to standardize traits to mean 0 and SD 1.
  Default is `TRUE`.

- group:

  Optional string specifying a grouping variable (e.g., "year", "site").

- bootstrap:

  Logical; add percentile intervals for the eigenvalues by resampling
  individuals, within groups when `group` is given. Default is `FALSE`.

- n_boot:

  Number of bootstrap resamples. Default is 200.

- conf:

  Confidence level of the bootstrap interval. Default is 0.95.

## Value

An object of class `"canonical_analysis"`: `gamma`, the gamma matrix;
`beta`; `M`, the loadings with one column per axis; `axes`, a data frame
with `lambda`, its standard error and p-value, `theta`, and with the
bootstrap `CI_lower`, `CI_upper` and `N_Boot`; `scores`, the prepared
data with the canonical scores `m1`, `m2`, ... added; the fitness type
and `n`.

## Details

The gamma matrix \\\gamma\\ comes from
[`selection_coefficients()`](https://human-augment-analytics.github.io/Lande/reference/selection_coefficients.md),
with the quadratic gradients doubled and the correlational gradients as
they stand. Its eigenvectors are the columns of `M` and its eigenvalues
the curvatures \\\lambda\\; directional selection along the axes is
\\\theta = M^T \beta\\. Standard errors and p-values come from the
double regression of Bisgaard and Ankenman (1996): fitness is refitted
on the canonical scores and their squares, without cross-products, and
twice the coefficient of a squared score is that axis's \\\lambda\\. For
survival and counts the p-values come from the logistic or count model
on the same terms, as they do for the gradients.

The tests treat the axes as known when they were estimated from the same
data, so they are anticonservative (Reynolds et al. 2010). Sampling
error in \\\gamma\\ also spreads its eigenvalues, so the largest
curvatures are overestimated, especially with many traits and few
individuals (Morrissey 2014). Read the eigenvalues beside the fitness
surface along the same axes, which
[`plot_canonical_axes()`](https://human-augment-analytics.github.io/Lande/reference/plot_canonical_axes.md)
draws, and rely more on the bootstrap intervals. With the bootstrap,
each resample's axes are matched to the original ones by their largest
absolute cosine, since order and sign change between resamples.

## References

Bisgaard, S. and Ankenman, B. (1996) Standard errors for the eigenvalues
in second-order response surface models. Technometrics 38, 238-246.
Blows, M. W. and Brooks, R. (2003) Measuring nonlinear selection. The
American Naturalist 162, 815-820. Morrissey, M. B. (2014) In search of
the best methods for multivariate selection analysis. Methods in Ecology
and Evolution 5, 1095-1109. Phillips, P. C. and Arnold, S. J. (1989)
Visualizing multivariate selection. Evolution 43, 1209-1222. Reynolds,
R. J., Childers, D. K. and Pajewski, N. M. (2010) The distribution and
hypothesis testing of eigenvalues from the canonical analysis of the
gamma matrix of quadratic and correlational selection gradients.
Evolution 64, 1076-1085.

## Examples

``` r
ca <- canonical_analysis(bumpus, "survival", c("total_length", "weight", "humerus"))
#> Warning: High multicollinearity detected (VIF > 5) - standard errors may be inflated
ca
#> Canonical analysis of gamma for total_length, weight, humerus on 136 individuals
#> 
#> Loadings (columns are the canonical axes):
#>                  m1     m2    m3
#> total_length -0.256 -0.056 0.965
#> weight        0.590  0.782 0.202
#> humerus       0.766 -0.621 0.167
#> 
#> Curvature along each axis (negative stabilising, positive disruptive):
#>  axis  lambda     se p_value  theta
#>    m1  0.0331 0.1138   0.678  0.281
#>    m2 -0.0248 0.2509   0.661 -0.507
#>    m3 -0.1511 0.0988   0.125 -0.282
#> 
#> The axes come from these data, so the tests are anticonservative.
```
