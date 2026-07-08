# RforEvolution

Stroud lab's Lande-Arnold toolkit for measuring phenotypic selection.

The package estimates linear (beta), quadratic (gamma), and correlational
(gamma_ij) selection gradients, computes selection differentials, and draws
univariate fitness functions, correlated fitness surfaces, and adaptive
landscapes.

Traits are standardised to mean 0 and SD 1 and fitness is relativised (w / mean w)
before fitting. Quadratic gradients and their standard errors are doubled
(Stinchcombe et al. 2008); correlational gradients are not. For binary fitness the
gradients come from OLS on relative fitness and the p-values from a logistic GLM;
count fitness gets its p-values from a Poisson GLM, or a negative binomial one when
the counts are overdispersed.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Human-Augment-Analytics/R-for-Evolution")
```

## Example

```r
library(RforEvolution)

data(bumpus)                              # Bumpus (1898) sparrow survival data
traits <- c("weight", "total_length", "humerus")

# Standardise traits and add relative fitness
prepared <- prepare_selection_data(bumpus, "survival", traits)

# Selection gradients
coefs <- selection_coefficients(bumpus, "survival", traits, fitness_type = "binary")

# Summary table with differentials and gradients together
selection_report(bumpus, "survival", traits, fitness_type = "binary")

# Bootstrap standard errors and confidence intervals
bootstrap_selection(bumpus, "survival", traits, fitness_type = "binary")
```

Fitness surfaces:

```r
uni <- univariate_spline(prepared, "survival", "weight", fitness_type = "binary")
plot_univariate_fitness(uni, "weight")

surf <- correlated_fitness_surface(prepared, "survival", c("weight", "total_length"))
plot_correlated_fitness(surf, c("weight", "total_length"))
```

## Notes

`detect_family()` classifies fitness as binary, count, proportion, or continuous
and picks the model used for the p-values; the gradients themselves always come
from OLS on relative fitness. The univariate surface is a cubic spline with the
smoothing chosen by GCV, and its confidence ribbon is bootstrapped.

See `vignette("evolutionary-selection-analysis")` for details.
