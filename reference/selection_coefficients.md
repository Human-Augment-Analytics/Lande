# Calculate selection coefficients

Main wrapper function for selection analysis following Lande & Arnold
(1983). Extracts linear (beta), quadratic (gamma), and correlational
selection gradients.

## Usage

``` r
selection_coefficients(
  data,
  fitness_col,
  trait_cols,
  fitness_type = c("auto", "binary", "count", "continuous"),
  standardize = TRUE,
  group = NULL,
  use_relative_for_fit = TRUE,
  return_grouped = FALSE
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

- use_relative_for_fit:

  Logical; if `TRUE` (default) the gradients are estimated on relative
  fitness \\W / \bar{W}\\ for every fitness type, as in Lande & Arnold
  (1983). Set `FALSE` only to reproduce coefficients on the absolute
  fitness scale.

- return_grouped:

  Logical indicating whether to return results grouped if a `group` is
  specified.

## Value

A data frame containing selection coefficients (Term, Type,
Beta_Coefficient, Standard_Error, P_Value, Variance).

## Examples

``` r
selection_coefficients(bumpus, "survival", c("total_length", "weight"), fitness_type = "binary")
#> there are higher-order terms (interactions) in this model
#> consider setting type = 'predictor'; see ?vif
#> Warning: High multicollinearity detected (VIF > 5) - standard errors may be inflated
#>                    Term          Type Beta_Coefficient Standard_Error
#> 1          total_length        Linear      -0.17811220      0.0974800
#> 2                weight        Linear      -0.09992466      0.0974800
#> 3         total_length²     Quadratic      -0.23342057      0.2111590
#> 4               weight²     Quadratic       0.02096490      0.1628337
#> 5 total_length × weight Correlational      -0.06949785      0.1585041
#>      P_Value    Variance
#> 1 0.07276336 0.009502351
#> 2 0.30287778 0.009502351
#> 3 0.24402273 0.044588139
#> 4 0.93873420 0.026514809
#> 5 0.57236384 0.025123549

# one set of gradients per sex, each sex standardised on its own
selection_coefficients(bumpus, "survival", c("total_length", "weight"),
                       group = "sex", return_grouped = TRUE)
#> there are higher-order terms (interactions) in this model
#> consider setting type = 'predictor'; see ?vif
#> there are higher-order terms (interactions) in this model
#> consider setting type = 'predictor'; see ?vif
#>                           Term          Type Beta_Coefficient Standard_Error
#> male.1            total_length        Linear    -0.3407483627     0.09731904
#> male.2                  weight        Linear    -0.0545617939     0.09731904
#> male.3           total_length²     Quadratic    -0.0273816426     0.19241802
#> male.4                 weight²     Quadratic     0.0462210677     0.17271044
#> male.5   total_length × weight Correlational    -0.0762079316     0.15832913
#> female.1          total_length        Linear     0.0004585805     0.20517997
#> female.2                weight        Linear    -0.2876716237     0.20517997
#> female.3         total_length²     Quadratic    -0.2985821975     0.48947090
#> female.4               weight²     Quadratic     0.0333986785     0.32247945
#> female.5 total_length × weight Correlational    -0.0520198090     0.33611268
#>              P_Value    Variance  Group
#> male.1   0.001978417 0.009470995   male
#> male.2   0.606430455 0.009470995   male
#> male.3   0.832288575 0.037024693   male
#> male.4   0.751956548 0.029828895   male
#> male.5   0.544449522 0.025068112   male
#> female.1 0.962802935 0.042098822 female
#> female.2 0.159163810 0.042098822 female
#> female.3 0.478153622 0.239581761 female
#> female.4 0.599203577 0.103992995 female
#> female.5 0.968497077 0.112971735 female
```
