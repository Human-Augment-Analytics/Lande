# RforEvolution 0.1.0

First release.

* Selection differentials and linear, quadratic and correlational gradients on
  relative fitness, with p-values from the logistic model for survival and from
  a Poisson or negative binomial model for counts.
* Cubic spline fitness functions with a bootstrap band, two-trait fitness
  surfaces drawn only where there are data, and adaptive landscapes of mean
  fitness for one or two traits.
* Several groups on one surface, with each group's mean and highest point
  marked, and interior peaks told from maxima at the edge of the data.
* Fitness functions and landscapes fitted period by period.
* Bootstrap standard errors and percentile intervals for the gradients.
* Assumption checks in one table: multivariate normality of the traits, VIF,
  rows per term, residuals or dispersion of the gradient model.
* Bumpus's sparrows bundled as data; pupfish and finch data in `extdata`.
* A Shiny app for the whole workflow, opened with `run_app()`.
* Standard errors and a band on the GAM fitness surface, drawn with the plots'
  `uncertainty` option, and `peak_difference()` to compare two points of a
  surface, or each with the dip between them.
* The adaptive landscape reports how much of each simulated population falls
  outside the data, and the plot can mark where it does.
* A `clamp` option on the landscape and the thin-plate surface holds predicted
  fitness inside the range of the fitness type; on by default.
