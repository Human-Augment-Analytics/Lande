# RforEvolution GUI

Shiny front-end for the package. Pick a dataset or upload a CSV (one row per
individual), choose the fitness column, the traits and optionally a group
column, then Run.

Tabs:

- **Data**: n, survival rate or mean fitness, group sizes, summary table,
  histograms, and the settings used (including the file name).
- **Selection gradients**: S, β ± SE and γ ± SE with p-values, correlational
  terms, a one-line summary per trait, a gradient plot, optional bootstrap
  intervals, CSV and PNG downloads.
- **Fitness functions**: cubic-spline fitness function per trait with a
  bootstrap band, and optionally the adaptive landscape for that trait.
- **Fitness surface**: two-trait surface with the individuals overlaid, drawn
  only where there are data. With a group, each group's mean and the highest
  point of the surface within its own range are marked. Untick "Standardise
  within group" to treat the groups as one population, for several species
  on one surface.
- **Adaptive landscape**: 2D contour and rotatable 3D view, with the optimum
  and the current population mean.
- **Groups**: gradients per group, table and plot.

Bootstrap, spline band and landscape simulation all run from the seed in
Advanced settings, so repeated runs match. Advanced settings also hold the
fitness type, the spline and surface basis (cubic regression, thin plate or
P-spline) and smoothing criterion (GCV, REML or ML), the surface basis size
k, a switch to draw the surface over the whole grid rather than only where
there are data, and a distance rule that blanks cells far from any
individual.

## Run locally

```r
install.packages(".", repos = NULL, type = "source")  # or devtools::load_all()
shiny::runApp("app")
```

Bundled datasets: Bumpus sparrows, Crescent Pond pupfish, Little Lake pupfish
(shipped in `inst/extdata`). Uploaded CSVs need numeric fitness and trait
columns.

## Optional packages

- `fields`: thin-plate spline surfaces, static 3D landscape
- `plotly`: rotatable 3D landscape
- `car`: Type III ANOVA tables

## GitHub Pages via Shinylive

```r
shinylive::export("app", "docs")
# commit docs/ and set Pages to deploy from /docs
```

Everything runs in the visitor's browser (first load fetches WebR and the
packages, roughly 20 to 40 MB). Keep the landscape grid and simulation sliders
modest. RforEvolution has no compiled code, so under WebR you load the CRAN
dependencies and source `R/*.R` directly.
