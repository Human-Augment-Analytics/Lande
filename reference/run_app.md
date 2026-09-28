# Open the Shiny app

Starts the app that ships with the package: pick a bundled dataset or
upload a CSV, choose the fitness column, the traits and an optional
group, and read the gradients, fitness functions, surfaces and
landscapes from the tabs. Needs the `shiny` package; `plotly` adds the
rotatable 3D landscape and `fields` the static one.

## Usage

``` r
run_app(...)
```

## Arguments

- ...:

  Passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html), for
  example `port` or `launch.browser`.

## Value

Whatever [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html)
returns; called to open the app.

## Examples

``` r
if (interactive()) run_app()
```
