#' Open the Shiny app
#'
#' Starts the app that ships with the package: pick a bundled dataset or upload
#' a CSV, choose the fitness column, the traits and an optional group, and read
#' the gradients, fitness functions, surfaces and landscapes from the tabs.
#' Needs the \code{shiny} package; \code{plotly} adds the rotatable 3D
#' landscape and \code{fields} the static one.
#'
#' @param ... Passed to \code{shiny::runApp()}, for example \code{port} or
#'   \code{launch.browser}.
#' @return Whatever \code{shiny::runApp()} returns; called to open the app.
#' @examples
#' if (interactive()) run_app()
#' @export
run_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("run_app() needs the shiny package: install.packages(\"shiny\")", call. = FALSE)
  }
  dir <- system.file("app", package = "RforEvolution")
  if (!nzchar(dir)) stop("The app folder is missing from this installation", call. = FALSE)
  shiny::runApp(dir, ...)
}
