# run a call without its messages, warnings or printed output
quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))
  out
}
