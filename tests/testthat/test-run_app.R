library(testthat)

test_that("the app ships with the package and run_app points at it", {
  dir <- system.file("app", package = "RforEvolution")
  expect_true(nzchar(dir))
  expect_true(file.exists(file.path(dir, "app.R")))
  expect_true(is.function(run_app))
  # the app file parses and refers to the package
  expr <- parse(file.path(dir, "app.R"))
  expect_gt(length(expr), 10)
  calls <- vapply(expr, function(e) paste(deparse(e), collapse = ""), "")
  expect_true(any(grepl("library(RforEvolution)", calls, fixed = TRUE)))
})
