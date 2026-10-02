test_that("HC3 standard errors are the leave-one-out sum and change nothing else", {
  set.seed(4)
  n <- 150
  d <- data.frame(z1 = rnorm(n), z2 = rexp(n))
  d$w <- 2 + 0.3 * d$z1 - 0.2 * d$z1^2 + rnorm(n, 0, 0.3 + 0.3 * abs(d$z1))
  ols <- suppressMessages(selection_coefficients(d, "w", c("z1", "z2"), fitness_type = "continuous"))
  hc3 <- suppressMessages(selection_coefficients(d, "w", c("z1", "z2"), fitness_type = "continuous", se_type = "HC3"))
  expect_equal(attr(ols, "se_type"), "ols")
  expect_equal(attr(hc3, "se_type"), "hc3")
  expect_equal(hc3$Beta_Coefficient, ols$Beta_Coefficient)
  expect_false(isTRUE(all.equal(hc3$Standard_Error, ols$Standard_Error)))

  # HC3 is the sum over rows of (b - b without row i)(b - b without row i)'
  z <- as.data.frame(scale(d[c("z1", "z2")]))
  z$w <- d$w / mean(d$w)
  loo_se <- function(f) {
    full <- coef(lm(f, data = z))
    drops <- t(sapply(seq_len(n), function(i) full - coef(lm(f, data = z[-i, ]))))
    sqrt(diag(crossprod(drops)))
  }
  lin <- loo_se(w ~ z1 + z2)
  quad <- loo_se(w ~ z1 + z2 + I(z1^2) + I(z2^2) + z1:z2)
  expect_equal(hc3$Standard_Error, unname(c(lin[c("z1", "z2")], 2 * quad[c("I(z1^2)", "I(z2^2)")], quad["z1:z2"])),
               tolerance = 1e-10)
  # continuous p-values follow the HC3 errors
  expect_equal(hc3$P_Value[1], 2 * pt(-abs(hc3$Beta_Coefficient[1] / hc3$Standard_Error[1]), df = n - 3))
})

test_that("with survival HC3 moves the errors but not the logistic p-values", {
  b0 <- suppressWarnings(suppressMessages(selection_coefficients(bumpus, "survival", c("total_length", "weight"))))
  b3 <- suppressWarnings(suppressMessages(selection_coefficients(bumpus, "survival", c("total_length", "weight"), se_type = "hc3")))
  expect_equal(b3$P_Value, b0$P_Value)
  expect_false(isTRUE(all.equal(b3$Standard_Error, b0$Standard_Error)))
  g3 <- suppressWarnings(suppressMessages(selection_coefficients(bumpus, "survival", c("total_length", "weight"),
                                                                 group = "sex", return_grouped = TRUE, se_type = "hc3")))
  expect_equal(attr(g3, "se_type"), "hc3")
  r3 <- suppressWarnings(suppressMessages(selection_report(bumpus, "survival", c("total_length", "weight"), se_type = "hc3")))
  expect_equal(attr(r3, "se_type"), "hc3")
  expect_output(print(r3), "heteroscedasticity-consistent \\(HC3\\)")
  expect_error(selection_coefficients(bumpus, "survival", "weight", se_type = "hc9"))
})
