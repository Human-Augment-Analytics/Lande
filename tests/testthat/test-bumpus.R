# Validation against the bundled Bumpus sparrow dataset: confirm standardisation
# and the factor-of-2 quadratic convention on real, published data.

test_that("bumpus dataset is available and well-formed", {
  expect_s3_class(bumpus, "data.frame")
  expect_equal(nrow(bumpus), 136)
  expect_true(all(c("survival", "weight", "total_length") %in% names(bumpus)))
  expect_true(all(bumpus$survival %in% c(0, 1)))
})

test_that("prepare_selection_data standardises bumpus traits", {
  out <- suppressWarnings(suppressMessages(
    prepare_selection_data(bumpus, "survival", c("weight", "total_length"),
                           name_relative = "survival_relative")
  ))
  expect_equal(mean(out$weight), 0, tolerance = 1e-8)
  expect_equal(sd(out$weight), 1, tolerance = 1e-8)
})

test_that("quadratic gradient on bumpus is twice the OLS squared-term coefficient", {
  traits <- c("weight", "total_length")
  res <- suppressWarnings(suppressMessages(
    selection_coefficients(bumpus, "survival", traits, fitness_type = "binary")
  ))

  # Manual Lande-Arnold reference: standardised traits, relative fitness,
  # full second-order OLS, gamma = 2 * quadratic coefficient.
  z <- as.data.frame(scale(bumpus[, traits]))
  names(z) <- traits
  w <- bumpus$survival / mean(bumpus$survival)
  fit <- lm(w ~ weight + total_length + I(weight^2) + I(total_length^2) +
              weight:total_length, data = cbind(w = w, z))
  gamma_weight <- 2 * coef(fit)["I(weight^2)"]

  got <- res$Beta_Coefficient[res$Type == "Quadratic" & startsWith(res$Term, "weight")]
  expect_equal(unname(got), unname(gamma_weight), tolerance = 1e-6)
})
