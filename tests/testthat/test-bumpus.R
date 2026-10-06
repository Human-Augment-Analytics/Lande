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

test_that("on logged traits the gradients within each sex match Janzen and Stern (1998)", {
  traits9 <- c("total_length", "wingspread", "weight", "head_length", "humerus",
               "femur", "tibiotarsus", "skull_width", "sternum")
  # their Tables 1 and 2: least-squares gradients and standard errors
  js <- list(male = list(beta = c(-0.516, 0.097, -0.272, 0.051, 0.164, 0.116, -0.018, 0.137, 0.164),
                         se = c(0.100, 0.128, 0.093, 0.093, 0.168, 0.169, 0.128, 0.085, 0.087)),
             female = list(beta = c(-0.213, -0.143, -0.531, 0.061, 0.328, -0.174, 0.427, -0.019, 0.169),
                           se = c(0.284, 0.318, 0.253, 0.322, 0.369, 0.360, 0.349, 0.260, 0.237)))
  d <- bumpus
  d[traits9] <- log(d[traits9])
  for (sx in names(js)) {
    r <- suppressWarnings(suppressMessages(
      selection_coefficients(d[d$sex == sx, ], "survival", traits9, fitness_type = "binary")))
    lin <- r[r$Type == "Linear", ]
    lin <- lin[match(traits9, lin$Term), ]
    expect_lt(max(abs(lin$Beta_Coefficient - js[[sx]]$beta)), 0.001)
    expect_lt(max(abs(lin$Standard_Error - js[[sx]]$se)), 0.001)
  }
})
