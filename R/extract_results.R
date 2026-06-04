# ============================================================================
# extract_results.R
# Extract selection coefficients from model results
# ============================================================================

#' @noRd
# internal utility: get coefficient-level p-value column name from summary()
.p_col_from_summary <- function(coef_mat) {
  pcols <- intersect(colnames(coef_mat), c("Pr(>|t|)", "Pr(>|z|)"))
  if (length(pcols)) pcols[1] else NA_character_
}

#' @noRd
# internal utility: safely extract ANOVA p-value (car::Anova)
.ano_p <- function(anova_obj, term) {
  if (is.null(anova_obj)) {
    return(NA_real_)
  }
  rn <- rownames(anova_obj)
  if (is.null(rn) || !term %in% rn) {
    return(NA_real_)
  }
  pcols <- intersect(colnames(anova_obj), c("Pr(>F)", "Pr(>Chisq)"))
  if (!length(pcols)) {
    return(NA_real_)
  }
  as.numeric(anova_obj[term, pcols[1]])
}
#' @noRd
# Helper to get appropriate summary and p-value column
.get_summary_and_pcol <- function(results) {
  if (results$fitness_type == "binary") {
    # cat("DEBUG: In .get_summary_and_pcol (binary branch)\n")
    # cat("DEBUG: names(results$summary):", paste(names(results$summary), collapse = ", "), "\n")

    # Check if summary$glm exists
    if (is.null(results$summary$glm)) {
      cat("DEBUG: summary$glm is NULL, falling back to OLS\n")
      sm <- results$summary$ols
      pcol <- .p_col_from_summary(coef(sm))
      if (is.na(pcol)) pcol <- "Pr(>|t|)"

      return(list(
        coef_mat = coef(sm),
        pcol = pcol,
        anova = results$anova,
        is_binary = FALSE
      ))
    }

    # cat("DEBUG: summary$glm exists!\n")

    sm_ols <- results$summary$ols
    sm_glm <- results$summary$glm

    coef_ols <- coef(sm_ols)
    coef_glm <- coef(sm_glm)

    if (is.null(coef_glm)) {
      stop("GLM coefficients are NULL. GLM model may have failed.")
    }

    pcol_glm <- .p_col_from_summary(coef(sm_glm))
    if (is.na(pcol_glm)) pcol_glm <- "Pr(>|z|)"

    if (!is.matrix(coef_glm)) {
      if (is.vector(coef_glm) && length(coef_glm) > 0) {
        coef_glm <- matrix(coef_glm, nrow = 1)
        rownames(coef_glm) <- names(coef_glm)
      }
    }

    return(list(
      coef_mat_ols = coef_ols,
      coef_mat_glm = coef_glm,
      pcol_glm = pcol_glm,
      anova = results$anova,
      is_binary = TRUE
    ))
  } else {
    # Continuous case
    sm <- results$summary
    pcol <- .p_col_from_summary(coef(sm))
    if (is.na(pcol)) pcol <- "Pr(>|t|)"

    return(list(
      coef_mat = coef(sm),
      pcol = pcol,
      anova = results$anova,
      is_binary = FALSE
    ))
  }
}

#' Extract linear selection coefficients
#'
#' @param trait_cols A character vector of trait column names.
#' @param results A model results object returned by \code{analyze_linear_selection()}.
#'
#' @return A data frame with linear selection coefficients and statistics.
#' @export
extract_linear_coefficients <- function(trait_cols, results) {
  obj <- .get_summary_and_pcol(results)

  if (obj$is_binary) {
    # Binary case: gradients from OLS, p-values from GLM
    coef_ols <- obj$coef_mat_ols
    coef_glm <- obj$coef_mat_glm
    pcol_glm <- obj$pcol_glm
    anova_obj <- obj$anova

    keep <- intersect(rownames(coef_ols), trait_cols)

    rows <- lapply(keep, function(t) {
      # p-value from GLM
      if (t %in% rownames(coef_glm)) {
        p_val <- coef_glm[t, pcol_glm]
      } else {
        p_val <- .ano_p(anova_obj, t)
      }

      data.frame(
        Term = t,
        Type = "Linear",
        Beta_Coefficient = coef_ols[t, "Estimate"],
        Standard_Error = coef_ols[t, "Std. Error"],
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  } else {
    # Continuous case: everything from OLS
    coef_mat <- obj$coef_mat
    pcol <- obj$pcol
    anova_obj <- obj$anova

    keep <- intersect(rownames(coef_mat), trait_cols)

    rows <- lapply(keep, function(t) {
      p_val <- if (!is.na(pcol) && t %in% rownames(coef_mat)) {
        coef_mat[t, pcol]
      } else {
        .ano_p(anova_obj, t)
      }

      data.frame(
        Term = t,
        Type = "Linear",
        Beta_Coefficient = coef_mat[t, "Estimate"],
        Standard_Error = coef_mat[t, "Std. Error"],
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  }
}

#' Extract quadratic selection coefficients
#'
#' @param trait_cols A character vector of trait column names.
#' @param results A model results object returned by \code{analyze_nonlinear_selection()}.
#'
#' @return A data frame with quadratic selection coefficients (doubled estimates) and statistics.
#' @export
extract_quadratic_coefficients <- function(trait_cols, results) {
  obj <- .get_summary_and_pcol(results)

  if (obj$is_binary) {
    # Binary case: gradients from OLS, p-values from GLM
    coef_ols <- obj$coef_mat_ols
    coef_glm <- obj$coef_mat_glm
    pcol_glm <- obj$pcol_glm
    anova_obj <- obj$anova

    rows <- lapply(trait_cols, function(t) {
      term <- paste0("I(", t, "^2)")
      if (!term %in% rownames(coef_ols)) {
        return(NULL)
      }

      # Multiply estimate and SE by 2 for gamma_ii
      est <- 2 * coef_ols[term, "Estimate"]
      se <- 2 * coef_ols[term, "Std. Error"]

      # p-value from GLM
      if (term %in% rownames(coef_glm)) {
        p_val <- coef_glm[term, pcol_glm]
      } else {
        p_val <- .ano_p(anova_obj, term)
      }

      data.frame(
        Term = paste0(t, "²"),
        Type = "Quadratic",
        Beta_Coefficient = as.numeric(est),
        Standard_Error = as.numeric(se),
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  } else {
    # Continuous case: everything from OLS
    coef_mat <- obj$coef_mat
    pcol <- obj$pcol
    anova_obj <- obj$anova

    rows <- lapply(trait_cols, function(t) {
      term <- paste0("I(", t, "^2)")
      if (!term %in% rownames(coef_mat)) {
        return(NULL)
      }

      est <- 2 * coef_mat[term, "Estimate"]
      se <- 2 * coef_mat[term, "Std. Error"]
      p_val <- if (!is.na(pcol) && term %in% rownames(coef_mat)) {
        coef_mat[term, pcol]
      } else {
        .ano_p(anova_obj, term)
      }

      data.frame(
        Term = paste0(t, "²"),
        Type = "Quadratic",
        Beta_Coefficient = as.numeric(est),
        Standard_Error = as.numeric(se),
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  }
}


#' Extract correlational selection coefficients
#'
#' @param trait_cols A character vector of trait column names.
#' @param results A model results object returned by \code{analyze_nonlinear_selection()}.
#'
#' @return A data frame with correlational (interaction) selection coefficients and statistics.
#' @export
extract_interaction_coefficients <- function(trait_cols, results) {
  if (length(trait_cols) < 2) {
    return(data.frame(
      Term = character(), Type = character(),
      Beta_Coefficient = numeric(), Standard_Error = numeric(),
      P_Value = numeric(), check.names = FALSE
    ))
  }

  obj <- .get_summary_and_pcol(results)
  pairs <- utils::combn(trait_cols, 2, simplify = FALSE)

  if (obj$is_binary) {
    # Binary case: gradients from OLS, p-values from GLM
    coef_ols <- obj$coef_mat_ols
    coef_glm <- obj$coef_mat_glm
    pcol_glm <- obj$pcol_glm
    anova_obj <- obj$anova

    rows <- lapply(pairs, function(p) {
      raw_term <- paste0(p[1], ":", p[2])
      if (!raw_term %in% rownames(coef_ols)) {
        return(NULL)
      }

      est <- coef_ols[raw_term, "Estimate"]
      se <- coef_ols[raw_term, "Std. Error"]

      # p-value from GLM
      if (raw_term %in% rownames(coef_glm)) {
        p_val <- coef_glm[raw_term, pcol_glm]
      } else {
        p_val <- .ano_p(anova_obj, raw_term)
      }

      data.frame(
        Term = paste(p[1], "×", p[2]),
        Type = "Correlational",
        Beta_Coefficient = as.numeric(est),
        Standard_Error = as.numeric(se),
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  } else {
    # Continuous case: everything from OLS
    coef_mat <- obj$coef_mat
    pcol <- obj$pcol
    anova_obj <- obj$anova

    rows <- lapply(pairs, function(p) {
      raw_term <- paste0(p[1], ":", p[2])
      if (!raw_term %in% rownames(coef_mat)) {
        return(NULL)
      }

      est <- coef_mat[raw_term, "Estimate"]
      se <- coef_mat[raw_term, "Std. Error"]
      p_val <- if (!is.na(pcol) && raw_term %in% rownames(coef_mat)) {
        coef_mat[raw_term, pcol]
      } else {
        .ano_p(anova_obj, raw_term)
      }

      data.frame(
        Term = paste(p[1], "×", p[2]),
        Type = "Correlational",
        Beta_Coefficient = as.numeric(est),
        Standard_Error = as.numeric(se),
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  }
}
