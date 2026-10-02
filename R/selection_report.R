# ============================================================================
# selection_report
#
# One table of selection differentials and linear, quadratic and
# correlational gradients, with the same columns in every analysis.
# ============================================================================

#' Standardised selection analysis report
#'
#' Runs a Lande-Arnold analysis and returns one table of selection
#' differentials (S) and linear (beta), quadratic (gamma) and correlational
#' (gamma_ij) gradients. All estimates use traits standardised to mean 0, SD 1
#' and relative fitness, so tables from different studies can be compared.
#'
#' @inheritParams selection_coefficients
#' @param include_differentials Logical; if \code{TRUE} (default) selection
#'   differentials S are included alongside the gradients.
#' @param digits Integer number of digits used when printing. Default is 4.
#'
#' @details Differentials and gradients are computed on the same individuals
#'   (those with complete fitness and trait values) and on the same trait and
#'   fitness scales. S is the
#'   population covariance (divides by n) while the OLS gradient on
#'   sd-standardised traits corresponds to the sample covariance (n - 1), so
#'   for a single trait the Differential and Linear rows differ by the factor
#'   (n - 1) / n.
#'
#' @return A data frame of class \code{"selection_report"} with columns
#'   \code{Term}, \code{Type}, \code{Estimate}, \code{Std_Error}, and
#'   \code{P_Value}.
#' @export
#'
#' @examples
#' selection_report(bumpus, "survival", c("total_length", "weight", "humerus"))
selection_report <- function(data,
                             fitness_col,
                             trait_cols,
                             fitness_type = c("auto", "binary", "count", "continuous"),
                             standardize = TRUE,
                             group = NULL,
                             use_relative_for_fit = TRUE,
                             include_differentials = TRUE,
                             digits = 4,
                             se_type = c("ols", "hc3")) {
  fitness_type <- match.arg(fitness_type)
  se_type <- .se_type_arg(se_type)

  gradients <- suppressMessages(selection_coefficients(
    data, fitness_col, trait_cols,
    fitness_type = fitness_type,
    standardize = standardize,
    group = group,
    use_relative_for_fit = use_relative_for_fit,
    se_type = se_type
  ))

  tab <- data.frame(
    Term = gradients$Term,
    Type = gradients$Type,
    Estimate = gradients$Beta_Coefficient,
    Std_Error = gradients$Standard_Error,
    P_Value = gradients$P_Value,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (include_differentials) {
    # Prepare the data exactly as selection_coefficients() does and keep only
    # the rows that reach the gradient models, so S and beta share a sample
    # and a scale (standardised traits, relative fitness unless the caller
    # asked for absolute fitness).
    rel_col <- paste0(fitness_col, "_relative")
    prep <- suppressWarnings(suppressMessages(prepare_selection_data(
      data, fitness_col, trait_cols,
      standardize = standardize, group = group, add_relative = TRUE,
      na_action = "none", name_relative = rel_col
    )))
    # unlabelled rows count as one more group
    if (!is.null(group)) {
      lab <- as.character(prep[[group]])
      lab[is.na(lab)] <- "(no group)"
      prep[[group]] <- lab
    }
    prep <- prep[stats::complete.cases(prep[, c(fitness_col, trait_cols), drop = FALSE]), , drop = FALSE]
    w_col <- if (use_relative_for_fit && rel_col %in% names(prep)) rel_col else fitness_col
    S <- vapply(trait_cols, function(t) {
      suppressWarnings(suppressMessages(
        selection_differential(prep, w_col, t, standardized = TRUE, use_relative = FALSE, group = group)
      ))
    }, numeric(1))

    diff_tab <- data.frame(
      Term = trait_cols,
      Type = "Differential",
      Estimate = as.numeric(S),
      Std_Error = NA_real_,
      P_Value = NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    tab <- rbind(diff_tab, tab)
  }

  rownames(tab) <- NULL
  attr(tab, "digits") <- digits
  attr(tab, "fitness_type") <- attr(gradients, "fitness_type_used")
  attr(tab, "p_model") <- c(linear = attr(gradients, "model_family_used"),
                            quadratic = attr(gradients, "model_family_quadratic"))
  attr(tab, "se_type") <- se_type
  attr(tab, "scale") <- paste0(
    if (standardize) "standardised traits" else "unstandardised traits", ", ",
    if (use_relative_for_fit) "relative fitness" else "absolute fitness"
  )
  class(tab) <- c("selection_report", "data.frame")
  tab
}

#' Print a selection report
#'
#' @param x An object of class \code{"selection_report"}.
#' @param ... Additional arguments (ignored).
#' @return The input object \code{x}, invisibly.
#' @export
print.selection_report <- function(x, ...) {
  digits <- attr(x, "digits") %||% 4
  stars <- function(p) {
    ifelse(is.na(p), "",
      ifelse(p < 0.001, "***",
        ifelse(p < 0.01, "**",
          ifelse(p < 0.05, "*",
            ifelse(p < 0.1, ".", "")))))
  }

  out <- data.frame(
    Term = x$Term,
    Type = x$Type,
    Estimate = round(x$Estimate, digits),
    Std_Error = round(x$Std_Error, digits),
    P_Value = round(x$P_Value, digits),
    Sig = stars(x$P_Value),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cat("Selection analysis (", attr(x, "scale") %||% "standardised traits, relative fitness", ")\n", sep = "")
  if (!is.null(attr(x, "fitness_type"))) {
    cat("Fitness type:", attr(x, "fitness_type"), "\n")
  }
  p_model <- attr(x, "p_model")
  if (!is.null(p_model) && any(p_model != "gaussian")) {
    names_of <- c("binomial(logit)" = "logistic", "poisson(log)" = "Poisson", "negative binomial" = "negative binomial")
    p_model <- ifelse(p_model %in% names(names_of), names_of[p_model], p_model)
    which_model <- if (length(unique(p_model)) == 1) paste("a", p_model[1], "model on the same terms") else
      paste("a", p_model[1], "model for the linear gradients and a", p_model[2], "model for the rest")
    cat("p-values are from ", which_model, "\n", sep = "")
  }
  if (identical(attr(x, "se_type"), "hc3")) {
    cat("Standard errors: heteroscedasticity-consistent (HC3)\n")
  }
  cat("\n")
  print(out, row.names = FALSE)
  cat("\nSignif: *** 0.001  ** 0.01  * 0.05  . 0.1\n")
  invisible(x)
}
