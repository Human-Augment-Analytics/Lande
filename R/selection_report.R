# ============================================================================
# selection_report
#
# Assemble a single standardised results table (selection differentials and
# linear/quadratic/correlational gradients) with identical columns across
# analyses, so results are directly comparable between studies.
# ============================================================================

#' Standardised selection analysis report
#'
#' Runs the standard Lande-Arnold workflow and returns one tidy table combining
#' selection differentials (S) and linear (beta), quadratic (gamma), and
#' correlational (gamma_ij) gradients. All estimates use traits standardised to
#' mean 0, SD 1 and relative fitness, so the table is directly comparable across
#' studies.
#'
#' @inheritParams selection_coefficients
#' @param include_differentials Logical; if \code{TRUE} (default) selection
#'   differentials S are included alongside the gradients.
#' @param digits Integer number of digits used when printing. Default is 4.
#'
#' @return A data frame of class \code{"selection_report"} with columns
#'   \code{Term}, \code{Type}, \code{Estimate}, \code{Std_Error}, and
#'   \code{P_Value}.
#' @export
#'
#' @examples
#' \dontrun{
#' selection_report(my_data, "fitness", c("trait1", "trait2"))
#' }
selection_report <- function(data,
                             fitness_col,
                             trait_cols,
                             fitness_type = c("auto", "binary", "continuous"),
                             standardize = TRUE,
                             group = NULL,
                             use_relative_for_fit = TRUE,
                             include_differentials = TRUE,
                             digits = 4) {
  fitness_type <- match.arg(fitness_type)

  gradients <- suppressMessages(selection_coefficients(
    data, fitness_col, trait_cols,
    fitness_type = fitness_type,
    standardize = standardize,
    group = group,
    use_relative_for_fit = use_relative_for_fit
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
    S <- vapply(trait_cols, function(t) {
      suppressWarnings(suppressMessages(
        selection_differential(data, fitness_col, t, standardized = !standardize, group = group)
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
  cat("\n")
  print(out, row.names = FALSE)
  cat("\nSignif: *** 0.001  ** 0.01  * 0.05  . 0.1\n")
  invisible(x)
}
