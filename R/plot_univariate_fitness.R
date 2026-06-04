# ======================================================
# plot_univariate_fitness.R
# Visualize univariate correlated fitness function
# ======================================================


#' Plot Univariate Correlated Fitness Function
#'
#' Visualizes the univariate fitness function estimated via spline.
#'
#' @param uni Output object of class \code{"univariate_fitness"} from \code{univariate_spline()}.
#' @param trait_col Character string specifying the trait column name.
#' @param title Optional character string for the plot title.
#' @param show_points Logical indicating whether to show the original data points. Default is \code{TRUE}.
#' @param point_alpha Numeric value for point transparency. Default is 0.25.
#' @param ribbon_fill Color for the confidence interval ribbon. Default is \code{NA} (transparent).
#' @param ribbon_linetype Linetype for the confidence interval bounds. Default is \code{"dashed"}.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#'
#' @return A \code{ggplot} object representing the univariate fitness function.
#' @export
plot_univariate_fitness <- function(uni,
                                    trait_col,
                                    title = NULL,
                                    show_points = TRUE,
                                    point_alpha = 0.25,
                                    ribbon_fill = NA,
                                    ribbon_linetype = "dashed",
                                    ...) {
  # Input validation
  if (!inherits(uni, "univariate_fitness")) {
    warning("Object is not of class 'univariate_fitness'")
  }

  stopifnot(is.list(uni), "grid" %in% names(uni))
  grid <- uni$grid

  if (!all(c(trait_col, "fit", "lwr", "upr") %in% names(grid))) {
    stop("`uni$grid` must contain columns: trait, fit, lwr, upr")
  }

  # Set default title if not provided
  if (is.null(title)) {
    title <- "Univariate Correlated Fitness Function"
  }

  # Determine y-axis limits based on fitness type
  if (uni$fitness_type == "binary") {
    y_limits <- c(0, 1)
    y_label <- "Predicted survival probability"
  } else {
    # For continuous, use data range
    y_limits <- c(
      min(grid$lwr, na.rm = TRUE),
      max(grid$upr, na.rm = TRUE)
    )
    y_label <- "Predicted fitness"
  }


  p <- ggplot2::ggplot()

  # Add confidence ribbon (if fill is specified)
  if (!is.na(ribbon_fill)) {
    p <- p + ggplot2::geom_ribbon(
      data = grid,
      ggplot2::aes(
        x = .data[[trait_col]],
        ymin = .data[["lwr"]],
        ymax = .data[["upr"]]
      ),
      fill = ribbon_fill,
      alpha = 0.2
    )
  }

  # Add confidence bounds (as lines)
  p <- p + ggplot2::geom_line(
    data = grid,
    ggplot2::aes(
      x = .data[[trait_col]],
      y = .data[["lwr"]]
    ),
    color = "black",
    linetype = ribbon_linetype,
    linewidth = 0.8
  ) +
    ggplot2::geom_line(
      data = grid,
      ggplot2::aes(
        x = .data[[trait_col]],
        y = .data[["upr"]]
      ),
      color = "black",
      linetype = ribbon_linetype,
      linewidth = 0.8
    )

  # Add main fit line
  p <- p + ggplot2::geom_line(
    data = grid,
    ggplot2::aes(
      x = .data[[trait_col]],
      y = .data[["fit"]]
    ),
    color = "black",
    linewidth = 1
  )

  # Add data points if requested
  if (show_points && "data" %in% names(uni)) {
    # If original data is stored in the object
    data_points <- uni$data
  } else if (show_points && exists("original_data", where = uni)) {
    # Alternative storage location
    data_points <- uni$original_data
  } else {
    data_points <- NULL
  }

  if (show_points && !is.null(data_points) &&
    trait_col %in% names(data_points) &&
    ".y" %in% names(data_points)) {
    # For binary fitness, jitter points slightly to avoid overlap
    if (uni$fitness_type == "binary") {
      p <- p + ggplot2::geom_point(
        data = data_points,
        ggplot2::aes(
          x = .data[[trait_col]],
          y = .data[[".y"]]
        ),
        size = 2,
        alpha = point_alpha,
        position = ggplot2::position_jitter(width = 0, height = 0.02)
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = data_points,
        ggplot2::aes(
          x = .data[[trait_col]],
          y = .data[[".y"]]
        ),
        size = 2,
        alpha = point_alpha
      )
    }
  }

  p <- p + ggplot2::labs(
    x = trait_col,
    y = y_label,
    title = title,
    subtitle = "Individual fitness (correlated fitness function)",
    ...
  ) +
    ggplot2::ylim(y_limits) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = ggplot2::element_text(color = "black", size = 10),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40")
    )

  return(p)
}
