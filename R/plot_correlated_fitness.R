# ======================================================
# plot_correlated_fitness.R
# Visualize correlated fitness surface (individual-level)
#
# IMPORTANT CONCEPT:
# - Correlated Fitness Surface: individual fitness (w ~ z₁, z₂)
# - This is DIFFERENT from adaptive landscape (population-level)
#
# KEY PRINCIPLE:
# - Traits should already be standardized
# - DO NOT standardize again within this function
# ======================================================

#' Plot Correlated Fitness Surface
#'
#' @param tps Output list from \code{correlated_fitness_surface()}.
#' @param trait_cols Character vector of length 2 specifying the trait column names.
#' @param bins Integer specifying the number of contour bins. Default is 12.
#' @param point_alpha Numeric value for point transparency. Default is 0.7.
#' @param show_points Logical indicating whether to show original data points. Default is \code{FALSE}.
#' @param show_optimum Logical indicating whether to mark the optimum point. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#'
#' @return A \code{ggplot} object representing the correlated fitness surface.
#' @export
plot_correlated_fitness <- function(
  tps,
  trait_cols,
  bins = 12,
  point_alpha = 0.7,
  show_points = FALSE,
  show_optimum = TRUE,
  ...
) {
  # Input validation
  stopifnot(is.list(tps), "grid" %in% names(tps))

  df <- tps$grid

  if (length(trait_cols) != 2) {
    stop("trait_cols must be length 2.")
  }
  if (!all(trait_cols %in% names(df))) {
    stop("Trait columns not found in tps$grid.")
  }

  # Find fitness column
  if (".fit" %in% names(df)) {
    fitness_col <- ".fit"
  } else if ("fitness" %in% names(df)) {
    fitness_col <- "fitness"
  } else if ("pred" %in% names(df)) {
    fitness_col <- "pred"
  } else {
    stop("No fitness column found in grid.")
  }

  p <- ggplot2::ggplot(df) +
    # Filled contours
    ggplot2::geom_contour_filled(
      ggplot2::aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]],
        z = .data[[fitness_col]]
      ),
      bins = bins,
      inherit.aes = FALSE
    ) +
    # Contour lines
    ggplot2::geom_contour(
      ggplot2::aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]],
        z = .data[[fitness_col]]
      ),
      color = "black",
      alpha = 0.3,
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    # Labels
    ggplot2::labs(
      x = trait_cols[1],
      y = trait_cols[2],
      fill = "Fitness",
      title = "Correlated Fitness Surface",
      subtitle = "Individual-level fitness",
      ...
    ) +
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

  # Add data points if requested
  if (show_points && "original_data" %in% names(tps)) {
    orig_data <- tps$original_data
    if (!is.null(orig_data) && all(trait_cols %in% names(orig_data))) {
      p <- p + ggplot2::geom_point(
        data = orig_data,
        ggplot2::aes(
          x = .data[[trait_cols[1]]],
          y = .data[[trait_cols[2]]]
        ),
        size = 1.5,
        alpha = point_alpha,
        color = "gray40"
      )
    }
  }

  if (show_optimum) {
    opt <- df[which.max(df[[fitness_col]]), ]
    p <- p + ggplot2::geom_point(
      data = opt,
      ggplot2::aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]]
      ),
      color = "gold",
      size = 4,
      shape = 18
    )
  }

  return(p)
}


# ======================================================
# plot_correlated_fitness_enhanced
# Enhanced version with original data points colored by fitness
# ======================================================

#' Plot Enhanced Correlated Fitness Surface
#'
#' @param tps Output list from \code{correlated_fitness_surface()}.
#' @param trait_cols Character vector of length 2. If \code{NULL}, inferred from grid.
#' @param original_data Optional data frame of original data.
#' @param fitness_col Optional character string specifying the fitness column for coloring points.
#' @param bins Integer specifying the number of contour bins. Default is 12.
#' @param point_alpha Numeric value for point transparency. Default is 0.7.
#' @param show_optimum Logical indicating whether to mark the optimum point. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#'
#' @return A \code{ggplot} object with enhanced visualizations.
#' @export
plot_correlated_fitness_enhanced <- function(
  tps,
  trait_cols = NULL,
  original_data = NULL,
  fitness_col = NULL,
  bins = 12,
  point_alpha = 0.7,
  show_optimum = TRUE,
  ...
) {
  # Input validation
  df <- tps$grid

  # Determine trait columns
  if (!is.null(trait_cols) && length(trait_cols) == 2) {
    trait1 <- trait_cols[1]
    trait2 <- trait_cols[2]
    cat("Using provided traits:", trait1, trait2, "\n")
  } else {
    # Infer from grid columns
    possible_traits <- setdiff(
      names(df),
      c(".fit", "fitness", "pred", "fit", "lwr", "upr", "type", "surface_type")
    )
    if (length(possible_traits) >= 2) {
      trait1 <- possible_traits[1]
      trait2 <- possible_traits[2]
      cat("Inferred traits:", trait1, trait2, "\n")
    } else {
      stop("Cannot determine trait columns. Please provide trait_cols.")
    }
  }

  # Determine fitness column in grid
  if (".fit" %in% names(df)) {
    fit_col <- ".fit"
  } else if ("fitness" %in% names(df)) {
    fit_col <- "fitness"
  } else {
    stop("No fitness column found in grid")
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_contour_filled(
      data = df,
      ggplot2::aes(
        x = .data[[trait1]],
        y = .data[[trait2]],
        z = .data[[fit_col]]
      ),
      bins = bins
    ) +
    ggplot2::geom_contour(
      data = df,
      ggplot2::aes(
        x = .data[[trait1]],
        y = .data[[trait2]],
        z = .data[[fit_col]]
      ),
      color = "black",
      alpha = 0.3,
      linewidth = 0.3
    ) +
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

  if (!is.null(original_data) && !is.null(fitness_col)) {
    if (all(c(trait1, trait2, fitness_col) %in% names(original_data))) {
      unique_vals <- unique(original_data[[fitness_col]])
      is_binary <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))

      if (is_binary) {
        p <- p +
          ggplot2::geom_point(
            data = original_data,
            ggplot2::aes(
              x = .data[[trait1]],
              y = .data[[trait2]],
              color = as.factor(.data[[fitness_col]])
            ),
            size = 2,
            alpha = point_alpha
          ) +
          ggplot2::scale_color_manual(
            values = c("0" = "red", "1" = "green"),
            labels = c("0" = "Perished", "1" = "Survived"),
            name = "Outcome"
          )
      } else {
        p <- p +
          ggplot2::geom_point(
            data = original_data,
            ggplot2::aes(
              x = .data[[trait1]],
              y = .data[[trait2]],
              color = .data[[fitness_col]]
            ),
            size = 2,
            alpha = point_alpha
          ) +
          ggplot2::scale_color_viridis_c(name = "Fitness")
      }
    } else {
      warning("Required columns not found in original_data")
    }
  }

  if (show_optimum) {
    opt <- df[which.max(df[[fit_col]]), ]
    p <- p +
      ggplot2::geom_point(
        data = opt,
        ggplot2::aes(
          x = .data[[trait1]],
          y = .data[[trait2]]
        ),
        color = "gold",
        size = 4,
        shape = 18
      ) +
      # Add label for optimum (optional)
      ggplot2::annotate(
        "text",
        x = opt[[trait1]],
        y = opt[[trait2]],
        label = "Optimum",
        vjust = -1,
        size = 3,
        color = "gray30"
      )
  }

  p <- p + ggplot2::labs(
    x = trait1,
    y = trait2,
    title = "Correlated Fitness Surface",
    subtitle = "Individual-level fitness",
    ...
  )

  return(p)
}
