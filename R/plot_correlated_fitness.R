# ======================================================
# plot_correlated_fitness.R
# Visualize correlated fitness surface (individual-level)
#
# IMPORTANT CONCEPT:
# - Correlated Fitness Surface: individual fitness (w ~ z1, z2)
# - This is DIFFERENT from adaptive landscape (population-level)
#
# KEY PRINCIPLE:
# - Traits should already be standardized
# - DO NOT standardize again within this function
# ======================================================

#' @noRd
# internal utility: what to draw for a surface. When the grid carries the
# unmasked predictions the full surface is drawn and the blanked region is
# covered: the outside of the hull by a polygon, cells too far from any
# individual by a white contour band on the distance field. Both give clean
# edges. Otherwise masked (NA) cells are simply dropped.
.surface_layers <- function(tps, df, trait1, trait2, fit_col) {
  if (!".fit_all" %in% names(df)) {
    return(list(df = df[!is.na(df[[fit_col]]), , drop = FALSE], fit_col = fit_col, cover = NULL))
  }
  cover <- list()
  hull <- tps$hull
  if (!is.null(hull)) {
    frame <- .outside_frame(hull, trait1, trait2, range(df[[trait1]]), range(df[[trait2]]))
    cover <- c(cover, list(ggplot2::geom_polygon(
      data = frame,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = "white", colour = "white", linewidth = 0.5, inherit.aes = FALSE
    )))
  }
  too_far <- tps$too_far
  if (!is.null(too_far) && ".dist" %in% names(df) && any(df$.dist > too_far)) {
    # two bands of the distance field, near and far; only the far one is drawn
    cover <- c(cover, list(
      ggplot2::layer(
        stat = .StatDistanceBand, geom = "polygon", position = "identity",
        data = df,
        mapping = ggplot2::aes(
          x = .data[[trait1]], y = .data[[trait2]], z = .data$.dist,
          alpha = ggplot2::after_stat(.data$level)
        ),
        params = list(breaks = c(-1, too_far, max(df$.dist) + 1), fill = "white", na.rm = FALSE),
        inherit.aes = FALSE
      ),
      ggplot2::scale_alpha_manual(values = c(0, 1), guide = "none")
    ))
  }
  list(df = df, fit_col = ".fit_all", cover = cover)
}

#' @noRd
# internal utility: the region outside the hull as one ring: round the grid
# rectangle, along a bridge to the nearest hull vertex, once round the hull
# the opposite way, and back. Renders as a frame with a hull-shaped window.
.outside_frame <- function(hull, trait1, trait2, xr, yr) {
  h <- cbind(hull[[trait1]], hull[[trait2]])
  h <- h[-nrow(h), , drop = FALSE]
  corner <- c(xr[1], yr[1])
  i <- which.min((h[, 1] - corner[1])^2 + (h[, 2] - corner[2])^2)
  ring <- h[c(i:nrow(h), seq_len(i)), , drop = FALSE]
  frame <- rbind(
    c(xr[1], yr[1]), c(xr[2], yr[1]), c(xr[2], yr[2]), c(xr[1], yr[2]), c(xr[1], yr[1]),
    ring, c(xr[1], yr[1])
  )
  data.frame(x = frame[, 1], y = frame[, 2])
}

#' @noRd
# the filled-contour stat without its fill mapping, so the white band drawn
# over cells far from the data adds no keys to the fitness legend
.StatDistanceBand <- ggplot2::ggproto(
  "StatDistanceBand", ggplot2::StatContourFilled,
  default_aes = ggplot2::aes(order = ggplot2::after_stat(level))
)

#' @noRd
# internal utility: group means (open circles, labelled) and the highest
# point of the surface within each group's hull (filled triangles), joined by
# a dashed line. Only present when the surface was fitted with a group.
.group_layers <- function(tps, trait1, trait2) {
  g <- tps$groups
  if (is.null(g) || !nrow(g)) return(NULL)
  m1 <- paste0("mean_", trait1)
  m2 <- paste0("mean_", trait2)
  p1 <- paste0("peak_", trait1)
  p2 <- paste0("peak_", trait2)
  if (!all(c(m1, m2, p1, p2) %in% names(g))) return(NULL)
  with_peak <- g[!is.na(g[[p1]]), , drop = FALSE]
  list(
    ggplot2::geom_segment(
      data = with_peak,
      ggplot2::aes(x = .data[[m1]], y = .data[[m2]], xend = .data[[p1]], yend = .data[[p2]]),
      colour = "black", linewidth = 0.4, linetype = "dashed", inherit.aes = FALSE
    ),
    ggplot2::geom_point(
      data = g, ggplot2::aes(x = .data[[m1]], y = .data[[m2]], shape = "Group mean"),
      fill = "white", colour = "black", size = 3, inherit.aes = FALSE
    ),
    ggplot2::geom_point(
      data = with_peak, ggplot2::aes(x = .data[[p1]], y = .data[[p2]], shape = "Group peak"),
      fill = "black", colour = "white", size = 2.8, inherit.aes = FALSE
    ),
    ggplot2::geom_text(
      data = g, ggplot2::aes(x = .data[[m1]], y = .data[[m2]], label = .data$group),
      vjust = -1, size = 3, colour = "black", inherit.aes = FALSE
    )
  )
}

#' @noRd
# one key for the marks drawn on a surface: the optimum and, with a group,
# the group means and peaks. Only the marks present appear in it.
.mark_key <- function() {
  ggplot2::scale_shape_manual(
    name = NULL,
    values = c("Optimum" = 18, "Group mean" = 21, "Group peak" = 24),
    guide = ggplot2::guide_legend(order = 1)
  )
}

#' Plot Correlated Fitness Surface
#'
#' @param tps Output list from \code{correlated_fitness_surface()}.
#' @param trait_cols Character vector of length 2 specifying the trait column names.
#' @param bins Integer specifying the number of contour bins. Default is 12.
#' @param point_alpha Numeric value for point transparency. Default is 0.7.
#' @param show_points Logical indicating whether to show original data points. Default is \code{FALSE}.
#' @param show_optimum Logical indicating whether to mark the optimum point. Default is \code{TRUE}.
#' @param show_groups Logical; when the surface was fitted with a \code{group},
#'   draw each group's mean (open circle, labelled) and the highest point of
#'   the surface within that group's hull (filled triangle). Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#'
#' @return A \code{ggplot} object representing the correlated fitness surface.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#' plot_correlated_fitness(surf, c("total_length", "weight"), show_points = TRUE)
#' @export
plot_correlated_fitness <- function(
  tps,
  trait_cols,
  bins = 12,
  point_alpha = 0.7,
  show_points = FALSE,
  show_optimum = TRUE,
  show_groups = TRUE,
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

  # Masked surfaces: draw the full surface, then cover the outside of the data hull
  layers <- .surface_layers(tps, df, trait_cols[1], trait_cols[2], fitness_col)
  draw_df <- layers$df
  z_col <- layers$fit_col

  p <- ggplot2::ggplot(draw_df) +
    # Filled contours
    ggplot2::geom_contour_filled(
      ggplot2::aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]],
        z = .data[[z_col]]
      ),
      bins = bins,
      inherit.aes = FALSE
    ) +
    # Contour lines
    ggplot2::geom_contour(
      ggplot2::aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]],
        z = .data[[z_col]]
      ),
      color = "black",
      alpha = 0.3,
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    layers$cover +
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
        y = .data[[trait_cols[2]]],
        shape = "Optimum"
      ),
      color = "gold",
      size = 4
    )
  }

  if (show_groups) {
    p <- p + .group_layers(tps, trait_cols[1], trait_cols[2])
  }
  p <- p + .mark_key()

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
#' @param show_groups Logical; when the surface was fitted with a \code{group},
#'   draw each group's mean and the highest point of the surface within that
#'   group's hull. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#'
#' @return A \code{ggplot} object with enhanced visualizations.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#' plot_correlated_fitness_enhanced(surf, c("total_length", "weight"),
#'                                  original_data = prep, fitness_col = "survival")
#' @export
plot_correlated_fitness_enhanced <- function(
  tps,
  trait_cols = NULL,
  original_data = NULL,
  fitness_col = NULL,
  bins = 12,
  point_alpha = 0.7,
  show_optimum = TRUE,
  show_groups = TRUE,
  ...
) {
  # Input validation
  df <- tps$grid

  # Determine trait columns
  if (!is.null(trait_cols) && length(trait_cols) == 2) {
    trait1 <- trait_cols[1]
    trait2 <- trait_cols[2]
    message("Using provided traits: ", trait1, ", ", trait2)
  } else if (!is.null(tps$trait_cols) && length(tps$trait_cols) == 2) {
    trait1 <- tps$trait_cols[1]
    trait2 <- tps$trait_cols[2]
    message("Using traits from surface object: ", trait1, ", ", trait2)
  } else {
    # Infer from grid columns, excluding any grouping column carried along
    possible_traits <- setdiff(
      names(df),
      c(".fit", ".fit_all", ".inside", ".dist", "fitness", "pred", "fit", "lwr", "upr", "type", "surface_type", tps$group_used)
    )
    if (length(possible_traits) >= 2) {
      trait1 <- possible_traits[1]
      trait2 <- possible_traits[2]
      message("Inferred traits: ", trait1, ", ", trait2)
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

  # Masked surfaces: draw the full surface, then cover the outside of the data hull
  layers <- .surface_layers(tps, df, trait1, trait2, fit_col)
  draw_df <- layers$df
  z_col <- layers$fit_col

  p <- ggplot2::ggplot() +
    ggplot2::geom_contour_filled(
      data = draw_df,
      ggplot2::aes(
        x = .data[[trait1]],
        y = .data[[trait2]],
        z = .data[[z_col]]
      ),
      bins = bins
    ) +
    ggplot2::geom_contour(
      data = draw_df,
      ggplot2::aes(
        x = .data[[trait1]],
        y = .data[[trait2]],
        z = .data[[z_col]]
      ),
      color = "black",
      alpha = 0.3,
      linewidth = 0.3
    ) +
    layers$cover +
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
            values = c("0" = "#D55E00", "1" = "#0072B2"),
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
          y = .data[[trait2]],
          shape = "Optimum"
        ),
        color = "gold",
        size = 4
      )
  }

  if (show_groups) {
    p <- p + .group_layers(tps, trait1, trait2)
  }
  p <- p + .mark_key()

  p <- p + ggplot2::labs(
    x = trait1,
    y = trait2,
    title = "Correlated Fitness Surface",
    subtitle = "Individual-level fitness",
    ...
  )

  return(p)
}
