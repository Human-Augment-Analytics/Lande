# ======================================================
# plot_adaptive_landscape.R
# Visualize Adaptive Landscape (population-level)
#
# IMPORTANT CONCEPT:
# - Adaptive Landscape: mean fitness (W̄ ~ z̄₁, z̄₂)
# - This is DIFFERENT from correlated fitness surface (individual-level)
#
# KEY PRINCIPLE:
# - Traits should already be standardized
# - DO NOT standardize again within this function
# ======================================================

#' Plot Adaptive Landscape (2D Contour)
#'
#' @param landscape Output object of class \code{"adaptive_landscape"}.
#' @param trait_cols Character vector of length 2 specifying the trait column names.
#' @param original_data Optional data frame of original data points. Default is \code{NULL}.
#' @param group_col Optional character string specifying a grouping variable for labels.
#' @param bins Integer specifying the number of contour bins. Default is 12.
#' @param show_optimum Logical indicating whether to display the optimum point. Default is \code{TRUE}.
#' @param show_actual_means Logical indicating whether to display actual population means. Default is \code{TRUE}.
#' @param point_alpha Numeric value for point transparency. Default is 0.8.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#' @return A \code{ggplot} object representing the adaptive landscape.
#' @export
plot_adaptive_landscape <- function(
  landscape,
  trait_cols,
  original_data = NULL,
  group_col = NULL,
  bins = 12,
  show_optimum = TRUE,
  show_actual_means = TRUE,
  point_alpha = 0.8,
  ...
) {
    # Input validation
    if (!inherits(landscape, "adaptive_landscape")) {
        warning("Object is not of class 'adaptive_landscape'")
    }

    df <- landscape$grid

    # Check if trait columns exist
    if (!all(trait_cols %in% names(df))) {
        stop(
            "Trait columns not found in landscape$grid. ",
            "Expected: ", paste(trait_cols, collapse = ", "), "\n",
            "Found: ", paste(names(df), collapse = ", ")
        )
    }

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = .data[[trait_cols[1]]],
            y = .data[[trait_cols[2]]],
            z = .data[[".mean_fit"]]
        )
    ) +
        # Filled contours
        ggplot2::geom_contour_filled(bins = bins) +
        # Contour lines
        ggplot2::geom_contour(color = "black", alpha = 0.3, linewidth = 0.3) +
        # Color scale (continuous)
        ggplot2::scale_fill_viridis_d(name = "Mean Fitness", option = "plasma") +
        # Labels
        ggplot2::labs(
            x = trait_cols[1],
            y = trait_cols[2],
            title = "Adaptive Landscape",
            subtitle = "Population-level mean fitness",
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
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40"),
            aspect.ratio = 0.8
        )

    # Add actual population means
    if (show_actual_means && !is.null(landscape$actual_population_means)) {
        actual_df <- landscape$actual_population_means

        actual_df <- actual_df[, c(trait_cols, group_col), drop = FALSE]

        p <- p +
            ggplot2::geom_point(
                data = actual_df,
                ggplot2::aes(
                    x = .data[[trait_cols[1]]],
                    y = .data[[trait_cols[2]]]
                ),
                color = "red",
                size = 4,
                shape = 19,
                alpha = point_alpha,
                inherit.aes = FALSE
            )

        # Add labels if group_col provided
        if (!is.null(group_col) && group_col %in% names(actual_df)) {
            if (requireNamespace("ggrepel", quietly = TRUE)) {
                p <- p +
                    ggrepel::geom_text_repel(
                        data = actual_df,
                        ggplot2::aes(
                            x = .data[[trait_cols[1]]],
                            y = .data[[trait_cols[2]]],
                            label = .data[[group_col]]
                        ),
                        size = 3,
                        box.padding = 0.5,
                        point.padding = 0.3,
                        inherit.aes = FALSE
                    )
            } else {
                p <- p +
                    ggplot2::geom_text(
                        data = actual_df,
                        ggplot2::aes(
                            x = .data[[trait_cols[1]]],
                            y = .data[[trait_cols[2]]],
                            label = .data[[group_col]]
                        ),
                        size = 3,
                        vjust = -1,
                        hjust = 0.5,
                        inherit.aes = FALSE
                    )
            }
        }
    }

    # Add optimum point
    if (show_optimum && !is.null(landscape$optimum)) {
        opt_df <- landscape$optimum

        opt_df <- opt_df[, trait_cols, drop = FALSE]

        p <- p +
            ggplot2::geom_point(
                data = opt_df,
                ggplot2::aes(
                    x = .data[[trait_cols[1]]],
                    y = .data[[trait_cols[2]]]
                ),
                color = "gold",
                size = 5,
                shape = 18,
                inherit.aes = FALSE
            ) +
            ggplot2::annotate(
                "text",
                x = opt_df[[trait_cols[1]]],
                y = opt_df[[trait_cols[2]]],
                label = "Optimum",
                vjust = -1,
                size = 3,
                color = "gray30"
            )
    }

    return(p)
}


# ======================================================
# plot_adaptive_landscape_3d
# 3D perspective view of adaptive landscape
# ======================================================

#' Plot Adaptive Landscape (3D Perspective)
#'
#' @param landscape Output object of class \code{"adaptive_landscape"}.
#' @param trait_cols Character vector of length 2 specifying the trait column names.
#' @param theta Numeric azimuthal viewing angle. Default is -30.
#' @param phi Numeric colatitude viewing angle. Default is 30.
#' @param grid_n Integer specifying the resolution of the interpolation grid. Default is 200.
#' @param color_palette Optional vector of colors for the surface. Defaults to viridis plasma.
#' @param ... Additional arguments passed to \code{fields::drape.plot()}.
#'
#' @return A 3D plot produced by \code{fields::drape.plot()}.
#' @export
plot_adaptive_landscape_3d <- function(
  landscape,
  trait_cols,
  theta = -30,
  phi = 30,
  grid_n = 200,
  color_palette = NULL,
  ...
) {
    # Input validation
    if (!inherits(landscape, "adaptive_landscape")) {
        warning("Object is not of class 'adaptive_landscape'")
    }

    df <- landscape$grid

    if (!all(trait_cols %in% names(df))) {
        stop("Trait columns not found in landscape$grid")
    }

    x <- df[[trait_cols[1]]]
    y <- df[[trait_cols[2]]]
    z <- df$.mean_fit

    # Add small jitter to avoid collinear warnings
    set.seed(42)
    x <- x + rnorm(length(x), 0, 1e-8)
    y <- y + rnorm(length(y), 0, 1e-8)

    interp <- akima::interp(
        x = x,
        y = y,
        z = z,
        xo = seq(min(x), max(x), length = grid_n),
        yo = seq(min(y), max(y), length = grid_n),
        linear = FALSE, # Spline interpolation
        extrap = TRUE,
        duplicate = "mean"
    )

    if (is.null(color_palette)) {
        if (requireNamespace("viridis", quietly = TRUE)) {
            color_palette <- viridis::plasma(100)
        } else {
            color_palette <- grDevices::heat.colors(100)
        }
    }

    fields::drape.plot(
        x = interp$x,
        y = interp$y,
        z = interp$z,
        theta = theta,
        phi = phi,
        xlab = trait_cols[1],
        ylab = trait_cols[2],
        zlab = "Mean Fitness",
        main = "Adaptive Landscape (3D View)",
        col = color_palette,
        border = NA,
        shade = 0.5,
        ...
    )
}
