# ======================================================
# plot_adaptive_landscape.R
# Visualize Adaptive Landscape (population-level)
#
# IMPORTANT CONCEPT:
# - Adaptive Landscape: mean fitness (Wbar ~ zbar1, zbar2)
# - This is DIFFERENT from correlated fitness surface (individual-level)
#
# KEY PRINCIPLE:
# - Traits should already be standardized
# - DO NOT standardize again within this function
# ======================================================

#' Plot Adaptive Landscape
#'
#' Two traits give a contour map of mean fitness against the two population
#' means. One trait gives a curve of mean fitness against the population mean,
#' with the individual fitness function drawn alongside it for comparison.
#'
#' @param landscape Output object of class \code{"adaptive_landscape"}.
#' @param trait_cols One or two trait column names, matching the landscape.
#' @param original_data Optional data frame of original data points. Default is \code{NULL}.
#' @param group_col Optional character string specifying a grouping variable for labels.
#' @param bins Integer specifying the number of contour bins. Default is 12.
#' @param show_optimum Logical indicating whether to display the optimum point. Default is \code{TRUE}.
#' @param show_actual_means Logical indicating whether to display actual population means. Default is \code{TRUE}.
#' @param show_individual Logical; for a single trait, also draw the
#'   individual fitness function (dashed). Default is \code{TRUE}.
#' @param point_alpha Numeric value for point transparency. Default is 0.8.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#' @return A \code{ggplot} object representing the adaptive landscape.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#' land <- adaptive_landscape(prep, surf$model, c("total_length", "weight"),
#'                            simulation_n = 100, grid_n = 15)
#' plot_adaptive_landscape(land, c("total_length", "weight"))
#' @export
plot_adaptive_landscape <- function(
  landscape,
  trait_cols,
  original_data = NULL,
  group_col = NULL,
  bins = 12,
  show_optimum = TRUE,
  show_actual_means = TRUE,
  show_individual = TRUE,
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

    if (length(trait_cols) == 1L) {
        return(.plot_landscape_curve(
            landscape, trait_cols, group_col, show_optimum,
            show_actual_means, show_individual, ...
        ))
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


#' @noRd
# One-trait landscape: mean fitness against the population mean, with the
# individual fitness function alongside so the smoothing by within-population
# variance is visible.
.plot_landscape_curve <- function(landscape, trait, group_col, show_optimum,
                                  show_actual_means, show_individual, ...) {
    df <- landscape$grid
    curves <- data.frame(
        x = df[[trait]], y = df$.mean_fit,
        curve = "Adaptive landscape (mean fitness)"
    )
    if (show_individual && ".ind_fit" %in% names(df)) {
        curves <- rbind(curves, data.frame(
            x = df[[trait]], y = df$.ind_fit,
            curve = "Fitness function (individual)"
        ))
    }
    curves$curve <- factor(curves$curve, levels = unique(curves$curve))

    p <- ggplot2::ggplot(curves, ggplot2::aes(x = .data$x, y = .data$y, linetype = .data$curve)) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
        ggplot2::labs(
            x = trait,
            y = "Mean fitness",
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
            legend.position = "bottom"
        )

    if (show_actual_means && !is.null(landscape$actual_population_means)) {
        means <- landscape$actual_population_means
        p <- p + ggplot2::geom_vline(
            xintercept = means[[trait]], color = "red", alpha = 0.6, linewidth = 0.5
        )
        if (!is.null(group_col) && group_col %in% names(means)) {
            p <- p + ggplot2::annotate(
                "text",
                x = means[[trait]], y = max(curves$y, na.rm = TRUE),
                label = means[[group_col]], angle = 90, vjust = -0.4, hjust = 1,
                size = 3, color = "red"
            )
        }
    }

    if (show_optimum && !is.null(landscape$optimum)) {
        opt <- landscape$optimum
        p <- p +
            ggplot2::annotate(
                "point",
                x = opt[[trait]], y = opt$.mean_fit,
                color = "gold", size = 5, shape = 18
            ) +
            ggplot2::annotate(
                "text",
                x = opt[[trait]], y = opt$.mean_fit,
                label = "Optimum", vjust = -1, size = 3, color = "gray30"
            )
    }

    p
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
#' @examples
#' if (requireNamespace("fields", quietly = TRUE)) {
#'   prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#'   surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#'   land <- adaptive_landscape(prep, surf$model, c("total_length", "weight"),
#'                              simulation_n = 100, grid_n = 15)
#'   plot_adaptive_landscape_3d(land, c("total_length", "weight"))
#' }
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

    if (length(trait_cols) != 2L) {
        stop("The 3D view needs two traits; a single-trait landscape is a curve, ",
             "use plot_adaptive_landscape()")
    }

    if (!requireNamespace("fields", quietly = TRUE)) {
        stop("Package 'fields' is required for 3D landscape plots. ",
             "Install it with install.packages('fields').")
    }

    df <- landscape$grid

    if (!all(trait_cols %in% names(df))) {
        stop("Trait columns not found in landscape$grid")
    }

    x <- df[[trait_cols[1]]]
    y <- df[[trait_cols[2]]]
    z <- df$.mean_fit

    # adaptive_landscape() evaluates mean fitness on a regular grid, so the
    # surface can be reshaped directly into a matrix. Interpolating it again
    # (as earlier versions did with akima) either overshoots at the edges or
    # leaves the border undefined, and is not needed. `grid_n` is kept for
    # backwards compatibility; the landscape's own resolution is used.
    xu <- sort(unique(x))
    yu <- sort(unique(y))
    zmat <- matrix(NA_real_, nrow = length(xu), ncol = length(yu))
    zmat[cbind(match(x, xu), match(y, yu))] <- z
    if (anyNA(zmat)) {
        stop("landscape$grid is not a complete regular grid; cannot draw a 3D surface")
    }

    if (is.null(color_palette)) {
        if (requireNamespace("viridis", quietly = TRUE)) {
            color_palette <- viridis::plasma(100)
        } else {
            color_palette <- grDevices::heat.colors(100)
        }
    }

    fields::drape.plot(
        x = xu,
        y = yu,
        z = zmat,
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
