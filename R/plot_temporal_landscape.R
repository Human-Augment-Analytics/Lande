# ======================================================
# plot_temporal_landscape.R
# Draw the fitness functions or surfaces of each period side by side, or,
# for one trait, as a heat map of fitness against trait and time.
# ======================================================

#' @noRd
# the highest fitness of each period: a gold diamond when it lies inside the
# data, an open one when it sits at the edge of the range.
# Carries the shape key for the period marks.
.period_optima <- function(s, xcol, ycol) {
  edge <- if ("optimum_edge" %in% names(s)) s$optimum_edge %in% TRUE else rep(FALSE, nrow(s))
  list(
    if (any(!edge)) ggplot2::geom_point(data = s[!edge, , drop = FALSE],
                                        ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], shape = "Optimum"),
                                        colour = "gold", size = 3.5, inherit.aes = FALSE),
    if (any(edge)) ggplot2::geom_point(data = s[edge, , drop = FALSE],
                                       ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], shape = "Edge maximum"),
                                       colour = "black", size = 3, stroke = 0.9, inherit.aes = FALSE),
    ggplot2::scale_shape_manual(name = NULL, values = c("Optimum" = 18, "Edge maximum" = 5, "Period mean" = 21))
  )
}

#' Plot fitness functions and landscapes over time
#'
#' @param tl Output of \code{temporal_landscape()}.
#' @param type \code{"panels"} draws one panel per period: the fitness
#'   function with its band and the individuals for one trait, the fitness
#'   surface for two. \code{"heatmap"}, for one trait only, draws fitness
#'   against trait and period in one panel, with the position of the highest
#'   fitness in each period marked.
#' @param show_points Logical; draw the individuals in each panel.
#' @param show_landscape Logical; for one trait, also draw each period's
#'   adaptive landscape as a dashed curve.
#' @param show_optimum Logical; mark the highest fitted fitness in each period,
#'   as a gold diamond when it lies inside the data and an open one when it
#'   sits at the edge of the range.
#' @param connect Logical; in the heat map, join the highest fitness of
#'   successive periods with a line. Default is \code{FALSE}, since a maximum at
#'   the edge of the range is not a peak.
#' @param bins Contour bins for two-trait panels.
#' @param ncol Number of panel columns.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#'
#' @return A \code{ggplot} object.
#' @examples
#' finch <- read.csv(system.file("extdata", "finch_yearly.csv", package = "RforEvolution"))
#' prep <- prepare_selection_data(finch, "survived", "beak_pc1")
#' years <- temporal_landscape(prep, "survived", "beak_pc1", "year", landscape = FALSE)
#' plot_temporal_landscape(years, ncol = 4)
#' plot_temporal_landscape(years, type = "heatmap")
#' @export
plot_temporal_landscape <- function(
  tl,
  type = c("panels", "heatmap"),
  show_points = TRUE,
  show_landscape = TRUE,
  show_optimum = TRUE,
  connect = FALSE,
  bins = 10,
  ncol = NULL,
  ...
) {
  type <- match.arg(type)
  stopifnot(inherits(tl, "temporal_landscape"))
  one <- length(tl$trait_cols) == 1L
  tr <- tl$trait_cols
  lev <- as.character(tl$times)
  as_time <- function(x) factor(as.character(x), levels = lev)
  theme_plain <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = ggplot2::element_text(color = "black", size = 10),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40"),
      strip.background = ggplot2::element_rect(fill = "grey95")
    )
  s <- tl$summary
  s$time <- as_time(s$time)

  if (type == "heatmap") {
    if (!one) stop("The heat map is for one trait; use type = 'panels' for two")
    h <- tl$heat
    h$time <- as_time(h$time)
    p <- ggplot2::ggplot(h, ggplot2::aes(x = .data[[tr]], y = .data$time, fill = .data$fit)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(na.value = "white", name = "Fitness") +
      ggplot2::labs(x = tr, y = tl$time_col, title = "Fitness function over time",
                    subtitle = "Blank where a period has no individuals", ...) +
      theme_plain
    if (show_optimum) {
      if (connect) {
        p <- p + ggplot2::geom_path(data = s, ggplot2::aes(x = .data[[paste0("optimum_", tr)]], y = .data$time, group = 1),
                                    colour = "white", linewidth = 0.5, inherit.aes = FALSE)
      }
      p <- p + .period_optima(s, paste0("optimum_", tr), "time")
    }
    return(p)
  }

  if (one) {
    g <- tl$grid
    g$time <- as_time(g$time)
    p <- ggplot2::ggplot(g, ggplot2::aes(x = .data[[tr]], y = .data$fit))
    if (show_points) {
      pts <- tl$data
      pts$time <- as_time(pts[[tl$time_col]])
      p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = .data[[tr]], y = .data[[tl$fitness_col]]),
                                   colour = "grey55", alpha = 0.25, size = 1, inherit.aes = FALSE)
    }
    if (all(c("lwr", "upr") %in% names(g))) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lwr, ymax = .data$upr), fill = "steelblue", alpha = 0.25)
    }
    p <- p + ggplot2::geom_line(linewidth = 0.9)
    if (show_landscape && !is.null(tl$landscape_grid)) {
      lg <- tl$landscape_grid
      lg$time <- as_time(lg$time)
      p <- p + ggplot2::geom_line(data = lg, ggplot2::aes(x = .data[[tr]], y = .data$.mean_fit),
                                  linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE)
    }
    if (show_optimum) {
      p <- p + .period_optima(s, paste0("optimum_", tr), "optimum_fit")
    }
    p <- p +
      ggplot2::facet_wrap(~ time, ncol = ncol) +
      ggplot2::labs(x = tr, y = tl$fitness_col, title = "Fitness function by period",
                    subtitle = if (show_landscape && !is.null(tl$landscape_grid)) "Solid: fitness function. Dashed: adaptive landscape" else NULL, ...) +
      theme_plain
    return(p)
  }

  # two traits: one surface per panel, the outside of each period's hull covered
  g <- tl$grid
  g$time <- as_time(g$time)
  z_col <- if (".fit_all" %in% names(g)) ".fit_all" else ".fit"
  if (!is.null(tl$too_far) && ".dist" %in% names(g)) g <- g[!(g$.dist > tl$too_far), , drop = FALSE]
  p <- ggplot2::ggplot(g[!is.na(g[[z_col]]), , drop = FALSE]) +
    ggplot2::geom_contour_filled(ggplot2::aes(x = .data[[tr[1]]], y = .data[[tr[2]]], z = .data[[z_col]]), bins = bins)
  if (!is.null(tl$hulls)) {
    frames <- lapply(names(tl$hulls), function(key) {
      hull <- tl$hulls[[key]]
      if (is.null(hull)) return(NULL)
      gk <- g[g$time == key, , drop = FALSE]
      f <- .outside_frame(hull, tr[1], tr[2], range(gk[[tr[1]]]), range(gk[[tr[2]]]))
      f$time <- key
      f
    })
    frames <- do.call(rbind, frames)
    if (!is.null(frames) && nrow(frames)) {
      frames$time <- as_time(frames$time)
      p <- p + ggplot2::geom_polygon(data = frames, ggplot2::aes(x = .data$x, y = .data$y, group = .data$time),
                                     fill = "white", colour = "white", linewidth = 0.5, inherit.aes = FALSE)
    }
  }
  if (show_points) {
    pts <- tl$data
    pts$time <- as_time(pts[[tl$time_col]])
    p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = .data[[tr[1]]], y = .data[[tr[2]]]),
                                 colour = "grey30", alpha = 0.35, size = 0.8, inherit.aes = FALSE)
  }
  p <- p + ggplot2::geom_point(data = s, ggplot2::aes(x = .data[[paste0("mean_", tr[1])]], y = .data[[paste0("mean_", tr[2])]], shape = "Period mean"),
                               fill = "white", colour = "black", size = 2.6, inherit.aes = FALSE)
  p <- p + if (show_optimum) {
    .period_optima(s, paste0("optimum_", tr[1]), paste0("optimum_", tr[2]))
  } else {
    ggplot2::scale_shape_manual(name = NULL, values = c("Period mean" = 21))
  }
  p +
    ggplot2::facet_wrap(~ time, ncol = ncol) +
    ggplot2::labs(x = tr[1], y = tr[2], fill = "Fitness", title = "Fitness surface by period", ...) +
    theme_plain
}
