# ======================================================
# peak_difference.R
# Difference in fitted fitness between two points of a GAM surface, with its
# standard error from the covariance of the coefficients.
# ======================================================

#' Compare the fitted fitness at two points of a surface
#'
#' Takes two points on a GAM surface, given as group names or trait values,
#' and returns the difference in fitted fitness on the scale of the link,
#' with its standard error. Two peaks are separate only if each is higher
#' than the pass between them, the lowest point on the highest route from one
#' to the other, which \code{valley = TRUE} compares.
#'
#' @param surface Output of \code{correlated_fitness_surface()} fitted with
#'   \code{method = "gam"}.
#' @param from,to A group name, standing for that group's highest point in
#'   \code{surface$groups}, or a numeric vector of the two trait values.
#' @param valley Logical; also find the valley between the two points and
#'   compare each end with it. Default is \code{FALSE}.
#' @param route How the valley is found: \code{"pass"} (the default), the
#'   pass over the kept cells of the surface's grid, or \code{"line"}, the
#'   lowest point on the straight line between the two points, which is how
#'   Beausoleil et al. (2023) measured valley depths. A straight line can
#'   cross a trough that a curving ridge goes round, so it can make two peaks
#'   look more separate than they are.
#' @param n_path Number of points along the straight line with
#'   \code{route = "line"}. Default is 50.
#'
#' @details Differences are on the scale of the link, log fitness for counts
#'   and log odds for survival. The standard error comes from the covariance
#'   of the model's coefficients, corrected for smoothing parameter
#'   uncertainty when the fit was by REML or ML. The points are taken as
#'   fixed, although the peaks and the valley were found on the same fitted
#'   surface, so the comparisons describe the fitted surface and are not
#'   planned tests. The valley is the lowest point on the route, so its
#'   differences are biased upwards and their z values too large, which is why
#'   the valley rows have no p-value. The pass is found cell by cell over the
#'   grid, joining each kept cell to its eight neighbours, from the two kept
#'   cells nearest the points in grid steps; its precision follows
#'   \code{grid_n}, and a point far from any kept cell gets a message.
#'   If the route never drops below the cell a point starts from, or the lower
#'   point is no higher than the pass, there is no dip to compare. The
#'   straight line ignores the mask. For overdispersed counts fit the surface with
#'   \code{count_family = "quasipoisson"} first, since Poisson standard errors
#'   are then too small.
#'
#' @return A data frame with one row per comparison: the fitted fitness at the
#'   two points (\code{fit_a}, \code{fit_b}), their difference on the link
#'   scale, its standard error, \code{z} and the two-sided \code{p_value}
#'   (\code{NA} for the valley rows).
#'   With \code{valley = TRUE} the attribute \code{"valley"} holds the trait
#'   values of the pass or of the lowest point on the line.
#' @export
#'
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#' peak_difference(surf, c(-1, -1), c(1, 1))
peak_difference <- function(surface, from, to, valley = FALSE, route = c("pass", "line"), n_path = 50) {
  route <- match.arg(route)
  stopifnot(is.list(surface), "model" %in% names(surface))
  model <- surface$model
  if (!inherits(model, "gam")) {
    stop("peak_difference() needs a surface fitted with method = \"gam\"")
  }
  if (identical(surface$count_family, "poisson") && isTRUE(surface$dispersion > 1.5)) {
    warning("Counts are overdispersed (dispersion ", round(surface$dispersion, 2),
            "); these Poisson standard errors are too small, and refitting with count_family = \"quasipoisson\" corrects them")
  }
  tr <- surface$trait_cols

  point <- function(p) {
    if (is.character(p) && length(p) == 1L) {
      g <- surface$groups
      if (is.null(g) || !p %in% g$group) stop("No group '", p, "' on this surface")
      xy <- as.numeric(g[g$group == p, paste0("peak_", tr)])
      if (anyNA(xy)) stop("Group '", p, "' has no highest point on the kept surface")
      return(xy)
    }
    if (is.numeric(p) && length(p) == 2L && !anyNA(p)) return(as.numeric(p))
    stop("`from` and `to` must each be a group name or the two trait values")
  }
  label <- function(p) if (is.character(p)) p else paste0("(", paste(signif(p, 3), collapse = ", "), ")")

  a <- point(from)
  b <- point(to)
  pts <- rbind(a, b)
  labels <- c(label(from), label(to))

  # a surface fitted with the group as a fixed effect is compared at the
  # reference level, as it is drawn
  frame <- function(m) {
    d <- stats::setNames(as.data.frame(m), tr)
    if (isTRUE(surface$group_effect)) {
      d[[surface$group_used]] <- .reference_group(model$model[[surface$group_used]])
    }
    d
  }
  beta <- stats::coef(model)
  V <- if (!is.null(model$Vc)) model$Vc else model$Vp
  linkinv <- model$family$linkinv

  low <- NULL
  if (valley && route == "line") {
    t <- seq(0, 1, length.out = n_path)
    path <- cbind(a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]))
    eta_path <- drop(stats::predict(model, newdata = frame(path), type = "lpmatrix") %*% beta)
    i <- which.min(eta_path)
    if (i == 1L || i == n_path) {
      message("No dip between the two points: fitted fitness falls all the way from one to the other")
    } else {
      low <- path[i, ]
    }
  } else if (valley) {
    pass <- .surface_pass(surface$grid, tr, a, b)
    if (!is.null(pass) && any(pass$far)) {
      message("A point lies away from the kept surface; the route starts from the nearest kept cell")
    }
    if (is.null(pass)) {
      message("The two points are not joined within the kept surface; there is no pass to compare")
    } else {
      # no dip if the route never drops below the cell a point starts from, or,
      # judged at the points themselves (which need not sit on the grid), if
      # the lower point is no higher than the pass
      at <- rbind(a, b, as.numeric(surface$grid[pass$cell, tr]))
      eta <- drop(stats::predict(model, newdata = frame(at), type = "lpmatrix") %*% beta)
      if (pass$cell %in% pass$ends || eta[3] >= min(eta[1:2]) - 1e-10) {
        message("No dip between the two points: the lower one is no higher than the pass")
      } else {
        low <- at[3, ]
      }
    }
  }
  if (!is.null(low)) {
    pts <- rbind(pts, low)
    labels <- c(labels, "valley")
  }

  Xp <- stats::predict(model, newdata = frame(pts), type = "lpmatrix")
  fitted <- linkinv(drop(Xp %*% beta))
  compare <- function(i, j) {
    d <- Xp[i, ] - Xp[j, ]
    est <- sum(d * beta)
    se <- sqrt(drop(d %*% V %*% d))
    data.frame(comparison = paste(labels[i], "-", labels[j]), fit_a = fitted[i], fit_b = fitted[j],
               difference = est, se = se, z = est / se, p_value = 2 * stats::pnorm(-abs(est / se)),
               stringsAsFactors = FALSE)
  }
  out <- compare(1, 2)
  if (!is.null(low)) {
    out <- rbind(out, compare(1, 3), compare(2, 3))
    out$p_value[2:3] <- NA_real_
  }
  rownames(out) <- NULL
  if (!is.null(low)) attr(out, "valley") <- stats::setNames(as.numeric(low), tr)
  out
}

#' @noRd
# internal utility: the pass between two points of the kept surface. Kept
# cells are added from the highest down, each joined to the kept neighbours
# in its 3 x 3 block already added, until the cells nearest the two points
# fall in one piece; the cell that joins them is the pass. NULL if they never
# do.
.surface_pass <- function(grid, trait_cols, a, b) {
  gx <- grid[[trait_cols[1]]]
  gy <- grid[[trait_cols[2]]]
  ux <- sort(unique(gx))
  uy <- sort(unique(gy))
  ix <- match(gx, ux)
  iy <- match(gy, uy)
  cell <- matrix(NA_integer_, length(ux), length(uy))
  cell[cbind(ix, iy)] <- seq_along(gx)
  kept <- !is.na(grid$.fit)
  # distances in grid steps, so the two traits count alike whatever their units
  sx <- min(diff(ux))
  sy <- min(diff(uy))
  steps <- function(p, i) sqrt(((gx[i] - p[1]) / sx)^2 + ((gy[i] - p[2]) / sy)^2)
  nearest <- function(p) which.min(ifelse(kept, steps(p, seq_along(gx)), Inf))
  ends <- c(nearest(a), nearest(b))
  far <- c(steps(a, ends[1]), steps(b, ends[2])) > 1.5

  parent <- seq_along(gx)
  root <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  added <- rep(FALSE, length(gx))
  for (i in order(grid$.fit, decreasing = TRUE, na.last = NA)) {
    added[i] <- TRUE
    for (dx in -1:1) for (dy in -1:1) {
      jx <- ix[i] + dx
      jy <- iy[i] + dy
      if (jx < 1 || jy < 1 || jx > length(ux) || jy > length(uy)) next
      j <- cell[jx, jy]
      if (is.na(j) || !added[j]) next
      ri <- root(i)
      rj <- root(j)
      if (ri != rj) parent[ri] <- rj
    }
    if (all(added[ends]) && root(ends[1]) == root(ends[2])) return(list(cell = i, ends = ends, far = far))
  }
  NULL
}
