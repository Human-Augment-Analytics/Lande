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
#' than the dip between them, which \code{valley = TRUE} checks.
#'
#' @param surface Output of \code{correlated_fitness_surface()} fitted with
#'   \code{method = "gam"}.
#' @param from,to A group name, standing for that group's highest point in
#'   \code{surface$groups}, or a numeric vector of the two trait values.
#' @param valley Logical; also find the lowest fitted fitness on the straight
#'   line between the two points and compare each end with it. Default is
#'   \code{FALSE}.
#' @param n_path Number of points along that line. Default is 50.
#'
#' @details Differences are on the scale of the link, log fitness for counts
#'   and log odds for survival. The standard error comes from the covariance
#'   of the model's coefficients, corrected for smoothing parameter
#'   uncertainty when the fit was by REML or ML. The valley is picked as the
#'   lowest fitted point on the line, so its comparisons describe the fitted
#'   surface and are not planned tests.
#'
#' @return A data frame with one row per comparison: the fitted fitness at the
#'   two points (\code{fit_a}, \code{fit_b}), their difference on the link
#'   scale, its standard error, \code{z} and the two-sided \code{p_value}.
#'   With \code{valley = TRUE} the attribute \code{"valley"} holds the trait
#'   values of the lowest point on the line.
#' @export
#'
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#' peak_difference(surf, c(-1, -1), c(1, 1))
peak_difference <- function(surface, from, to, valley = FALSE, n_path = 50) {
  stopifnot(is.list(surface), "model" %in% names(surface))
  model <- surface$model
  if (!inherits(model, "gam")) {
    stop("peak_difference() needs a surface fitted with method = \"gam\"")
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
  if (valley) {
    t <- seq(0, 1, length.out = n_path)
    path <- cbind(a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]))
    eta_path <- drop(stats::predict(model, newdata = frame(path), type = "lpmatrix") %*% beta)
    i <- which.min(eta_path)
    if (i == 1L || i == n_path) {
      message("No dip between the two points: fitted fitness falls all the way from one to the other")
    } else {
      low <- path[i, ]
      pts <- rbind(pts, low)
      labels <- c(labels, "valley")
    }
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
  if (!is.null(low)) out <- rbind(out, compare(1, 3), compare(2, 3))
  rownames(out) <- NULL
  if (!is.null(low)) attr(out, "valley") <- stats::setNames(as.numeric(low), tr)
  out
}
