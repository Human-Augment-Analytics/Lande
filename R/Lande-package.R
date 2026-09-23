#' @section Data:
#' \code{bumpus} has its own help page. The files in \code{inst/extdata} come
#' from the data archived with three papers:
#' \itemize{
#'   \item \code{crescent_pond_pupfish.csv} and \code{little_lake_pupfish.csv}:
#'     Martin, C. H. (2016) Context dependence in complex adaptive landscapes:
#'     frequency and trait-dependent selection surfaces within an adaptive
#'     radiation of Caribbean pupfishes. Evolution 70, 1265-1282. Every
#'     enclosure treatment is there; his analyses use the high-density fish,
#'     \code{density == "H"}, as do the examples and the app.
#'   \item \code{finch_yearly.csv}: Beausoleil, M.-O. et al. (2019) Temporally
#'     varying disruptive selection in the medium ground finch (Geospiza
#'     fortis). Proceedings of the Royal Society B 286, 20192290. Built by
#'     \code{data-raw/finch_yearly.R}: one row per bird per year, survival to the
#'     next year, beak size as the first principal component of the beak medians.
#'   \item \code{finch_community.csv}: Beausoleil, M.-O. et al. (2023) The
#'     fitness landscape of a community of Darwin's finches. Evolution 77,
#'     2533-2546. Dataset \doi{10.5683/SP3/0YIWSE}; the source file is the one in
#'     the authors' code repository (GPL-3). Built by
#'     \code{data-raw/finch_community.R}: one row per bird with its mean beak
#'     measurements and the number of later years it was seen.
#' }
#' @keywords internal
"_PACKAGE"

#' @importFrom dplyr %>% .data
#' @importFrom stats aggregate as.formula binomial coef complete.cases cor cov gaussian glm lm median na.omit poisson predict sd var
#' @importFrom utils combn
#' @importFrom grDevices colorRampPalette heat.colors
NULL
