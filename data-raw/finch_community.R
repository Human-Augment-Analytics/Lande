# Builds inst/extdata/finch_community.csv: the five finch species of El
# Garrapatero from the public data of Beausoleil et al. (2023), one row per bird
# with its mean beak measurements and the number of later years it was seen
# again (their fitness measure). Source file kept in validation/data.
# Run from the package root.
load("validation/data/beausoleil2023_bird.data.RData")
d <- bird.data[complete.cases(bird.data[, c("avg.mbl", "avg.mbd", "avg.mbw", "mxcpois", "sp2")]), ]
out <- data.frame(
  band = d$BANDFINAL,
  species = d$sp2,
  beak_length = d$avg.mbl,
  beak_depth = d$avg.mbd,
  beak_width = d$avg.mbw,
  recaptures = as.integer(d$mxcpois)
)
rownames(out) <- NULL
write.csv(out, "inst/extdata/finch_community.csv", row.names = FALSE)
