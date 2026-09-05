# Builds inst/extdata/finch_yearly.csv: one row per medium ground finch per
# year it was seen at El Garrapatero, with survival to the next year and beak
# size, from the public data of Beausoleil et al. (2019) kept in
# validation/data/bird_data.csv. Beak size is the first principal component of
# the three beak medians, signed so that larger beaks score higher.
# Run from the package root.
raw <- read.csv("validation/data/bird_data.csv")
raw <- raw[raw$Species1 == "fortis", ]
beak <- c("MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth")
ok <- complete.cases(raw[, beak])
pca <- prcomp(raw[ok, beak], scale. = TRUE)
pc1 <- rep(NA_real_, nrow(raw))
pc1[ok] <- pca$x[, 1]
if (cor(pc1[ok], raw$MedianBeakDepth[ok]) < 0) pc1 <- -pc1
raw$beak_pc1 <- pc1
rows <- list()
for (yr in 2004:2010) {
  d <- raw[raw[[paste0("y.", yr)]] == 1, ]
  rows[[as.character(yr)]] <- data.frame(
    band = d$BANDFINAL, year = yr,
    survived = as.integer(d[[paste0("y.", yr + 1)]] == 1),
    beak_pc1 = round(d$beak_pc1, 4),
    beak_length = d$MedianBeakLength, beak_width = d$MedianBeakWidth, beak_depth = d$MedianBeakDepth
  )
}
out <- do.call(rbind, rows)
out <- out[complete.cases(out[, c("survived", "beak_pc1")]), ]
rownames(out) <- NULL
write.csv(out, "inst/extdata/finch_yearly.csv", row.names = FALSE)
