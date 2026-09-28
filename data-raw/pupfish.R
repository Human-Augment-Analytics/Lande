# Builds data/crescent_pond_pupfish.rda and data/little_lake_pupfish.rda from
# the size-corrected files of Martin (2016) on Dryad (doi:10.5061/dryad.n3mj3,
# CC0), kept as inst/extdata/crescent_pond_pupfish.csv and
# little_lake_pupfish.csv. Their first column, X, is the row of each fish in
# his raw data (validation/data/rawdata.csv) and is dropped.
# Run from the package root.
for (name in c("crescent_pond_pupfish", "little_lake_pupfish")) {
  d <- read.csv(file.path("inst/extdata", paste0(name, ".csv")))
  d$X <- NULL
  assign(name, d)
  save(list = name, file = file.path("data", paste0(name, ".rda")), compress = "xz")
}
