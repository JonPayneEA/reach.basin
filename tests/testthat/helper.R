# Shared test fixtures for reach.basin tests
# Loaded automatically by testthat before each test file.

library(sf)

# A small square catchment in BNG (EPSG:27700)
# 10 km × 10 km centred around Easting 350000, Northing 280000
make_test_catchment <- function() {
  coords <- matrix(c(
    345000, 275000,
    355000, 275000,
    355000, 285000,
    345000, 285000,
    345000, 275000
  ), ncol = 2, byrow = TRUE)
  poly <- st_sfc(st_polygon(list(coords)), crs = 27700)
  st_as_sf(data.frame(id = 1L), geometry = poly)
}

# Three gauges well inside the test catchment
make_test_gauges <- function() {
  data.frame(
    notation = c("G1", "G2", "G3"),
    label    = c("Gauge 1", "Gauge 2", "Gauge 3"),
    easting  = c(347000, 351000, 353000),
    northing = c(277000, 281000, 283000),
    stringsAsFactors = FALSE
  )
}

# Minimal hypsometric data
make_test_hyps <- function() {
  data.frame(
    elevation     = c(50, 100, 200, 300, 400),
    area_fraction = c(1.0, 0.85, 0.55, 0.25, 0.0)
  )
}
