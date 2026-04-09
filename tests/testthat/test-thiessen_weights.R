test_that("thiessen_weights returns one row per gauge", {
  catchment <- make_test_catchment()
  gauges    <- make_test_gauges()

  result <- thiessen_weights(gauges, catchment)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(gauges))
})

test_that("thiessen_weights sums to 1 (within floating-point tolerance)", {
  catchment <- make_test_catchment()
  gauges    <- make_test_gauges()

  result <- thiessen_weights(gauges, catchment)

  expect_equal(sum(result$thiessen_weight), 1, tolerance = 1e-6)
})

test_that("thiessen_weights errors when gauge is outside catchment bbox", {
  catchment <- make_test_catchment()

  gauges_outside <- data.frame(
    notation = c("G1", "G_OUT"),
    easting  = c(347000, 999999),   # G_OUT is far outside
    northing = c(277000, 999999),
    stringsAsFactors = FALSE
  )

  expect_error(
    thiessen_weights(gauges_outside, catchment),
    "outside the catchment bounding box"
  )
})

test_that("thiessen_weights accepts sf POINT input", {
  catchment <- make_test_catchment()
  gauges    <- make_test_gauges()

  pts <- sf::st_as_sf(gauges, coords = c("easting", "northing"), crs = 27700)
  result <- thiessen_weights(pts, catchment)

  expect_equal(nrow(result), 3L)
  expect_equal(sum(result$thiessen_weight), 1, tolerance = 1e-6)
})

test_that("thiessen_weights accepts long/lat input and reprojects", {
  catchment <- make_test_catchment()
  gauges_ll <- make_test_gauges()

  # Convert BNG to WGS84
  pts_bng <- sf::st_as_sf(gauges_ll, coords = c("easting", "northing"), crs = 27700)
  pts_wgs <- sf::st_transform(pts_bng, crs = 4326)
  coords_wgs <- sf::st_coordinates(pts_wgs)

  gauges_latlong <- data.frame(
    notation = gauges_ll$notation,
    long     = coords_wgs[, "X"],
    lat      = coords_wgs[, "Y"],
    stringsAsFactors = FALSE
  )

  result <- thiessen_weights(gauges_latlong, catchment)
  expect_equal(sum(result$thiessen_weight), 1, tolerance = 1e-6)
})
