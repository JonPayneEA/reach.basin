test_that("load_catchment errors on missing file", {
  expect_error(
    load_catchment("/nonexistent/path/catchment.shp"),
    "File not found"
  )
})

test_that("load_catchment returns a single-row sf with correct CRS", {
  tmp <- tempfile(fileext = ".gpkg")
  on.exit(unlink(tmp))

  catchment <- make_test_catchment()
  sf::st_write(catchment, tmp, quiet = TRUE)

  result <- load_catchment(tmp)

  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 1L)
  expect_equal(sf::st_crs(result)$epsg, 27700L)
})

test_that("load_catchment warns on non-BNG CRS", {
  tmp <- tempfile(fileext = ".gpkg")
  on.exit(unlink(tmp))

  # Write in WGS84
  wgs <- sf::st_transform(make_test_catchment(), crs = 4326)
  sf::st_write(wgs, tmp, quiet = TRUE)

  expect_warning(
    load_catchment(tmp),
    "EPSG:4326"
  )
})

test_that("load_catchment dissolves multi-part polygon to one row", {
  tmp <- tempfile(fileext = ".gpkg")
  on.exit(unlink(tmp))

  # Two separate squares
  sq1 <- sf::st_polygon(list(matrix(c(
    0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE)))
  sq2 <- sf::st_polygon(list(matrix(c(
    2, 0, 3, 0, 3, 1, 2, 1, 2, 0), ncol = 2, byrow = TRUE)))

  multi <- sf::st_as_sf(data.frame(id = 1:2),
                         geometry = sf::st_sfc(sq1, sq2, crs = 27700))
  sf::st_write(multi, tmp, quiet = TRUE)

  # Suppress the dissolve message
  result <- suppressMessages(load_catchment(tmp))

  expect_equal(nrow(result), 1L)
})
