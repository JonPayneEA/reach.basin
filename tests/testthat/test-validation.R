test_that("compare_flood_extent returns correct structure", {
  skip_if_not_installed("terra")

  # Create a 10×10 raster (1 km resolution) with some flooded cells
  r <- terra::rast(
    xmin = 345000, xmax = 355000,
    ymin = 275000, ymax = 285000,
    resolution = 1000, crs = "EPSG:27700"
  )
  # Set the central 4×4 cells to depth 0.5 m, rest to 0
  terra::values(r) <- 0
  terra::values(r)[c(34, 35, 44, 45)] <- 0.5

  # Observed extent covers the central 3×3 area (slightly offset)
  obs_coords <- matrix(c(
    348000, 278000,
    352000, 278000,
    352000, 282000,
    348000, 282000,
    348000, 278000
  ), ncol = 2, byrow = TRUE)
  observed <- sf::st_as_sf(
    data.frame(id = 1L),
    geometry = sf::st_sfc(sf::st_polygon(list(obs_coords)), crs = 27700)
  )

  result <- compare_flood_extent(r, observed, threshold = 0)

  expect_type(result, "list")
  expect_true(all(c("csi", "hit_rate", "false_alarm_ratio", "areas", "geometry") %in% names(result)))
  expect_true(result$csi >= 0 && result$csi <= 1)
})

test_that("compare_flood_extent: perfect match gives CSI = 1", {
  skip_if_not_installed("terra")

  r <- terra::rast(
    xmin = 0, xmax = 10, ymin = 0, ymax = 10,
    resolution = 1, crs = "EPSG:27700"
  )
  terra::values(r) <- 1  # all flooded

  # Observed exactly matches the raster extent
  obs <- sf::st_as_sf(
    data.frame(id = 1L),
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE))),
      crs = 27700
    )
  )

  result <- compare_flood_extent(r, obs, threshold = 0)
  expect_equal(result$csi, 1, tolerance = 0.01)
})

test_that("compare_flood_extent: no overlap gives CSI = 0", {
  skip_if_not_installed("terra")

  r <- terra::rast(
    xmin = 0, xmax = 5, ymin = 0, ymax = 5,
    resolution = 1, crs = "EPSG:27700"
  )
  terra::values(r) <- 1  # modelled flooded area 0–5, 0–5

  # Observed is completely separate (10–15, 10–15)
  obs <- sf::st_as_sf(
    data.frame(id = 1L),
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(c(10,10, 15,10, 15,15, 10,15, 10,10), ncol=2, byrow=TRUE))),
      crs = 27700
    )
  )

  result <- compare_flood_extent(r, obs, threshold = 0)
  expect_equal(result$csi, 0, tolerance = 0.01)
})
