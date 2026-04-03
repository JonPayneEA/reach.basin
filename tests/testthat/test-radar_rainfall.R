# ---- ARF tests (unchanged) --------------------------------------------------

test_that("areal_reduction_factor returns values in (0, 1]", {
  arf <- areal_reduction_factor(area_km2 = 100, duration_hr = 24)
  expect_length(arf, 1L)
  expect_true(arf > 0 && arf <= 1)
})

test_that("areal_reduction_factor is monotonically decreasing with area", {
  areas <- c(1, 10, 100, 1000, 5000)
  arfs  <- areal_reduction_factor(area_km2 = areas, duration_hr = 24)
  expect_true(all(diff(arfs) < 0))
})

test_that("areal_reduction_factor is monotonically increasing with duration", {
  durs <- c(1, 2, 6, 12, 24, 48, 96)
  arfs <- areal_reduction_factor(area_km2 = 100, duration_hr = durs)
  expect_true(all(diff(arfs) > 0))
})

test_that("areal_reduction_factor known value: 1 km², 24 hr ≈ 0.996", {
  arf <- areal_reduction_factor(area_km2 = 1, duration_hr = 24)
  expect_equal(arf, 0.996, tolerance = 0.001)
})

test_that("areal_reduction_factor known value: 10000 km², 1 hr ≈ 0.604", {
  arf <- areal_reduction_factor(area_km2 = 10000, duration_hr = 1)
  expect_equal(arf, 0.604, tolerance = 0.001)
})

test_that("areal_reduction_factor warns on area > 10000 km²", {
  expect_warning(
    areal_reduction_factor(area_km2 = 15000, duration_hr = 24),
    "10 000 km"
  )
})

test_that("areal_reduction_factor is vectorised over area", {
  arfs <- areal_reduction_factor(area_km2 = c(10, 100, 500), duration_hr = 12)
  expect_length(arfs, 3L)
  expect_true(all(arfs > 0 & arfs <= 1))
})


# ---- QMED tests (also in this file per original layout) --------------------

test_that("estimate_qmed returns positive numeric", {
  qmed <- estimate_qmed(list(area_km2 = 150, saar = 1050, farl = 0.98, bfihost = 0.44))
  expect_length(qmed, 1L)
  expect_true(qmed > 0)
})

test_that("estimate_qmed errors on missing descriptors", {
  expect_error(
    estimate_qmed(list(area_km2 = 150, saar = 1050)),
    "missing elements"
  )
})

test_that("estimate_qmed: larger area gives larger QMED", {
  q_small <- estimate_qmed(list(area_km2 = 50,  saar = 1000, farl = 1, bfihost = 0.4))
  q_large <- estimate_qmed(list(area_km2 = 500, saar = 1000, farl = 1, bfihost = 0.4))
  expect_true(q_large > q_small)
})


# ---- Helpers for NetCDF-based tests -----------------------------------------

# Create a tiny 3×3 SpatRaster with n_layers layers in EPSG:27700 and write it
# to a temporary NetCDF file.  Returns the file path.
make_test_nc <- function(n_layers, var_name = "rainfall_amount",
                          xmin = 345000, xmax = 355000,
                          ymin = 275000, ymax = 285000,
                          seed = 42L) {
  skip_if_not_installed("terra")

  r <- terra::rast(
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    resolution = 5000, crs = "EPSG:27700",
    nlyrs = n_layers
  )

  set.seed(seed)
  terra::values(r) <- runif(terra::ncell(r) * n_layers, 0, 5)

  # Assign simple time metadata (hourly from a fixed origin)
  origin <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  terra::time(r) <- origin + (seq_len(n_layers) - 1L) * 3600

  tmp <- tempfile(fileext = ".nc")
  terra::writeCDF(r, tmp, varname = var_name, overwrite = TRUE)
  tmp
}


# ---- extract_radar_rainfall tests -------------------------------------------

test_that("extract_radar_rainfall: backward-compatible SpatRaster input", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  r <- terra::rast(
    xmin = 345000, xmax = 355000, ymin = 275000, ymax = 285000,
    resolution = 5000, crs = "EPSG:27700", nlyrs = 4L
  )
  terra::values(r) <- 1
  terra::time(r) <- as.POSIXct("2024-01-01", tz = "UTC") + (0:3) * 900

  result <- extract_radar_rainfall(catchment, r)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4L)
  expect_true(all(c("timestep", "datetime", "areal_rainfall_mm") %in% names(result)))
})

test_that("extract_radar_rainfall: single NetCDF file path", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  f <- make_test_nc(n_layers = 4L)
  on.exit(unlink(f))

  result <- extract_radar_rainfall(catchment, f)

  expect_equal(nrow(result), 4L)
  expect_true(all(c("timestep", "datetime", "areal_rainfall_mm") %in% names(result)))
  expect_true(all(!is.na(result$areal_rainfall_mm)))
})

test_that("extract_radar_rainfall: vector of two files concatenated correctly", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  f1 <- make_test_nc(n_layers = 4L, seed = 1L)
  f2 <- make_test_nc(n_layers = 4L, seed = 2L)
  on.exit({ unlink(f1); unlink(f2) })

  result <- extract_radar_rainfall(catchment, c(f1, f2))

  expect_equal(nrow(result), 8L)
  expect_equal(result$timestep, 1:8)
})

test_that("extract_radar_rainfall: directory path loads all .nc files", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  # Write two nc files with date-based names so lex sort = chronological
  f1 <- make_test_nc(n_layers = 3L, seed = 1L)
  f2 <- make_test_nc(n_layers = 3L, seed = 2L)
  file.copy(f1, file.path(tmp_dir, "h19_20240101.nc"))
  file.copy(f2, file.path(tmp_dir, "h19_20240102.nc"))
  unlink(f1); unlink(f2)

  result <- extract_radar_rainfall(catchment, tmp_dir)
  expect_equal(nrow(result), 6L)
})

test_that("extract_radar_rainfall: errors clearly on empty directory", {
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  expect_error(
    extract_radar_rainfall(make_test_catchment(), tmp_dir),
    "No .nc files found"
  )
})

test_that("extract_radar_rainfall: warns when file layer count not multiple of 96", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  # 5 layers is not a multiple of 96
  f <- make_test_nc(n_layers = 5L)
  on.exit(unlink(f))

  expect_warning(
    extract_radar_rainfall(catchment, f),
    "multiple of 96"
  )
})


# ---- extract_moses_pe tests -------------------------------------------------

test_that("extract_moses_pe: single NetCDF file returns correct structure", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  f <- make_test_nc(n_layers = 4L, var_name = "peti")
  on.exit(unlink(f))

  result <- extract_moses_pe(catchment, f, var_name = "peti")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("timestep", "datetime", "areal_pe_mm") %in% names(result)))
  expect_equal(nrow(result), 4L)
  expect_true(all(!is.na(result$areal_pe_mm)))
})

test_that("extract_moses_pe: two files stacked correctly", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  f1 <- make_test_nc(n_layers = 3L, var_name = "peti", seed = 10L)
  f2 <- make_test_nc(n_layers = 3L, var_name = "peti", seed = 11L)
  on.exit({ unlink(f1); unlink(f2) })

  result <- extract_moses_pe(catchment, c(f1, f2), var_name = "peti")
  expect_equal(nrow(result), 6L)
})

test_that("extract_moses_pe: directory path works", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  tmp_dir   <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  f1 <- make_test_nc(n_layers = 2L, var_name = "peti", seed = 20L)
  f2 <- make_test_nc(n_layers = 2L, var_name = "peti", seed = 21L)
  file.copy(f1, file.path(tmp_dir, "moses_pe_20240101.nc"))
  file.copy(f2, file.path(tmp_dir, "moses_pe_20240102.nc"))
  unlink(f1); unlink(f2)

  result <- extract_moses_pe(catchment, tmp_dir, var_name = "peti")
  expect_equal(nrow(result), 4L)
})

test_that("extract_moses_pe: warns when file layer count != 24", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  # 5 layers is not 24
  f <- make_test_nc(n_layers = 5L, var_name = "peti")
  on.exit(unlink(f))

  expect_warning(
    extract_moses_pe(catchment, f, var_name = "peti"),
    "24"
  )
})

test_that("extract_moses_pe: accepts pre-loaded SpatRaster", {
  skip_if_not_installed("terra")
  skip_if_not_installed("exactextractr")

  catchment <- make_test_catchment()
  r <- terra::rast(
    xmin = 345000, xmax = 355000, ymin = 275000, ymax = 285000,
    resolution = 5000, crs = "EPSG:27700", nlyrs = 3L
  )
  terra::values(r) <- 0.5
  result <- extract_moses_pe(catchment, r)

  expect_equal(nrow(result), 3L)
  expect_true("areal_pe_mm" %in% names(result))
})
