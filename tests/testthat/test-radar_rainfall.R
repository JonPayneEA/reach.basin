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

test_that("areal_reduction_factor is monotonically decreasing with duration (for fixed area)", {
  durs <- c(1, 2, 6, 12, 24, 48, 96)
  arfs <- areal_reduction_factor(area_km2 = 100, duration_hr = durs)
  expect_true(all(diff(arfs) > 0))  # longer duration → higher ARF
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
