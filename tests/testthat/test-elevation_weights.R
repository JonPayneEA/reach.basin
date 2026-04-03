test_that("elevation_weights returns one row per gauge", {
  hyps   <- make_test_hyps()
  gauges <- data.frame(notation = c("G1", "G2"), elevation = c(120, 280))
  result <- elevation_weights(gauges, hyps)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2L)
  expect_true("elevation_weight" %in% names(result))
})

test_that("elevation_weights sums to 1", {
  hyps   <- make_test_hyps()
  gauges <- data.frame(notation = c("G1", "G2", "G3"), elevation = c(100, 200, 350))
  result <- elevation_weights(gauges, hyps)

  expect_equal(sum(result$elevation_weight), 1, tolerance = 1e-6)
})

test_that("elevation_weights errors on missing columns", {
  hyps   <- make_test_hyps()
  gauges <- data.frame(label = "G1", elev = 100)

  expect_error(elevation_weights(gauges, hyps), "missing columns")
})

test_that("gauge_to_catchment_weights dispatches to thiessen", {
  catchment <- make_test_catchment()
  gauges    <- make_test_gauges()
  w <- gauge_to_catchment_weights(gauges, catchment, method = "thiessen")

  expect_type(w, "double")
  expect_equal(sum(w), 1, tolerance = 1e-6)
  expect_named(w)
})

test_that("gauge_to_catchment_weights dispatches to distance_decay", {
  catchment <- make_test_catchment()
  gauges    <- make_test_gauges()
  w <- gauge_to_catchment_weights(gauges, catchment, method = "distance_decay")

  expect_equal(sum(w), 1, tolerance = 1e-6)
})

test_that("gauge_to_catchment_weights errors when elevation missing hypsometric_data", {
  catchment <- make_test_catchment()
  gauges    <- cbind(make_test_gauges(), elevation = c(100, 200, 300))

  expect_error(
    gauge_to_catchment_weights(gauges, catchment, method = "elevation"),
    "hypsometric_data"
  )
})
