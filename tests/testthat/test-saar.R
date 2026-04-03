test_that("saar_adjusted_weights returns three weight columns", {
  thiessen_w <- data.frame(
    notation        = c("G1", "G2", "G3"),
    thiessen_weight = c(0.4, 0.35, 0.25),
    stringsAsFactors = FALSE
  )
  gauge_saar <- data.frame(
    notation = c("G1", "G2", "G3"),
    saar     = c(1100, 1200, 1050),
    stringsAsFactors = FALSE
  )

  result <- saar_adjusted_weights(thiessen_w, gauge_saar, catchment_saar = 1120)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3L)
  expect_true(all(c("notation", "thiessen_weight", "saar_adjusted",
                    "saar_adjusted_rescaled") %in% names(result)))
})

test_that("saar_adjusted_weights: both adjusted columns sum to 1", {
  thiessen_w <- data.frame(
    notation        = c("G1", "G2", "G3"),
    thiessen_weight = c(0.4, 0.35, 0.25)
  )
  gauge_saar <- data.frame(
    notation = c("G1", "G2", "G3"),
    saar     = c(1100, 1200, 1050)
  )
  result <- saar_adjusted_weights(thiessen_w, gauge_saar, catchment_saar = 1120)

  expect_equal(sum(result$saar_adjusted),          1, tolerance = 1e-6)
  expect_equal(sum(result$saar_adjusted_rescaled), 1, tolerance = 1e-9)
})

test_that("saar_adjusted_weights errors on missing columns", {
  expect_error(
    saar_adjusted_weights(
      data.frame(notation = "G1", weight = 1),
      data.frame(notation = "G1", saar = 1000),
      catchment_saar = 1000
    ),
    "thiessen_w is missing columns"
  )
})

test_that("saar_adjusted_weights: thiessen_weight column is unchanged", {
  tw <- data.frame(notation = c("A", "B"), thiessen_weight = c(0.6, 0.4))
  gs <- data.frame(notation = c("A", "B"), saar = c(900, 1100))
  result <- saar_adjusted_weights(tw, gs, catchment_saar = 1000)

  expect_equal(result$thiessen_weight, c(0.6, 0.4))
})
