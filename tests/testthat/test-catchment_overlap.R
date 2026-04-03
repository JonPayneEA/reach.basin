test_that("catchment_overlap_matrix returns list with matrix and flagged", {
  # Two overlapping squares
  sq1 <- sf::st_polygon(list(matrix(c(0,0, 4,0, 4,4, 0,4, 0,0), ncol=2, byrow=TRUE)))
  sq2 <- sf::st_polygon(list(matrix(c(2,0, 6,0, 6,4, 2,4, 2,0), ncol=2, byrow=TRUE)))

  reg <- sf::st_as_sf(
    data.frame(notation = c("A", "B")),
    geometry = sf::st_sfc(sq1, sq2, crs = 27700)
  )

  result <- catchment_overlap_matrix(reg, threshold = 0.4)

  expect_type(result, "list")
  expect_true(all(c("matrix", "flagged") %in% names(result)))
  expect_equal(dim(result$matrix), c(2L, 2L))
  expect_equal(diag(result$matrix), c(1, 1))
})

test_that("catchment_overlap_matrix flags overlapping pair", {
  sq1 <- sf::st_polygon(list(matrix(c(0,0, 4,0, 4,4, 0,4, 0,0), ncol=2, byrow=TRUE)))
  sq2 <- sf::st_polygon(list(matrix(c(2,0, 6,0, 6,4, 2,4, 2,0), ncol=2, byrow=TRUE)))

  reg <- sf::st_as_sf(
    data.frame(notation = c("A", "B")),
    geometry = sf::st_sfc(sq1, sq2, crs = 27700)
  )

  result <- catchment_overlap_matrix(reg, threshold = 0.4)
  expect_equal(nrow(result$flagged), 1L)
})

test_that("catchment_overlap_matrix: non-overlapping catchments produce no flags", {
  sq1 <- sf::st_polygon(list(matrix(c(0,0, 2,0, 2,2, 0,2, 0,0), ncol=2, byrow=TRUE)))
  sq2 <- sf::st_polygon(list(matrix(c(5,0, 7,0, 7,2, 5,2, 5,0), ncol=2, byrow=TRUE)))

  reg <- sf::st_as_sf(
    data.frame(notation = c("A", "B")),
    geometry = sf::st_sfc(sq1, sq2, crs = 27700)
  )

  result <- catchment_overlap_matrix(reg, threshold = 0.1)
  expect_equal(nrow(result$flagged), 0L)
})
