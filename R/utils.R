#' @keywords internal
#' @noRd
#'
#' Coerce a data.frame / data.table with coordinate columns to an sf POINT object.
#'
#' Accepted column combinations (checked in order):
#'   1. `easting` + `northing`  → EPSG:27700 (OSGB / BNG)
#'   2. `long` + `lat`          → EPSG:4326  (WGS84)
#'
#' If `x` is already an sf object it is returned as-is (after a CRS sanity
#' check when `crs` is supplied).
#'
#' @param x  sf object, data.frame, or data.table with coordinate columns.
#' @param crs  Integer EPSG code to transform the result to.  When `NULL`
#'   (default) the native CRS is kept.
#' @return An sf POINT object.
as_sf_points <- function(x, crs = NULL) {
  if (inherits(x, "sf")) {
    pts <- x
  } else if (all(c("easting", "northing") %in% names(x))) {
    pts <- sf::st_as_sf(
      as.data.frame(x),
      coords = c("easting", "northing"),
      crs    = 27700
    )
  } else if (all(c("long", "lat") %in% names(x))) {
    pts <- sf::st_as_sf(
      as.data.frame(x),
      coords = c("long", "lat"),
      crs    = 4326
    )
  } else {
    stop(
      "coords must be an sf object or a data.frame with either ",
      "'easting'/'northing' (EPSG:27700) or 'long'/'lat' (EPSG:4326) columns.",
      call. = FALSE
    )
  }

  if (!is.null(crs)) {
    pts <- sf::st_transform(pts, crs = crs)
  }

  pts
}


#' @keywords internal
#' @noRd
#'
#' Warn if the CRS of an sf object does not match the expected EPSG code.
#'
#' @param x        sf object.
#' @param expected Integer EPSG code (default 27700 for OSGB).
#' @param name     Name used in the warning message (e.g. "catchment").
check_crs <- function(x, expected = 27700L, name = deparse(substitute(x))) {
  actual <- sf::st_crs(x)$epsg
  if (!isTRUE(actual == expected)) {
    warning(
      sprintf(
        "'%s' has CRS EPSG:%s; expected EPSG:%s (OSGB / British National Grid). ",
        name, actual %||% "unknown", expected
      ),
      "Results may be incorrect if units differ.",
      call. = FALSE
    )
  }
  invisible(x)
}


# Minimal null-coalescing operator used internally.
`%||%` <- function(a, b) if (!is.null(a)) a else b
