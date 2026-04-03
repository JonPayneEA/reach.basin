#' Load and validate a catchment shapefile
#'
#' Reads a vector file (typically a shapefile or GeoPackage) and validates that
#' it represents a single catchment polygon suitable for use in the reach.basin
#' workflow.  Multi-part (MULTIPOLYGON) geometries are split via
#' [sf::st_cast()] and then dissolved back to a single POLYGON with
#' [sf::st_union()], so shapefiles that were saved as MULTIPOLYGON by GIS tools
#' are handled correctly.
#'
#' A warning is raised (not an error) when the CRS is not EPSG:27700 (OSGB /
#' British National Grid), because all Environment Agency hydrometric data
#' uses that projection.  The function still returns the object so the caller
#' can decide whether to reproject.
#'
#' @param filepath Character. Path to the spatial file to read.  Any format
#'   supported by GDAL / [sf::st_read()] is accepted (`.shp`, `.gpkg`,
#'   `.geojson`, etc.).
#' @param quiet Logical. Passed to [sf::st_read()] to suppress driver/layer
#'   messages.  Default `TRUE`.
#'
#' @return An sf object with exactly one row and a (MULTI)POLYGON geometry
#'   column, in the CRS of the source file.
#'
#' @seealso [thiessen_weights()], [gauge_to_catchment_weights()]
#'
#' @examples
#' \dontrun{
#' catchment <- load_catchment("data/wye_catchment.shp")
#' }
#'
#' @export
load_catchment <- function(filepath, quiet = TRUE) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath, call. = FALSE)
  }

  raw <- sf::st_read(filepath, quiet = quiet)
  raw <- sf::st_make_valid(raw)

  # Explode to individual parts so we can count real polygons
  parts <- sf::st_cast(raw, "POLYGON", warn = FALSE)

  n_parts <- nrow(parts)

  if (n_parts == 0L) {
    stop("No polygon features found in '", filepath, "'.", call. = FALSE)
  }

  if (n_parts > 1L) {
    message(
      sprintf(
        "load_catchment: %d polygon parts found; dissolving to a single catchment boundary.",
        n_parts
      )
    )
    geom <- sf::st_union(parts)
    catchment <- sf::st_as_sf(data.frame(id = 1L), geometry = geom)
  } else {
    catchment <- parts[1L, ]
  }

  check_crs(catchment, expected = 27700L, name = "catchment")

  catchment
}
