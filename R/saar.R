#' Valid SAAR dataset keys
#' @keywords internal
#' @noRd
SAAR_DATASETS <- c("HadUK_1991_2020", "HadUK_1961_1990", "FEH_1961_1990")


#' Extract spatially averaged SAAR at gauge locations
#'
#' For each gauge a 1 km² bounding box is constructed around the gauge point
#' and intersected with a gridded SAAR (Standard Average Annual Rainfall)
#' shapefile.  The area-weighted mean SAAR within that box is returned.
#'
#' Three gridded SAAR datasets are referenced in current EA operational
#' workflows.  The `dataset` argument documents which period the caller is
#' using; the function itself does not fetch the data — you must supply the
#' pre-loaded `saar_sf` object.
#'
#' @param coords  An sf POINT object or data.frame with `easting`/`northing`
#'   (EPSG:27700) or `long`/`lat` (EPSG:4326) columns.  Typically the output
#'   of `reach.io::find_stations()`.
#' @param saar_sf  An sf polygon object with a `saar` column containing the
#'   gridded SAAR values (mm/year).
#' @param dataset  Character. One of `"HadUK_1991_2020"` (default),
#'   `"HadUK_1961_1990"`, or `"FEH_1961_1990"`.  Used only for validation and
#'   to annotate the output; does not change the computation.
#'
#' @return A data.frame with columns `notation`, `easting`, `northing`, and
#'   `saar` (area-weighted mean SAAR in mm/year).
#'
#' @seealso [saar_adjusted_weights()], [thiessen_weights()]
#'
#' @examples
#' \dontrun{
#' saar_grid <- sf::st_read("data/haduk_saar_1991_2020.gpkg")
#' gauges    <- reach.io::find_stations(names = "wye")
#' gauge_saar <- saar_at_gauges(gauges, saar_grid, dataset = "HadUK_1991_2020")
#' }
#'
#' @export
saar_at_gauges <- function(coords, saar_sf, dataset = "HadUK_1991_2020") {
  dataset <- match.arg(dataset, choices = SAAR_DATASETS)

  if (!"saar" %in% names(saar_sf)) {
    stop(
      "saar_sf must contain a column named 'saar' with SAAR values (mm/year).",
      call. = FALSE
    )
  }

  # Work in BNG
  pts     <- as_sf_points(coords, crs = 27700)
  saar_bng <- sf::st_transform(saar_sf, crs = 27700)

  xy <- sf::st_coordinates(pts)

  saar_vals <- vapply(seq_len(nrow(pts)), function(i) {
    # 1 km² box centred on the gauge (500 m each side)
    box <- sf::st_as_sf(data.frame(x = xy[i, "X"], y = xy[i, "Y"]),
                        coords = c("x", "y"), crs = 27700)
    box <- sf::st_buffer(box, dist = 500, endCapStyle = "SQUARE",
                         joinStyle = "MITRE")

    clipped <- suppressWarnings(sf::st_intersection(saar_bng, box))

    if (nrow(clipped) == 0L) {
      warning(
        sprintf("No SAAR grid cells found within 1 km\u00b2 of gauge %d.", i),
        call. = FALSE
      )
      return(NA_real_)
    }

    # Area-weighted mean
    areas <- as.numeric(sf::st_area(clipped))
    sum(clipped$saar * areas) / sum(areas)
  }, numeric(1L))

  data.frame(
    notation  = if ("notation" %in% names(pts)) pts$notation
                else seq_len(nrow(pts)),
    easting   = xy[, "X"],
    northing  = xy[, "Y"],
    saar      = saar_vals,
    dataset   = dataset,
    stringsAsFactors = FALSE
  )
}


#' Build a three-method gauge weight table using SAAR adjustment
#'
#' Takes Thiessen area weights and per-gauge SAAR values and returns a table
#' with three weighting schemes side-by-side so they can be compared directly:
#'
#' | Column | Method |
#' |---|---|
#' | `thiessen_weight` | Unmodified Voronoi area proportion |
#' | `saar_adjusted` | Thiessen weight × gauge SAAR, normalised to sum to 1 |
#' | `saar_adjusted_rescaled` | `saar_adjusted` re-normalised to exactly 1 (floating-point safe) |
#'
#' The two SAAR-adjusted columns will normally be identical; `saar_adjusted_rescaled`
#' guards against floating-point residuals accumulating in downstream code.
#'
#' @param thiessen_w  data.frame from [thiessen_weights()]; must contain
#'   `notation` and `thiessen_weight` columns.
#' @param gauge_saar  data.frame from [saar_at_gauges()]; must contain
#'   `notation` and `saar` columns.
#' @param catchment_saar  Numeric.  Mean SAAR (mm/year) for the whole
#'   catchment.  Required for the SAAR-adjusted method but not used in the
#'   current calculation — retained for parity with the operational script §7
#'   interface and for future rescaling variants.
#'
#' @return A data.frame with columns `notation`, `thiessen_weight`,
#'   `saar_adjusted`, and `saar_adjusted_rescaled`.
#'
#' @seealso [thiessen_weights()], [saar_at_gauges()]
#'
#' @examples
#' \dontrun{
#' catchment   <- load_catchment("wye.shp")
#' gauges      <- reach.io::find_stations(names = "wye")
#' tw          <- thiessen_weights(gauges, catchment)
#' saar_grid   <- sf::st_read("haduk_saar.gpkg")
#' gs          <- saar_at_gauges(gauges, saar_grid)
#' weight_tbl  <- saar_adjusted_weights(tw, gs, catchment_saar = 1150)
#' }
#'
#' @export
saar_adjusted_weights <- function(thiessen_w, gauge_saar, catchment_saar) {
  if (!is.numeric(catchment_saar) || length(catchment_saar) != 1L) {
    stop("catchment_saar must be a single numeric value (mm/year).", call. = FALSE)
  }

  required_t <- c("notation", "thiessen_weight")
  required_s <- c("notation", "saar")
  miss_t <- setdiff(required_t, names(thiessen_w))
  miss_s <- setdiff(required_s, names(gauge_saar))
  if (length(miss_t)) stop("thiessen_w is missing columns: ", paste(miss_t, collapse = ", "), call. = FALSE)
  if (length(miss_s)) stop("gauge_saar is missing columns: ", paste(miss_s, collapse = ", "), call. = FALSE)

  merged <- merge(
    thiessen_w[, c("notation", "thiessen_weight")],
    gauge_saar[, c("notation", "saar")],
    by = "notation",
    all.x = TRUE,
    sort  = FALSE
  )

  if (anyNA(merged$saar)) {
    warning(
      "Some gauges have no matching SAAR value; their SAAR-adjusted weights will be NA.",
      call. = FALSE
    )
  }

  raw_adj <- merged$thiessen_weight * merged$saar
  total   <- sum(raw_adj, na.rm = TRUE)

  if (total == 0) {
    stop("Sum of (thiessen_weight * saar) is zero; cannot compute SAAR-adjusted weights.", call. = FALSE)
  }

  merged$saar_adjusted          <- raw_adj / total
  merged$saar_adjusted_rescaled <- merged$saar_adjusted / sum(merged$saar_adjusted, na.rm = TRUE)

  merged[, c("notation", "thiessen_weight", "saar_adjusted", "saar_adjusted_rescaled")]
}
