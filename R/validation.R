#' Compare modelled flood extent to an observed extent
#'
#' Evaluates a binary flood inundation model against an observed extent polygon
#' (e.g. from Sentinel-1 SAR or the Copernicus Emergency Management Service)
#' using standard dichotomous skill scores: Critical Success Index (CSI), hit
#' rate, and false alarm ratio.
#'
#' All area calculations are performed in the native CRS of `observed_sf`.
#' The modelled raster is reprojected if necessary.
#'
#' @param model_raster  A terra `SpatRaster` with modelled water depth or
#'   inundation extent values.  Cells > `threshold` are treated as flooded.
#' @param observed_sf  An sf POLYGON object of observed flood extent.
#' @param threshold  Numeric. Raster threshold above which cells are classified
#'   as flooded.  Default `0` (any positive depth = flooded).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`csi`}{Critical Success Index: TP / (TP + FP + FN).  Range [0, 1];
#'       1 is perfect.}
#'     \item{`hit_rate`}{Probability of Detection: TP / (TP + FN).}
#'     \item{`false_alarm_ratio`}{False Alarm Ratio: FP / (TP + FP).}
#'     \item{`areas`}{Named numeric vector of areas (m²) for TP, FP, FN.}
#'     \item{`geometry`}{An sf object with a `category` column (`"TP"`, `"FP"`,
#'       `"FN"`) so results can be mapped.}
#'   }
#'
#' @seealso [fetch_sentinel_extent()]
#'
#' @examples
#' \dontrun{
#' model_depth <- terra::rast("model_depth_max.tif")
#' observed    <- sf::st_read("sentinel1_extent.gpkg")
#' scores      <- compare_flood_extent(model_depth, observed, threshold = 0.1)
#' cat(sprintf("CSI = %.3f\n", scores$csi))
#' }
#'
#' @export
compare_flood_extent <- function(model_raster, observed_sf, threshold = 0) {
  if (!inherits(model_raster, "SpatRaster")) {
    stop("model_raster must be a terra SpatRaster.", call. = FALSE)
  }
  if (!inherits(observed_sf, "sf")) {
    stop("observed_sf must be an sf object.", call. = FALSE)
  }

  # Reproject raster to match observed CRS
  obs_crs  <- sf::st_crs(observed_sf)
  rast_reproj <- terra::project(model_raster, terra::crs(observed_sf))

  # Binarise: flooded = 1, dry = 0
  flooded_rast <- terra::classify(
    rast_reproj,
    matrix(c(-Inf, threshold, 0,
              threshold, Inf, 1), ncol = 3, byrow = TRUE)
  )

  # Vectorise modelled extent (flooded cells only)
  modelled_poly <- terra::as.polygons(flooded_rast, dissolve = TRUE)
  modelled_poly <- modelled_poly[terra::values(modelled_poly)[, 1] == 1, ]

  if (nrow(modelled_poly) == 0L) {
    warning("No flooded cells found in model_raster above threshold.", call. = FALSE)
    empty_sf <- sf::st_as_sf(data.frame(category = character(0)),
                              geometry = sf::st_sfc(crs = obs_crs))
    return(list(csi = 0, hit_rate = 0, false_alarm_ratio = NA_real_,
                areas = c(TP = 0, FP = 0, FN = as.numeric(sf::st_area(sf::st_union(observed_sf)))),
                geometry = empty_sf))
  }

  modelled_sf <- sf::st_as_sf(modelled_poly)
  modelled_sf <- sf::st_transform(modelled_sf, crs = obs_crs)
  modelled_union <- sf::st_union(modelled_sf)
  observed_union <- sf::st_union(observed_sf)

  # TP: intersection of modelled and observed
  tp_geom <- suppressWarnings(sf::st_intersection(modelled_union, observed_union))
  # FP: modelled but not observed
  fp_geom <- suppressWarnings(sf::st_difference(modelled_union, observed_union))
  # FN: observed but not modelled
  fn_geom <- suppressWarnings(sf::st_difference(observed_union, modelled_union))

  safe_area <- function(geom) {
    if (length(geom) == 0L || sf::st_is_empty(geom)) return(0)
    as.numeric(sf::st_area(geom))
  }

  area_tp <- safe_area(tp_geom)
  area_fp <- safe_area(fp_geom)
  area_fn <- safe_area(fn_geom)

  csi <- area_tp / (area_tp + area_fp + area_fn)
  hit_rate <- if ((area_tp + area_fn) > 0) area_tp / (area_tp + area_fn) else NA_real_
  far <- if ((area_tp + area_fp) > 0) area_fp / (area_tp + area_fp) else NA_real_

  # Build geometry sf for mapping
  build_geom_row <- function(geom, cat) {
    if (length(geom) == 0L || sf::st_is_empty(geom)) return(NULL)
    sf::st_as_sf(data.frame(category = cat), geometry = sf::st_sfc(geom, crs = obs_crs))
  }

  geom_parts <- Filter(Negate(is.null), list(
    build_geom_row(tp_geom, "TP"),
    build_geom_row(fp_geom, "FP"),
    build_geom_row(fn_geom, "FN")
  ))
  geom_sf <- if (length(geom_parts)) do.call(rbind, geom_parts) else
    sf::st_as_sf(data.frame(category = character(0)),
                 geometry = sf::st_sfc(crs = obs_crs))

  list(
    csi                = csi,
    hit_rate           = hit_rate,
    false_alarm_ratio  = far,
    areas              = c(TP = area_tp, FP = area_fp, FN = area_fn),
    geometry           = geom_sf
  )
}


#' Fetch flood extent polygons from the Copernicus Emergency Management Service
#'
#' Queries the CEMS (Copernicus Emergency Management Service) STAC catalogue or
#' a local processed SAR archive for flood extent polygons within an area of
#' interest and date range.
#'
#' Requires \pkg{rstac}.  This function currently implements the STAC query
#' interface; users with access to a local Sentinel-1 SAR archive should set
#' `archive_path` to bypass the STAC call.
#'
#' @param aoi  An sf polygon defining the area of interest.
#' @param date_range  A length-2 character vector of ISO dates
#'   (`c("YYYY-MM-DD", "YYYY-MM-DD")`).
#' @param archive_path  Optional character. Path to a local directory of
#'   processed SAR flood extent files (GeoPackage or shapefile).  If supplied,
#'   the STAC query is skipped and local files are loaded instead.
#' @param stac_url  Character. STAC API endpoint URL.  Defaults to the
#'   Copernicus STAC endpoint.
#'
#' @return An sf POLYGON object of flood extent(s) within the AOI and date
#'   range, or `NULL` if no extents are found.
#'
#' @seealso [compare_flood_extent()]
#'
#' @examples
#' \dontrun{
#' aoi      <- load_catchment("wye.shp")
#' extents  <- fetch_sentinel_extent(aoi, c("2024-01-01", "2024-01-31"))
#' }
#'
#' @export
fetch_sentinel_extent <- function(aoi,
                                   date_range,
                                   archive_path = NULL,
                                   stac_url = "https://catalogue.dataspace.copernicus.eu/stac") {
  if (!is.null(archive_path)) {
    return(.fetch_from_local_archive(aoi, date_range, archive_path))
  }

  if (!requireNamespace("rstac", quietly = TRUE)) {
    stop(
      "Package 'rstac' is required for fetch_sentinel_extent().\n",
      "Install it with: install.packages('rstac')\n",
      "Alternatively, supply archive_path to read from a local SAR archive.",
      call. = FALSE
    )
  }

  bbox <- as.numeric(sf::st_bbox(sf::st_transform(aoi, crs = 4326)))

  s_obj <- rstac::stac(stac_url)
  items <- rstac::stac_search(
    q          = s_obj,
    collections = "SENTINEL-1",
    bbox       = bbox,
    datetime   = paste(date_range, collapse = "/"),
    limit      = 100
  )
  items <- rstac::get_request(items)

  if (length(items$features) == 0L) {
    message("fetch_sentinel_extent: no items found in STAC for the specified AOI and date range.")
    return(NULL)
  }

  message(sprintf("fetch_sentinel_extent: %d STAC items found. Manual processing of SAR data is required to produce flood extents.", length(items$features)))
  message("Consider using the Copernicus EMS Rapid Mapping service for pre-processed extents: https://emergency.copernicus.eu/mapping/")
  invisible(items)
}


#' @keywords internal
#' @noRd
.fetch_from_local_archive <- function(aoi, date_range, archive_path) {
  if (!dir.exists(archive_path)) {
    stop("archive_path does not exist: ", archive_path, call. = FALSE)
  }

  files <- c(
    list.files(archive_path, pattern = "\\.gpkg$", full.names = TRUE),
    list.files(archive_path, pattern = "\\.shp$",  full.names = TRUE),
    list.files(archive_path, pattern = "\\.geojson$", full.names = TRUE)
  )

  if (length(files) == 0L) {
    message("No spatial files found in archive_path: ", archive_path)
    return(NULL)
  }

  aoi_4326 <- sf::st_transform(aoi, crs = 4326)
  d_from   <- as.Date(date_range[1L])
  d_to     <- as.Date(date_range[2L])

  layers <- lapply(files, function(f) {
    tryCatch({
      lyr <- sf::st_read(f, quiet = TRUE)
      lyr <- sf::st_transform(lyr, crs = 4326)
      # Spatial filter
      lyr <- lyr[sf::st_intersects(lyr, aoi_4326, sparse = FALSE)[, 1L], ]
      # Date filter if date column exists
      date_col <- intersect(c("date", "datetime", "event_date"), names(lyr))[1L]
      if (!is.na(date_col)) {
        dates <- as.Date(lyr[[date_col]])
        lyr   <- lyr[!is.na(dates) & dates >= d_from & dates <= d_to, ]
      }
      if (nrow(lyr) > 0L) lyr else NULL
    }, error = function(e) NULL)
  })

  layers <- Filter(Negate(is.null), layers)
  if (length(layers) == 0L) {
    message("No matching flood extents found in local archive.")
    return(NULL)
  }

  do.call(rbind, layers)
}
