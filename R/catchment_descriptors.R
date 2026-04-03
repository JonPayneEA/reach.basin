#' Extract physical catchment descriptors from national datasets
#'
#' Attempts to extract key FEH catchment descriptors by spatially querying
#' national gridded datasets within the catchment boundary.  Because these
#' datasets are large and must be obtained from their respective sources, this
#' function requires that the caller supplies the relevant sf/terra objects
#' rather than fetching them automatically.
#'
#' Required datasets and their sources:
#' \describe{
#'   \item{HOST (Hydrology of Soil Types)}{Available from NERC/BGS.
#'     `host_sf` should be an sf polygon grid with a `bfihost` column.}
#'   \item{Land Cover Map (LCM)}{Available from UKCEH.
#'     `lcm_rast` should be a terra SpatRaster of land-cover class codes.}
#'   \item{OS Terrain 50 / 5 DEM}{Available from Ordnance Survey.
#'     `dem_rast` should be a terra SpatRaster of elevation in metres.}
#' }
#'
#' Any dataset not supplied is silently skipped and the corresponding output
#' column is set to `NA`.
#'
#' @param catchment  An sf polygon from [load_catchment()].
#' @param host_sf    Optional sf object with a `bfihost` column (HOST dataset).
#' @param lcm_rast   Optional terra SpatRaster (LCM land-cover codes).
#' @param dem_rast   Optional terra SpatRaster (DEM, metres).
#'
#' @return A named list (one row data.frame) with elements:
#'   \describe{
#'     \item{`area_km2`}{Catchment plan area in km².}
#'     \item{`bfihost`}{Area-weighted mean BFI from HOST; `NA` if not supplied.}
#'     \item{`mean_slope_pct`}{Mean slope in percent from DEM; `NA` if not supplied.}
#'     \item{`woodland_fraction`}{Fraction of catchment classified as woodland
#'       (LCM class 1); `NA` if not supplied.}
#'     \item{`dominant_soil`}{Modal HOST soil class code; `NA` if not supplied.}
#'   }
#'
#' @seealso [estimate_qmed()]
#'
#' @examples
#' \dontrun{
#' catchment <- load_catchment("wye.shp")
#' host      <- sf::st_read("host_grid.gpkg")
#' dem       <- terra::rast("os_terrain_50.tif")
#' attrs     <- extract_catchment_attrs(catchment, host_sf = host, dem_rast = dem)
#' }
#'
#' @export
extract_catchment_attrs <- function(catchment,
                                    host_sf  = NULL,
                                    lcm_rast = NULL,
                                    dem_rast = NULL) {
  if (!inherits(catchment, "sf")) {
    stop("catchment must be an sf object.", call. = FALSE)
  }

  # Area in km²
  area_km2 <- as.numeric(sf::st_area(catchment)) / 1e6

  # ---- BFI from HOST --------------------------------------------------------
  bfihost        <- NA_real_
  dominant_soil  <- NA_character_

  if (!is.null(host_sf)) {
    if (!"bfihost" %in% names(host_sf)) {
      warning("host_sf does not have a 'bfihost' column; skipping HOST extraction.", call. = FALSE)
    } else {
      host_proj <- sf::st_transform(host_sf, crs = sf::st_crs(catchment))
      clipped   <- suppressWarnings(sf::st_intersection(host_proj, catchment))
      if (nrow(clipped) > 0L) {
        areas   <- as.numeric(sf::st_area(clipped))
        bfihost <- sum(clipped$bfihost * areas, na.rm = TRUE) / sum(areas, na.rm = TRUE)
        if ("soil_class" %in% names(clipped)) {
          dominant_soil <- names(which.max(tapply(areas, clipped$soil_class, sum)))
        }
      }
    }
  }

  # ---- Mean slope from DEM --------------------------------------------------
  mean_slope_pct <- NA_real_

  if (!is.null(dem_rast)) {
    if (!inherits(dem_rast, "SpatRaster")) {
      warning("dem_rast must be a terra SpatRaster; skipping slope extraction.", call. = FALSE)
    } else if (!requireNamespace("exactextractr", quietly = TRUE)) {
      warning(
        "Package 'exactextractr' is needed for DEM extraction; skipping slope.\n",
        "Install with: install.packages('exactextractr')",
        call. = FALSE
      )
    } else {
      dem_reproj <- terra::project(dem_rast, terra::crs(catchment))
      slope_rast <- terra::terrain(dem_reproj, v = "slope", unit = "degrees")
      catch_reproj <- sf::st_transform(catchment, crs = terra::crs(slope_rast))
      slope_vals   <- exactextractr::exact_extract(slope_rast, catch_reproj, "mean",
                                                   progress = FALSE)
      mean_slope_pct <- tan(slope_vals[[1L]][1L] * pi / 180) * 100
    }
  }

  # ---- Woodland fraction from LCM ------------------------------------------
  woodland_fraction <- NA_real_

  if (!is.null(lcm_rast)) {
    if (!inherits(lcm_rast, "SpatRaster")) {
      warning("lcm_rast must be a terra SpatRaster; skipping LCM extraction.", call. = FALSE)
    } else if (!requireNamespace("exactextractr", quietly = TRUE)) {
      warning(
        "Package 'exactextractr' is needed for LCM extraction; skipping woodland fraction.\n",
        "Install with: install.packages('exactextractr')",
        call. = FALSE
      )
    } else {
      catch_reproj <- sf::st_transform(catchment, crs = terra::crs(lcm_rast))
      # Extract with coverage fractions to compute area-weighted class proportions
      vals <- exactextractr::exact_extract(lcm_rast, catch_reproj, progress = FALSE)[[1L]]
      if (!is.null(vals) && nrow(vals) > 0L) {
        # LCM class 1 = broadleaved woodland; class 2 = coniferous woodland
        woodland_fraction <- sum(vals$coverage_fraction[vals[[1L]] %in% c(1L, 2L)], na.rm = TRUE) /
          sum(vals$coverage_fraction, na.rm = TRUE)
      }
    }
  }

  list(
    area_km2          = area_km2,
    bfihost           = bfihost,
    mean_slope_pct    = mean_slope_pct,
    woodland_fraction = woodland_fraction,
    dominant_soil     = dominant_soil
  )
}


#' Estimate QMED using FEH statistical method
#'
#' Computes the median annual maximum flood (QMED, m³/s) from catchment
#' descriptors using the FEH statistical regression equation:
#'
#' ```
#' QMED = 8.3062 × AREA^0.8510 × 0.1536^(SAAR/1000) × FARL^3.4451 × 0.0460^(BFIHOST²)
#' ```
#'
#' Reference: Kjeldsen, T.R. (2008) *The Revitalised FSR/FEH Rainfall-Runoff
#' Method*, FEH Supplementary Report No. 1, CEH Wallingford.
#'
#' @param descriptors  A named list or one-row data.frame with elements:
#'   \describe{
#'     \item{`area_km2`}{Catchment area (km²).}
#'     \item{`saar`}{Standard Average Annual Rainfall, 1961–1990 (mm/year).}
#'     \item{`farl`}{Flood Attenuation by Reservoirs and Lakes index [0, 1].
#'       Use 1 for catchments with no significant reservoir influence.}
#'     \item{`bfihost`}{Baseflow index from HOST soils [0, 1].}
#'   }
#'
#' @return Numeric. QMED estimate in m³/s.
#'
#' @references
#' Kjeldsen, T.R. (2008) *The Revitalised FSR/FEH Rainfall-Runoff Method*.
#' FEH Supplementary Report No. 1.  Centre for Ecology and Hydrology, Wallingford.
#'
#' @examples
#' estimate_qmed(list(area_km2 = 150, saar = 1050, farl = 0.98, bfihost = 0.44))
#'
#' @export
estimate_qmed <- function(descriptors) {
  required <- c("area_km2", "saar", "farl", "bfihost")
  miss <- setdiff(required, names(descriptors))
  if (length(miss)) {
    stop("descriptors is missing elements: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  area    <- descriptors[["area_km2"]]
  saar    <- descriptors[["saar"]]
  farl    <- descriptors[["farl"]]
  bfihost <- descriptors[["bfihost"]]

  if (any(c(area, saar, farl, bfihost) <= 0)) {
    stop("All descriptor values must be positive.", call. = FALSE)
  }
  if (farl > 1 || bfihost > 1) {
    stop("farl and bfihost must be in [0, 1].", call. = FALSE)
  }

  qmed <- 8.3062 *
    area^0.8510 *
    0.1536^(saar / 1000) *
    farl^3.4451 *
    0.0460^(bfihost^2)

  unname(qmed)
}
