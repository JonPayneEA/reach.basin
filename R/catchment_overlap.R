#' Identify nested or overlapping catchments in a gauge registry
#'
#' Given a collection of catchment polygons (e.g. one per gauge in a registry),
#' this function computes pairwise overlap proportions and flags pairs where
#' one catchment substantially contains another.  Overlapping catchments
#' introduce double-counting risk when aggregating areal rainfall estimates
#' across sites.
#'
#' @param registry  An sf object with one row per catchment.  Must have a
#'   column that uniquely identifies each catchment; by default the function
#'   looks for `notation`, then `id`, then falls back to row numbers.
#' @param threshold  Numeric [0, 1].  Pairs where the smaller catchment's
#'   overlap with the larger exceeds this proportion are flagged.  Default
#'   `0.5` (50 %).
#' @param id_col  Character.  Name of the column in `registry` used as the
#'   catchment identifier.  If `NULL` (default) the function auto-detects as
#'   described above.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`matrix`}{A symmetric numeric matrix of overlap proportions.
#'     Entry `[i, j]` is the area of the intersection of catchments *i* and *j*
#'     divided by the area of catchment *i* (i.e. how much of *i* is covered
#'     by *j*).  Diagonal is 1.}
#'   \item{`flagged`}{A data.frame of pairs that exceed `threshold`, with
#'     columns `id_a`, `id_b`, `overlap_a` (fraction of a covered by b), and
#'     `overlap_b` (fraction of b covered by a).  Zero rows if no pairs
#'     exceed the threshold.}
#' }
#'
#' @examples
#' \dontrun{
#' catchments <- sf::st_read("data/gauge_catchments.gpkg")
#' result <- catchment_overlap_matrix(catchments, threshold = 0.5)
#' print(result$flagged)
#' }
#'
#' @export
catchment_overlap_matrix <- function(registry, threshold = 0.5, id_col = NULL) {
  if (!inherits(registry, "sf")) {
    stop("registry must be an sf object.", call. = FALSE)
  }
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("threshold must be a numeric value in [0, 1].", call. = FALSE)
  }

  # Determine ID column
  if (is.null(id_col)) {
    if ("notation" %in% names(registry)) {
      id_col <- "notation"
    } else if ("id" %in% names(registry)) {
      id_col <- "id"
    } else {
      registry$..row_id <- seq_len(nrow(registry))
      id_col <- "..row_id"
    }
  } else {
    if (!id_col %in% names(registry)) {
      stop("id_col '", id_col, "' not found in registry.", call. = FALSE)
    }
  }

  ids   <- registry[[id_col]]
  n     <- nrow(registry)
  areas <- as.numeric(sf::st_area(registry))

  registry <- sf::st_make_valid(registry)

  # Build overlap matrix
  overlap_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(ids, ids))
  diag(overlap_mat) <- 1

  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      inter <- suppressWarnings(sf::st_intersection(
        registry[i, ],
        registry[j, ]
      ))
      if (nrow(inter) == 0L) next

      inter_area <- as.numeric(sum(sf::st_area(inter)))
      if (inter_area == 0) next

      overlap_mat[i, j] <- inter_area / areas[i]
      overlap_mat[j, i] <- inter_area / areas[j]
    }
  }

  # Identify flagged pairs (upper triangle only to avoid duplicates)
  flagged_rows <- list()
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      oa <- overlap_mat[i, j]
      ob <- overlap_mat[j, i]
      if (oa > threshold || ob > threshold) {
        flagged_rows[[length(flagged_rows) + 1L]] <- data.frame(
          id_a      = ids[i],
          id_b      = ids[j],
          overlap_a = oa,
          overlap_b = ob,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  flagged <- if (length(flagged_rows)) {
    do.call(rbind, flagged_rows)
  } else {
    data.frame(
      id_a = character(0), id_b = character(0),
      overlap_a = numeric(0), overlap_b = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  list(matrix = overlap_mat, flagged = flagged)
}
