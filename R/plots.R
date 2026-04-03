#' Plot a hypsometric curve with optional gauge elevation overlay
#'
#' Draws the hypsometric (area–elevation) curve for a catchment, showing the
#' cumulative fraction of catchment area at or above each elevation.  If gauge
#' elevations are provided they are overlaid as labelled points so the user
#' can see how the gauge network samples the altitudinal range.
#'
#' Requires \pkg{ggplot2}.  The returned object is a standard `ggplot` so
#' callers can add further layers (themes, titles, etc.) as normal.
#'
#' @param hyps_data  A data.frame with columns:
#'   \describe{
#'     \item{`elevation`}{Elevation above OD (metres).}
#'     \item{`area_fraction`}{Cumulative fraction of catchment area at or above
#'       this elevation, in [0, 1].}
#'   }
#' @param gauges  Optional data.frame with columns `elevation` and (optionally)
#'   `notation` (used as point labels).  If `NULL` (default) no gauge points
#'   are plotted.
#'
#' @return A `ggplot` object.
#'
#' @seealso [elevation_weights()]
#'
#' @examples
#' \dontrun{
#' hyps <- data.frame(
#'   elevation     = c(50, 100, 200, 300, 400, 500),
#'   area_fraction = c(1.0, 0.90, 0.68, 0.42, 0.18, 0.0)
#' )
#' gauges <- data.frame(notation = c("G1", "G2"), elevation = c(120, 360))
#' plot_hypsometric_curve(hyps, gauges)
#' }
#'
#' @export
plot_hypsometric_curve <- function(hyps_data, gauges = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plot_hypsometric_curve().\n",
      "Install it with: install.packages('ggplot2')",
      call. = FALSE
    )
  }

  required_h <- c("elevation", "area_fraction")
  miss_h <- setdiff(required_h, names(hyps_data))
  if (length(miss_h)) {
    stop("hyps_data is missing columns: ", paste(miss_h, collapse = ", "), call. = FALSE)
  }

  p <- ggplot2::ggplot(hyps_data, ggplot2::aes(x = .data$elevation, y = .data$area_fraction)) +
    ggplot2::geom_line(colour = "#1f77b4", linewidth = 1) +
    ggplot2::labs(
      x = "Elevation (m OD)",
      y = "Fraction of catchment area at or above elevation",
      title = "Hypsometric curve"
    ) +
    ggplot2::scale_y_continuous(labels = scales_percent_if_available, limits = c(0, 1)) +
    ggplot2::theme_bw()

  if (!is.null(gauges)) {
    miss_g <- setdiff("elevation", names(gauges))
    if (length(miss_g)) {
      stop("gauges is missing column: elevation", call. = FALSE)
    }

    # Interpolate area_fraction for each gauge elevation
    gauge_fracs <- stats::approx(
      x    = hyps_data$elevation,
      y    = hyps_data$area_fraction,
      xout = gauges$elevation,
      rule = 2L
    )$y

    gauge_df <- data.frame(
      elevation     = gauges$elevation,
      area_fraction = gauge_fracs,
      label         = if ("notation" %in% names(gauges)) gauges$notation
                      else seq_len(nrow(gauges)),
      stringsAsFactors = FALSE
    )

    p <- p +
      ggplot2::geom_point(
        data    = gauge_df,
        mapping = ggplot2::aes(x = .data$elevation, y = .data$area_fraction),
        colour  = "#d62728",
        size    = 3
      ) +
      ggplot2::geom_text(
        data    = gauge_df,
        mapping = ggplot2::aes(x = .data$elevation, y = .data$area_fraction,
                               label = .data$label),
        vjust   = -0.8,
        size    = 3,
        colour  = "#d62728"
      )
  }

  p
}


# Internal helper: use scales::percent if available, else identity formatter
scales_percent_if_available <- function(x) {
  if (requireNamespace("scales", quietly = TRUE)) {
    scales::percent(x)
  } else {
    as.character(x)
  }
}
