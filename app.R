#' AltimeterPlots Shiny Application
#'
#' This app provides a minimal UI to submit peptide fragmentation prediction queries
#' to the Altimeter spline model hosted on Koina. It currently supports submitting
#' a single peptide sequence/charge pair and displays the raw predictions table.

library(shiny)
library(httr)
library(jsonlite)
library(plotly)
library(splines)
library(ggplot2)
library(OrgMassSpecR)

koina_infer_url <- "https://koina.wilhelmlab.org:443/v2/models/Altimeter_2024_splines/infer"
koina_isotope_infer_url <- "https://koina.wilhelmlab.org:443/v2/models/Altimeter_2024_isotopes/infer"
proton_mass <- 1.007276466812

`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs)) lhs else rhs
}

# Helper that issues an HTTP inference request for a single peptide/charge pair
# and parses the spline metadata.
fetch_predictions <- function(sequence, charge) {
  stopifnot(is.character(sequence), length(sequence) == 1)
  stopifnot(is.numeric(charge), length(charge) == 1)

  payload <- toJSON(
    list(
      id = "0",
      inputs = list(
        list(
          name = "peptide_sequences",
          shape = c(1, 1),
          datatype = "BYTES",
          data = list(sequence)
        ),
        list(
          name = "precursor_charges",
          shape = c(1, 1),
          datatype = "INT32",
          data = list(as.integer(charge))
        )
      )
    ),
    auto_unbox = TRUE
  )

  response <- httr::POST(koina_infer_url, body = payload, encode = "json")
  status <- httr::http_status(response)

  if (!identical(status$category, "Success")) {
    stop("Koina request failed (", status$reason, ")")
  }

  parsed <- httr::content(response)
  parse_prediction_outputs(parsed$outputs)
}

parse_prediction_outputs <- function(outputs) {
  if (length(outputs) == 0) {
    stop("Koina response did not include any outputs")
  }

  output_map <- create_output_entry_map(outputs)

  fragments <- unlist(resolve_output_data(output_map, c("annotations")))
  coefficients <- as.numeric(unlist(resolve_output_data(output_map, c("coefficients"))))
  knots <- as.numeric(unlist(resolve_output_data(output_map, c("knots"))))
  mzs <- as.numeric(unlist(resolve_output_data(output_map, c("mz"))))

  if (is.null(fragments) || is.null(coefficients) || is.null(knots) || is.null(mzs)){
    stop("Koina response missing fragments, coefficients, mzs, or knot vector")
  }

  num_fragments <- length(fragments)
  coeffs_per_fragment <- length(coefficients) / num_fragments

  if (coeffs_per_fragment != 4) {
    stop("Unexpected number of coefficients per fragment: ", coeffs_per_fragment)
  }

  # Koina returns coefficient 1 for all fragments followed by coefficient 2, etc.
  # Reshape so each row corresponds to a fragment in the same order as the
  # annotation vector.
  coeff_matrix <- t(matrix(coefficients, nrow = coeffs_per_fragment, byrow = TRUE))

  valid_fragments <- !is.na(fragments) & nzchar(trimws(fragments))

  list(
    fragments = fragments[valid_fragments],
    mz = mzs[valid_fragments],
    coefficients = coeff_matrix[valid_fragments, , drop = FALSE],
    knots = knots
  )
}

create_output_entry_map <- function(outputs) {
  if (length(outputs) == 0) {
    return(list())
  }

  setNames(outputs, vapply(outputs, function(entry) entry$name, character(1)))
}

resolve_output_data <- function(output_map, candidates) {
  for (name in candidates) {
    if (name %in% names(output_map)) {
      return(output_map[[name]]$data)
    }
  }
  NULL
}


compute_isotope_profile <- function(sequence, charge, max_isotopes = 8, min_relative = 1e-3) {
  stopifnot(is.character(sequence), length(sequence) == 1)
  stopifnot(is.numeric(charge), length(charge) == 1, charge > 0)

  composition <- OrgMassSpecR::ConvertPeptide(sequence, output = "elements", IAA = TRUE)

  if (is.null(composition) || length(composition) == 0) {
    stop("Unable to determine elemental composition for peptide")
  }

  distribution <- OrgMassSpecR::IsotopicDistribution(formula = composition, charge = charge)

  if (!is.data.frame(distribution) || nrow(distribution) == 0) {
    stop("Unable to compute isotopic distribution for peptide")
  }

  distribution <- distribution[seq_len(min(nrow(distribution), max_isotopes)), , drop = FALSE]

  list(
    distribution = distribution,
    monoisotopic_mz = distribution$mz[[1]],
    charge = charge
  )
}

generate_isolation_offsets <- function(limit = 2, step = 0.005) {
  seq(-abs(limit), abs(limit), by = abs(step))
}

isolation_curve_value <- function(mz, center, half_width, scale = 100, slope = 10) {
  if (!is.finite(mz) || !is.finite(center) || !is.finite(half_width) || half_width <= 0) {
    return(0)
  }

  relative <- abs((mz - center) / half_width)
  scale / (1 + (relative)^(2 * slope))
}

compute_isolation_efficiency_matrix <- function(profile, width, offsets, max_isotopes = 5) {
  if (is.null(profile) || is.null(profile$distribution) || length(offsets) == 0) {
    return(matrix(0, nrow = 0, ncol = max_isotopes))
  }

  distribution <- profile$distribution

  if (is.null(distribution$mz)) {
    return(matrix(0, nrow = 0, ncol = max_isotopes))
  }

  isotope_mz <- distribution$mz[seq_len(min(length(distribution$mz), max_isotopes))]

  if (length(isotope_mz) < max_isotopes) {
    isotope_mz <- c(isotope_mz, rep(NA_real_, max_isotopes - length(isotope_mz)))
  }

  half_width <- max(0, width %||% 0) / 2
  offsets <- as.numeric(offsets)

  efficiency_matrix <- matrix(0, nrow = length(offsets), ncol = max_isotopes)

  if (half_width <= 0 || !is.finite(profile$monoisotopic_mz)) {
    return(efficiency_matrix)
  }

  for (i in seq_along(offsets)) {
    center <- profile$monoisotopic_mz + offsets[[i]]
    efficiency_matrix[i, ] <- vapply(
      isotope_mz,
      function(mz) isolation_curve_value(mz, center, half_width) / 100,
      numeric(1)
    )
  }

  efficiency_matrix
}

closest_offset_index <- function(offsets, target, tolerance = 1e-6) {
  if (length(offsets) == 0 || is.na(target)) {
    return(NA_integer_)
  }

  diffs <- abs(offsets - target)
  idx <- which.min(diffs)

  if (length(idx) == 0 || diffs[[idx]] > tolerance) {
    return(NA_integer_)
  }

  idx
}

fetch_isotope_intensities <- function(sequence, charge, nce, efficiency_matrix) {
  stopifnot(is.character(sequence), length(sequence) == 1)
  stopifnot(is.numeric(charge), length(charge) == 1)

  if (is.null(efficiency_matrix) || !is.matrix(efficiency_matrix) || nrow(efficiency_matrix) == 0) {
    stop("No isolation efficiencies available for Koina isotope request")
  }

  num_requests <- nrow(efficiency_matrix)
  num_isotopes <- ncol(efficiency_matrix)

  payload <- toJSON(
    list(
      id = "0",
      inputs = list(
        list(
          name = "peptide_sequences",
          shape = c(num_requests, 1),
          datatype = "BYTES",
          data = rep(sequence, num_requests)
        ),
        list(
          name = "precursor_charges",
          shape = c(num_requests, 1),
          datatype = "INT32",
          data = rep(as.integer(charge), num_requests)
        ),
        list(
          name = "collision_energies",
          shape = c(num_requests, 1),
          datatype = "FP32",
          data = rep(as.numeric(nce), num_requests)
        ),
        list(
          name = "isotope_isolation_efficiencies",
          shape = c(num_requests, num_isotopes),
          datatype = "FP32",
          data = as.numeric(t(efficiency_matrix))
        )
      )
    ),
    auto_unbox = TRUE
  )

  response <- httr::POST(koina_isotope_infer_url, body = payload, encode = "json")
  status <- httr::http_status(response)

  if (!identical(status$category, "Success")) {
    stop("Koina isotope request failed (", status$reason, ")")
  }

  parsed <- httr::content(response)
  parse_isotope_prediction_outputs(parsed$outputs, num_requests)
}

parse_isotope_prediction_outputs <- function(outputs, num_requests) {
  if (length(outputs) == 0) {
    stop("Koina isotope response did not include any outputs")
  }

  output_map <- create_output_entry_map(outputs)

  fragments <- unlist(resolve_output_data(output_map, c("annotations", "fragments")))
  mzs <- as.numeric(unlist(resolve_output_data(output_map, c("mz", "mzs"))))
  intensity_data <- resolve_output_data(output_map, c("intensities", "intensity", "predictions"))

  if (is.null(fragments) || length(fragments) == 0) {
    stop("Koina isotope response missing fragment annotations")
  }

  if (is.null(intensity_data)) {
    stop("Koina isotope response missing intensity predictions")
  }

  intensities <- as.numeric(unlist(intensity_data))
  num_fragments <- length(fragments)

  if (length(intensities) == 0 || (num_requests * num_fragments) != length(intensities)) {
    stop("Koina isotope response has unexpected number of intensities")
  }

  intensity_matrix <- matrix(intensities, nrow = num_requests, byrow = TRUE)

  list(
    fragments = fragments,
    mz = mzs,
    intensities = intensity_matrix
  )
}

expand_fragment_coefficients <- function(coeff_matrix, knots, degree = 3) {
  coeffs_per_fragment <- ncol(coeff_matrix)
  total_coeffs <- length(knots) - degree - 1
  padding <- total_coeffs - coeffs_per_fragment

  if (padding < 0 || padding %% 2 != 0) {
    stop("Unable to align fragment coefficients with knot vector")
  }

  pad_each_side <- padding / 2

  left_pad <- if (pad_each_side > 0) matrix(0, nrow = nrow(coeff_matrix), ncol = pad_each_side) else NULL
  right_pad <- if (pad_each_side > 0) matrix(0, nrow = nrow(coeff_matrix), ncol = pad_each_side) else NULL

  cbind(left_pad, coeff_matrix, right_pad)
}

evaluate_fragment_curves <- function(predictions, degree = 3, x_range = c(20, 40), points = 200, max_fragments = 24) {
  all_fragments <- predictions$fragments
  all_mz_values <- predictions$mz
  all_coeff_matrix <- predictions$coefficients
  knots <- predictions$knots

  if (length(all_fragments) == 0 || nrow(all_coeff_matrix) == 0) {
    return(data.frame())
  }

  fragments_to_use <- seq_len(min(length(all_fragments), max_fragments))
  fragments <- all_fragments[fragments_to_use]
  coeff_matrix <- all_coeff_matrix[fragments_to_use, , drop = FALSE]

  if (length(x_range) != 2 || any(is.na(x_range))) {
    return(data.frame())
  }

  x_start <- x_range[1]
  x_end <- x_range[2]

  if (!is.finite(x_start) || !is.finite(x_end) || x_start >= x_end) {
    return(data.frame())
  }

  x_vals <- seq(x_start, x_end, length.out = points)
  basis <- splines::splineDesign(knots = knots, x = x_vals, ord = degree + 1, outer.ok = TRUE)
  full_coeffs <- expand_fragment_coefficients(coeff_matrix, knots, degree)
  value_matrix <- basis %*% t(full_coeffs)

  curve_data <- data.frame(
    Fragment = factor(rep(fragments, each = length(x_vals)), levels = fragments),
    Position = rep(x_vals, times = length(fragments)),
    Intensity = as.vector(value_matrix),
    stringsAsFactors = FALSE
  )

  if (nrow(curve_data) == 0) {
    return(curve_data)
  }

  attr(curve_data, "fragment_metadata") <- list(
    fragments = fragments,
    mz = if (!is.null(all_mz_values)) all_mz_values[fragments_to_use] else NULL,
    coefficients = coeff_matrix,
    knots = knots,
    degree = degree,
    full_fragments = all_fragments,
    full_mz = all_mz_values,
    full_coefficients = all_coeff_matrix
  )

  curve_data
}

fragment_intensity_ranges <- function(curve_df) {
  if (nrow(curve_df) == 0) {
    return(list())
  }

  fragments <- levels(curve_df$Fragment)

  lapply(fragments, function(fragment) {
    subset_df <- curve_df[curve_df$Fragment == fragment, , drop = FALSE]
    range_vals <- range(subset_df$Intensity, na.rm = TRUE)

    if (!all(is.finite(range_vals)) || range_vals[1] == range_vals[2]) {
      # Expand to avoid zero height guides
      midpoint <- ifelse(is.finite(range_vals[1]), range_vals[1], 0)
      range_vals <- c(midpoint - 0.5, midpoint + 0.5)
    }

    range_vals
  }) |>
    setNames(fragments)
}

fragment_intensity_at_nce <- function(curve_df, nce_value, use_full = FALSE) {
  if (nrow(curve_df) == 0 || is.null(nce_value) || !is.finite(nce_value)) {
    return(numeric())
  }

  metadata <- attr(curve_df, "fragment_metadata")

  if (is.null(metadata)) {
    return(numeric())
  }

  fragments <- if (use_full && !is.null(metadata$full_fragments)) metadata$full_fragments else metadata$fragments
  coeff_matrix <- if (use_full && !is.null(metadata$full_coefficients)) metadata$full_coefficients else metadata$coefficients
  knots <- metadata$knots
  degree <- metadata$degree %||% 3

  if (is.null(fragments) || length(fragments) == 0 || nrow(coeff_matrix) == 0) {
    return(numeric())
  }

  full_coeffs <- expand_fragment_coefficients(coeff_matrix, knots, degree)
  basis <- splines::splineDesign(knots = knots, x = nce_value, ord = degree + 1, outer.ok = TRUE)
  basis_vec <- as.numeric(basis)
  values <- as.numeric(full_coeffs %*% basis_vec)

  setNames(values, fragments)
}

axis_ref_from_name <- function(axis_name, axis_letter) {
  if (is.null(axis_name) || axis_name == paste0(axis_letter, "axis")) {
    return(axis_letter)
  }

  sub("axis", "", axis_name)
}

extract_fragment_axis_map <- function(plotly_traces, fragments) {
  axis_map <- list()

  if (length(plotly_traces) == 0 || length(fragments) == 0) {
    return(axis_map)
  }

  valid_types <- c("scatter", "scattergl")
  line_traces <- Filter(
    function(trace) {
      trace_type <- trace$type %||% ""
      trace_mode <- trace$mode

      trace_type %in% valid_types && (is.null(trace_mode) || grepl("lines", trace_mode, fixed = TRUE))
    },
    plotly_traces
  )

  identify_fragment <- function(trace_name) {
    if (is.null(trace_name) || !nzchar(trace_name)) {
      return(NULL)
    }

    if (grepl("Fragment\\s*=", trace_name)) {
      candidate <- trimws(sub(".*Fragment\\s*=\\s*([^,]+).*", "\\1", trace_name))
      return(candidate)
    }

    if (grepl("Fragment\\s*:", trace_name)) {
      candidate <- trimws(sub(".*Fragment\\s*:\\s*([^,]+).*", "\\1", trace_name))
      return(candidate)
    }

    trimws(trace_name)
  }

  remaining_fragments <- fragments

  for (trace in line_traces) {
    fragment <- identify_fragment(trace$name)

    if (is.null(fragment) || !(fragment %in% fragments)) {
      if (length(remaining_fragments) == 0) {
        next
      }

      fragment <- remaining_fragments[[1]]
    }

    if (!is.null(axis_map[[fragment]])) {
      next
    }

    axis_map[[fragment]] <- list(
      xref = axis_ref_from_name(trace$xaxis %||% "xaxis", "x"),
      yref = axis_ref_from_name(trace$yaxis %||% "yaxis", "y"),
      xaxis = trace$xaxis %||% "xaxis",
      yaxis = trace$yaxis %||% "yaxis"
    )

    remaining_fragments <- remaining_fragments[remaining_fragments != fragment]

    if (length(axis_map) == length(fragments)) {
      break
    }
  }

  if (length(axis_map) < length(fragments) && length(remaining_fragments) > 0) {
    leftover_traces <- Filter(
      function(trace) {
        trace_type <- trace$type %||% ""
        trace_type %in% valid_types
      },
      plotly_traces
    )

    unused_traces <- Filter(
      function(trace) {
        !any(vapply(axis_map, function(entry) {
          identical(entry$xaxis, trace$xaxis %||% "xaxis") &&
            identical(entry$yaxis, trace$yaxis %||% "yaxis")
        }, logical(1)))
      },
      leftover_traces
    )

    for (fragment in remaining_fragments) {
      if (length(unused_traces) == 0) {
        break
      }

      trace <- unused_traces[[1]]
      axis_map[[fragment]] <- list(
        xref = axis_ref_from_name(trace$xaxis %||% "xaxis", "x"),
        yref = axis_ref_from_name(trace$yaxis %||% "yaxis", "y"),
        xaxis = trace$xaxis %||% "xaxis",
        yaxis = trace$yaxis %||% "yaxis"
      )

      unused_traces <- unused_traces[-1]
    }
  }

  axis_map
}

build_vertical_shapes <- function(nce_value, fragments, fragment_ranges, fragment_values, axis_map) {
  shapes <- lapply(fragments, function(fragment) {
    axes <- axis_map[[fragment]]
    range_vals <- fragment_ranges[[fragment]]
    target_value <- fragment_values[[fragment]]

    if (is.null(axes) || is.null(range_vals) || length(range_vals) < 2 || any(!is.finite(range_vals))) {
      return(NULL)
    }

    list(
      type = "line",
      x0 = nce_value,
      x1 = nce_value,
      xref = axes$xref,
      y0 = range_vals[1],
      y1 = range_vals[2],
      yref = axes$yref,
      line = list(color = "#7a7a7a", dash = "dash", width = 1)
    )
  })

  Filter(Negate(is.null), shapes)
}

build_fragment_point_traces <- function(nce_value, fragments, fragment_values, axis_map) {
  traces <- lapply(fragments, function(fragment) {
    axes <- axis_map[[fragment]]
    target_value <- fragment_values[[fragment]]

    if (is.null(axes) || is.na(target_value)) {
      return(NULL)
    }

    list(
      x = c(nce_value),
      y = c(target_value),
      type = "scatter",
      mode = "markers",
      marker = list(color = "#ab1f25", size = 6),
      hoverinfo = "x+y",
      showlegend = FALSE,
      xaxis = axes$xaxis,
      yaxis = axes$yaxis
    )
  })

  names(traces) <- fragments
  Filter(Negate(is.null), traces)
}

build_isotope_segment_vectors <- function(peaks) {
  if (is.null(peaks) || nrow(peaks) == 0) {
    return(list(x = numeric(0), y = numeric(0)))
  }

  percent <- peaks$percent

  if (is.null(percent)) {
    abundance <- peaks$abundance
    percent <- if (is.null(abundance)) rep(0, nrow(peaks)) else abundance * 100
  }

  percent[!is.finite(percent)] <- 0

  x <- rep(peaks$mz, each = 3)
  y <- numeric(length(x))

  indices <- seq_along(x)
  start_points <- indices[(indices - 1) %% 3 == 0]
  end_points <- indices[(indices - 2) %% 3 == 0]
  separators <- indices[indices %% 3 == 0]

  y[start_points] <- 0
  y[end_points] <- percent
  x[separators] <- NA_real_
  y[separators] <- NA_real_

  list(x = x, y = y)
}

normalize_intensity_matrix <- function(intensity_matrix) {
  if (is.null(intensity_matrix) || !is.matrix(intensity_matrix) || nrow(intensity_matrix) == 0) {
    return(matrix(0, nrow = 0, ncol = 0))
  }

  normalized <- intensity_matrix

  for (i in seq_len(nrow(intensity_matrix))) {
    row_values <- as.numeric(intensity_matrix[i, , drop = TRUE])
    max_val <- suppressWarnings(max(row_values, na.rm = TRUE))

    if (!is.finite(max_val) || max_val <= 0) {
      normalized[i, ] <- 0
    } else {
      normalized[i, ] <- row_values / max_val
    }
  }

  normalized[!is.finite(normalized)] <- 0
  normalized
}

build_stem_segment_x <- function(x_values) {
  if (length(x_values) == 0) {
    return(numeric(0))
  }

  segment_x <- rep(x_values, each = 3)
  segment_x[seq(3, length(segment_x), by = 3)] <- NA_real_
  segment_x
}

build_stem_segment_y <- function(heights) {
  if (length(heights) == 0) {
    return(numeric(0))
  }

  heights <- as.numeric(heights)
  heights[!is.finite(heights)] <- 0

  segment_y <- numeric(length(heights) * 3)
  start_points <- seq(1, length(segment_y), by = 3)
  end_points <- start_points + 1
  separators <- start_points + 2

  segment_y[start_points] <- 0
  segment_y[end_points] <- heights
  segment_y[separators] <- NA_real_

  segment_y
}

build_stem_segment_vectors <- function(x_values, heights) {
  if (length(x_values) == 0 || length(heights) == 0) {
    return(list(x = numeric(0), y = numeric(0)))
  }

  list(x = build_stem_segment_x(x_values), y = build_stem_segment_y(heights))
}

build_spectrum_tooltip_prefix <- function(fragments, mz_values) {
  if (length(fragments) == 0) {
    return(character(0))
  }

  if (length(mz_values) != length(fragments)) {
    mz_values <- rep_len(mz_values %||% NA_real_, length(fragments))
  }

  mz_display <- ifelse(is.finite(mz_values), sprintf("%.4f", mz_values), "NA")
  paste0("Fragment: ", fragments, "<br>m/z: ", mz_display, "<br>Normalized intensity: ")
}

build_spectrum_tooltip_vector <- function(prefixes, intensities) {
  if (length(prefixes) == 0 || length(intensities) == 0) {
    return(character(0))
  }

  normalized <- as.numeric(intensities)
  valid <- is.finite(normalized)
  normalized[!valid] <- 0
  formatted <- sprintf("%.3f", normalized)
  tooltip_text <- paste0(prefixes, formatted)
  tooltip_text[!valid] <- NA_character_

  text <- rep(NA_character_, length(prefixes) * 3)
  text_positions <- seq(2, length(text), by = 3)
  text[text_positions] <- tooltip_text
  text
}

build_spectrum_segment_template <- function(fragments, mz_values) {
  list(
    x = build_stem_segment_x(mz_values),
    tooltip_prefix = build_spectrum_tooltip_prefix(fragments, mz_values)
  )
}

is_isotope_annotation <- function(annotations) {
  if (length(annotations) == 0) {
    return(rep(FALSE, 0))
  }

  grepl("\\+i", annotations, ignore.case = TRUE)
}

build_isotope_plot_components <- function(profile, width, center_offset) {
  distribution <- profile$distribution
  monoisotopic_mz <- profile$monoisotopic_mz

  width <- width %||% 0
  center_offset <- center_offset %||% 0

  window_center <- monoisotopic_mz + center_offset
  half_width <- max(0, width) / 2
  x_limits <- if (is.finite(monoisotopic_mz)) monoisotopic_mz + c(-4, 4) else range(distribution$mz, na.rm = TRUE)

  if (any(!is.finite(x_limits))) {
    x_limits <- c(min(distribution$mz, na.rm = TRUE), max(distribution$mz, na.rm = TRUE))
  }

  window_bounds <- c(window_center - half_width, window_center + half_width)
  inside_window <- if (half_width > 0 && all(is.finite(window_bounds))) {
    distribution$mz >= window_bounds[1] & distribution$mz <= window_bounds[2]
  } else {
    rep(FALSE, nrow(distribution))
  }

  distribution$InsideWindow <- inside_window

  mz_range <- range(distribution$mz, na.rm = TRUE)

  if (!all(is.finite(mz_range))) {
    mz_range <- c(window_center - half_width, window_center + half_width)
  }

  span <- max(diff(mz_range), half_width * 2)

  if (!is.finite(span) || span <= 0) {
    span <- max(half_width * 2, 1)
  }

  isolation_curve <- list(x = numeric(0), y = numeric(0), fillcolor = "rgba(58,128,185,0)")

  if (half_width > 0 && is.finite(window_center)) {
    x_vals <- seq(window_center - span, window_center + span, length.out = 400)
    a <- half_width
    b <- 10
    d <- 100
    relative <- if (a > 0) abs((x_vals - window_center) / a) else 0
    y_vals <- d / (1 + (relative)^(2 * b))

    isolation_curve <- list(
      x = x_vals,
      y = y_vals,
      fillcolor = "rgba(58,128,185,0.25)"
    )
  }

  list(
    distribution = distribution,
    x_limits = x_limits,
    window_center = window_center,
    half_width = half_width,
    isolation_curve = isolation_curve,
    peak_segments = build_isotope_segment_vectors(distribution)
  )
}

ensure_plotly_vector <- function(values) {
  if (length(values) == 0) {
    return(NA_real_)
  }

  values
}

wrap_restyle_vector <- function(values) {
  list(ensure_plotly_vector(values))
}

wrap_restyle_value <- function(value) {
  list(value)
}

ui <- fluidPage(
  titlePanel("Altimeter Peptide Fragmentation Predictions"),

  sidebarLayout(
    sidebarPanel(
      textInput(
        inputId = "peptide",
        label = "Peptide sequence",
        placeholder = "e.g. AAAAAKAK"
      ),
      sliderInput(
        inputId = "charge",
        label = "Precursor charge state",
        min = 1,
        max = 7,
        value = 2,
        step = 1
      ),
      sliderInput(
        inputId = "nce",
        label = "Normalized collision energy (NCE)",
        min = 20,
        max = 40,
        value = 30,
        step = 0.1
      ),
      numericInput(
        inputId = "isolation_width",
        label = "Isolation window width (m/z)",
        min = 0.1,
        max = 10,
        value = 2.0,
        step = 0.1
      ),
      div(
        class = "d-flex gap-2",
        actionButton("isolation_width_decrement", "-0.1 width"),
        actionButton("isolation_width_increment", "+0.1 width")
      ),
      div(
        class = "d-flex gap-2 mt-2",
        actionButton("isolation_center_left", "Center -0.005 m/z"),
        actionButton("isolation_center_right", "Center +0.005 m/z")
      ),
      div(
        class = "d-flex gap-2 mt-2",
        actionButton("isolation_center_play", "Center Play"),
        actionButton("isolation_center_stop", "Center Stop")
      ),
      div(
        class = "d-flex gap-2", # rely on bootstrap utility classes bundled with shiny
        actionButton("nce_decrement", "-0.1 NCE"),
        actionButton("nce_increment", "+0.1 NCE")
      ),
      div(
        class = "d-flex gap-2 mt-2",
        actionButton("nce_play", "Play"),
        actionButton("nce_stop", "Stop")
      ),
      actionButton(
        inputId = "submit",
        label = "Predict"
      ),
      helpText(
        "Enter a single peptide sequence and choose a charge state",
        "between 1 and 7. Click Predict to submit the query to Koina."
      )
    ),
    mainPanel(
      h3("Prediction output"),
      uiOutput("status"),
      h4("Precursor isotope envelope"),
      plotlyOutput("isotope_plot", height = "300px"),
      fluidRow(
        column(
          width = 6,
          h4("Fragment spline curves"),
          plotlyOutput("fragment_plot", height = "600px")
        ),
        column(
          width = 6,
          h4("Predicted spectrum"),
          plotlyOutput("spectrum_plot", height = "400px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  nce_playing <- reactiveVal(FALSE)
  nce_direction <- reactiveVal(1)
  nce_timer <- reactiveTimer(200)
  isolation_center_offset <- reactiveVal(0)
  isolation_center_playing <- reactiveVal(FALSE)
  isolation_center_direction <- reactiveVal(1)
  isolation_center_timer <- reactiveTimer(100)
  isotope_prediction_table <- reactiveVal(NULL)
  isolation_offset_grid <- generate_isolation_offsets()
  session$userData$spectrum_plot_ready <- FALSE
  session$userData$spectrum_trace_indices <- NULL
  session$userData$spectrum_segment_x <- NULL

  predictions <- eventReactive(input$submit, {
    req(input$peptide)

    withProgress(message = "Submitting to Koina", value = 0, {
      incProgress(0.5)
      tryCatch(
        {
          result <- fetch_predictions(input$peptide, input$charge)
          incProgress(0.5)
          list(status = "success", data = result)
        },
        error = function(e) {
          list(status = "error", data = NULL, message = conditionMessage(e))
        }
      )
    })
  }, ignoreNULL = TRUE)

  precursor_profile <- eventReactive(input$submit, {
    req(input$peptide, input$charge)

    profile <- compute_isotope_profile(input$peptide, input$charge)
    profile$message <- NULL
    isolation_center_offset(0)
    isolation_center_playing(FALSE)
    isolation_center_direction(1)
    profile
  }, ignoreNULL = TRUE)

  observeEvent({
    list(precursor_profile(), input$isolation_width, input$nce)
  }, {
    profile <- precursor_profile()
    req(profile)
    req(input$peptide, input$charge)

    offsets <- isolation_offset_grid
    efficiency_matrix <- compute_isolation_efficiency_matrix(profile, input$isolation_width, offsets)

    if (is.null(efficiency_matrix) || nrow(efficiency_matrix) == 0) {
      isotope_prediction_table(NULL)
      session$userData$spectrum_segment_x <- NULL
      return()
    }

    nce_value <- input$nce %||% 30
    isotope_prediction_table(NULL)
    session$userData$spectrum_segment_x <- NULL

    tryCatch({
      withProgress(message = "Fetching isolation-aware spectrum", value = 0, {
        result <- fetch_isotope_intensities(input$peptide, input$charge, nce_value, efficiency_matrix)
        normalized_intensities <- normalize_intensity_matrix(result$intensities)
        isotope_mask <- is_isotope_annotation(result$fragments)
        segment_template <- build_spectrum_segment_template(result$fragments, result$mz)
        isotope_segment_template <- build_spectrum_segment_template(result$fragments[isotope_mask], result$mz[isotope_mask])

        isotope_prediction_table(list(
          fragments = result$fragments,
          mz = result$mz,
          offsets = offsets,
          intensities = result$intensities,
          normalized_intensities = normalized_intensities,
          efficiencies = efficiency_matrix,
          nce = nce_value,
          segment_template = segment_template,
          isotope_mask = isotope_mask,
          isotope_segment_template = isotope_segment_template
        ))
        session$userData$spectrum_segment_x <- NULL
      })
    }, error = function(e) {
      showNotification(
        paste("Unable to fetch isolation-aware spectrum:", conditionMessage(e)),
        type = "error",
        duration = 5
      )
      isotope_prediction_table(NULL)
      session$userData$spectrum_segment_x <- NULL
    })
  }, ignoreNULL = TRUE)

  output$status <- renderUI({
    preds <- predictions()
    req(preds)

    if (identical(preds$status, "error")) {
      div(class = "text-danger", paste("Request failed:", preds$message))
    } else {
      div(class = "text-success", "Prediction successful.")
    }
  })

  curve_data <- reactive({
    preds <- predictions()
    req(preds)

    if (identical(preds$status, "error")) {
      return(data.frame())
    }

    evaluate_fragment_curves(preds$data)
  })

  output$isotope_plot <- renderPlotly({
    profile <- precursor_profile()

    if (is.null(profile)) {
      return(NULL)
    }

    session$userData$isotope_plot_ready <- FALSE

    width <- isolate(input$isolation_width %||% 0)
    center_offset <- isolate(isolation_center_offset() %||% 0)
    components <- build_isotope_plot_components(profile, width, center_offset)

    isolation_curve <- components$isolation_curve
    peak_segments <- components$peak_segments

    plot <- plot_ly()

    plot <- plot %>% add_trace(
      x = ensure_plotly_vector(isolation_curve$x),
      y = ensure_plotly_vector(isolation_curve$y),
      type = "scatter",
      mode = "lines",
      line = list(color = "#3a80b9", width = 1),
      fill = "tozeroy",
      fillcolor = isolation_curve$fillcolor,
      hoverinfo = "none",
      name = "Isolation window"
    )

    plot <- plot %>% add_trace(
      x = ensure_plotly_vector(peak_segments$x),
      y = ensure_plotly_vector(peak_segments$y),
      type = "scatter",
      mode = "lines",
      line = list(color = "#193c55", width = 4),
      hoverinfo = "x+y",
      name = "Isotope peaks"
    )

    plot <- plot %>% layout(
      title = paste0("Theoretical isotope distribution (z = ", profile$charge, "+)"),
      xaxis = list(title = "m/z", range = components$x_limits, zeroline = FALSE),
      yaxis = list(title = "Relative intensity", range = c(0, 100.05), ticksuffix = "%"),
      hovermode = "x unified",
      showlegend = FALSE,
      margin = list(t = 40)
    )

    session$userData$isotope_profile_cache <- profile
    session$userData$isotope_trace_indices <- list(isolation = 0, peaks = 1)
    session$userData$isotope_plot_ready <- TRUE

    plot
  })

  observeEvent({
    list(input$isolation_width, isolation_center_offset())
  }, {
    profile <- session$userData$isotope_profile_cache
    traces <- session$userData$isotope_trace_indices

    if (is.null(profile) || is.null(traces) || !isTRUE(session$userData$isotope_plot_ready)) {
      return()
    }

    width <- input$isolation_width %||% 0
    center_offset <- isolation_center_offset() %||% 0
    components <- build_isotope_plot_components(profile, width, center_offset)

    proxy <- plotlyProxy("isotope_plot", session)

    proxy <- plotlyProxyInvoke(
      proxy,
      "restyle",
      list(
        x = wrap_restyle_vector(components$isolation_curve$x),
        y = wrap_restyle_vector(components$isolation_curve$y),
        fillcolor = wrap_restyle_value(components$isolation_curve$fillcolor)
      ),
      list(traces$isolation)
    )

    proxy <- plotlyProxyInvoke(
      proxy,
      "restyle",
      list(
        x = wrap_restyle_vector(components$peak_segments$x),
        y = wrap_restyle_vector(components$peak_segments$y)
      ),
      list(traces$peaks)
    )
  }, ignoreNULL = FALSE)

  output$fragment_plot <- renderPlotly({
    curve_data_df <- curve_data()

    if (nrow(curve_data_df) == 0) {
      session$userData$fragment_axis_map <- NULL
      session$userData$fragment_point_traces <- NULL
      session$userData$fragment_order <- NULL
      return(NULL)
    }

    fragments <- levels(curve_data_df$Fragment)
    fragment_ranges <- fragment_intensity_ranges(curve_data_df)
    initial_nce <- isolate(if (!is.null(input$nce)) input$nce else 30)
    fragment_values <- fragment_intensity_at_nce(curve_data_df, initial_nce)

    plot <- ggplot(curve_data_df, aes(x = Position, y = Intensity, text = paste("Fragment:", Fragment))) +
      geom_line(color = "#3a80b9") +
      facet_wrap(~Fragment, scales = "free_y", nrow = 4, ncol = 6) +
      scale_x_continuous(breaks = c(20, 30, 40)) +
      labs(
        x = "Normalized Collision Energy % (NCE)",
        y = "Predicted abundance",
        title = "Altimeter fragment spline curves"
      ) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "#193c55"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(color = "#193c55"),
        axis.title.y = element_text(color = "#193c55"),
        axis.ticks = element_line(color = "#193c55", size = 0.5),
        axis.line = element_line(color = "#193c55", size = 0.5),
        panel.spacing = grid::unit(0.8, "lines"),
        strip.text = element_text(face = "bold")
      )

    built_plot <- ggplotly(plot, tooltip = c("x", "y", "text")) |>
      plotly_build()

    axis_map <- extract_fragment_axis_map(built_plot$x$data, fragments)
    shapes <- build_vertical_shapes(initial_nce, fragments, fragment_ranges, fragment_values, axis_map)
    built_plot$x$layout$shapes <- shapes

    point_traces <- build_fragment_point_traces(initial_nce, fragments, fragment_values, axis_map)
    existing_traces <- length(built_plot$x$data)
    built_plot$x$data <- c(built_plot$x$data, point_traces)

    if (length(point_traces) > 0) {
      point_indexes <- existing_traces + seq_along(point_traces) - 1
      names(point_indexes) <- names(point_traces)
      session$userData$fragment_point_traces <- point_indexes
    } else {
      session$userData$fragment_point_traces <- NULL
    }

    session$userData$fragment_axis_map <- axis_map
    session$userData$fragment_order <- fragments

    built_plot
  })

  output$spectrum_plot <- renderPlotly({
    session$userData$spectrum_plot_ready <- FALSE
    session$userData$spectrum_trace_indices <- NULL

    plot <- plot_ly()

    plot <- plot %>% add_trace(
      x = ensure_plotly_vector(NA_real_),
      y = ensure_plotly_vector(NA_real_),
      type = "scatter",
      mode = "lines",
      line = list(color = "#193c55", width = 2),
      hoverinfo = "text",
      text = "",
      name = "Predicted peaks"
    )

    plot <- plot %>% add_trace(
      x = ensure_plotly_vector(NA_real_),
      y = ensure_plotly_vector(NA_real_),
      type = "scatter",
      mode = "lines",
      line = list(color = "rgba(255,127,14,0.5)", width = 4),
      hoverinfo = "text",
      text = "",
      name = "Isotope peaks"
    )

    plot <- plot %>% layout(
      title = "Predicted spectrum",
      xaxis = list(title = "m/z"),
      yaxis = list(title = "Normalized intensity", range = c(0, 1.05)),
      hovermode = "x",
      showlegend = FALSE,
      margin = list(t = 40)
    )

    session$userData$spectrum_trace_indices <- list(peaks = 0, isotopes = 1)
    session$userData$spectrum_plot_ready <- TRUE

    plot
  })

  observeEvent({
    list(isotope_prediction_table(), isolation_center_offset())
  }, {
    traces <- session$userData$spectrum_trace_indices

    if (is.null(traces) || is.null(traces$peaks) || !isTRUE(session$userData$spectrum_plot_ready)) {
      return()
    }

    proxy <- plotlyProxy("spectrum_plot", session)
    predictions <- isotope_prediction_table()

    reset_spectrum_traces <- function(message = "Predicted spectrum") {
      session$userData$spectrum_segment_x <- NULL

      plotlyProxyInvoke(
        proxy,
        "restyle",
        list(
          x = wrap_restyle_vector(NA_real_),
          y = wrap_restyle_vector(NA_real_),
          text = wrap_restyle_vector(NA_character_)
        ),
        list(traces$peaks)
      )

      if (!is.null(traces$isotopes)) {
        plotlyProxyInvoke(
          proxy,
          "restyle",
          list(
            x = wrap_restyle_vector(NA_real_),
            y = wrap_restyle_vector(NA_real_),
            text = wrap_restyle_vector(NA_character_)
          ),
          list(traces$isotopes)
        )
      }

      plotlyProxyInvoke(proxy, "relayout", list(title = list(text = message)))
    }

    if (is.null(predictions) || is.null(predictions$offsets) || is.null(predictions$normalized_intensities)) {
      reset_spectrum_traces()
      return()
    }

    offsets <- predictions$offsets
    idx <- closest_offset_index(offsets, isolation_center_offset() %||% 0)

    normalized_matrix <- predictions$normalized_intensities

    if (is.na(idx) || idx > nrow(normalized_matrix)) {
      reset_spectrum_traces()
      return()
    }

    normalized_row <- as.numeric(normalized_matrix[idx, , drop = TRUE])

    if (length(normalized_row) == 0) {
      reset_spectrum_traces()
      return()
    }

    segment_template <- predictions$segment_template %||%
      build_spectrum_segment_template(predictions$fragments, predictions$mz)

    isotope_mask <- predictions$isotope_mask %||% is_isotope_annotation(predictions$fragments)
    if (length(isotope_mask) != length(normalized_row)) {
      isotope_mask <- rep(FALSE, length(normalized_row))
    }
    isotope_template <- predictions$isotope_segment_template %||%
      build_spectrum_segment_template(predictions$fragments[isotope_mask], predictions$mz[isotope_mask])

    segment_cache <- session$userData$spectrum_segment_x

    if (is.null(segment_cache) || !is.list(segment_cache)) {
      segment_cache <- list(peaks = NULL, isotopes = NULL)
    }

    primary_heights <- normalized_row
    primary_heights[isotope_mask] <- 0

    isotope_heights <- normalized_row[isotope_mask]

    y_values <- build_stem_segment_y(primary_heights)
    text_values <- build_spectrum_tooltip_vector(segment_template$tooltip_prefix, ifelse(isotope_mask, NA, normalized_row))

    restyle_payload <- list(
      y = wrap_restyle_vector(if (length(y_values) > 0) y_values else NA_real_),
      text = wrap_restyle_vector(if (length(text_values) > 0) text_values else NA_character_)
    )

    if (is.null(segment_cache$peaks)) {
      restyle_payload$x <- wrap_restyle_vector(if (length(segment_template$x) > 0) segment_template$x else NA_real_)
      segment_cache$peaks <- segment_template$x
    }

    plotlyProxyInvoke(
      proxy,
      "restyle",
      restyle_payload,
      list(traces$peaks)
    )

    isotope_payload <- list(
      y = wrap_restyle_vector(if (length(isotope_heights) > 0) build_stem_segment_y(isotope_heights) else NA_real_),
      text = wrap_restyle_vector(
        if (length(isotope_heights) > 0) {
          build_spectrum_tooltip_vector(isotope_template$tooltip_prefix, isotope_heights)
        } else {
          NA_character_
        }
      )
    )

    if (is.null(segment_cache$isotopes)) {
      isotope_payload$x <- wrap_restyle_vector(if (length(isotope_template$x) > 0) isotope_template$x else NA_real_)
      segment_cache$isotopes <- isotope_template$x
    }

    if (!is.null(traces$isotopes)) {
      plotlyProxyInvoke(
        proxy,
        "restyle",
        isotope_payload,
        list(traces$isotopes)
      )
    }

    session$userData$spectrum_segment_x <- segment_cache

    title_text <- sprintf(
      "Predicted spectrum at %s NCE (offset %+0.3f m/z)",
      predictions$nce %||% input$nce %||% 30,
      offsets[[idx]] %||% 0
    )

    plotlyProxyInvoke(
      proxy,
      "relayout",
      list(title = list(text = title_text))
    )
  }, ignoreNULL = FALSE)

  observeEvent(input$nce_increment, {
    current <- isolate(input$nce)
    new_value <- min(40, ifelse(is.null(current), 30, current + 0.1))
    updateSliderInput(session, "nce", value = new_value)
  })

  observeEvent(input$nce_decrement, {
    current <- isolate(input$nce)
    new_value <- max(20, ifelse(is.null(current), 30, current - 0.1))
    updateSliderInput(session, "nce", value = new_value)
  })

  observeEvent(input$isolation_width_increment, {
    current <- isolate(input$isolation_width %||% 0)
    new_value <- min(10, current + 0.1)
    updateNumericInput(session, "isolation_width", value = new_value)
  })

  observeEvent(input$isolation_width_decrement, {
    current <- isolate(input$isolation_width %||% 0)
    new_value <- max(0.1, current - 0.1)
    updateNumericInput(session, "isolation_width", value = new_value)
  })

  observeEvent(input$isolation_center_left, {
    offset <- isolation_center_offset() %||% 0
    isolation_center_offset(offset - 0.005)
  })

  observeEvent(input$isolation_center_right, {
    offset <- isolation_center_offset() %||% 0
    isolation_center_offset(offset + 0.005)
  })

  observeEvent(input$nce_play, {
    nce_playing(TRUE)
  })

  observeEvent(input$nce_stop, {
    nce_playing(FALSE)
  })

  observeEvent(input$isolation_center_play, {
    isolation_center_playing(TRUE)
  })

  observeEvent(input$isolation_center_stop, {
    isolation_center_playing(FALSE)
  })

  observe({
    nce_timer()

    if (!isTRUE(nce_playing())) {
      return()
    }

    current <- isolate(input$nce %||% 30)
    direction <- nce_direction()
    step <- 0.1
    upper <- 40
    lower <- 20

    next_value <- current + (step * direction)

    if (next_value >= upper) {
      next_value <- upper
      direction <- -1
    } else if (next_value <= lower) {
      next_value <- lower
      direction <- 1
    }

    nce_direction(direction)
    updateSliderInput(session, "nce", value = next_value)
  })

  observe({
    req(isTRUE(isolation_center_playing()))

    isolation_center_timer()

    profile <- isolate(precursor_profile())

    if (is.null(profile) || is.null(profile$monoisotopic_mz) || !is.finite(profile$monoisotopic_mz)) {
      return()
    }

    direction <- isolation_center_direction()
    step <- 0.005
    limit <- 2
    offset <- isolation_center_offset() %||% 0
    next_offset <- offset + (step * direction)

    if (next_offset >= limit) {
      next_offset <- limit
      direction <- -1
    } else if (next_offset <= -limit) {
      next_offset <- -limit
      direction <- 1
    }

    isolation_center_direction(direction)
    isolation_center_offset(next_offset)
  })

  observeEvent({
    list(input$nce, curve_data())
  }, {
    curve_data_df <- curve_data()

    if (nrow(curve_data_df) == 0) {
      return(NULL)
    }

    fragments <- levels(curve_data_df$Fragment)
    fragment_ranges <- fragment_intensity_ranges(curve_data_df)
    fragment_values <- fragment_intensity_at_nce(curve_data_df, input$nce)
    axis_map <- session$userData$fragment_axis_map %||% list()

    proxy <- plotlyProxy("fragment_plot", session)
    proxy <- plotlyProxyInvoke(
      proxy,
      "relayout",
      list(shapes = build_vertical_shapes(input$nce, fragments, fragment_ranges, fragment_values, axis_map))
    )

    point_indexes <- session$userData$fragment_point_traces
    if (!is.null(point_indexes) && length(point_indexes) > 0) {
      tracked_fragments <- intersect(fragments, names(point_indexes))
      tracked_fragments <- tracked_fragments[!is.na(fragment_values[tracked_fragments])]

      if (length(tracked_fragments) > 0) {
        trace_indices <- as.list(unname(unlist(point_indexes[tracked_fragments])))

        if (length(trace_indices) > 0) {
          x_updates <- lapply(tracked_fragments, function(fragment) list(input$nce))
          y_updates <- lapply(tracked_fragments, function(fragment) list(fragment_values[[fragment]]))

          proxy <- plotlyProxyInvoke(
            proxy,
            "restyle",
            list(
              x = x_updates,
              y = y_updates,
              marker = list(color = "#ab1f25", size = 6)
            ),
            trace_indices
          )
        }
      }
    }
  }, ignoreNULL = FALSE)
}

shinyApp(ui = ui, server = server)
