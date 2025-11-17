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

koina_infer_url <- "https://koina.wilhelmlab.org:443/v2/models/Altimeter_2024_splines/infer"

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

  output_map <- setNames(
    lapply(outputs, function(entry) entry$data),
    vapply(outputs, function(entry) entry$name, character(1))
  )

  fragments <- unlist(resolve_output(output_map, c("annotations")))
  coefficients <- as.numeric(unlist(resolve_output(output_map, c("coefficients"))))
  knots <- as.numeric(unlist(resolve_output(output_map, c("knots"))))
  mzs <- as.numeric(unlist(resolve_output(output_map, c("mz"))))

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

resolve_output <- function(output_map, candidates) {
  for (name in candidates) {
    if (name %in% names(output_map)) {
      return(output_map[[name]])
    }
  }
  NULL
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
  fragments <- predictions$fragments
  coeff_matrix <- predictions$coefficients
  knots <- predictions$knots

  if (length(fragments) == 0 || nrow(coeff_matrix) == 0) {
    return(data.frame())
  }

  fragments_to_use <- seq_len(min(length(fragments), max_fragments))
  fragments <- fragments[fragments_to_use]
  coeff_matrix <- coeff_matrix[fragments_to_use, , drop = FALSE]

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
    coefficients = coeff_matrix,
    knots = knots,
    degree = degree
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

fragment_intensity_at_nce <- function(curve_df, nce_value) {
  if (nrow(curve_df) == 0 || is.null(nce_value) || !is.finite(nce_value)) {
    return(numeric())
  }

  metadata <- attr(curve_df, "fragment_metadata")

  if (is.null(metadata)) {
    return(numeric())
  }

  fragments <- metadata$fragments
  coeff_matrix <- metadata$coefficients
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

  if (length(plotly_traces) == 0) {
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

    if (grepl("Fragment=", trace_name, fixed = TRUE)) {
      candidate <- trimws(sub(".*Fragment=([^,]+).*", "\\1", trace_name))
      return(candidate)
    }

    if (grepl("Fragment:", trace_name, fixed = TRUE)) {
      candidate <- trimws(sub(".*Fragment:\\s*([^,]+).*", "\\1", trace_name))
      return(candidate)
    }

    trimws(trace_name)
  }

  for (trace in line_traces) {
    fragment <- identify_fragment(trace$name)

    if (is.null(fragment) || !(fragment %in% fragments)) {
      next
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
      marker = list(color = "#7a7a7a", size = 6),
      hoverinfo = "x+y",
      showlegend = FALSE,
      xaxis = axes$xaxis,
      yaxis = axes$yaxis
    )
  })

  names(traces) <- fragments
  Filter(Negate(is.null), traces)
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
        step = 0.5
      ),
      div(
        class = "d-flex gap-2", # rely on bootstrap utility classes bundled with shiny
        actionButton("nce_decrement", "-0.5 NCE"),
        actionButton("nce_increment", "+0.5 NCE")
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
      plotlyOutput("fragment_plot", height = "600px")
    )
  )
)

server <- function(input, output, session) {
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

  observeEvent(input$nce_increment, {
    current <- isolate(input$nce)
    new_value <- min(40, ifelse(is.null(current), 30, current + 0.5))
    updateSliderInput(session, "nce", value = new_value)
  })

  observeEvent(input$nce_decrement, {
    current <- isolate(input$nce)
    new_value <- max(20, ifelse(is.null(current), 30, current - 0.5))
    updateSliderInput(session, "nce", value = new_value)
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

      if (length(tracked_fragments) > 0) {
        new_x <- lapply(tracked_fragments, function(fragment) {
          c(input$nce)
        })

        new_y <- lapply(tracked_fragments, function(fragment) {
          c(fragment_values[[fragment]])
        })

        proxy <- plotlyProxyInvoke(
          proxy,
          "restyle",
          list(x = new_x, y = new_y),
          as.list(unname(point_indexes[tracked_fragments]))
        )
      }
    }
  }, ignoreNULL = FALSE)
}

shinyApp(ui = ui, server = server)
