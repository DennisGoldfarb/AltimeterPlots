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
  total_coeffs <- ncol(basis)
  coeffs_per_fragment <- ncol(coeff_matrix)
  padding <- total_coeffs - coeffs_per_fragment

  if (padding < 0 || padding %% 2 != 0) {
    stop("Unable to align fragment coefficients with knot vector")
  }

  pad_each_side <- padding / 2

  curve_dfs <- lapply(seq_along(fragments), function(idx) {
    fragment_coeffs <- coeff_matrix[idx, ]
    full_coeffs <- c(rep(0, pad_each_side), fragment_coeffs, rep(0, pad_each_side))
    values <- as.numeric(basis %*% full_coeffs)

    data.frame(
      Fragment = fragments[[idx]],
      Position = x_vals,
      Intensity = values,
      stringsAsFactors = FALSE
    )
  })

  curve_data <- do.call(rbind, curve_dfs)

  if (nrow(curve_data) == 0) {
    return(curve_data)
  }

  ordered_facets <- unique(fragments)
  curve_data$Fragment <- factor(curve_data$Fragment, levels = ordered_facets)

  curve_data
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
    req(input$nce)

    curve_data_df <- curve_data()

    if (nrow(curve_data_df) == 0) {
      return(NULL)
    }

    plot <- ggplot(curve_data_df, aes(x = Position, y = Intensity, text = paste("Fragment:", Fragment))) +
      geom_line(color = "#3a80b9") +
      geom_vline(xintercept = input$nce, color = "#d95f02", linetype = "dashed") +
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

    ggplotly(plot, tooltip = c("x", "y", "text"))
  })
}

shinyApp(ui = ui, server = server)
