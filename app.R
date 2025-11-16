#' AltimeterPlots Shiny Application
#'
#' This app provides a minimal UI to submit peptide fragmentation prediction queries
#' to the Altimeter spline model hosted on Koina. It currently supports submitting
#' a single peptide sequence/charge pair and displays the raw predictions table.

library(shiny)
library(httr)
library(jsonlite)

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

  coeff_matrix <- matrix(coefficients, ncol = coeffs_per_fragment, byrow = TRUE)
  coefficient_strings <- apply(coeff_matrix, 1, function(row) paste(formatC(row, format = "f", digits = 6), collapse = ", "))

  list(
    table = data.frame(
      Fragment = fragments,
      Mz = mzs,
      Coefficients = coefficient_strings,
      stringsAsFactors = FALSE
    ),
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
      tableOutput("prediction_table"),
      h4("Shared knot vector"),
      uiOutput("knot_list")
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

  output$prediction_table <- renderTable({
    preds <- predictions()
    req(preds)

    validate(need(identical(preds$status, "success"), preds$message))
    preds$data$table
  }, striped = TRUE, spacing = "s")

  output$knot_list <- renderUI({
    preds <- predictions()
    req(preds)

    validate(need(identical(preds$status, "success"), preds$message))
    knots <- preds$data$knots

    if (length(knots) == 0) {
      return(div(class = "text-muted", "No knot vector returned."))
    }

    tags$ul(
      lapply(knots, function(value) tags$li(formatC(value, format = "f", digits = 6)))
    )
  })
}

shinyApp(ui = ui, server = server)
