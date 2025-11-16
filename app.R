#' AltimeterPlots Shiny Application
#'
#' This app provides a minimal UI to submit peptide fragmentation prediction queries
#' to the Altimeter spline model hosted on Koina. It currently supports submitting
#' a single peptide sequence/charge pair and displays the raw predictions table.

library(shiny)
library(koinar)

# Helper that creates a Koina model client and requests predictions for
# a single peptide/charge pair. Returns a data.frame with the response.
fetch_predictions <- function(sequence, charge) {
  stopifnot(is.character(sequence), length(sequence) == 1)
  stopifnot(is.numeric(charge), length(charge) == 1)

  koina_model <- koinar::Koina$new(
    model_name = "Altimeter_2024_splines",
    server_url = "koina.wilhelmlab.org:443",
    ssl = TRUE
  )

  payload <- data.frame(
    peptide_sequences = sequence,
    precursor_charges = charge
  )

  koina_model$predict(payload)
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
      h4("Shared spline knots"),
      uiOutput("knot_vector"),
      h4("Fragment spline coefficients"),
      tableOutput("prediction_table")
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

  parse_prediction_response <- function(raw_prediction) {
    stopifnot(is.data.frame(raw_prediction))

    locate_column <- function(candidates) {
      matched <- intersect(candidates, names(raw_prediction))
      validate(need(length(matched) > 0, paste(
        "Prediction response is missing expected columns:",
        paste(candidates, collapse = ", ")
      )))
      matched[[1]]
    }

    knot_col <- locate_column(c("knot_vector", "knots", "spline_knots"))
    coef_col <- locate_column(c("coefficients", "spline_coefficients"))

    knot_vector <- raw_prediction[[knot_col]][[1]]
    validate(need(length(knot_vector) > 0, "Koina response did not include knot values."))

    coefficient_strings <- vapply(raw_prediction[[coef_col]], function(coeffs) {
      validate(need(length(coeffs) == 4, "Each fragment should provide 4 coefficients."))
      paste(formatC(coeffs, digits = 6, format = "fg", flag = "-"), collapse = ", ")
    }, character(1))

    fragments <- raw_prediction
    fragments[[knot_col]] <- NULL
    fragments[[coef_col]] <- coefficient_strings

    list(knots = knot_vector, fragments = fragments, coefficient_column = coef_col)
  }

  parsed_predictions <- reactive({
    preds <- predictions()
    req(preds)

    validate(need(identical(preds$status, "success"), preds$message))
    parse_prediction_response(preds$data)
  })

  output$status <- renderUI({
    preds <- predictions()
    req(preds)

    if (identical(preds$status, "error")) {
      div(class = "text-danger", paste("Request failed:", preds$message))
    } else {
      div(class = "text-success", "Prediction successful.")
    }
  })

  output$knot_vector <- renderUI({
    parsed <- parsed_predictions()
    req(parsed$knots)

    tags$ul(
      lapply(parsed$knots, function(value) {
        tags$li(format(value, digits = 6))
      })
    )
  })

  output$prediction_table <- renderTable({
    parsed <- parsed_predictions()
    parsed$fragments
  })
}

shinyApp(ui = ui, server = server)
