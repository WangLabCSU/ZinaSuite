# General Analysis Module UI
mod_general_analysis_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shinydashboard::box(
      title = "Analysis Parameters", width = 4, status = "primary", solidHeader = TRUE,
      shiny::selectInput(ns("analysis_type"), "Analysis Type:",
        choices = c("Correlation" = "correlation", "Survival" = "survival")
      ),
      shiny::textInput(ns("gene1"), "Gene 1:", value = "TP53"),
      shiny::textInput(ns("gene2"), "Gene 2:", value = "BRCA1"),
      shiny::actionButton(ns("run_btn"), "Run Analysis", icon = shiny::icon("play"), class = "btn-primary")
    ),
    shinydashboard::box(
      title = "Results", width = 8, status = "info", solidHeader = TRUE,
      shiny::verbatimTextOutput(ns("results"))
    )
  )
}

# General Analysis Module Server
mod_general_analysis_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$run_btn, {
      shiny::withProgress(message = "Running analysis...", {
        tryCatch({
          shiny::incProgress(0.3, detail = "Querying data")
          g1 <- query_gene_expression(input$gene1)
          g2 <- query_gene_expression(input$gene2)
          
          shiny::incProgress(0.7, detail = "Computing correlation")
          common <- intersect(names(g1), names(g2))
          result <- analyze_correlation(g1[common], g2[common])
          
          app_state$analysis_results <- result
          shiny::showNotification("Analysis complete!", type = "message")
        }, error = function(e) {
          shiny::showNotification(paste("Error:", e$message), type = "error")
        })
      })
    })
    
    output$results <- shiny::renderPrint({
      shiny::req(app_state$analysis_results)
      print(app_state$analysis_results)
    })
  })
}
