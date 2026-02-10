#' Batch Analysis Module UI
mod_batch_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shinydashboard::box(
      title = "Batch Parameters", width = 4, status = "primary", solidHeader = TRUE,
      shiny::textAreaInput(ns("genes"), "Gene List (one per line):", 
        value = "TP53\nBRCA1\nEGFR", height = "150px"),
      shiny::selectInput(ns("batch_type"), "Analysis Type:",
        choices = c("Correlation" = "correlation", "Cox" = "cox")
      ),
      shiny::numericInput(ns("n_workers"), "Workers:", value = 2, min = 1, max = 4),
      shiny::actionButton(ns("run_btn"), "Run Batch", icon = shiny::icon("tasks"), class = "btn-success")
    ),
    shinydashboard::box(
      title = "Progress & Results", width = 8, status = "info", solidHeader = TRUE,
      shiny::verbatimTextOutput(ns("progress")),
      shiny::hr(),
      DT::dataTableOutput(ns("results_table"))
    )
  )
}

#' Batch Analysis Module Server
mod_batch_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    progress_text <- shiny::reactiveVal("Ready")
    
    shiny::observeEvent(input$run_btn, {
      genes <- strsplit(input$genes, "\n")[[1]]
      genes <- trimws(genes[genes != ""])
      
      if (length(genes) < 2) {
        shiny::showNotification("Please enter at least 2 genes", type = "error")
        return()
      }
      
      progress_text(paste("Running batch analysis on", length(genes), "genes..."))
      
      shiny::withProgress(message = "Batch analysis...", value = 0, {
        tryCatch({
          if (input$batch_type == "correlation") {
            result <- analyze_correlation_batch(
              target_gene = genes[1],
              candidate_genes = genes[-1],
              n_workers = input$n_workers,
              .progress = FALSE
            )
          } else {
            result <- analyze_unicox_batch(
              genes = genes,
              n_workers = input$n_workers,
              .progress = FALSE
            )
          }
          
          app_state$batch_results <- result
          progress_text(paste("Complete! Analyzed", nrow(result), "genes"))
          shiny::showNotification("Batch analysis complete!", type = "message")
        }, error = function(e) {
          progress_text(paste("Error:", e$message))
          shiny::showNotification(paste("Error:", e$message), type = "error")
        })
      })
    })
    
    output$progress <- shiny::renderText({ progress_text() })
    
    output$results_table <- DT::renderDataTable({
      shiny::req(app_state$batch_results)
      DT::datatable(app_state$batch_results, options = list(pageLength = 10))
    })
  })
}
