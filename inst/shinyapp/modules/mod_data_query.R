# Data Query Module UI
mod_data_query_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shinydashboard::box(
      title = "Query Parameters", width = 4, status = "primary", solidHeader = TRUE,
      shiny::selectInput(ns("data_source"), "Data Source:",
        choices = c("TCGA" = "tcga", "PCAWG" = "pcawg", "CCLE" = "ccle")
      ),
      shiny::selectInput(ns("data_type"), "Data Type:",
        choices = c("mRNA" = "mRNA", "Mutation" = "mutation", "CNV" = "cnv")
      ),
      shiny::textInput(ns("gene_symbol"), "Gene Symbol:", value = "TP53"),
      shiny::actionButton(ns("query_btn"), "Query", icon = shiny::icon("search"), class = "btn-primary")
    ),
    shinydashboard::box(
      title = "Results", width = 8, status = "info", solidHeader = TRUE,
      shiny::verbatimTextOutput(ns("query_result"))
    )
  )
}

# Data Query Module Server
mod_data_query_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$query_btn, {
      shiny::withProgress(message = "Querying data...", {
        tryCatch({
          shiny::incProgress(0.3, detail = "Fetching from Xena")
          result <- query_molecule(input$gene_symbol, data_type = input$data_type, source = input$data_source)
          
          app_state$data <- result
          shiny::showNotification("Query successful!", type = "message")
        }, error = function(e) {
          shiny::showNotification(paste("Error:", e$message), type = "error")
        })
      })
    })
    
    output$query_result <- shiny::renderPrint({
      shiny::req(app_state$data)
      summary(app_state$data)
    })
  })
}
