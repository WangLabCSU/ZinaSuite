#' Data Query Module UI
#'
#' @param id Module ID
#' @return UI elements
mod_data_query_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::fluidRow(
    shinydashboard::box(
      title = "Query Parameters",
      width = 4,
      status = "primary",
      solidHeader = TRUE,

      shiny::selectInput(
        ns("data_source"),
        "Data Source:",
        choices = c("TCGA" = "tcga", "PCAWG" = "pcawg", "CCLE" = "ccle"),
        selected = "tcga"
      ),

      shiny::selectInput(
        ns("data_type"),
        "Data Type:",
        choices = c(
          "Gene Expression (TPM)" = "mRNA",
          "Gene Expression (FPKM)" = "fpkm",
          "Mutation Status" = "mutation",
          "Copy Number Variation" = "cnv",
          "Methylation" = "methylation",
          "miRNA Expression" = "miRNA",
          "Protein Expression" = "protein"
        ),
        selected = "mRNA"
      ),

      shiny::textInput(
        ns("gene_symbol"),
        "Gene Symbol:",
        value = "TP53",
        placeholder = "e.g., TP53, BRCA1, EGFR"
      ),

      shiny::actionButton(
        ns("query_btn"),
        "Query Data",
        icon = shiny::icon("search"),
        class = "btn-primary"
      ),

      shiny::hr(),

      shiny::h5("Query Status:"),
      shiny::verbatimTextOutput(ns("query_status"))
    ),

    shinydashboard::box(
      title = "Query Results",
      width = 8,
      status = "info",
      solidHeader = TRUE,

      shiny::tabsetPanel(
        id = ns("results_tabs"),

        shiny::tabPanel(
          title = "Summary",
          shiny::br(),
          shiny::verbatimTextOutput(ns("data_summary"))
        ),

        shiny::tabPanel(
          title = "Data Table",
          shiny::br(),
          DT::dataTableOutput(ns("data_table"))
        ),

        shiny::tabPanel(
          title = "Distribution Plot",
          shiny::br(),
          shiny::plotOutput(ns("dist_plot"), height = "400px")
        )
      )
    )
  )
}

#' Data Query Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
mod_data_query_server <- function(id, app_state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Query status
    query_status <- shiny::reactiveVal("Ready to query")

    # Perform query when button is clicked
    shiny::observeEvent(input$query_btn, {
      gene <- input$gene_symbol
      data_type <- input$data_type
      source <- input$data_source

      if (gene == "") {
        shiny::showNotification("Please enter a gene symbol", type = "error")
        return()
      }

      query_status(paste("Querying", gene, "from", toupper(source), "..."))

      # Show progress
      shiny::withProgress(message = "Querying data...", value = 0, {
        tryCatch({
          shiny::incProgress(0.3, detail = "Fetching from UCSC Xena")

          # Query data based on source
          result <- switch(source,
            "tcga" = query_molecule(gene, data_type = data_type, source = "tcga"),
            "pcawg" = query_molecule(gene, data_type = data_type, source = "pcawg"),
            "ccle" = query_molecule(gene, data_type = data_type, source = "ccle"),
            stop("Unknown data source")
          )

          shiny::incProgress(0.7, detail = "Processing results")

          if (is.null(result) || length(result) == 0) {
            query_status(paste("No data found for", gene))
            shiny::showNotification(paste("No data found for", gene), type = "warning")
          } else {
            # Store in app state
            app_state$data <- list(
              gene = gene,
              data_type = data_type,
              source = source,
              values = result
            )

            query_status(paste("Successfully queried", gene, "-", length(result), "samples"))
            shiny::showNotification(paste("Query successful:", length(result), "samples"), type = "message")
          }

          shiny::incProgress(1, detail = "Complete")

        }, error = function(e) {
          query_status(paste("Error:", conditionMessage(e)))
          shiny::showNotification(paste("Query failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    # Output query status
    output$query_status <- shiny::renderText({
      query_status()
    })

    # Output data summary
    output$data_summary <- shiny::renderPrint({
      shiny::req(app_state$data)

      data <- app_state$data$values

      cat("Gene:", app_state$data$gene, "\n")
      cat("Data Type:", app_state$data$data_type, "\n")
      cat("Source:", toupper(app_state$data$source), "\n")
      cat("Samples:", length(data), "\n")
      cat("\n")
      cat("Summary Statistics:\n")
      print(summary(data))
    })

    # Output data table
    output$data_table <- DT::renderDataTable({
      shiny::req(app_state$data)

      df <- data.frame(
        Sample = names(app_state$data$values),
        Value = as.numeric(app_state$data$values),
        stringsAsFactors = FALSE
      )

      DT::datatable(
        df,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    })

    # Output distribution plot
    output$dist_plot <- shiny::renderPlot({
      shiny::req(app_state$data)

      values <- as.numeric(app_state$data$values)
      gene <- app_state$data$gene

      viz <- Visualization$new()
      viz$plot_distribution(
        values,
        type = "histogram",
        title = paste("Distribution of", gene, "Expression")
      )
    })
  })
}
