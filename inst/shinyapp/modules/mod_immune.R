#' Immune Analysis Module UI
#'
#' @param id Module ID
#' @return UI elements
mod_immune_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::fluidRow(
    shinydashboard::box(
      title = "Immune Analysis Parameters",
      width = 4,
      status = "primary",
      solidHeader = TRUE,

      shiny::textInput(
        ns("gene_symbol"),
        "Gene Symbol:",
        value = "TP53",
        placeholder = "e.g., TP53, BRCA1, EGFR"
      ),

      shiny::selectInput(
        ns("analysis_type"),
        "Analysis Type:",
        choices = c(
          "Immune Correlation" = "immune_cor",
          "TIL Analysis" = "til",
          "MSI Correlation" = "msi",
          "TMB Correlation" = "tmb"
        ),
        selected = "immune_cor"
      ),

      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'til'", ns("analysis_type")),
        shiny::selectInput(
          ns("cell_type"),
          "Cell Type:",
          choices = c(
            "CD8 T Cells" = "CD8_T_cells",
            "CD4 T Cells" = "CD4_T_cells",
            "Macrophages" = "Macrophages",
            "NK Cells" = "NK_cells",
            "B Cells" = "B_cells",
            "Tregs" = "Tregs"
          ),
          selected = "CD8_T_cells"
        )
      ),

      shiny::selectInput(
        ns("cancer_type"),
        "Cancer Type:",
        choices = c(
          "All Cancers" = "all",
          "BRCA" = "BRCA",
          "LUAD" = "LUAD",
          "LUSC" = "LUSC",
          "COAD" = "COAD",
          "READ" = "READ",
          "HNSC" = "HNSC",
          "STAD" = "STAD"
        ),
        selected = "all"
      ),

      shiny::actionButton(
        ns("analyze_btn"),
        "Run Analysis",
        icon = shiny::icon("play"),
        class = "btn-primary"
      )
    ),

    shinydashboard::box(
      title = "Analysis Results",
      width = 8,
      status = "info",
      solidHeader = TRUE,

      shiny::tabsetPanel(
        id = ns("results_tabs"),

        shiny::tabPanel(
          title = "Plot",
          shiny::br(),
          shiny::plotOutput(ns("analysis_plot"), height = "500px")
        ),

        shiny::tabPanel(
          title = "Statistics",
          shiny::br(),
          shiny::verbatimTextOutput(ns("stats_output"))
        ),

        shiny::tabPanel(
          title = "Data Table",
          shiny::br(),
          DT::dataTableOutput(ns("data_table"))
        )
      )
    )
  )
}

#' Immune Analysis Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
#' @param async_compute Async compute engine (optional)
mod_immune_server <- function(id, app_state, async_compute = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Analysis results
    analysis_results <- shiny::reactiveVal(NULL)

    # Run analysis when button is clicked
    shiny::observeEvent(input$analyze_btn, {
      gene <- input$gene_symbol
      analysis_type <- input$analysis_type
      cancer_type <- input$cancer_type

      if (gene == "") {
        shiny::showNotification("Please enter a gene symbol", type = "error")
        return()
      }

      # Show progress
      shiny::withProgress(message = "Running immune analysis...", value = 0, {
        tryCatch({
          cancers <- if (cancer_type == "all") NULL else cancer_type

          result <- switch(analysis_type,
            "immune_cor" = {
              shiny::incProgress(0.3, detail = "Querying immune data")
              vis_gene_immune_cor(gene, cancers = cancers)
            },
            "til" = {
              shiny::incProgress(0.3, detail = "Querying TIL data")
              cell_type <- input$cell_type
              vis_gene_TIL_cor(gene, cell_type = cell_type, cancers = cancers)
            },
            "msi" = {
              shiny::incProgress(0.3, detail = "Querying MSI data")
              vis_gene_msi_cor(gene, cancers = cancers)
            },
            "tmb" = {
              shiny::incProgress(0.3, detail = "Querying TMB data")
              vis_gene_tmb_cor(gene, cancers = cancers)
            }
          )

          analysis_results(result)
          shiny::incProgress(1, detail = "Complete")
          shiny::showNotification("Analysis complete!", type = "message")

        }, error = function(e) {
          shiny::showNotification(paste("Analysis failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    # Output plot
    output$analysis_plot <- shiny::renderPlot({
      shiny::req(analysis_results())
      result <- analysis_results()

      # Handle different return types
      if (is.ggplot(result)) {
        print(result)
      } else if (is.list(result) && !is.null(result$scatter)) {
        # For TMB correlation which returns list
        print(result$scatter)
      } else {
        # Try to extract plot from result
        plot(result)
      }
    })

    # Output statistics
    output$stats_output <- shiny::renderPrint({
      shiny::req(analysis_results())
      result <- analysis_results()

      cat("Immune Analysis Results\n")
      cat("=======================\n\n")
      cat("Gene:", input$gene_symbol, "\n")
      cat("Analysis Type:", input$analysis_type, "\n")
      cat("Cancer Type:", input$cancer_type, "\n\n")

      # Print result summary based on type
      if (is.list(result) && !is.null(result$correlation)) {
        print(result$correlation)
      } else {
        print(summary(result))
      }
    })

    # Output data table
    output$data_table <- DT::renderDataTable({
      shiny::req(analysis_results())
      result <- analysis_results()

      # Extract data frame from result if possible
      if (is.list(result) && !is.null(result$correlation)) {
        DT::datatable(result$correlation, options = list(pageLength = 10))
      } else {
        DT::datatable(data.frame(Message = "Data table not available for this analysis type"))
      }
    })
  })
}

#' Check if object is a ggplot
#'
#' @param x Object to check
#' @return TRUE if ggplot object
is.ggplot <- function(x) {
  inherits(x, "ggplot")
}
