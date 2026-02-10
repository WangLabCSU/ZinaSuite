#' Mutation Analysis Module UI
#'
#' @param id Module ID
#' @return UI elements
mod_mutation_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::fluidRow(
    shinydashboard::box(
      title = "Mutation Analysis Parameters",
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
          "Mutation Frequency" = "frequency",
          "Mutation vs Expression" = "mut_expr",
          "Mutation Survival" = "mut_survival"
        ),
        selected = "frequency"
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
          "READ" = "READ"
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
      title = "Mutation Analysis Results",
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

#' Mutation Analysis Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
mod_mutation_server <- function(id, app_state) {
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
      shiny::withProgress(message = "Running mutation analysis...", value = 0, {
        tryCatch({
          cancers <- if (cancer_type == "all") NULL else cancer_type

          result <- switch(analysis_type,
            "frequency" = {
              shiny::incProgress(0.3, detail = "Querying mutation data")
              vis_mutation_frequency(gene, cancers = cancers)
            },
            "mut_expr" = {
              shiny::incProgress(0.3, detail = "Comparing mutation vs expression")
              # Query mutation and expression data
              mut_data <- query_mutation(gene, source = "tcga")
              expr_data <- query_gene_expression(gene, source = "tcga")

              # Create comparison
              common_samples <- intersect(names(mut_data), names(expr_data))
              plot_data <- data.frame(
                Mutation = mut_data[common_samples],
                Expression = as.numeric(expr_data[common_samples])
              )

              # Create boxplot
              ggplot2::ggplot(plot_data, ggplot2::aes(x = Mutation, y = Expression, fill = Mutation)) +
                ggplot2::geom_boxplot() +
                ggplot2::theme_minimal() +
                ggplot2::labs(
                  title = paste(gene, "Expression by Mutation Status"),
                  x = "Mutation Status",
                  y = "Expression"
                )
            },
            "mut_survival" = {
              shiny::incProgress(0.3, detail = "Analyzing survival by mutation")
              # This would need mutation-based survival analysis
              stop("Mutation survival analysis not yet implemented")
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

      if (is.ggplot(result)) {
        print(result)
      } else {
        plot(result)
      }
    })

    # Output statistics
    output$stats_output <- shiny::renderPrint({
      shiny::req(analysis_results())

      cat("Mutation Analysis Results\n")
      cat("=========================\n\n")
      cat("Gene:", input$gene_symbol, "\n")
      cat("Analysis Type:", input$analysis_type, "\n")
      cat("Cancer Type:", input$cancer_type, "\n\n")

      result <- analysis_results()
      if (is.list(result) && !is.null(result$data)) {
        print(summary(result$data))
      } else {
        cat("Plot object generated successfully\n")
      }
    })

    # Output data table
    output$data_table <- DT::renderDataTable({
      shiny::req(analysis_results())
      result <- analysis_results()

      if (is.list(result) && !is.null(result$data)) {
        DT::datatable(result$data, options = list(pageLength = 10))
      } else {
        DT::datatable(data.frame(Message = "Data table not available for this analysis type"))
      }
    })
  })
}
