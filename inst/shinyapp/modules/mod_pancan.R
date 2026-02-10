#' Pan-Cancer Analysis Module UI
#'
#' @param id Module ID
#' @return UI elements
mod_pancan_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::fluidRow(
    shinydashboard::box(
      title = "Pan-Cancer Analysis Parameters",
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
          "Expression Distribution" = "distribution",
          "Cancer Comparison" = "comparison",
          "Correlation" = "correlation",
          "Survival" = "survival"
        ),
        selected = "distribution"
      ),

      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'correlation'", ns("analysis_type")),
        shiny::textInput(
          ns("gene2"),
          "Second Gene (for correlation):",
          value = "BRCA1",
          placeholder = "e.g., BRCA1, EGFR"
        )
      ),

      shiny::selectInput(
        ns("data_type"),
        "Data Type:",
        choices = c(
          "Gene Expression (TPM)" = "mRNA",
          "Mutation Status" = "mutation",
          "Copy Number Variation" = "cnv"
        ),
        selected = "mRNA"
      ),

      shiny::actionButton(
        ns("analyze_btn"),
        "Run Analysis",
        icon = shiny::icon("play"),
        class = "btn-primary"
      )
    ),

    shinydashboard::box(
      title = "Pan-Cancer Results",
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

#' Pan-Cancer Analysis Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
mod_pancan_server <- function(id, app_state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Analysis results
    analysis_results <- shiny::reactiveVal(NULL)

    # Run analysis when button is clicked
    shiny::observeEvent(input$analyze_btn, {
      gene <- input$gene_symbol
      analysis_type <- input$analysis_type
      data_type <- input$data_type

      if (gene == "") {
        shiny::showNotification("Please enter a gene symbol", type = "error")
        return()
      }

      # Show progress
      shiny::withProgress(message = "Running pan-cancer analysis...", value = 0, {
        tryCatch({
          result <- switch(analysis_type,
            "distribution" = {
              shiny::incProgress(0.3, detail = "Querying pan-cancer data")
              vis_pancan_distribution(gene, data_type = data_type)
            },
            "comparison" = {
              shiny::incProgress(0.3, detail = "Comparing across cancers")
              # Query data for all cancers
              expr <- query_gene_expression(gene, source = "tcga")
              sample_info <- load_data("tcga_gtex")

              # Create comparison data
              plot_data <- data.frame(
                Expression = as.numeric(expr),
                Cancer = sample_info$tissue[match(names(expr), sample_info$sample)]
              )
              plot_data <- plot_data[!is.na(plot_data$Cancer), ]

              # Create boxplot
              ggplot2::ggplot(plot_data, ggplot2::aes(x = Cancer, y = Expression, fill = Cancer)) +
                ggplot2::geom_boxplot() +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
                ggplot2::labs(title = paste(gene, "Expression Across Cancer Types"))
            },
            "correlation" = {
              shiny::incProgress(0.3, detail = "Computing correlation")
              gene2 <- input$gene2
              if (gene2 == "") {
                stop("Please enter a second gene for correlation analysis")
              }
              vis_pancan_correlation(gene, gene2, data_type = data_type)
            },
            "survival" = {
              shiny::incProgress(0.3, detail = "Analyzing survival")
              # This would need to be implemented for pan-cancer survival
              stop("Pan-cancer survival analysis not yet implemented")
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

      cat("Pan-Cancer Analysis Results\n")
      cat("===========================\n\n")
      cat("Gene:", input$gene_symbol, "\n")
      cat("Analysis Type:", input$analysis_type, "\n")
      cat("Data Type:", input$data_type, "\n\n")

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
