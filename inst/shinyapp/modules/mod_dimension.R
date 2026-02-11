#' Dimensionality Reduction Module UI
#'
#' @param id Module ID
#' @return UI elements
mod_dimension_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::fluidRow(
    shinydashboard::box(
      title = "Dimensionality Reduction Parameters",
      width = 4,
      status = "primary",
      solidHeader = TRUE,

      shiny::textAreaInput(
        ns("gene_list"),
        "Gene List (one per line):",
        value = "TP53\nBRCA1\nEGFR\nKRAS\nMYC",
        height = "150px",
        placeholder = "Enter gene symbols, one per line"
      ),

      shiny::selectInput(
        ns("method"),
        "Method:",
        choices = c(
          "PCA" = "pca",
          "t-SNE" = "tsne",
          "UMAP" = "umap"
        ),
        selected = "pca"
      ),

      shiny::sliderInput(
        ns("n_components"),
        "Number of Components:",
        min = 2,
        max = 10,
        value = 2,
        step = 1
      ),

      shiny::selectInput(
        ns("color_by"),
        "Color By:",
        choices = c(
          "Cancer Type" = "cancer",
          "Sample Type" = "sample_type",
          "None" = "none"
        ),
        selected = "cancer"
      ),

      shiny::actionButton(
        ns("analyze_btn"),
        "Run Analysis",
        icon = shiny::icon("play"),
        class = "btn-primary"
      )
    ),

    shinydashboard::box(
      title = "Dimensionality Reduction Results",
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
          title = "Variance Explained",
          shiny::br(),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'pca'", ns("method")),
            shiny::plotOutput(ns("variance_plot"), height = "300px")
          ),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] != 'pca'", ns("method")),
            shiny::p("Variance explained plot only available for PCA")
          )
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

#' Dimensionality Reduction Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
#' @param async_compute Async compute engine (optional)
mod_dimension_server <- function(id, app_state, async_compute = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Analysis results
    analysis_results <- shiny::reactiveVal(NULL)

    # Run analysis when button is clicked
    shiny::observeEvent(input$analyze_btn, {
      genes <- strsplit(input$gene_list, "\n")[[1]]
      genes <- trimws(genes[genes != ""])
      method <- input$method
      n_components <- input$n_components
      color_by <- input$color_by

      if (length(genes) < 3) {
        shiny::showNotification("Please enter at least 3 genes", type = "error")
        return()
      }

      # Show progress
      shiny::withProgress(message = "Running dimensionality reduction...", value = 0, {
        tryCatch({
          shiny::incProgress(0.2, detail = "Querying gene expression data")

          # Query expression data for all genes
          expr_list <- query_molecules(
            identifiers = genes,
            data_type = "mRNA",
            source = "tcga",
            n_workers = 4,
            .progress = FALSE
          )

          shiny::incProgress(0.5, detail = "Preparing data matrix")

          # Create expression matrix
          expr_matrix <- do.call(rbind, lapply(expr_list, function(x) {
            values <- as.numeric(x)
            names(values) <- names(x)
            values
          }))

          # Remove rows with all NA
          expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) < ncol(expr_matrix), ]

          # Get color variable
          color_var <- NULL
          if (color_by != "none") {
            sample_info <- load_data("tcga_gtex")
            if (color_by == "cancer") {
              color_var <- stats::setNames(
                sample_info$tissue[match(colnames(expr_matrix), sample_info$sample)],
                colnames(expr_matrix)
              )
            }
          }

          shiny::incProgress(0.8, detail = paste("Running", toupper(method)))

          # Run dimensionality reduction
          result <- vis_identifier_dim_dist(
            data = t(expr_matrix),
            method = method,
            n_components = n_components,
            color_by = color_var
          )

          analysis_results(list(
            plot = result,
            method = method,
            data = expr_matrix
          ))

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

      if (!is.null(result$plot)) {
        print(result$plot)
      }
    })

    # Output variance plot (for PCA)
    output$variance_plot <- shiny::renderPlot({
      shiny::req(analysis_results())
      result <- analysis_results()

      if (result$method == "pca") {
        # Compute PCA to get variance
        pca_result <- stats::prcomp(t(result$data), scale. = TRUE)
        var_explained <- summary(pca_result)$importance[2, ]

        # Plot variance explained
        var_df <- data.frame(
          PC = paste0("PC", 1:length(var_explained)),
          Variance = var_explained
        )

        ggplot2::ggplot(var_df, ggplot2::aes(x = PC, y = Variance)) +
          ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            title = "Variance Explained by Principal Components",
            x = "Principal Component",
            y = "Proportion of Variance Explained"
          ) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      }
    })

    # Output data table
    output$data_table <- DT::renderDataTable({
      shiny::req(analysis_results())
      result <- analysis_results()

      if (!is.null(result$data)) {
        # Show summary of data
        summary_df <- data.frame(
          Metric = c("Number of Genes", "Number of Samples", "Method"),
          Value = c(nrow(result$data), ncol(result$data), toupper(result$method))
        )
        DT::datatable(summary_df, options = list(pageLength = 10), rownames = FALSE)
      } else {
        DT::datatable(data.frame(Message = "No data available"))
      }
    })
  })
}
