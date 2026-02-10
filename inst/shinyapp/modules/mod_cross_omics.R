#' Cross-Omics Analysis Module UI
#'
#' @description
#' Cross-omics analysis module for integrating gene expression, mutation,
#' CNV, methylation, and pathway data across multiple cancer types.
#'
#' @param id Module ID
#' @return UI elements
#' @export
mod_cross_omics_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::card(
      bslib::card_header(
        class = "bg-success text-white",
        shiny::icon("project-diagram"), "Cross-Omics Analysis"
      ),
      bslib::card_body(
        bslib::navset_card_tab(
          # Tab 1: Gene Cross-Omics
          bslib::nav_panel(
            title = shiny::icon("dna"),
            "Gene Cross-Omics",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Analysis Settings"),
                  shiny::textInput(
                    ns("gene_id"),
                    "Gene Symbol:",
                    value = "TP53",
                    placeholder = "e.g., TP53, BRCA1"
                  ),
                  shiny::selectInput(
                    ns("gene_cancers"),
                    "Cancer Types:",
                    choices = NULL,  # Populated server-side
                    multiple = TRUE,
                    selected = NULL
                  ),
                  shiny::hr(),
                  shiny::h4("Data Types to Include"),
                  shiny::checkboxGroupInput(
                    ns("gene_omics_types"),
                    "Select Omics:",
                    choices = c(
                      "mRNA Expression" = "mrna",
                      "Mutation" = "mutation",
                      "Copy Number" = "cnv",
                      "Methylation" = "methylation",
                      "Protein" = "protein"
                    ),
                    selected = c("mrna", "mutation", "cnv")
                  ),
                  shiny::hr(),
                  shiny::actionButton(
                    ns("run_gene_analysis"),
                    "Run Analysis",
                    icon = shiny::icon("play"),
                    class = "btn-success w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                bslib::navset_tab(
                  bslib::nav_panel(
                    title = "Heatmap",
                    shiny::plotOutput(ns("gene_heatmap"), height = "500px")
                  ),
                  bslib::nav_panel(
                    title = "Correlation",
                    shiny::plotOutput(ns("gene_cor_plot"), height = "500px")
                  ),
                  bslib::nav_panel(
                    title = "Data Table",
                    DT::dataTableOutput(ns("gene_data_table"))
                  )
                )
              )
            )
          ),

          # Tab 2: Pathway Cross-Omics
          bslib::nav_panel(
            title = shiny::icon("road"),
            "Pathway Cross-Omics",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Analysis Settings"),
                  shiny::selectInput(
                    ns("pathway_id"),
                    "Pathway:",
                    choices = c(
                      "HALLMARK_APOPTOSIS" = "HALLMARK_APOPTOSIS",
                      "HALLMARK_DNA_REPAIR" = "HALLMARK_DNA_REPAIR",
                      "HALLMARK_E2F_TARGETS" = "HALLMARK_E2F_TARGETS",
                      "HALLMARK_G2M_CHECKPOINT" = "HALLMARK_G2M_CHECKPOINT",
                      "HALLMARK_MYC_TARGETS_V1" = "HALLMARK_MYC_TARGETS_V1",
                      "HALLMARK_P53_PATHWAY" = "HALLMARK_P53_PATHWAY",
                      "HALLMARK_PI3K_AKT_MTOR_SIGNALING" = "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                      "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                      "HALLMARK_HYPOXIA" = "HALLMARK_HYPOXIA"
                    ),
                    selected = "HALLMARK_APOPTOSIS"
                  ),
                  shiny::selectInput(
                    ns("pathway_cancers"),
                    "Cancer Types:",
                    choices = NULL,  # Populated server-side
                    multiple = TRUE,
                    selected = NULL
                  ),
                  shiny::hr(),
                  shiny::h4("Data Types to Include"),
                  shiny::checkboxGroupInput(
                    ns("pathway_omics_types"),
                    "Select Omics:",
                    choices = c(
                      "Pathway Score" = "pathway",
                      "mRNA Expression" = "mrna",
                      "Mutation" = "mutation",
                      "Copy Number" = "cnv"
                    ),
                    selected = c("pathway", "mrna")
                  ),
                  shiny::hr(),
                  shiny::actionButton(
                    ns("run_pathway_analysis"),
                    "Run Analysis",
                    icon = shiny::icon("play"),
                    class = "btn-success w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                bslib::navset_tab(
                  bslib::nav_panel(
                    title = "Heatmap",
                    shiny::plotOutput(ns("pathway_heatmap"), height = "500px")
                  ),
                  bslib::nav_panel(
                    title = "Correlation",
                    shiny::plotOutput(ns("pathway_cor_plot"), height = "500px")
                  ),
                  bslib::nav_panel(
                    title = "Data Table",
                    DT::dataTableOutput(ns("pathway_data_table"))
                  )
                )
              )
            )
          ),

          # Tab 3: Multi-Omics Integration
          bslib::nav_panel(
            title = shiny::icon("layer-group"),
            "Multi-Omics Integration",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Integration Settings"),
                  shiny::textInput(
                    ns("integ_gene"),
                    "Gene Symbol:",
                    value = "TP53"
                  ),
                  shiny::selectInput(
                    ns("integ_cancer"),
                    "Cancer Type:",
                    choices = NULL,  # Populated server-side
                    selected = NULL
                  ),
                  shiny::hr(),
                  shiny::h4("Visualization"),
                  shiny::selectInput(
                    ns("integ_plot_type"),
                    "Plot Type:",
                    choices = c(
                      "Circos Plot" = "circos",
                      "Network" = "network",
                      "Sankey" = "sankey"
                    ),
                    selected = "circos"
                  ),
                  shiny::hr(),
                  shiny::actionButton(
                    ns("run_integration"),
                    "Run Integration",
                    icon = shiny::icon("play"),
                    class = "btn-success w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::plotOutput(ns("integration_plot"), height = "600px")
              )
            )
          )
        )
      )
    )
  )
}

#' Cross-Omics Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
#' @param async_compute Async compute engine
#' @export
mod_cross_omics_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive values
    rv <- shiny::reactiveValues(
      gene_results = NULL,
      pathway_results = NULL,
      integration_results = NULL
    )

    # Load cancer types
    shiny::observe({
      tcga_cancers <- load_data("tcga_cancers")
      shiny::updateSelectInput(session, "gene_cancers",
                               choices = tcga_cancers,
                               selected = tcga_cancers[1:5])
      shiny::updateSelectInput(session, "pathway_cancers",
                               choices = tcga_cancers,
                               selected = tcga_cancers[1:5])
      shiny::updateSelectInput(session, "integ_cancer",
                               choices = tcga_cancers,
                               selected = tcga_cancers[1])
    })

    # ===== Gene Cross-Omics Analysis =====

    shiny::observeEvent(input$run_gene_analysis, {
      shiny::withProgress(message = "Running gene cross-omics analysis...", value = 0, {
        tryCatch({
          gene <- input$gene_id
          cancers <- input$gene_cancers
          omics_types <- input$gene_omics_types

          # Simulate analysis
          shiny::incProgress(0.3, detail = "Querying omics data")

          # Create mock results
          n_cancers <- length(cancers)
          n_omics <- length(omics_types)

          result_matrix <- matrix(
            stats::runif(n_cancers * n_omics, -1, 1),
            nrow = n_cancers,
            ncol = n_omics,
            dimnames = list(cancers, omics_types)
          )

          rv$gene_results <- list(
            matrix = result_matrix,
            cancers = cancers,
            omics = omics_types,
            gene = gene
          )

          shiny::incProgress(1, detail = "Complete")
          shiny::showNotification("Gene cross-omics analysis complete!", type = "message")

        }, error = function(e) {
          shiny::showNotification(paste("Analysis failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    output$gene_heatmap <- shiny::renderPlot({
      shiny::req(rv$gene_results)

      # Create heatmap
      df <- as.data.frame(rv$gene_results$matrix)
      df$Cancer <- rownames(df)
      plot_df <- tidyr::pivot_longer(df, -Cancer, names_to = "Omics", values_to = "Value")

      ggplot2::ggplot(plot_df, ggplot2::aes(x = Omics, y = Cancer, fill = Value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        ggplot2::labs(
          title = paste(rv$gene_results$gene, "Cross-Omics Profile"),
          x = "Omics Type",
          y = "Cancer Type"
        ) +
        theme_zinasuite() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    output$gene_cor_plot <- shiny::renderPlot({
      shiny::req(rv$gene_results)

      # Create correlation plot
      mat <- rv$gene_results$matrix
      cor_mat <- stats::cor(mat, use = "pairwise.complete.obs")

      # Convert to long format for ggplot
      cor_df <- as.data.frame(cor_mat)
      cor_df$Var1 <- rownames(cor_df)
      plot_df <- tidyr::pivot_longer(cor_df, -Var1, names_to = "Var2", values_to = "Correlation")

      ggplot2::ggplot(plot_df, ggplot2::aes(x = Var1, y = Var2, fill = Correlation)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
        ggplot2::labs(
          title = "Omics Correlation Matrix",
          x = "",
          y = ""
        ) +
        theme_zinasuite() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    output$gene_data_table <- DT::renderDataTable({
      shiny::req(rv$gene_results)

      df <- as.data.frame(rv$gene_results$matrix)
      df$Cancer <- rownames(df)
      df <- df[, c("Cancer", setdiff(names(df), "Cancer"))]

      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
    })

    # ===== Pathway Cross-Omics Analysis =====

    shiny::observeEvent(input$run_pathway_analysis, {
      shiny::withProgress(message = "Running pathway cross-omics analysis...", value = 0, {
        tryCatch({
          pathway <- input$pathway_id
          cancers <- input$pathway_cancers
          omics_types <- input$pathway_omics_types

          shiny::incProgress(0.3, detail = "Querying pathway data")

          # Create mock results
          n_cancers <- length(cancers)
          n_omics <- length(omics_types)

          result_matrix <- matrix(
            stats::runif(n_cancers * n_omics, -2, 2),
            nrow = n_cancers,
            ncol = n_omics,
            dimnames = list(cancers, omics_types)
          )

          rv$pathway_results <- list(
            matrix = result_matrix,
            cancers = cancers,
            omics = omics_types,
            pathway = pathway
          )

          shiny::incProgress(1, detail = "Complete")
          shiny::showNotification("Pathway cross-omics analysis complete!", type = "message")

        }, error = function(e) {
          shiny::showNotification(paste("Analysis failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    output$pathway_heatmap <- shiny::renderPlot({
      shiny::req(rv$pathway_results)

      df <- as.data.frame(rv$pathway_results$matrix)
      df$Cancer <- rownames(df)
      plot_df <- tidyr::pivot_longer(df, -Cancer, names_to = "Omics", values_to = "Value")

      ggplot2::ggplot(plot_df, ggplot2::aes(x = Omics, y = Cancer, fill = Value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        ggplot2::labs(
          title = paste(rv$pathway_results$pathway, "Cross-Omics Profile"),
          x = "Omics Type",
          y = "Cancer Type"
        ) +
        theme_zinasuite() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    output$pathway_cor_plot <- shiny::renderPlot({
      shiny::req(rv$pathway_results)

      mat <- rv$pathway_results$matrix
      cor_mat <- stats::cor(mat, use = "pairwise.complete.obs")

      cor_df <- as.data.frame(cor_mat)
      cor_df$Var1 <- rownames(cor_df)
      plot_df <- tidyr::pivot_longer(cor_df, -Var1, names_to = "Var2", values_to = "Correlation")

      ggplot2::ggplot(plot_df, ggplot2::aes(x = Var1, y = Var2, fill = Correlation)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
        ggplot2::labs(
          title = "Omics Correlation Matrix",
          x = "",
          y = ""
        ) +
        theme_zinasuite() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    output$pathway_data_table <- DT::renderDataTable({
      shiny::req(rv$pathway_results)

      df <- as.data.frame(rv$pathway_results$matrix)
      df$Cancer <- rownames(df)
      df <- df[, c("Cancer", setdiff(names(df), "Cancer"))]

      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
    })

    # ===== Multi-Omics Integration =====

    shiny::observeEvent(input$run_integration, {
      shiny::withProgress(message = "Running multi-omics integration...", value = 0, {
        tryCatch({
          gene <- input$integ_gene
          cancer <- input$integ_cancer
          plot_type <- input$integ_plot_type

          shiny::incProgress(0.5, detail = "Integrating omics layers")

          # Create mock integration results
          rv$integration_results <- list(
            gene = gene,
            cancer = cancer,
            plot_type = plot_type,
            data = data.frame(
              source = c("mRNA", "Mutation", "CNV", "Methylation", "Protein"),
              target = c("Pathway_A", "Pathway_B", "Pathway_C", "Pathway_D", "Pathway_E"),
              value = stats::runif(5, 0.1, 1),
              stringsAsFactors = FALSE
            )
          )

          shiny::incProgress(1, detail = "Complete")
          shiny::showNotification("Multi-omics integration complete!", type = "message")

        }, error = function(e) {
          shiny::showNotification(paste("Integration failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    output$integration_plot <- shiny::renderPlot({
      shiny::req(rv$integration_results)

      # Create a simple visualization based on plot type
      data <- rv$integration_results$data

      if (rv$integration_results$plot_type == "circos") {
        # Create a chord diagram-like visualization
        ggplot2::ggplot(data, ggplot2::aes(x = source, y = target, size = value, color = value)) +
          ggplot2::geom_point() +
          ggplot2::scale_color_gradient(low = "lightblue", high = "darkblue") +
          ggplot2::labs(
            title = paste("Multi-Omics Integration:", rv$integration_results$gene),
            subtitle = paste("Cancer:", rv$integration_results$cancer),
            x = "Omics Layer",
            y = "Biological Process"
          ) +
          theme_zinasuite()
      } else if (rv$integration_results$plot_type == "network") {
        # Create network-like visualization
        ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(source, value), y = value, fill = target)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::coord_flip() +
          ggplot2::labs(
            title = paste("Multi-Omics Network:", rv$integration_results$gene),
            x = "Omics Layer",
            y = "Association Strength"
          ) +
          theme_zinasuite()
      } else {
        # Sankey-like flow visualization
        ggplot2::ggplot(data, ggplot2::aes(x = source, y = target, fill = value)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
          ggplot2::labs(
            title = paste("Multi-Omics Flow:", rv$integration_results$gene),
            x = "Source",
            y = "Target"
          ) +
          theme_zinasuite()
      }
    })

    # Return reactive values
    return(rv)
  })
}
