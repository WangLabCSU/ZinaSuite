#' PharmacoGenomics Module
#'
#' @description
#' Module for drug-omics correlation analysis and pharmacogenomics visualization.
#' Following UCSCXenaShiny PharmacoGenomics functionality.
#'
#' @param id Module ID
#' @return Shiny module UI
#' @export
mod_pharmacogenomics_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::card(
      bslib::card_header(
        class = "bg-primary text-white",
        shiny::icon("pills"), "PharmacoGenomics Analysis"
      ),
      bslib::card_body(
        bslib::navset_card_tab(
          # Drug-Omics Correlation Tab
          bslib::nav_panel(
            title = shiny::tagList(shiny::icon("pills"), "Drug-Omics Correlation"),

            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                open = TRUE,
                width = 300,

                shiny::selectInput(
                  ns("data_source"),
                  "Data Source:",
                  choices = c(
                    "CCLE" = "ccle",
                    "GDSC" = "gdsc"
                  ),
                  selected = "ccle"
                ),

                shiny::textInput(
                  ns("gene_id"),
                  "Gene Symbol:",
                  value = "TP53",
                  placeholder = "e.g., TP53, BRCA1"
                ),

                shiny::selectInput(
                  ns("drug_id"),
                  "Drug:",
                  choices = c("Select drug..." = ""),
                  selected = ""
                ),

                shiny::selectInput(
                  ns("sensitivity_metric"),
                  "Sensitivity Metric:",
                  choices = c(
                    "AUC" = "auc",
                    "IC50" = "ic50"
                  ),
                  selected = "auc"
                ),

                shiny::selectInput(
                  ns("cor_method"),
                  "Correlation Method:",
                  choices = c(
                    "Spearman" = "spearman",
                    "Pearson" = "pearson"
                  ),
                  selected = "spearman"
                ),

                shiny::hr(),

                shiny::actionButton(
                  ns("run_correlation"),
                  "Run Analysis",
                  icon = shiny::icon("play"),
                  class = "btn-primary w-100"
                )
              ),

              # Main content
              bslib::navset_card_tab(
                bslib::nav_panel(
                  title = shiny::icon("chart-scatter"),
                  shiny::plotOutput(ns("correlation_plot"), height = "500px")
                ),
                bslib::nav_panel(
                  title = shiny::icon("table"),
                  DT::dataTableOutput(ns("correlation_table"))
                ),
                bslib::nav_panel(
                  title = shiny::icon("info-circle"),
                  shiny::verbatimTextOutput(ns("correlation_stats"))
                )
              )
            )
          ),

          # Drug Sensitivity by Mutation Tab
          bslib::nav_panel(
            title = shiny::tagList(shiny::icon("dna"), "Drug-Mutation Analysis"),

            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                open = TRUE,
                width = 300,

                shiny::textInput(
                  ns("mut_gene"),
                  "Gene (Mutation):",
                  value = "BRCA1",
                  placeholder = "e.g., BRCA1, TP53"
                ),

                shiny::selectInput(
                  ns("mut_drug"),
                  "Drug:",
                  choices = c("Select drug..." = ""),
                  selected = ""
                ),

                shiny::selectInput(
                  ns("mut_metric"),
                  "Sensitivity Metric:",
                  choices = c(
                    "AUC" = "auc",
                    "IC50" = "ic50"
                  ),
                  selected = "auc"
                ),

                shiny::hr(),

                shiny::actionButton(
                  ns("run_mutation_analysis"),
                  "Run Analysis",
                  icon = shiny::icon("play"),
                  class = "btn-primary w-100"
                )
              ),

              bslib::navset_card_tab(
                bslib::nav_panel(
                  title = shiny::icon("chart-bar"),
                  shiny::plotOutput(ns("mutation_plot"), height = "500px")
                ),
                bslib::nav_panel(
                  title = shiny::icon("table"),
                  DT::dataTableOutput(ns("mutation_table"))
                ),
                bslib::nav_panel(
                  title = shiny::icon("info-circle"),
                  shiny::verbatimTextOutput(ns("mutation_stats"))
                )
              )
            )
          ),

          # Batch Gene-Drug Analysis Tab
          bslib::nav_panel(
            title = shiny::tagList(shiny::icon("list"), "Batch Analysis"),

            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                open = TRUE,
                width = 300,

                shiny::textAreaInput(
                  ns("batch_genes"),
                  "Genes (one per line):",
                  value = "TP53\nBRCA1\nEGFR\nMYC",
                  rows = 6
                ),

                shiny::selectInput(
                  ns("batch_drug"),
                  "Drug:",
                  choices = c("Select drug..." = ""),
                  selected = ""
                ),

                shiny::selectInput(
                  ns("batch_metric"),
                  "Sensitivity Metric:",
                  choices = c(
                    "AUC" = "auc",
                    "IC50" = "ic50"
                  ),
                  selected = "auc"
                ),

                shiny::hr(),

                shiny::actionButton(
                  ns("run_batch"),
                  "Run Batch Analysis",
                  icon = shiny::icon("play"),
                  class = "btn-primary w-100"
                )
              ),

              bslib::navset_card_tab(
                bslib::nav_panel(
                  title = shiny::icon("chart-bar"),
                  shiny::plotOutput(ns("batch_plot"), height = "500px")
                ),
                bslib::nav_panel(
                  title = shiny::icon("table"),
                  DT::dataTableOutput(ns("batch_table"))
                )
              )
            )
          )
        )
      )
    )
  )
}

#' PharmacoGenomics Module Server
#'
#' @param id Module ID
#' @param app_state Reactive values for shared state
#' @param async_compute Async compute engine
#' @export
mod_pharmacogenomics_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Update drug choices based on data source
    shiny::observeEvent(input$data_source, {
      drugs <- get_available_drugs(source = input$data_source)
      shiny::updateSelectInput(session, "drug_id", choices = c("Select drug..." = "", drugs))
    })

    shiny::observeEvent(input$data_source, {
      drugs <- get_available_drugs(source = input$data_source)
      shiny::updateSelectInput(session, "mut_drug", choices = c("Select drug..." = "", drugs))
      shiny::updateSelectInput(session, "batch_drug", choices = c("Select drug..." = "", drugs))
    })

    # Drug-Gene Correlation Analysis
    correlation_results <- shiny::eventReactive(input$run_correlation, {
      shiny::req(input$gene_id, input$drug_id)

      shiny::withProgress(message = 'Analyzing drug-gene correlation...', value = 0, {
        shiny::incProgress(0.3, detail = "Querying data")

        result <- analyze_drug_gene_cor(
          gene = input$gene_id,
          drug = input$drug_id,
          metric = input$sensitivity_metric,
          cor_method = input$cor_method
        )

        shiny::incProgress(1, detail = "Complete")
        result
      })
    })

    output$correlation_plot <- shiny::renderPlot({
      shiny::req(correlation_results())

      result <- correlation_results()
      data <- result$data

      cor_label <- sprintf(
        "%s = %.3f, p = %.2e, n = %d",
        ifelse(input$cor_method == "spearman", "rho", "r"),
        result$correlation$estimate,
        result$correlation$pvalue,
        result$n_samples
      )

      ggplot2::ggplot(data, ggplot2::aes(x = Gene, y = Drug)) +
        ggplot2::geom_point(alpha = 0.6, color = "steelblue", size = 2) +
        ggplot2::geom_smooth(method = "lm", color = "red", se = TRUE) +
        ggplot2::labs(
          title = paste(input$gene_id, "Expression vs", input$drug_id, "Sensitivity"),
          subtitle = cor_label,
          x = paste(input$gene_id, "Expression"),
          y = paste(input$drug_id, toupper(input$sensitivity_metric))
        ) +
        theme_zinasuite()
    })

    output$correlation_table <- DT::renderDataTable({
      shiny::req(correlation_results())

      result <- correlation_results()
      data <- result$data
      data$Gene <- round(data$Gene, 4)
      data$Drug <- round(data$Drug, 4)

      DT::datatable(
        data,
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$correlation_stats <- shiny::renderPrint({
      shiny::req(correlation_results())

      result <- correlation_results()
      cor <- result$correlation

      cat("Drug-Gene Correlation Analysis\n")
      cat("==============================\n\n")
      cat("Gene:", result$gene, "\n")
      cat("Drug:", result$drug, "\n")
      cat("Metric:", toupper(result$metric), "\n")
      cat("Method:", cor$method, "\n\n")
      cat("Results:\n")
      cat("  Correlation coefficient:", round(cor$estimate, 4), "\n")
      cat("  P-value:", format(cor$pvalue, digits = 4), "\n")
      cat("  Sample size:", result$n_samples, "\n")
    })

    # Drug-Mutation Analysis
    mutation_results <- shiny::eventReactive(input$run_mutation_analysis, {
      shiny::req(input$mut_gene, input$mut_drug)

      shiny::withProgress(message = 'Analyzing drug-mutation association...', value = 0, {
        shiny::incProgress(0.3, detail = "Querying data")

        result <- analyze_drug_mutation(
          drug = input$mut_drug,
          gene = input$mut_gene,
          metric = input$mut_metric
        )

        shiny::incProgress(1, detail = "Complete")
        result
      })
    })

    output$mutation_plot <- shiny::renderPlot({
      shiny::req(mutation_results())

      result <- mutation_results()
      data <- result$data

      ggplot2::ggplot(data, ggplot2::aes(x = Mutated, y = Sensitivity, fill = Mutated)) +
        ggplot2::geom_boxplot(alpha = 0.7) +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
        ggplot2::scale_fill_manual(values = c("Wild-type" = "steelblue", "Mutated" = "coral")) +
        ggplot2::labs(
          title = paste(input$mut_drug, "Sensitivity by", input$mut_gene, "Status"),
          subtitle = paste("Wilcoxon p-value:", format(result$test$p.value, digits = 4)),
          x = paste(input$mut_gene, "Status"),
          y = paste(toupper(input$mut_metric))
        ) +
        theme_zinasuite() +
        ggplot2::theme(legend.position = "none")
    })

    output$mutation_table <- DT::renderDataTable({
      shiny::req(mutation_results())

      result <- mutation_results()

      summary_df <- data.frame(
        Group = c("Wild-type", "Mutated"),
        N = c(result$summary$wt_n, result$summary$mut_n),
        Median_Sensitivity = c(
          round(result$summary$wt_median, 4),
          round(result$summary$mut_median, 4)
        ),
        P_value = c(format(result$test$p.value, digits = 4), NA)
      )

      DT::datatable(
        summary_df,
        options = list(pageLength = 10, dom = 't'),
        rownames = FALSE
      )
    })

    output$mutation_stats <- shiny::renderPrint({
      shiny::req(mutation_results())

      result <- mutation_results()

      cat("Drug-Mutation Association Analysis\n")
      cat("==================================\n\n")
      cat("Drug:", result$drug, "\n")
      cat("Gene:", result$gene, "\n")
      cat("Metric:", toupper(result$metric), "\n\n")
      cat("Summary:\n")
      cat("  Wild-type median:", round(result$summary$wt_median, 4), "\n")
      cat("  Mutated median:", round(result$summary$mut_median, 4), "\n")
      cat("  Wild-type n:", result$summary$wt_n, "\n")
      cat("  Mutated n:", result$summary$mut_n, "\n\n")
      cat("Wilcoxon test p-value:", format(result$test$p.value, digits = 4), "\n")
    })

    # Batch Analysis
    batch_results <- shiny::eventReactive(input$run_batch, {
      shiny::req(input$batch_genes, input$batch_drug)

      genes <- strsplit(input$batch_genes, "\n")[[1]]
      genes <- trimws(genes)
      genes <- genes[genes != ""]

      if (length(genes) < 2) {
        shiny::showNotification("Please enter at least 2 genes", type = "error")
        return(NULL)
      }

      shiny::withProgress(message = 'Running batch analysis...', value = 0, {
        result <- analyze_drug_gene_batch(
          genes = genes,
          drug = input$batch_drug,
          metric = input$batch_metric,
          cor_method = "spearman",
          n_workers = 2,
          .progress = TRUE
        )

        shiny::incProgress(1, detail = "Complete")
        result
      })
    })

    output$batch_plot <- shiny::renderPlot({
      shiny::req(batch_results())

      result <- batch_results()
      result$significant <- result$padj < 0.05

      ggplot2::ggplot(result, ggplot2::aes(x = stats::reorder(gene, cor), y = cor, fill = significant)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
        ggplot2::labs(
          title = paste("Gene-Drug Correlations:", input$batch_drug),
          subtitle = paste("Metric:", toupper(input$batch_metric)),
          x = "Gene",
          y = "Correlation Coefficient",
          fill = "Significant (FDR < 0.05)"
        ) +
        theme_zinasuite()
    })

    output$batch_table <- DT::renderDataTable({
      shiny::req(batch_results())

      result <- batch_results()
      result$cor <- round(result$cor, 4)
      result$pvalue <- format(result$pvalue, digits = 4)
      result$padj <- format(result$padj, digits = 4)

      DT::datatable(
        result,
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

  })
}
