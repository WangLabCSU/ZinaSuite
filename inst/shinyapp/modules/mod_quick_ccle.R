' CCLE Quick Analysis Module UI
'
' @description
' Quick analysis module for CCLE (Cancer Cell Line Encyclopedia) data,
' providing analysis of cell line data including distribution, correlation,
' and drug sensitivity analysis.
'
' @param id Module ID
' @return UI elements
' @export
mod_quick_ccle_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::card(
      bslib::card_header(
        class = "bg-warning text-dark",
        shiny::icon("bolt"), "CCLE Quick Analysis"
      ),
      bslib::card_body(
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            open = TRUE,
            width = 350,

            ' Analysis Type Selection
            shiny::selectInput(
              ns("analysis_type"),
              "Analysis Type:",
              choices = c(
                "Expression Distribution" = "distribution",
                "Gene Correlation" = "correlation",
                "Drug Sensitivity" = "drug_sensitivity",
                "Drug Response" = "drug_response"
              ),
              selected = "distribution"
            ),

            shiny::hr(),

            ' Gene Input
            shiny::textInput(
              ns("gene"),
              "Gene Symbol:",
              value = "TP53",
              placeholder = "e.g., TP53, BRCA1, EGFR"
            ),

            ' Data Type Selection
            shiny::selectInput(
              ns("data_type"),
              "Data Type:",
              choices = c(
                "mRNA" = "mRNA",
                "Protein" = "protein",
                "CNV" = "cnv",
                "Mutation" = "mutation"
              ),
              selected = "mRNA"
            ),

            ' Conditional panels for different analysis types
            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'correlation'", ns("analysis_type")),
              shiny::textInput(
                ns("gene2"),
                "Second Gene:",
                value = "BRCA1",
                placeholder = "e.g., BRCA1, EGFR"
              )
            ),

            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'drug_sensitivity' || input['%s'] == 'drug_response'", ns("analysis_type"), ns("analysis_type")),
              shiny::selectInput(
                ns("drug"),
                "Drug:",
                choices = c(
                  "Cisplatin" = "cisplatin",
                  "Paclitaxel" = "paclitaxel",
                  "Doxorubicin" = "doxorubicin",
                  "Gemcitabine" = "gemcitabine"
                ),
                selected = "cisplatin"
              )
            ),

            shiny::hr(),

            ' Action Button
            shiny::actionButton(
              ns("run_btn"),
              "Run Analysis",
              icon = shiny::icon("play"),
              class = "btn-warning w-100"
            ),

            shiny::hr(),

            ' Download Options
            shiny::h5("Download Results"),
            shiny::downloadButton(
              ns("download_plot"),
              "Download Plot",
              class = "btn-outline-warning btn-sm w-100 mb-2"
            ),
            shiny::downloadButton(
              ns("download_data"),
              "Download Data",
              class = "btn-outline-secondary btn-sm w-100"
            )
          ),

          ' Main Content Area
          bslib::navset_card_tab(
            bslib::nav_panel(
              title = shiny::icon("chart-bar"),
              "Plot",
              shiny::plotOutput(ns("plot_output"), height = "600px")
            ),
            bslib::nav_panel(
              title = shiny::icon("table"),
              "Data Table",
              DT::dataTableOutput(ns("data_table"))
            ),
            bslib::nav_panel(
              title = shiny::icon("info-circle"),
              "Statistics",
              shiny::verbatimTextOutput(ns("stats_output"))
            )
          )
        )
      )
    )
  )
}

' CCLE Quick Analysis Module Server
'
' @param id Module ID
' @param app_state Shared reactive state
' @param async_compute Async compute engine
' @export
mod_quick_ccle_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ' Reactive values to store results
    results <- shiny::reactiveValues(
      plot = NULL,
      data = NULL,
      stats = NULL
    )

    ' Run analysis when button is clicked
    shiny::observeEvent(input$run_btn, {
      gene <- input$gene
      analysis_type <- input$analysis_type
      data_type <- input$data_type

      ' Validate input
      if (gene == "") {
        shiny::showNotification("Please enter a gene symbol", type = "error")
        return()
      }

      ' Show progress
      shiny::withProgress(message = "Running CCLE analysis...", value = 0, {
        tryCatch({
          result <- switch(analysis_type,
            "distribution" = {
              shiny::incProgress(0.3, detail = "Querying CCLE expression data")
              run_ccle_distribution(gene, data_type)
            },
            "correlation" = {
              shiny::incProgress(0.3, detail = "Querying gene data")
              run_ccle_correlation(gene, input$gene2, data_type)
            },
            "drug_sensitivity" = {
              shiny::incProgress(0.3, detail = "Analyzing drug sensitivity")
              run_ccle_drug_sensitivity(gene, data_type, input$drug)
            },
            "drug_response" = {
              shiny::incProgress(0.3, detail = "Analyzing drug response")
              run_ccle_drug_response(gene, data_type, input$drug)
            },
            stop("Unknown analysis type: ", analysis_type)
          )

          results$plot <- result$plot
          results$data <- result$data
          results$stats <- result$stats

          shiny::incProgress(1, detail = "Complete")
          shiny::showNotification("CCLE analysis complete!", type = "message")

        }, error = function(e) {
          shiny::showNotification(paste("Analysis failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    ' Output plot
    output$plot_output <- shiny::renderPlot({
      shiny::req(results$plot)
      print(results$plot)
    }, res = 100)

    ' Output data table
    output$data_table <- DT::renderDataTable({
      shiny::req(results$data)
      DT::datatable(results$data, options = list(pageLength = 10, scrollX = TRUE))
    })

    ' Output statistics
    output$stats_output <- shiny::renderPrint({
      shiny::req(results$stats)
      print(results$stats)
    })

    ' Download handlers
    output$download_plot <- shiny::downloadHandler(
      filename = function() {
        paste0(input$gene, "_ccle_", input$analysis_type, ".pdf")
      },
      content = function(file) {
        shiny::req(results$plot)
        ggplot2::ggsave(file, plot = results$plot, width = 12, height = 8, dpi = 300)
      }
    )

    output$download_data <- shiny::downloadHandler(
      filename = function() {
        paste0(input$gene, "_ccle_", input$analysis_type, "_data.csv")
      },
      content = function(file) {
        shiny::req(results$data)
        utils::write.csv(results$data, file, row.names = FALSE)
      }
    )
  })
}

' CCLE Analysis Functions ----------------------------------------------------

' Run CCLE Distribution Analysis
' @keywords internal
run_ccle_distribution <- function(gene, data_type) {
  ' Query gene expression from CCLE
  gene_expr <- query_molecule(gene, data_type = data_type, source = "ccle")

  ' Load CCLE sample information
  sample_info <- load_data("ccle_info")

  ' Prepare data
  common_samples <- intersect(names(gene_expr), sample_info$sample)

  plot_data <- data.frame(
    Sample = common_samples,
    Expression = as.numeric(gene_expr[common_samples]),
    Line = sample_info$cell_line[match(common_samples, sample_info$sample)],
    Primary = sample_info$primary_site[match(common_samples, sample_info$sample)],
    Subtype = sample_info$subtype[match(common_samples, sample_info$sample)],
    stringsAsFactors = FALSE
  )

  ' Create plot by primary site
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Primary, y = .data$Expression)) +
    ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.7) +
    ggplot2::labs(
      title = paste(gene, "Expression in CCLE Cell Lines"),
      x = "Primary Site",
      y = paste(data_type, "Expression")
    ) +
    theme_zinasuite() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  stats <- summary(plot_data$Expression)

  list(plot = plot, data = plot_data, stats = stats)
}

' Run CCLE Correlation Analysis
' @keywords internal
run_ccle_correlation <- function(gene1, gene2, data_type) {
  ' Query gene expressions from CCLE
  expr1 <- query_molecule(gene1, data_type = data_type, source = "ccle")
  expr2 <- query_molecule(gene2, data_type = data_type, source = "ccle")

  ' Find common samples
  common_samples <- intersect(names(expr1), names(expr2))

  if (length(common_samples) < 10) {
    stop("Insufficient common samples for correlation analysis")
  }

  ' Prepare data
  plot_data <- data.frame(
    Gene1 = as.numeric(expr1[common_samples]),
    Gene2 = as.numeric(expr2[common_samples]),
    stringsAsFactors = FALSE
  )

  ' Calculate correlation
  cor_test <- stats::cor.test(plot_data$Gene1, plot_data$Gene2)

  ' Create plot
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Gene1, y = .data$Gene2)) +
    ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
    ggplot2::geom_smooth(method = "lm", color = "red") +
    ggplot2::labs(
      title = paste(gene1, "vs", gene2, "(CCLE)"),
      subtitle = sprintf("r = %.3f, p = %.2e", cor_test$estimate, cor_test$p.value),
      x = paste(gene1, "Expression"),
      y = paste(gene2, "Expression")
    ) +
    theme_zinasuite()

  stats <- list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    n_samples = nrow(plot_data)
  )

  list(plot = plot, data = plot_data, stats = stats)
}

' Run CCLE Drug Sensitivity Analysis
' @keywords internal
run_ccle_drug_sensitivity <- function(gene, data_type, drug) {
  ' Use the new drug response analysis function
  ' This creates a volcano plot showing gene-drug associations

  plot <- vis_gene_drug_response_asso(
    Gene = gene,
    x_axis_type = "mean.diff",
    output_form = "ggplot2"
  )

  ' Get the underlying data
  data <- analyze_gene_drug_response_asso(gene)

  stats <- list(
    n_drugs = nrow(data),
    significant_drugs = sum(data$p.value < 0.05),
    top_correlation = data$cor[1],
    top_drug = data$drugs[1]
  )

  list(plot = plot, data = data, stats = stats)
}

' Run CCLE Drug Response Analysis
' @keywords internal
run_ccle_drug_response <- function(gene, data_type, drug) {
  ' Use the new drug response difference function
  ' This creates a dot plot comparing IC50 between high/low expression groups

  ' Get primary sites
  ccle_data <- load_data("ccle_expr_and_drug_response")
  available_sites <- unique(ccle_data$drug_info$`Site Primary`)
  selected_site <- available_sites[1]  ' Use first available site

  plot <- vis_gene_drug_response_diff(
    Gene = gene,
    tissue = selected_site,
    Show.P.label = TRUE,
    Method = "wilcox.test"
  )

  ' Get the underlying data
  data <- analyze_gene_drug_response_diff(gene, tissue = selected_site)

  stats <- list(
    n_comparisons = length(unique(data$drug_target)),
    n_cell_lines = length(unique(data$ccle_name)),
    selected_site = selected_site
  )

  list(plot = plot, data = data, stats = stats)
}
