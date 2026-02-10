#' TCGA Quick Analysis Module UI
#'
#' @description
#' Quick analysis module for TCGA data, providing fast access to common analyses
#' including tumor vs normal comparison, correlation, survival, and more.
#'
#' @param id Module ID
#' @return UI elements
#' @export
mod_quick_tcga_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::card(
      bslib::card_header(
        class = "bg-primary text-white",
        shiny::icon("bolt"), "TCGA Quick Analysis"
      ),
      bslib::card_body(
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            open = TRUE,
            width = 350,

            # Analysis Type Selection
            shiny::selectInput(
              ns("analysis_type"),
              "Analysis Type:",
              choices = c(
                "Tumor vs Normal" = "tumor_normal",
                "Anatomy Map" = "anatomy",
                "Gene Correlation" = "correlation",
                "TIL Correlation" = "til",
            "Immune Correlation" = "immune",
            "TMB Correlation" = "tmb",
            "MSI Correlation" = "msi",
            "Stemness Correlation" = "stemness",
            "Pathway Analysis" = "pathway",
            "Mutation Frequency" = "mutation",
                "Survival (KM)" = "survival_km",
                "Survival (Cox)" = "survival_cox",
                "Dimension Reduction" = "dimension"
              ),
              selected = "tumor_normal"
            ),

            shiny::hr(),

            # Gene Input
            shiny::textInput(
              ns("gene"),
              "Gene Symbol:",
              value = "TP53",
              placeholder = "e.g., TP53, BRCA1, EGFR"
            ),

            # Data Type Selection
            shiny::selectInput(
              ns("data_type"),
              "Data Type:",
              choices = c(
                "mRNA" = "mRNA",
                "Transcript" = "transcript",
                "Protein" = "protein",
                "miRNA" = "miRNA",
                "Methylation" = "methylation",
                "CNV" = "cnv"
              ),
              selected = "mRNA"
            ),

            # Conditional panels for different analysis types
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
              condition = sprintf("input['%s'] == 'tumor_normal'", ns("analysis_type")),
              shiny::radioButtons(
                ns("mode"),
                "Analysis Mode:",
                choices = c("Pan-cancer" = "pancan", "Single Cancer" = "single"),
                selected = "pancan",
                inline = TRUE
              ),
              shiny::conditionalPanel(
                condition = sprintf("input['%s'] == 'single'", ns("mode")),
                shiny::selectInput(
                  ns("cancer"),
                  "Select Cancer:",
                  choices = c("BRCA", "LUAD", "LUSC", "COAD", "READ", "HNSC", "STAD"),
                  selected = "BRCA"
                )
              ),
              shiny::checkboxInput(
                ns("include_gtex"),
                "Include GTEx Normal Samples",
                value = TRUE
              ),
              shiny::radioButtons(
                ns("plot_type"),
                "Plot Type:",
                choices = c("Boxplot" = "boxplot", "Violin" = "violin"),
                selected = "boxplot",
                inline = TRUE
              )
            ),

            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'til' || input['%s'] == 'immune' || input['%s'] == 'tmb' || input['%s'] == 'stemness'", ns("analysis_type"), ns("analysis_type"), ns("analysis_type"), ns("analysis_type")),
              shiny::selectInput(
                ns("cor_method"),
                "Correlation Method:",
                choices = c("Spearman" = "spearman", "Pearson" = "pearson"),
                selected = "spearman"
              )
            ),

            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'tmb' || input['%s'] == 'msi' || input['%s'] == 'stemness'", ns("analysis_type"), ns("analysis_type"), ns("analysis_type")),
              shiny::radioButtons(
                ns("index_plot_type"),
                "Plot Type:",
                choices = c("Scatter" = "scatter", "Summary" = "summary"),
                selected = "scatter",
                inline = TRUE
              )
            ),

            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'pathway'", ns("analysis_type")),
              shiny::selectInput(
                ns("pathway_name"),
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
                ns("pw_cor_method"),
                "Correlation Method:",
                choices = c("Spearman" = "spearman", "Pearson" = "pearson"),
                selected = "spearman"
              ),
              shiny::checkboxInput(
                ns("use_regline"),
                "Show Regression Line",
                value = TRUE
              )
            ),

            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'survival_km' || input['%s'] == 'survival_cox'", ns("analysis_type"), ns("analysis_type"))),
              shiny::selectInput(
                ns("surv_measure"),
                "Survival Measure:",
                choices = c("Overall Survival" = "OS", "Progression-Free Survival" = "PFI"),
                selected = "OS"
              )
            ),

            shiny::hr(),

            # Action Button
            shiny::actionButton(
              ns("run_btn"),
              "Run Analysis",
              icon = shiny::icon("play"),
              class = "btn-primary w-100"
            ),

            shiny::hr(),

            # Download Options
            shiny::h5("Download Results"),
            shiny::downloadButton(
              ns("download_plot"),
              "Download Plot",
              class = "btn-outline-primary btn-sm w-100 mb-2"
            ),
            shiny::downloadButton(
              ns("download_data"),
              "Download Data",
              class = "btn-outline-secondary btn-sm w-100"
            )
          ),

          # Main Content Area
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
}

#' TCGA Quick Analysis Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
#' @param async_compute Async compute engine
#' @export
mod_quick_tcga_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive values to store results
    results <- shiny::reactiveValues(
      plot = NULL,
      data = NULL,
      stats = NULL
    )

    # Run analysis when button is clicked
    shiny::observeEvent(input$run_btn, {
      gene <- input$gene
      analysis_type <- input$analysis_type
      data_type <- input$data_type

      # Validate input
      if (gene == "") {
        shiny::showNotification("Please enter a gene symbol", type = "error")
        return()
      }

      # Show progress
      shiny::withProgress(message = "Running analysis...", value = 0, {
        tryCatch({
          result <- switch(analysis_type,
            "tumor_normal" = {
              shiny::incProgress(0.3, detail = "Querying expression data")
              run_tumor_normal_analysis(gene, data_type, input$mode, input$cancer, input$include_gtex, input$plot_type)
            },
            "anatomy" = {
              shiny::incProgress(0.3, detail = "Generating anatomy map")
              run_anatomy_analysis(gene, data_type)
            },
            "correlation" = {
              shiny::incProgress(0.3, detail = "Querying gene data")
              run_correlation_analysis(gene, input$gene2, data_type)
            },
            "til" = {
              shiny::incProgress(0.3, detail = "Querying TIL data")
              run_til_analysis(gene, data_type, input$cor_method)
            },
            "immune" = {
              shiny::incProgress(0.3, detail = "Querying immune data")
              run_immune_analysis(gene, data_type, input$cor_method)
            },
            "tmb" = {
              shiny::incProgress(0.3, detail = "Querying TMB data")
              run_index_analysis(gene, data_type, "tmb", input$cor_method, input$index_plot_type)
            },
            "msi" = {
              shiny::incProgress(0.3, detail = "Querying MSI data")
              run_index_analysis(gene, data_type, "msi", input$cor_method, input$index_plot_type)
            },
            "stemness" = {
              shiny::incProgress(0.3, detail = "Querying stemness data")
              run_index_analysis(gene, data_type, "stemness", input$cor_method, input$index_plot_type)
            },
            "pathway" = {
              shiny::incProgress(0.3, detail = "Querying pathway data")
              run_pathway_analysis(gene, data_type, input$pathway_name, input$cancer, input$pw_cor_method, input$use_regline)
            },
            "mutation" = {
              shiny::incProgress(0.3, detail = "Querying mutation data")
              run_mutation_analysis(gene, data_type)
            },
            "survival_km" = {
              shiny::incProgress(0.3, detail = "Querying survival data")
              run_survival_km_analysis(gene, data_type, input$surv_measure)
            },
            "survival_cox" = {
              shiny::incProgress(0.3, detail = "Running Cox analysis")
              run_survival_cox_analysis(gene, data_type, input$surv_measure)
            },
            "dimension" = {
              shiny::incProgress(0.3, detail = "Running dimension reduction")
              run_dimension_analysis(gene, data_type)
            },
            stop("Unknown analysis type: ", analysis_type)
          )

          results$plot <- result$plot
          results$data <- result$data
          results$stats <- result$stats

          shiny::incProgress(1, detail = "Complete")
          shiny::showNotification("Analysis complete!", type = "message")

        }, error = function(e) {
          shiny::showNotification(paste("Analysis failed:", conditionMessage(e)), type = "error")
        })
      })
    })

    # Output plot
    output$plot_output <- shiny::renderPlot({
      shiny::req(results$plot)
      print(results$plot)
    }, res = 100)

    # Output data table
    output$data_table <- DT::renderDataTable({
      shiny::req(results$data)
      DT::datatable(results$data, options = list(pageLength = 10, scrollX = TRUE))
    })

    # Output statistics
    output$stats_output <- shiny::renderPrint({
      shiny::req(results$stats)
      print(results$stats)
    })

    # Download handlers
    output$download_plot <- shiny::downloadHandler(
      filename = function() {
        paste0(input$gene, "_", input$analysis_type, ".pdf")
      },
      content = function(file) {
        shiny::req(results$plot)
        ggplot2::ggsave(file, plot = results$plot, width = 12, height = 8, dpi = 300)
      }
    )

    output$download_data <- shiny::downloadHandler(
      filename = function() {
        paste0(input$gene, "_", input$analysis_type, "_data.csv")
      },
      content = function(file) {
        shiny::req(results$data)
        utils::write.csv(results$data, file, row.names = FALSE)
      }
    )
  })
}

# Analysis Functions ---------------------------------------------------------

#' Run Tumor vs Normal Analysis
#' @keywords internal
run_tumor_normal_analysis <- function(gene, data_type, mode, cancer, include_gtex, plot_type) {
  # Query gene expression
  gene_expr <- query_molecule(gene, data_type = data_type, source = "tcga")

  # Load sample information
  sample_info <- load_data("tcga_gtex")

  # Prepare data
  plot_data <- prepare_tumor_normal_data(gene_expr, sample_info, include_gtex)

  if (mode == "single" && !is.null(cancer)) {
    plot_data <- plot_data[plot_data$Cancer == cancer, ]
  }

  # Create plot
  plot_fn <- if (plot_type == "boxplot") create_tumor_normal_boxplot else create_tumor_normal_violin
  plot <- plot_fn(plot_data, gene, data_type)

  # Calculate statistics
  stats <- calculate_tumor_normal_stats(plot_data)

  list(plot = plot, data = plot_data, stats = stats)
}

#' Prepare Tumor vs Normal Data
#' @keywords internal
prepare_tumor_normal_data <- function(gene_expr, sample_info, include_gtex) {
  common_samples <- intersect(names(gene_expr), sample_info$sample)

  data.frame(
    Sample = common_samples,
    Expression = as.numeric(gene_expr[common_samples]),
    Type = sample_info$type2[match(common_samples, sample_info$sample)],
    Cancer = sample_info$tissue[match(common_samples, sample_info$sample)],
    Dataset = ifelse(substr(common_samples, 1, 4) == "TCGA", "TCGA", "GTEx"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(.data$Type %in% c("tumor", "normal")) %>%
    dplyr::filter(include_gtex | .data$Dataset == "TCGA")
}

#' Create Tumor vs Normal Boxplot
#' @keywords internal
create_tumor_normal_boxplot <- function(data, gene, data_type) {
  ggplot2::ggplot(data, ggplot2::aes(x = .data$Cancer, y = .data$Expression, fill = .data$Type)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    ggplot2::scale_fill_manual(values = c("tumor" = "#DF2020", "normal" = "#DDDF21")) +
    ggplot2::labs(
      title = paste(gene, "Expression: Tumor vs Normal"),
      x = "Cancer Type",
      y = paste(data_type, "Expression")
    ) +
    theme_zinasuite() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Create Tumor vs Normal Violin Plot
#' @keywords internal
create_tumor_normal_violin <- function(data, gene, data_type) {
  ggplot2::ggplot(data, ggplot2::aes(x = .data$Cancer, y = .data$Expression, fill = .data$Type)) +
    ggplot2::geom_violin(alpha = 0.7, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.1, alpha = 0.5) +
    ggplot2::scale_fill_manual(values = c("tumor" = "#DF2020", "normal" = "#DDDF21")) +
    ggplot2::labs(
      title = paste(gene, "Expression: Tumor vs Normal"),
      x = "Cancer Type",
      y = paste(data_type, "Expression")
    ) +
    theme_zinasuite() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Calculate Tumor vs Normal Statistics
#' @keywords internal
calculate_tumor_normal_stats <- function(data) {
  cancers <- unique(data$Cancer)
  stats <- lapply(cancers, function(cancer) {
    cancer_data <- data[data$Cancer == cancer, ]
    tumor_vals <- cancer_data$Expression[cancer_data$Type == "tumor"]
    normal_vals <- cancer_data$Expression[cancer_data$Type == "normal"]

    if (length(tumor_vals) > 0 && length(normal_vals) > 0) {
      test <- stats::wilcox.test(tumor_vals, normal_vals)
      data.frame(
        Cancer = cancer,
        Tumor_Median = stats::median(tumor_vals, na.rm = TRUE),
        Normal_Median = stats::median(normal_vals, na.rm = TRUE),
        P_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  })

  do.call(rbind, stats)
}

#' Run Anatomy Analysis
#' @keywords internal
run_anatomy_analysis <- function(gene, data_type) {
  plot <- vis_pancan_anatomy(gene, data_type = data_type, plot_type = "heatmap")
  data <- plot$data
  stats <- summary(data$Expression)

  list(plot = plot, data = data, stats = stats)
}

#' Run Correlation Analysis
#' @keywords internal
run_correlation_analysis <- function(gene1, gene2, data_type) {
  plot <- vis_gene_cor(gene1, gene2, source = "tcga")

  # Extract correlation statistics
  cor_test <- stats::cor.test(plot$data$Gene1, plot$data$Gene2)
  stats <- list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    n_samples = nrow(plot$data)
  )

  list(plot = plot, data = plot$data, stats = stats)
}

#' Run TIL Analysis
#' @keywords internal
run_til_analysis <- function(gene, data_type, cor_method) {
  plot <- vis_gene_TIL_cor(gene, data_type = data_type, method = cor_method)
  data <- plot$data
  stats <- summary(data$Correlation)

  list(plot = plot, data = data, stats = stats)
}

#' Run Immune Analysis
#' @keywords internal
run_immune_analysis <- function(gene, data_type, cor_method) {
  plot <- vis_gene_immune_cor(gene, data_type = data_type, method = cor_method)
  data <- plot$data
  stats <- summary(data$Correlation)

  list(plot = plot, data = data, stats = stats)
}

#' Run Index Analysis (TMB/MSI/Stemness)
#' @keywords internal
run_index_analysis <- function(gene, data_type, index_type, cor_method, plot_type) {
  # Map index type to visualization function
  plot <- switch(index_type,
    "tmb" = vis_gene_tmb_cor(gene, data_type = data_type, method = cor_method),
    "msi" = vis_gene_msi_cor(gene, data_type = data_type, method = cor_method),
    "stemness" = vis_gene_stemness_cor(gene, data_type = data_type, method = cor_method),
    stop("Unknown index type: ", index_type)
  )

  data <- plot$data

  # Calculate correlation statistics
  cor_test <- stats::cor.test(data$Gene, data$Index, method = cor_method)

  stats <- list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    n_samples = nrow(data),
    index_type = index_type
  )

  list(plot = plot, data = data, stats = stats)
}

#' Run Pathway Analysis
#' @keywords internal
run_pathway_analysis <- function(gene, data_type, pathway_name, cancer, cor_method, use_regline) {
  # Use vis_gene_pw_cor for visualization
  plot <- vis_gene_pw_cor(
    gene = gene,
    data_type = data_type,
    pw_name = pathway_name,
    cancer_choose = cancer,
    cor_method = cor_method,
    use_regline = use_regline
  )

  # Extract data from plot attribute
  data <- attr(plot, "data")

  # Calculate correlation statistics
  cor_test <- stats::cor.test(data$Gene, data$Pathway, method = cor_method)

  stats <- list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    method = cor_method,
    n_samples = nrow(data),
    pathway = pathway_name
  )

  list(plot = plot, data = data, stats = stats)
}

#' Run Mutation Analysis
#' @keywords internal
run_mutation_analysis <- function(gene, data_type) {
  # Get mutation data
  mut_data <- query_mutation(gene, source = "tcga")

  # Calculate frequencies by cancer
  sample_info <- load_data("tcga_gtex")
  plot_data <- data.frame(
    Sample = names(mut_data),
    Mutation = as.character(mut_data),
    Cancer = sample_info$tissue[match(names(mut_data), sample_info$sample)],
    stringsAsFactors = FALSE
  )

  # Calculate mutation frequency by cancer
  mut_freq <- plot_data %>%
    dplyr::group_by(.data$Cancer) %>%
    dplyr::summarise(
      Mutation_Rate = mean(.data$Mutation != "WT", na.rm = TRUE) * 100,
      N = dplyr::n(),
      .groups = "drop"
    )

  # Create plot
  plot <- ggplot2::ggplot(mut_freq, ggplot2::aes(x = stats::reorder(.data$Cancer, -.data$Mutation_Rate), y = .data$Mutation_Rate)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(
      title = paste(gene, "Mutation Frequency by Cancer"),
      x = "Cancer Type",
      y = "Mutation Rate (%)"
    ) +
    theme_zinasuite() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  list(plot = plot, data = mut_freq, stats = summary(mut_freq$Mutation_Rate))
}

#' Run Survival KM Analysis
#' @keywords internal
run_survival_km_analysis <- function(gene, data_type, measure) {
  # Query gene expression
  gene_expr <- query_molecule(gene, data_type = data_type, source = "tcga")

  # Get survival data
  surv_data <- load_data("tcga_surv")

  # Prepare data
  plot_data <- data.frame(
    Sample = names(gene_expr),
    Expression = as.numeric(gene_expr),
    stringsAsFactors = FALSE
  )

  plot_data$Time <- surv_data[[paste0(measure, ".time")]][match(plot_data$Sample, surv_data$sample)]
  plot_data$Status <- surv_data[[measure]][match(plot_data$Sample, surv_data$sample)]
  plot_data <- plot_data[stats::complete.cases(plot_data), ]

  # Create expression groups
  plot_data$Group <- ifelse(
    plot_data$Expression > stats::median(plot_data$Expression, na.rm = TRUE),
    "High", "Low"
  )

  # Create survival plot
  fit <- survival::survfit(survival::Surv(Time, Status) ~ Group, data = plot_data)
  plot <- survminer::ggsurvplot(
    fit,
    data = plot_data,
    pval = TRUE,
    risk.table = TRUE,
    title = paste(gene, "Survival Analysis"),
    xlab = "Time (days)",
    ylab = "Survival Probability"
  )

  list(plot = plot$plot, data = plot_data, stats = fit)
}

#' Run Survival Cox Analysis
#' @keywords internal
run_survival_cox_analysis <- function(gene, data_type, measure) {
  # Use the existing vis_unicox_tree function
  plot <- vis_unicox_tree(gene, measure = measure, data_type = data_type)

  list(plot = plot, data = NULL, stats = NULL)
}

#' Run Dimension Analysis
#' @keywords internal
run_dimension_analysis <- function(gene, data_type) {
  # Query gene expression
  gene_expr <- query_molecule(gene, data_type = data_type, source = "tcga")

  # Get sample info
  sample_info <- load_data("tcga_gtex")

  # Prepare data for PCA
  # For single gene, we'll show distribution across cancers
  plot_data <- data.frame(
    Sample = names(gene_expr),
    Expression = as.numeric(gene_expr),
    Cancer = sample_info$tissue[match(names(gene_expr), sample_info$sample)],
    Type = sample_info$type2[match(names(gene_expr), sample_info$sample)],
    stringsAsFactors = FALSE
  )

  plot_data <- plot_data[stats::complete.cases(plot_data), ]

  # Create dimension plot (t-SNE or PCA placeholder)
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Cancer, y = .data$Expression, color = .data$Type)) +
    ggplot2::geom_jitter(alpha = 0.3, width = 0.2) +
    ggplot2::labs(
      title = paste(gene, "Expression Distribution"),
      x = "Cancer Type",
      y = "Expression"
    ) +
    theme_zinasuite() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  list(plot = plot, data = plot_data, stats = summary(plot_data$Expression))
}
