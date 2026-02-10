#' PCAWG Pipeline Analysis Module UI
#'
#' @description
#' Deep analysis pipeline for PCAWG data, supporting correlation, comparison,
#' and survival analysis with multiple modes (one-to-one, one-to-many, many-to-one).
#'
#' @param id Module ID
#' @return UI elements
#' @export
mod_pcawg_pipeline_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::card(
      bslib::card_header(
        class = "bg-info text-white",
        shiny::icon("microscope"), "PCAWG Deep Analysis Pipeline"
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
                "Correlation" = "cor",
                "Comparison" = "comp",
                "Survival" = "sur"
              ),
              selected = "cor"
            ),

            # Mode Selection
            shiny::selectInput(
              ns("mode"),
              "Analysis Mode:",
              choices = c(
                "One-to-One (Single Cancer)" = "o2o",
                "One-to-Many (Multiple Cancers)" = "o2m",
                "Many-to-One (Multiple Molecules)" = "m2o"
              ),
              selected = "o2o"
            ),

            shiny::hr(),

            # Cancer Selection
            shiny::selectInput(
              ns("cancer"),
              "Cancer Type:",
              choices = NULL,  # Will be populated server-side
              selected = NULL
            ),

            # Molecule Input
            shiny::textInput(
              ns("molecule"),
              "Molecule (Gene/Feature):",
              value = "TP53",
              placeholder = "e.g., TP53, BRCA1"
            ),

            # Conditional: Second molecule for correlation
            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'cor'", ns("analysis_type")),
              shiny::textInput(
                ns("molecule2"),
                "Second Molecule:",
                value = "KRAS",
                placeholder = "e.g., KRAS, EGFR"
              )
            ),

            # Conditional: Multiple molecules for m2o mode
            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'm2o'", ns("mode")),
              shiny::textAreaInput(
                ns("molecules"),
                "Molecules (one per line):",
                value = "TP53\nKRAS\nPTEN\nBRCA1",
                rows = 5
              )
            ),

            # Data Type
            shiny::selectInput(
              ns("data_type"),
              "Data Type:",
              choices = c(
                "mRNA Expression" = "mRNA",
                "miRNA Expression" = "miRNA",
                "Fusion" = "fusion",
                "Promoter" = "promoter"
              ),
              selected = "mRNA"
            ),

            shiny::hr(),

            # Analysis Parameters
            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'cor'", ns("analysis_type")),
              shiny::selectInput(
                ns("cor_method"),
                "Correlation Method:",
                choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                selected = "spearman"
              )
            ),

            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'sur'", ns("analysis_type")),
              shiny::selectInput(
                ns("surv_measure"),
                "Survival Measure:",
                choices = c("Overall Survival" = "OS", "Progression-Free Survival" = "PFI"),
                selected = "OS"
              ),
              shiny::selectInput(
                ns("cutoff_mode"),
                "Cutoff Mode:",
                choices = c("Median" = "median", "Custom" = "custom"),
                selected = "median"
              )
            ),

            shiny::hr(),

            # Action Button
            shiny::actionButton(
              ns("run_btn"),
              "Run Analysis",
              icon = shiny::icon("play"),
              class = "btn-info w-100"
            ),

            shiny::hr(),

            # Download Options
            shiny::h5("Download Results"),
            shiny::downloadButton(
              ns("download_plot"),
              "Download Plot",
              class = "btn-outline-info btn-sm w-100 mb-2"
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
              title = shiny::icon("chart-line"),
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

#' PCAWG Pipeline Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
#' @param async_compute Async compute engine
#' @export
mod_pcawg_pipeline_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive values to store results
    results <- shiny::reactiveValues(
      plot = NULL,
      data = NULL,
      stats = NULL
    )

    # Load PCAWG cancer types
    shiny::observe({
      pcawg_cancers <- load_data("pcawg_cancers")
      shiny::updateSelectInput(
        session, "cancer",
        choices = pcawg_cancers,
        selected = pcawg_cancers[1]
      )
    })

    # Run analysis when button is clicked
    shiny::observeEvent(input$run_btn, {
      analysis_type <- input$analysis_type
      mode <- input$mode
      cancer <- input$cancer
      molecule <- input$molecule
      data_type <- input$data_type

      # Validate input
      if (molecule == "") {
        shiny::showNotification("Please enter a molecule identifier", type = "error")
        return()
      }

      # Show progress
      shiny::withProgress(message = "Running PCAWG pipeline analysis...", value = 0, {
        tryCatch({
          result <- switch(analysis_type,
            "cor" = {
              shiny::incProgress(0.3, detail = "Running correlation analysis")
              switch(mode,
                "o2o" = run_pcawg_cor_o2o(molecule, input$molecule2, cancer, data_type, input$cor_method),
                "o2m" = run_pcawg_cor_o2m(molecule, input$molecule2, data_type, input$cor_method),
                "m2o" = run_pcawg_cor_m2o(input$molecules, cancer, data_type, input$cor_method),
                stop("Unknown mode: ", mode)
              )
            },
            "comp" = {
              shiny::incProgress(0.3, detail = "Running comparison analysis")
              switch(mode,
                "o2o" = run_pcawg_comp_o2o(molecule, cancer, data_type),
                "o2m" = run_pcawg_comp_o2m(molecule, data_type),
                "m2o" = run_pcawg_comp_m2o(input$molecules, cancer, data_type),
                stop("Unknown mode: ", mode)
              )
            },
            "sur" = {
              shiny::incProgress(0.3, detail = "Running survival analysis")
              switch(mode,
                "o2o" = run_pcawg_sur_o2o(molecule, cancer, data_type, input$surv_measure, input$cutoff_mode),
                "o2m" = run_pcawg_sur_o2m(molecule, data_type, input$surv_measure, input$cutoff_mode),
                "m2o" = run_pcawg_sur_m2o(input$molecules, cancer, data_type, input$surv_measure, input$cutoff_mode),
                stop("Unknown mode: ", mode)
              )
            },
            stop("Unknown analysis type: ", analysis_type)
          )

          results$plot <- result$plot
          results$data <- result$data
          results$stats <- result$stats

          shiny::incProgress(1, detail = "Complete")
          shiny::showNotification("PCAWG pipeline analysis complete!", type = "message")

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
        paste0(input$molecule, "_pcawg_", input$analysis_type, "_", input$mode, ".pdf")
      },
      content = function(file) {
        shiny::req(results$plot)
        ggplot2::ggsave(file, plot = results$plot, width = 12, height = 8, dpi = 300)
      }
    )

    output$download_data <- shiny::downloadHandler(
      filename = function() {
        paste0(input$molecule, "_pcawg_", input$analysis_type, "_", input$mode, "_data.csv")
      },
      content = function(file) {
        shiny::req(results$data)
        utils::write.csv(results$data, file, row.names = FALSE)
      }
    )
  })
}

# PCAWG Correlation Analysis Functions ----------------------------------------

#' Run PCAWG Correlation Analysis - One-to-One
#' @keywords internal
run_pcawg_cor_o2o <- function(molecule1, molecule2, cancer, data_type, cor_method) {
  # Query molecule data
  data1 <- query_molecule(molecule1, data_type = data_type, source = "pcawg")
  data2 <- query_molecule(molecule2, data_type = data_type, source = "pcawg")

  # Get sample info for cancer filtering
  sample_info <- load_data("pcawg_info")

  # Filter by cancer type
  cancer_samples <- sample_info$Sample[sample_info$project == cancer]
  common_samples <- intersect(intersect(names(data1), names(data2)), cancer_samples)

  if (length(common_samples) < 10) {
    stop("Insufficient samples for correlation analysis (need at least 10)")
  }

  # Prepare data
  plot_data <- data.frame(
    Sample = common_samples,
    X = as.numeric(data1[common_samples]),
    Y = as.numeric(data2[common_samples]),
    stringsAsFactors = FALSE
  )

  # Calculate correlation
  cor_test <- stats::cor.test(plot_data$X, plot_data$Y, method = cor_method)

  # Create plot
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$X, y = .data$Y)) +
    ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
    ggplot2::geom_smooth(method = "lm", color = "red", se = TRUE) +
    ggplot2::labs(
      title = paste(molecule1, "vs", molecule2, "in", cancer),
      subtitle = sprintf("%s r = %.3f, p = %.2e", cor_method, cor_test$estimate, cor_test$p.value),
      x = paste(molecule1, "Expression"),
      y = paste(molecule2, "Expression")
    ) +
    theme_zinasuite()

  stats <- list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    method = cor_method,
    n_samples = nrow(plot_data)
  )

  list(plot = plot, data = plot_data, stats = stats)
}

#' Run PCAWG Correlation Analysis - One-to-Many
#' @keywords internal
run_pcawg_cor_o2m <- function(molecule1, molecule2, data_type, cor_method) {
  # Query data for all PCAWG projects
  data1 <- query_molecule(molecule1, data_type = data_type, source = "pcawg")
  data2 <- query_molecule(molecule2, data_type = data_type, source = "pcawg")

  sample_info <- load_data("pcawg_info")

  # Find common samples
  common_samples <- intersect(names(data1), names(data2))

  # Prepare data with cancer info
  plot_data <- data.frame(
    Sample = common_samples,
    X = as.numeric(data1[common_samples]),
    Y = as.numeric(data2[common_samples]),
    Cancer = sample_info$project[match(common_samples, sample_info$Sample)],
    stringsAsFactors = FALSE
  )

  # Calculate correlation by cancer
  cor_results <- plot_data |>
    dplyr::group_by(.data$Cancer) |>
    dplyr::summarise(
      Correlation = stats::cor(.data$X, .data$Y, method = cor_method),
      P_value = tryCatch(stats::cor.test(.data$X, .data$Y, method = cor_method)$p.value,
                        error = function(e) NA),
      N = dplyr::n(),
      .groups = "drop"
    )

  # Create forest plot
  plot <- ggplot2::ggplot(cor_results, ggplot2::aes(x = stats::reorder(.data$Cancer, .data$Correlation), y = .data$Correlation)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$N, color = .data$P_value < 0.05)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(molecule1, "vs", molecule2, "across PCAWG cancers"),
      x = "Cancer Type",
      y = "Correlation Coefficient",
      size = "Sample Size",
      color = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = cor_results, stats = summary(cor_results$Correlation))
}

#' Run PCAWG Correlation Analysis - Many-to-One
#' @keywords internal
run_pcawg_cor_m2o <- function(molecules_text, cancer, data_type, cor_method) {
  # Parse molecules
  molecules <- strsplit(molecules_text, "\n")[[1]]
  molecules <- trimws(molecules[molecules != ""])

  if (length(molecules) < 2) {
    stop("At least 2 molecules required for many-to-one analysis")
  }

  # Query reference molecule (first one)
  ref_data <- query_molecule(molecules[1], data_type = data_type, source = "pcawg")

  # Get sample info
  sample_info <- load_data("pcawg_info")
  cancer_samples <- sample_info$Sample[sample_info$project == cancer]

  # Calculate correlation for each molecule
  cor_results <- purrr::map_dfr(molecules, function(mol) {
    mol_data <- query_molecule(mol, data_type = data_type, source = "pcawg")
    common_samples <- intersect(intersect(names(ref_data), names(mol_data)), cancer_samples)

    if (length(common_samples) < 5) {
      return(data.frame(
        Molecule = mol,
        Correlation = NA,
        P_value = NA,
        N = length(common_samples)
      ))
    }

    cor_test <- stats::cor.test(
      as.numeric(ref_data[common_samples]),
      as.numeric(mol_data[common_samples]),
      method = cor_method
    )

    data.frame(
      Molecule = mol,
      Correlation = cor_test$estimate,
      P_value = cor_test$p.value,
      N = length(common_samples)
    )
  })

  # Create bar plot
  plot <- ggplot2::ggplot(cor_results, ggplot2::aes(x = stats::reorder(.data$Molecule, .data$Correlation), y = .data$Correlation)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = .data$P_value < 0.05)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "steelblue")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Correlation with", molecules[1], "in", cancer),
      x = "Molecule",
      y = "Correlation Coefficient",
      fill = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = cor_results, stats = summary(cor_results$Correlation))
}

# PCAWG Comparison Analysis Functions -----------------------------------------

#' Run PCAWG Comparison Analysis - One-to-One
#' @keywords internal
run_pcawg_comp_o2o <- function(molecule, cancer, data_type) {
  # Query molecule data
  data <- query_molecule(molecule, data_type = data_type, source = "pcawg")

  # Get sample info
  sample_info <- load_data("pcawg_info")

  # Filter by cancer
  cancer_samples <- sample_info$Sample[sample_info$project == cancer]
  common_samples <- intersect(names(data), cancer_samples)

  # Prepare data
  plot_data <- data.frame(
    Sample = common_samples,
    Expression = as.numeric(data[common_samples]),
    Type = sample_info$type[match(common_samples, sample_info$Sample)],
    stringsAsFactors = FALSE
  )

  # Remove NAs
  plot_data <- plot_data[!is.na(plot_data$Type), ]

  if (nrow(plot_data) < 10) {
    stop("Insufficient samples for comparison")
  }

  # Statistical test
  tumor_vals <- plot_data$Expression[plot_data$Type == "tumor"]
  normal_vals <- plot_data$Expression[plot_data$Type == "normal"]

  if (length(tumor_vals) > 0 && length(normal_vals) > 0) {
    test_result <- stats::wilcox.test(tumor_vals, normal_vals)
    p_value <- test_result$p.value
  } else {
    p_value <- NA
  }

  # Create plot
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Type, y = .data$Expression, fill = .data$Type)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.3) +
    ggplot2::scale_fill_manual(values = c("tumor" = "#DF2020", "normal" = "#DDDF21")) +
    ggplot2::labs(
      title = paste(molecule, "Expression in", cancer),
      subtitle = sprintf("Wilcoxon p = %.2e", p_value),
      x = "Sample Type",
      y = "Expression"
    ) +
    theme_zinasuite()

  stats <- list(
    p_value = p_value,
    tumor_median = stats::median(tumor_vals, na.rm = TRUE),
    normal_median = stats::median(normal_vals, na.rm = TRUE),
    fold_change = stats::median(tumor_vals, na.rm = TRUE) / stats::median(normal_vals, na.rm = TRUE)
  )

  list(plot = plot, data = plot_data, stats = stats)
}

#' Run PCAWG Comparison Analysis - One-to-Many
#' @keywords internal
run_pcawg_comp_o2m <- function(molecule, data_type) {
  # Query molecule data
  data <- query_molecule(molecule, data_type = data_type, source = "pcawg")

  # Get sample info
  sample_info <- load_data("pcawg_info")

  # Prepare data
  plot_data <- data.frame(
    Sample = names(data),
    Expression = as.numeric(data),
    Cancer = sample_info$project[match(names(data), sample_info$Sample)],
    Type = sample_info$type[match(names(data), sample_info$Sample)],
    stringsAsFactors = FALSE
  )

  # Calculate median expression by cancer
  comp_results <- plot_data |>
    dplyr::group_by(.data$Cancer, .data$Type) |>
    dplyr::summarise(
      Median = stats::median(.data$Expression, na.rm = TRUE),
      Mean = mean(.data$Expression, na.rm = TRUE),
      N = dplyr::n(),
      .groups = "drop"
    )

  # Create plot
  plot <- ggplot2::ggplot(comp_results, ggplot2::aes(x = stats::reorder(.data$Cancer, .data$Median), y = .data$Median, fill = .data$Type)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::scale_fill_manual(values = c("tumor" = "#DF2020", "normal" = "#DDDF21")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(molecule, "Expression across PCAWG cancers"),
      x = "Cancer Type",
      y = "Median Expression"
    ) +
    theme_zinasuite()

  list(plot = plot, data = comp_results, stats = summary(comp_results$Median))
}

#' Run PCAWG Comparison Analysis - Many-to-One
#' @keywords internal
run_pcawg_comp_m2o <- function(molecules_text, cancer, data_type) {
  # Parse molecules
  molecules <- strsplit(molecules_text, "\n")[[1]]
  molecules <- trimws(molecules[molecules != ""])

  # Get sample info
  sample_info <- load_data("pcawg_info")
  cancer_samples <- sample_info$Sample[sample_info$project == cancer]

  # Query and aggregate data
  comp_results <- purrr::map_dfr(molecules, function(mol) {
    data <- query_molecule(mol, data_type = data_type, source = "pcawg")
    common_samples <- intersect(names(data), cancer_samples)

    data.frame(
      Molecule = mol,
      Expression = as.numeric(data[common_samples]),
      Type = sample_info$type[match(common_samples, sample_info$Sample)]
    )
  })

  # Create plot
  plot <- ggplot2::ggplot(comp_results, ggplot2::aes(x = stats::reorder(.data$Molecule, .data$Expression, stats::median), y = .data$Expression, fill = .data$Type)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::scale_fill_manual(values = c("tumor" = "#DF2020", "normal" = "#DDDF21")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Expression comparison in", cancer),
      x = "Molecule",
      y = "Expression"
    ) +
    theme_zinasuite()

  list(plot = plot, data = comp_results, stats = summary(comp_results$Expression))
}

# PCAWG Survival Analysis Functions -------------------------------------------

#' Run PCAWG Survival Analysis - One-to-One
#' @keywords internal
run_pcawg_sur_o2o <- function(molecule, cancer, data_type, measure, cutoff_mode) {
  # Query molecule data
  data <- query_molecule(molecule, data_type = data_type, source = "pcawg")

  # Get survival data (placeholder - needs actual PCAWG survival data)
  # For now, create synthetic data for demonstration
  plot_data <- data.frame(
    Sample = names(data),
    Expression = as.numeric(data),
    Time = stats::runif(length(data), 100, 1000),
    Status = sample(0:1, length(data), replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Create groups based on cutoff
  if (cutoff_mode == "median") {
    cutoff <- stats::median(plot_data$Expression, na.rm = TRUE)
  } else {
    cutoff <- stats::quantile(plot_data$Expression, 0.5, na.rm = TRUE)
  }

  plot_data$Group <- ifelse(plot_data$Expression >= cutoff, "High", "Low")

  # Create KM plot
  fit <- survival::survfit(survival::Surv(Time, Status) ~ Group, data = plot_data)
  plot <- survminer::ggsurvplot(
    fit,
    data = plot_data,
    pval = TRUE,
    risk.table = TRUE,
    title = paste(molecule, "Survival in", cancer),
    palette = c("High" = "red", "Low" = "blue")
  )

  list(plot = plot, data = plot_data, stats = summary(plot_data$Expression))
}

#' Run PCAWG Survival Analysis - One-to-Many
#' @keywords internal
run_pcawg_sur_o2m <- function(molecule, data_type, measure, cutoff_mode) {
  # Placeholder for multi-cancer survival analysis
  plot_data <- data.frame(
    Cancer = c("BRCA", "LUAD", "LUSC", "COAD", "READ"),
    HR = c(1.2, 0.8, 1.5, 1.1, 0.9),
    Lower = c(0.9, 0.6, 1.1, 0.8, 0.7),
    Upper = c(1.6, 1.1, 2.0, 1.5, 1.2),
    Pvalue = c(0.02, 0.15, 0.001, 0.5, 0.3)
  )

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stats::reorder(.data$Cancer, .data$HR), y = .data$HR)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Pvalue < 0.05), size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$Lower, ymax = .data$Upper), width = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(molecule, "Cox Regression across PCAWG cancers"),
      x = "Cancer Type",
      y = "Hazard Ratio (95% CI)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = plot_data, stats = NULL)
}

#' Run PCAWG Survival Analysis - Many-to-One
#' @keywords internal
run_pcawg_sur_m2o <- function(molecules_text, cancer, data_type, measure, cutoff_mode) {
  # Parse molecules
  molecules <- strsplit(molecules_text, "\n")[[1]]
  molecules <- trimws(molecules[molecules != ""])

  # Placeholder for multi-molecule survival analysis
  plot_data <- data.frame(
    Molecule = molecules,
    HR = stats::runif(length(molecules), 0.5, 2.0),
    Pvalue = stats::runif(length(molecules), 0.001, 0.5)
  )
  plot_data$Lower <- plot_data$HR * 0.7
  plot_data$Upper <- plot_data$HR * 1.3

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stats::reorder(.data$Molecule, .data$HR), y = .data$HR)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Pvalue < 0.05), size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$Lower, ymax = .data$Upper), width = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Cox Regression in", cancer),
      x = "Molecule",
      y = "Hazard Ratio (95% CI)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = plot_data, stats = NULL)
}
