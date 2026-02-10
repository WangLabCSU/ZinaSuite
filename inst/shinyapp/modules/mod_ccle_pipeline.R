#' CCLE Pipeline Analysis Module UI
#'
#' @description
#' Deep analysis pipeline for CCLE (Cancer Cell Line Encyclopedia) data,
#' supporting correlation and comparison analysis with multiple modes
#' (one-to-one, one-to-many, many-to-one).
#'
#' @param id Module ID
#' @return UI elements
#' @export
mod_ccle_pipeline_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::card(
      bslib::card_header(
        class = "bg-warning text-dark",
        shiny::icon("flask"), "CCLE Deep Analysis Pipeline"
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
                "Comparison" = "comp"
              ),
              selected = "cor"
            ),

            # Mode Selection
            shiny::selectInput(
              ns("mode"),
              "Analysis Mode:",
              choices = c(
                "One-to-One (Single Site)" = "o2o",
                "One-to-Many (Multiple Sites)" = "o2m",
                "Many-to-One (Multiple Molecules)" = "m2o"
              ),
              selected = "o2o"
            ),

            shiny::hr(),

            # Site Selection
            shiny::selectInput(
              ns("site"),
              "Primary Site:",
              choices = NULL,  # Will be populated server-side
              selected = NULL,
              multiple = TRUE
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
                "Protein Expression" = "protein",
                "CNV" = "cnv",
                "Mutation" = "mutation"
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

            shiny::hr(),

            # Action Button
            shiny::actionButton(
              ns("run_btn"),
              "Run Analysis",
              icon = shiny::icon("play"),
              class = "btn-warning w-100"
            ),

            shiny::hr(),

            # Download Options
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

#' CCLE Pipeline Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
#' @param async_compute Async compute engine
#' @export
mod_ccle_pipeline_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive values to store results
    results <- shiny::reactiveValues(
      plot = NULL,
      data = NULL,
      stats = NULL
    )

    # Load CCLE primary sites
    shiny::observe({
      ccle_sites <- load_data("ccle_sites")
      shiny::updateSelectInput(
        session, "site",
        choices = ccle_sites,
        selected = ccle_sites[1:3]
      )
    })

    # Run analysis when button is clicked
    shiny::observeEvent(input$run_btn, {
      analysis_type <- input$analysis_type
      mode <- input$mode
      site <- input$site
      molecule <- input$molecule
      data_type <- input$data_type

      # Validate input
      if (molecule == "") {
        shiny::showNotification("Please enter a molecule identifier", type = "error")
        return()
      }

      if (is.null(site) || length(site) == 0) {
        shiny::showNotification("Please select at least one primary site", type = "error")
        return()
      }

      # Show progress
      shiny::withProgress(message = "Running CCLE pipeline analysis...", value = 0, {
        tryCatch({
          result <- switch(analysis_type,
            "cor" = {
              shiny::incProgress(0.3, detail = "Running correlation analysis")
              switch(mode,
                "o2o" = run_ccle_cor_o2o(molecule, input$molecule2, site, data_type, input$cor_method),
                "o2m" = run_ccle_cor_o2m(molecule, input$molecule2, site, data_type, input$cor_method),
                "m2o" = run_ccle_cor_m2o(input$molecules, site, data_type, input$cor_method),
                stop("Unknown mode: ", mode)
              )
            },
            "comp" = {
              shiny::incProgress(0.3, detail = "Running comparison analysis")
              switch(mode,
                "o2o" = run_ccle_comp_o2o(molecule, site, data_type),
                "o2m" = run_ccle_comp_o2m(molecule, site, data_type),
                "m2o" = run_ccle_comp_m2o(input$molecules, site, data_type),
                stop("Unknown mode: ", mode)
              )
            },
            stop("Unknown analysis type: ", analysis_type)
          )

          results$plot <- result$plot
          results$data <- result$data
          results$stats <- result$stats

          shiny::incProgress(1, detail = "Complete")
          shiny::showNotification("CCLE pipeline analysis complete!", type = "message")

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
        paste0(input$molecule, "_ccle_", input$analysis_type, "_", input$mode, ".pdf")
      },
      content = function(file) {
        shiny::req(results$plot)
        ggplot2::ggsave(file, plot = results$plot, width = 12, height = 8, dpi = 300)
      }
    )

    output$download_data <- shiny::downloadHandler(
      filename = function() {
        paste0(input$molecule, "_ccle_", input$analysis_type, "_", input$mode, "_data.csv")
      },
      content = function(file) {
        shiny::req(results$data)
        utils::write.csv(results$data, file, row.names = FALSE)
      }
    )
  })
}

# CCLE Correlation Analysis Functions -----------------------------------------

#' Run CCLE Correlation Analysis - One-to-One
#' @keywords internal
run_ccle_cor_o2o <- function(molecule1, molecule2, sites, data_type, cor_method) {
  # Query molecule data
  data1 <- query_molecule(molecule1, data_type = data_type, source = "ccle")
  data2 <- query_molecule(molecule2, data_type = data_type, source = "ccle")

  # Get sample info for site filtering
  sample_info <- load_data("ccle_info")

  # Filter by primary sites
  site_samples <- sample_info$sample[sample_info$primary_site %in% sites]
  common_samples <- intersect(intersect(names(data1), names(data2)), site_samples)

  if (length(common_samples) < 10) {
    stop("Insufficient samples for correlation analysis (need at least 10)")
  }

  # Prepare data
  plot_data <- data.frame(
    Sample = common_samples,
    X = as.numeric(data1[common_samples]),
    Y = as.numeric(data2[common_samples]),
    Site = sample_info$primary_site[match(common_samples, sample_info$sample)],
    CellLine = sample_info$cell_line[match(common_samples, sample_info$sample)],
    stringsAsFactors = FALSE
  )

  # Calculate correlation
  cor_test <- stats::cor.test(plot_data$X, plot_data$Y, method = cor_method)

  # Create plot with color by site
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$X, y = .data$Y, color = .data$Site)) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::geom_smooth(method = "lm", color = "black", se = TRUE) +
    ggplot2::labs(
      title = paste(molecule1, "vs", molecule2),
      subtitle = sprintf("%s r = %.3f, p = %.2e", cor_method, cor_test$estimate, cor_test$p.value),
      x = paste(molecule1, "Expression"),
      y = paste(molecule2, "Expression")
    ) +
    theme_zinasuite()

  stats <- list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    method = cor_method,
    n_samples = nrow(plot_data),
    n_sites = length(unique(plot_data$Site))
  )

  list(plot = plot, data = plot_data, stats = stats)
}

#' Run CCLE Correlation Analysis - One-to-Many
#' @keywords internal
run_ccle_cor_o2m <- function(molecule1, molecule2, sites, data_type, cor_method) {
  # Query data
  data1 <- query_molecule(molecule1, data_type = data_type, source = "ccle")
  data2 <- query_molecule(molecule2, data_type = data_type, source = "ccle")

  sample_info <- load_data("ccle_info")

  # Filter by sites
  site_samples <- sample_info$sample[sample_info$primary_site %in% sites]

  # Find common samples
  common_samples <- intersect(intersect(names(data1), names(data2)), site_samples)

  # Prepare data with site info
  plot_data <- data.frame(
    Sample = common_samples,
    X = as.numeric(data1[common_samples]),
    Y = as.numeric(data2[common_samples]),
    Site = sample_info$primary_site[match(common_samples, sample_info$sample)],
    stringsAsFactors = FALSE
  )

  # Calculate correlation by site
  cor_results <- plot_data |>
    dplyr::group_by(.data$Site) |>
    dplyr::summarise(
      Correlation = stats::cor(.data$X, .data$Y, method = cor_method),
      P_value = tryCatch(stats::cor.test(.data$X, .data$Y, method = cor_method)$p.value,
                        error = function(e) NA),
      N = dplyr::n(),
      .groups = "drop"
    )

  # Create forest plot
  plot <- ggplot2::ggplot(cor_results, ggplot2::aes(x = stats::reorder(.data$Site, .data$Correlation), y = .data$Correlation)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$N, color = .data$P_value < 0.05)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(molecule1, "vs", molecule2, "across CCLE sites"),
      x = "Primary Site",
      y = "Correlation Coefficient",
      size = "Sample Size",
      color = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = cor_results, stats = summary(cor_results$Correlation))
}

#' Run CCLE Correlation Analysis - Many-to-One
#' @keywords internal
run_ccle_cor_m2o <- function(molecules_text, sites, data_type, cor_method) {
  # Parse molecules
  molecules <- strsplit(molecules_text, "\n")[[1]]
  molecules <- trimws(molecules[molecules != ""])

  if (length(molecules) < 2) {
    stop("At least 2 molecules required for many-to-one analysis")
  }

  # Query reference molecule (first one)
  ref_data <- query_molecule(molecules[1], data_type = data_type, source = "ccle")

  # Get sample info
  sample_info <- load_data("ccle_info")
  site_samples <- sample_info$sample[sample_info$primary_site %in% sites]

  # Calculate correlation for each molecule
  cor_results <- purrr::map_dfr(molecules, function(mol) {
    mol_data <- query_molecule(mol, data_type = data_type, source = "ccle")
    common_samples <- intersect(intersect(names(ref_data), names(mol_data)), site_samples)

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
      title = paste("Correlation with", molecules[1]),
      x = "Molecule",
      y = "Correlation Coefficient",
      fill = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = cor_results, stats = summary(cor_results$Correlation))
}

# CCLE Comparison Analysis Functions ------------------------------------------

#' Run CCLE Comparison Analysis - One-to-One
#' @keywords internal
run_ccle_comp_o2o <- function(molecule, sites, data_type) {
  # Query molecule data
  data <- query_molecule(molecule, data_type = data_type, source = "ccle")

  # Get sample info
  sample_info <- load_data("ccle_info")

  # Filter by sites
  site_samples <- sample_info$sample[sample_info$primary_site %in% sites]
  common_samples <- intersect(names(data), site_samples)

  # Prepare data
  plot_data <- data.frame(
    Sample = common_samples,
    Expression = as.numeric(data[common_samples]),
    Site = sample_info$primary_site[match(common_samples, sample_info$sample)],
    CellLine = sample_info$cell_line[match(common_samples, sample_info$sample)],
    stringsAsFactors = FALSE
  )

  if (nrow(plot_data) < 10) {
    stop("Insufficient samples for comparison")
  }

  # Create plot
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stats::reorder(.data$Site, .data$Expression, stats::median), y = .data$Expression)) +
    ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(molecule, "Expression across CCLE sites"),
      x = "Primary Site",
      y = "Expression"
    ) +
    theme_zinasuite()

  # Calculate summary statistics by site
  comp_results <- plot_data |>
    dplyr::group_by(.data$Site) |>
    dplyr::summarise(
      Median = stats::median(.data$Expression, na.rm = TRUE),
      Mean = mean(.data$Expression, na.rm = TRUE),
      SD = stats::sd(.data$Expression, na.rm = TRUE),
      N = dplyr::n(),
      .groups = "drop"
    )

  list(plot = plot, data = comp_results, stats = summary(plot_data$Expression))
}

#' Run CCLE Comparison Analysis - One-to-Many
#' @keywords internal
run_ccle_comp_o2m <- function(molecule, sites, data_type) {
  # Query molecule data
  data <- query_molecule(molecule, data_type = data_type, source = "ccle")

  # Get sample info
  sample_info <- load_data("ccle_info")

  # Filter by sites
  site_samples <- sample_info$sample[sample_info$primary_site %in% sites]
  common_samples <- intersect(names(data), site_samples)

  # Prepare data
  plot_data <- data.frame(
    Sample = common_samples,
    Expression = as.numeric(data[common_samples]),
    Site = sample_info$primary_site[match(common_samples, sample_info$sample)],
    stringsAsFactors = FALSE
  )

  # Calculate statistics by site
  comp_results <- plot_data |>
    dplyr::group_by(.data$Site) |>
    dplyr::summarise(
      Median = stats::median(.data$Expression, na.rm = TRUE),
      Mean = mean(.data$Expression, na.rm = TRUE),
      Min = min(.data$Expression, na.rm = TRUE),
      Max = max(.data$Expression, na.rm = TRUE),
      N = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$Median))

  # Create bar plot
  plot <- ggplot2::ggplot(comp_results, ggplot2::aes(x = stats::reorder(.data$Site, .data$Median), y = .data$Median)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$Min, ymax = .data$Max), width = 0.2) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(molecule, "Expression across CCLE sites"),
      x = "Primary Site",
      y = "Median Expression (range)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = comp_results, stats = summary(comp_results$Median))
}

#' Run CCLE Comparison Analysis - Many-to-One
#' @keywords internal
run_ccle_comp_m2o <- function(molecules_text, sites, data_type) {
  # Parse molecules
  molecules <- strsplit(molecules_text, "\n")[[1]]
  molecules <- trimws(molecules[molecules != ""])

  # Get sample info
  sample_info <- load_data("ccle_info")
  site_samples <- sample_info$sample[sample_info$primary_site %in% sites]

  # Query and aggregate data
  comp_results <- purrr::map_dfr(molecules, function(mol) {
    data <- query_molecule(mol, data_type = data_type, source = "ccle")
    common_samples <- intersect(names(data), site_samples)

    data.frame(
      Molecule = mol,
      Expression = as.numeric(data[common_samples]),
      Sample = common_samples
    )
  })

  # Create plot
  plot <- ggplot2::ggplot(comp_results, ggplot2::aes(x = stats::reorder(.data$Molecule, .data$Expression, stats::median), y = .data$Expression)) +
    ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Expression comparison across molecules",
      x = "Molecule",
      y = "Expression"
    ) +
    theme_zinasuite()

  # Calculate summary by molecule
  summary_results <- comp_results |>
    dplyr::group_by(.data$Molecule) |>
    dplyr::summarise(
      Median = stats::median(.data$Expression, na.rm = TRUE),
      Mean = mean(.data$Expression, na.rm = TRUE),
      N = dplyr::n(),
      .groups = "drop"
    )

  list(plot = plot, data = summary_results, stats = summary(comp_results$Expression))
}
