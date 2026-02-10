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

  shinydashboard::box(
    title = "PharmacoGenomics Analysis",
    status = "primary",
    solidHeader = TRUE,
    width = 12,

    shiny::tabsetPanel(
      id = ns("pharma_tabs"),

      # Drug-Omics Correlation Tab
      shiny::tabPanel(
        title = "Drug-Omics Correlation",
        icon = shiny::icon("pills"),

        shiny::fluidRow(
          shiny::column(
            width = 3,
            shiny::wellPanel(
              shiny::h4("Analysis Parameters"),

              shiny::selectInput(
                ns("data_source"),
                "Data Source:",
                choices = c(
                  "CCLE" = "ccle",
                  "GDSC" = "gdsc",
                  "CTRP" = "ctrp"
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
                choices = c(
                  "Select drug..." = "",
                  "Paclitaxel" = "paclitaxel",
                  "Doxorubicin" = "doxorubicin",
                  "Cisplatin" = "cisplatin",
                  "Gemcitabine" = "gemcitabine",
                  "5-Fluorouracil" = "5-fu",
                  "Temozolomide" = "temozolomide"
                ),
                selected = ""
              ),

              shiny::selectInput(
                ns("omic_type"),
                "Omics Type:",
                choices = c(
                  "Gene Expression" = "expression",
                  "Copy Number" = "cnv",
                  "Mutation" = "mutation",
                  "Protein" = "protein"
                ),
                selected = "expression"
              ),

              shiny::hr(),

              shiny::actionButton(
                ns("run_correlation"),
                "Run Analysis",
                icon = shiny::icon("play"),
                class = "btn-primary"
              )
            )
          ),

          shiny::column(
            width = 9,
            shiny::wellPanel(
              shiny::h4("Drug-Omics Correlation Results"),
              shiny::plotOutput(ns("correlation_plot"), height = "500px"),
              shiny::hr(),
              shiny::verbatimTextOutput(ns("correlation_stats"))
            )
          )
        )
      ),

      # Drug Sensitivity Profile Tab
      shiny::tabPanel(
        title = "Drug Sensitivity Profile",
        icon = shiny::icon("chart-bar"),

        shiny::fluidRow(
          shiny::column(
            width = 3,
            shiny::wellPanel(
              shiny::h4("Profile Parameters"),

              shiny::selectInput(
                ns("profile_source"),
                "Database:",
                choices = c(
                  "CCLE" = "ccle",
                  "GDSC1" = "gdsc1",
                  "GDSC2" = "gdsc2"
                ),
                selected = "ccle"
              ),

              shiny::selectInput(
                ns("profile_metric"),
                "Sensitivity Metric:",
                choices = c(
                  "IC50" = "ic50",
                  "AUC" = "auc",
                  "EC50" = "ec50"
                ),
                selected = "ic50"
              ),

              shiny::checkboxGroupInput(
                ns("cancer_types"),
                "Cancer Types:",
                choices = c(
                  "BRCA" = "BRCA",
                  "LUAD" = "LUAD",
                  "LUSC" = "LUSC",
                  "COAD" = "COAD",
                  "READ" = "READ",
                  "HNSC" = "HNSC"
                ),
                selected = c("BRCA", "LUAD")
              ),

              shiny::hr(),

              shiny::actionButton(
                ns("show_profile"),
                "Show Profile",
                icon = shiny::icon("chart-bar"),
                class = "btn-primary"
              )
            )
          ),

          shiny::column(
            width = 9,
            shiny::wellPanel(
              shiny::h4("Drug Sensitivity Profile"),
              shiny::plotOutput(ns("sensitivity_plot"), height = "500px"),
              shiny::hr(),
              DT::dataTableOutput(ns("sensitivity_table"))
            )
          )
        )
      ),

      # Feature Association Tab
      shiny::tabPanel(
        title = "Feature Association",
        icon = shiny::icon("project-diagram"),

        shiny::fluidRow(
          shiny::column(
            width = 3,
            shiny::wellPanel(
              shiny::h4("Association Parameters"),

              shiny::textInput(
                ns("feature_gene"),
                "Feature Gene:",
                value = "TP53"
              ),

              shiny::selectInput(
                ns("feature_type"),
                "Feature Type:",
                choices = c(
                  "Mutation" = "mutation",
                  "CNV" = "cnv",
                  "Expression" = "expression"
                ),
                selected = "mutation"
              ),

              shiny::selectInput(
                ns("association_drug"),
                "Drug:",
                choices = c(
                  "All Drugs" = "all",
                  "Paclitaxel" = "paclitaxel",
                  "Doxorubicin" = "doxorubicin",
                  "Cisplatin" = "cisplatin"
                ),
                selected = "all"
              ),

              shiny::hr(),

              shiny::actionButton(
                ns("run_association"),
                "Run Association",
                icon = shiny::icon("calculator"),
                class = "btn-primary"
              )
            )
          ),

          shiny::column(
            width = 9,
            shiny::wellPanel(
              shiny::h4("Feature-Drug Association Results"),
              shiny::plotOutput(ns("association_plot"), height = "400px"),
              shiny::hr(),
              DT::dataTableOutput(ns("association_table"))
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
#' @export
mod_pharmacogenomics_server <- function(id, app_state) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Drug-Omics Correlation Analysis
    correlation_results <- shiny::eventReactive(input$run_correlation, {
      shiny::req(input$gene_id, input$drug_id)

      shiny::withProgress(message = 'Running drug-omics correlation...', {
        # Simulate analysis (in real implementation, query actual data)
        set.seed(42)
        n_samples <- 200

        omic_values <- rnorm(n_samples, mean = 5, sd = 2)
        drug_sensitivity <- -0.5 * omic_values + rnorm(n_samples, mean = 0, sd = 1)

        data.frame(
          omic_value = omic_values,
          drug_sensitivity = drug_sensitivity,
          cell_line = paste0("Cell_", 1:n_samples)
        )
      })
    })

    output$correlation_plot <- shiny::renderPlot({
      shiny::req(correlation_results())

      data <- correlation_results()

      ggplot2::ggplot(data, ggplot2::aes(x = omic_value, y = drug_sensitivity)) +
        ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
        ggplot2::geom_smooth(method = "lm", color = "red", se = TRUE) +
        ggplot2::labs(
          title = paste("Drug-Omics Correlation:", input$gene_id, "vs", input$drug_id),
          x = paste(input$omic_type, "level"),
          y = "Drug Sensitivity"
        ) +
        ggplot2::theme_minimal()
    })

    output$correlation_stats <- shiny::renderPrint({
      shiny::req(correlation_results())

      data <- correlation_results()
      cor_test <- cor.test(data$omic_value, data$drug_sensitivity)

      cat("Correlation Analysis Results:\n")
      cat("============================\n")
      cat("Pearson correlation:", round(cor_test$estimate, 4), "\n")
      cat("P-value:", format(cor_test$p.value, digits = 4), "\n")
      cat("95% CI:", paste(round(cor_test$conf.int, 4), collapse = " to "), "\n")
      cat("Sample size:", nrow(data), "\n")
    })

    # Drug Sensitivity Profile
    sensitivity_data <- shiny::eventReactive(input$show_profile, {
      shiny::req(input$cancer_types)

      shiny::withProgress(message = 'Loading drug sensitivity data...', {
        # Simulate data
        set.seed(123)
        data_list <- lapply(input$cancer_types, function(cancer) {
          data.frame(
            cancer_type = cancer,
            drug_name = c("Paclitaxel", "Doxorubicin", "Cisplatin", "Gemcitabine"),
            ic50 = runif(4, min = 0.1, max = 10),
            auc = runif(4, min = 0.2, max = 1.0),
            n_cell_lines = sample(20:100, 4)
          )
        })
        do.call(rbind, data_list)
      })
    })

    output$sensitivity_plot <- shiny::renderPlot({
      shiny::req(sensitivity_data())

      data <- sensitivity_data()

      ggplot2::ggplot(data, ggplot2::aes(x = cancer_type, y = .data[[input$profile_metric]], fill = drug_name)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::labs(
          title = paste("Drug Sensitivity Profile (", toupper(input$profile_metric), ")"),
          x = "Cancer Type",
          y = toupper(input$profile_metric)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    output$sensitivity_table <- DT::renderDataTable({
      shiny::req(sensitivity_data())
      DT::datatable(
        sensitivity_data(),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })

    # Feature Association
    association_results <- shiny::eventReactive(input$run_association, {
      shiny::req(input$feature_gene)

      shiny::withProgress(message = 'Running feature association...', {
        # Simulate association results
        set.seed(456)
        drugs <- c("Paclitaxel", "Doxorubicin", "Cisplatin", "Gemcitabine",
                   "5-Fluorouracil", "Temozolomide", "Erlotinib", "Lapatinib")

        data.frame(
          drug = drugs,
          estimate = runif(8, -1, 1),
          p_value = runif(8, 0.001, 0.1),
          n_samples = sample(50:200, 8)
        )
      })
    })

    output$association_plot <- shiny::renderPlot({
      shiny::req(association_results())

      data <- association_results()
      data$significance <- ifelse(data$p_value < 0.05, "Significant", "Not Significant")

      ggplot2::ggplot(data, ggplot2::aes(x = reorder(drug, estimate), y = estimate, fill = significance)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = paste("Feature Association:", input$feature_gene),
          x = "Drug",
          y = "Effect Size"
        ) +
        ggplot2::scale_fill_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
        ggplot2::theme_minimal()
    })

    output$association_table <- DT::renderDataTable({
      shiny::req(association_results())

      data <- association_results()
      data$p_value <- format(data$p_value, digits = 4)
      data$estimate <- round(data$estimate, 4)

      DT::datatable(
        data,
        options = list(pageLength = 10),
        rownames = FALSE
      )
    })

  })
}
