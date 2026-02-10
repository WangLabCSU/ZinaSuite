#' Other Page Modules UI
#'
#' @description
#' Additional page modules including daily gene, pan-cancer search,
#' file upload, and download functionality.
#'
#' @param id Module ID
#' @return UI elements
#' @export
mod_other_pages_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::card(
      bslib::card_header(
        class = "bg-secondary text-white",
        shiny::icon("th-large"), "Other Tools & Features"
      ),
      bslib::card_body(
        bslib::navset_card_tab(
          # Tab 1: Daily Gene
          bslib::nav_panel(
            title = shiny::icon("calendar-day"),
            "Daily Gene",
            shiny::fluidRow(
              shiny::column(
                width = 12,
                shiny::div(
                  class = "jumbotron",
                  style = "padding: 2rem; background-color: #f8f9fa; border-radius: 0.5rem;",
                  shiny::h1("Gene of the Day", class = "display-4"),
                  shiny::p(class = "lead", "Explore a featured gene with comprehensive pan-cancer analysis"),
                  shiny::hr(class = "my-4"),
                  shiny::h2(shiny::textOutput(ns("daily_gene_name"))),
                  shiny::p(shiny::textOutput(ns("daily_gene_desc"))),
                  shiny::actionButton(
                    ns("analyze_daily_gene"),
                    "Analyze This Gene",
                    icon = shiny::icon("chart-line"),
                    class = "btn-primary btn-lg"
                  )
                )
              )
            ),
            shiny::hr(),
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::h4("Expression Overview"),
                shiny::plotOutput(ns("daily_gene_expr"), height = "400px")
              ),
              shiny::column(
                width = 6,
                shiny::h4("Survival Analysis"),
                shiny::plotOutput(ns("daily_gene_survival"), height = "400px")
              )
            )
          ),

          # Tab 2: Pan-Cancer Search
          bslib::nav_panel(
            title = shiny::icon("search"),
            "Pan-Cancer Search",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Search Parameters"),
                  shiny::textInput(
                    ns("search_query"),
                    "Gene or Pathway:",
                    value = "TP53",
                    placeholder = "Enter gene symbol or pathway name"
                  ),
                  shiny::selectInput(
                    ns("search_type"),
                    "Search Type:",
                    choices = c(
                      "Gene Symbol" = "gene",
                      "Pathway" = "pathway",
                      "Drug" = "drug"
                    )
                  ),
                  shiny::selectInput(
                    ns("search_cancers"),
                    "Cancer Types:",
                    choices = NULL,  # Populated server-side
                    multiple = TRUE,
                    selected = NULL
                  ),
                  shiny::hr(),
                  shiny::actionButton(
                    ns("run_search"),
                    "Search",
                    icon = shiny::icon("search"),
                    class = "btn-secondary w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::h4("Search Results"),
                DT::dataTableOutput(ns("search_results_table")),
                shiny::hr(),
                shiny::h4("Visual Summary"),
                shiny::plotOutput(ns("search_summary_plot"), height = "300px")
              )
            )
          ),

          # Tab 3: File Upload
          bslib::nav_panel(
            title = shiny::icon("file-upload"),
            "File Upload",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Upload Data"),
                  shiny::fileInput(
                    ns("upload_file"),
                    "Choose File:",
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv",
                      ".txt",
                      ".tsv"
                    )
                  ),
                  shiny::helpText("Supported formats: CSV, TSV, TXT"),
                  shiny::radioButtons(
                    ns("upload_type"),
                    "Data Type:",
                    choices = c(
                      "Gene List" = "genes",
                      "Sample IDs" = "samples",
                      "Custom Data" = "custom"
                    )
                  ),
                  shiny::hr(),
                  shiny::actionButton(
                    ns("process_upload"),
                    "Process File",
                    icon = shiny::icon("cog"),
                    class = "btn-secondary w-100"
                  ),
                  shiny::hr(),
                  shiny::downloadButton(
                    ns("download_template"),
                    "Download Template",
                    class = "btn-outline-secondary w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::h4("Uploaded Data Preview"),
                DT::dataTableOutput(ns("upload_preview_table")),
                shiny::hr(),
                shiny::verbatimTextOutput(ns("upload_stats"))
              )
            )
          ),

          # Tab 4: Data Download
          bslib::nav_panel(
            title = shiny::icon("download"),
            "Data Download",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Download Options"),
                  shiny::selectInput(
                    ns("download_dataset"),
                    "Dataset:",
                    choices = c(
                      "TCGA Pan-Cancer" = "tcga",
                      "PCAWG" = "pcawg",
                      "CCLE" = "ccle"
                    )
                  ),
                  shiny::selectInput(
                    ns("download_data_type"),
                    "Data Type:",
                    choices = c(
                      "Expression Matrix" = "expression",
                      "Clinical Data" = "clinical",
                      "Mutation Data" = "mutation",
                      "Copy Number" = "cnv",
                      "Methylation" = "methylation"
                    )
                  ),
                  shiny::selectInput(
                    ns("download_format"),
                    "Format:",
                    choices = c(
                      "CSV" = "csv",
                      "TSV" = "tsv",
                      "RData" = "rdata"
                    )
                  ),
                  shiny::hr(),
                  shiny::downloadButton(
                    ns("download_data_btn"),
                    "Download",
                    class = "btn-secondary w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::h4("Dataset Information"),
                shiny::tableOutput(ns("dataset_info_table")),
                shiny::hr(),
                shiny::h4("Preview"),
                shiny::verbatimTextOutput(ns("download_preview"))
              )
            )
          ),

          # Tab 5: ID Query Help
          bslib::nav_panel(
            title = shiny::icon("question-circle"),
            "ID Query Help",
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::h4("Supported ID Types"),
                shiny::tags$ul(
                  shiny::tags$li(shiny::strong("Gene Symbols:"), " TP53, BRCA1, EGFR, etc."),
                  shiny::tags$li(shiny::strong("Ensembl IDs:"), " ENSG00000141510, etc."),
                  shiny::tags$li(shiny::strong("Entrez IDs:"), " 7157, 672, etc."),
                  shiny::tags$li(shiny::strong("Sample IDs:"), " TCGA-XX-XXXX-01A, etc."),
                  shiny::tags$li(shiny::strong("Pathway IDs:"), " HALLMARK_APOPTOSIS, etc.")
                ),
                shiny::hr(),
                shiny::h4("ID Conversion"),
                shiny::textInput(
                  ns("convert_id_input"),
                  "Enter ID:",
                  value = "TP53"
                ),
                shiny::selectInput(
                  ns("convert_id_from"),
                  "From:",
                  choices = c("Gene Symbol", "Ensembl ID", "Entrez ID")
                ),
                shiny::selectInput(
                  ns("convert_id_to"),
                  "To:",
                  choices = c("Ensembl ID", "Entrez ID", "Gene Symbol")
                ),
                shiny::actionButton(
                  ns("convert_id_btn"),
                  "Convert",
                  icon = shiny::icon("exchange-alt"),
                  class = "btn-secondary"
                ),
                shiny::hr(),
                shiny::verbatimTextOutput(ns("convert_result"))
              ),
              shiny::column(
                width = 6,
                shiny::h4("Quick Reference"),
                shiny::tags$table(
                  class = "table table-striped",
                  shiny::tags$thead(
                    shiny::tags$tr(
                      shiny::tags$th("ID Type"),
                      shiny::tags$th("Example"),
                      shiny::tags$th("Source")
                    )
                  ),
                  shiny::tags$tbody(
                    shiny::tags$tr(
                      shiny::tags$td("Gene Symbol"),
                      shiny::tags$td("TP53"),
                      shiny::tags$td("HGNC")
                    ),
                    shiny::tags$tr(
                      shiny::tags$td("Ensembl"),
                      shiny::tags$td("ENSG00000141510"),
                      shiny::tags$td("Ensembl")
                    ),
                    shiny::tags$tr(
                      shiny::tags$td("Entrez"),
                      shiny::tags$td("7157"),
                      shiny::tags$td("NCBI")
                    ),
                    shiny::tags$tr(
                      shiny::tags$td("TCGA Sample"),
                      shiny::tags$td("TCGA-BRCA-01"),
                      shiny::tags$td("TCGA")
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

#' Other Pages Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
#' @param async_compute Async compute engine
#' @export
mod_other_pages_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive values
    rv <- shiny::reactiveValues(
      daily_gene = NULL,
      search_results = NULL,
      uploaded_data = NULL
    )

    # Load cancer types
    shiny::observe({
      tcga_cancers <- load_data("tcga_cancers")
      shiny::updateSelectInput(session, "search_cancers",
                               choices = tcga_cancers,
                               selected = tcga_cancers[1:5])
    })

    # ===== Daily Gene Tab =====

    # Set daily gene based on date
    shiny::observe({
      # Use date to deterministically select a gene
      genes <- c("TP53", "BRCA1", "EGFR", "PTEN", "MYC", "KRAS", "BRCA2", "ATM")
      day_of_year <- as.numeric(format(Sys.Date(), "%j"))
      selected_gene <- genes[(day_of_year %% length(genes)) + 1]

      rv$daily_gene <- list(
        symbol = selected_gene,
        name = paste(selected_gene, "Gene"),
        description = paste("Tumor protein", selected_gene, "is a crucial gene involved in cancer development.")
      )
    })

    output$daily_gene_name <- shiny::renderText({
      req(rv$daily_gene)
      rv$daily_gene$symbol
    })

    output$daily_gene_desc <- shiny::renderText({
      req(rv$daily_gene)
      rv$daily_gene$description
    })

    output$daily_gene_expr <- shiny::renderPlot({
      req(rv$daily_gene)

      # Create mock expression data
      cancers <- load_data("tcga_cancers")[1:10]
      expr_data <- data.frame(
        Cancer = cancers,
        Expression = stats::rnorm(10, 5, 2),
        stringsAsFactors = FALSE
      )

      ggplot2::ggplot(expr_data, ggplot2::aes(x = stats::reorder(Cancer, Expression), y = Expression)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = paste(rv$daily_gene$symbol, "Expression"),
          x = "Cancer Type",
          y = "Expression Level"
        ) +
        theme_zinasuite()
    })

    output$daily_gene_survival <- shiny::renderPlot({
      req(rv$daily_gene)

      # Create mock survival data
      time <- seq(0, 120, by = 10)
      survival <- exp(-time / 100)

      surv_data <- data.frame(
        Time = time,
        Survival = survival,
        stringsAsFactors = FALSE
      )

      ggplot2::ggplot(surv_data, ggplot2::aes(x = Time, y = Survival)) +
        ggplot2::geom_line(color = "red", size = 1) +
        ggplot2::labs(
          title = paste(rv$daily_gene$symbol, "Survival Curve"),
          x = "Time (months)",
          y = "Survival Probability"
        ) +
        theme_zinasuite()
    })

    shiny::observeEvent(input$analyze_daily_gene, {
      shiny::showNotification(
        paste("Analyzing", rv$daily_gene$symbol, "- redirecting to analysis module..."),
        type = "message"
      )
    })

    # ===== Pan-Cancer Search Tab =====

    shiny::observeEvent(input$run_search, {
      shiny::withProgress(message = "Searching...", value = 0, {
        query <- input$search_query
        type <- input$search_type
        cancers <- input$search_cancers

        shiny::incProgress(0.5, detail = "Querying database")

        # Create mock search results
        rv$search_results <- data.frame(
          Cancer = cancers,
          Gene = query,
          Expression = stats::runif(length(cancers), 0, 10),
          P_Value = stats::runif(length(cancers), 0, 0.1),
          Significant = stats::runif(length(cancers), 0, 0.1) < 0.05,
          stringsAsFactors = FALSE
        )

        shiny::incProgress(1, detail = "Complete")
        shiny::showNotification(paste("Found results for", query), type = "message")
      })
    })

    output$search_results_table <- DT::renderDataTable({
      req(rv$search_results)
      DT::datatable(rv$search_results, options = list(pageLength = 10))
    })

    output$search_summary_plot <- shiny::renderPlot({
      req(rv$search_results)

      ggplot2::ggplot(rv$search_results, ggplot2::aes(x = Cancer, y = Expression, fill = Significant)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Search Results Summary",
          x = "Cancer Type",
          y = "Expression Level"
        ) +
        theme_zinasuite()
    })

    # ===== File Upload Tab =====

    shiny::observeEvent(input$process_upload, {
      req(input$upload_file)

      tryCatch({
        # Read uploaded file
        ext <- tools::file_ext(input$upload_file$name)
        data <- switch(ext,
                       "csv" = utils::read.csv(input$upload_file$datapath, stringsAsFactors = FALSE),
                       "tsv" = utils::read.delim(input$upload_file$datapath, stringsAsFactors = FALSE),
                       "txt" = utils::read.delim(input$upload_file$datapath, stringsAsFactors = FALSE),
                       utils::read.csv(input$upload_file$datapath, stringsAsFactors = FALSE)
        )

        rv$uploaded_data <- data
        shiny::showNotification(paste("Loaded", nrow(data), "rows"), type = "message")
      }, error = function(e) {
        shiny::showNotification(paste("Error:", conditionMessage(e)), type = "error")
      })
    })

    output$upload_preview_table <- DT::renderDataTable({
      if (is.null(rv$uploaded_data)) {
        return(data.frame(Message = "No data uploaded"))
      }
      DT::datatable(rv$uploaded_data, options = list(pageLength = 10, scrollX = TRUE))
    })

    output$upload_stats <- shiny::renderPrint({
      if (is.null(rv$uploaded_data)) {
        cat("No data uploaded yet.")
      } else {
        cat("Rows:", nrow(rv$uploaded_data), "\n")
        cat("Columns:", ncol(rv$uploaded_data), "\n")
        cat("Column names:", paste(names(rv$uploaded_data), collapse = ", "))
      }
    })

    output$download_template <- shiny::downloadHandler(
      filename = function() {
        paste0(input$upload_type, "_template.csv")
      },
      content = function(file) {
        template <- switch(input$upload_type,
                           "genes" = data.frame(gene = c("TP53", "BRCA1", "EGFR"), stringsAsFactors = FALSE),
                           "samples" = data.frame(sample = c("TCGA-XX-XXXX-01", "TCGA-XX-XXXX-02"), stringsAsFactors = FALSE),
                           "custom" = data.frame(id = 1:3, value = c("A", "B", "C"), stringsAsFactors = FALSE)
        )
        utils::write.csv(template, file, row.names = FALSE)
      }
    )

    # ===== Data Download Tab =====

    output$dataset_info_table <- shiny::renderTable({
      data.frame(
        Dataset = c("TCGA", "PCAWG", "CCLE"),
        Samples = c("~11,000", "~2,800", "~1,000"),
        DataTypes = c("Expression, Clinical, Mutation, CNV, Methylation",
                      "Expression, Fusion, Promoter, miRNA",
                      "Expression, Protein, Drug Response"),
        LastUpdated = c("2024-01", "2024-01", "2024-01"),
        stringsAsFactors = FALSE
      )
    })

    output$download_preview <- shiny::renderPrint({
      cat("Download Settings:\n")
      cat("Dataset:", input$download_dataset, "\n")
      cat("Data Type:", input$download_data_type, "\n")
      cat("Format:", input$download_format, "\n")
      cat("\nNote: Actual data download would be implemented here.")
    })

    output$download_data_btn <- shiny::downloadHandler(
      filename = function() {
        ext <- switch(input$download_format,
                      "csv" = ".csv",
                      "tsv" = ".tsv",
                      "rdata" = ".RData")
        paste0(input$download_dataset, "_", input$download_data_type, ext)
      },
      content = function(file) {
        # Create sample data
        data <- data.frame(
          sample = paste0("sample_", 1:100),
          value = stats::rnorm(100),
          stringsAsFactors = FALSE
        )

        switch(input$download_format,
               "csv" = utils::write.csv(data, file, row.names = FALSE),
               "tsv" = utils::write.table(data, file, sep = "\t", row.names = FALSE),
               "rdata" = save(data, file = file)
        )
      }
    )

    # ===== ID Query Help Tab =====

    shiny::observeEvent(input$convert_id_btn, {
      input_id <- input$convert_id_input
      from <- input$convert_id_from
      to <- input$convert_id_to

      # Mock conversion
      result <- switch(from,
                       "Gene Symbol" = paste0("Converted ", input_id, " to ", to, ": ENSG00000141510"),
                       "Ensembl ID" = paste0("Converted ", input_id, " to ", to, ": TP53"),
                       "Entrez ID" = paste0("Converted ", input_id, " to ", to, ": TP53")
      )

      output$convert_result <- shiny::renderPrint({
        cat(result)
      })
    })

    # Return reactive values
    return(rv)
  })
}
