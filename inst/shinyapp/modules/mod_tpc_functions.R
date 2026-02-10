#' TPC (TCGA/PCAWG/CCLE) Functional Modules UI
#'
#' @description
#' Shared functional modules for TPC analysis including custom signatures,
#' metadata upload, sample filtering, and data download functionality.
#'
#' @param id Module ID
#' @return UI elements
#' @export
mod_tpc_functions_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::card(
      bslib::card_header(
        class = "bg-info text-white",
        shiny::icon("tools"), "TPC Functional Modules"
      ),
      bslib::card_body(
        bslib::navset_card_tab(
          # Tab 1: Custom Signature
          bslib::nav_panel(
            title = shiny::icon("signature"),
            "Custom Signature",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Create Custom Signature"),
                  shiny::textInput(
                    ns("sig_name"),
                    "Signature Name:",
                    value = "My_Signature"
                  ),
                  shiny::selectInput(
                    ns("sig_data_type"),
                    "Data Type:",
                    choices = c(
                      "mRNA Expression" = "mRNA",
                      "Protein Expression" = "protein",
                      "Methylation" = "methylation"
                    ),
                    selected = "mRNA"
                  ),
                  shiny::textAreaInput(
                    ns("sig_genes"),
                    "Genes (one per line):",
                    value = "TP53\nKRAS\nEGFR",
                    rows = 5
                  ),
                  shiny::actionButton(
                    ns("create_sig_btn"),
                    "Create Signature",
                    icon = shiny::icon("plus"),
                    class = "btn-info w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::h4("Signature Preview"),
                shiny::verbatimTextOutput(ns("sig_preview")),
                shiny::hr(),
                shiny::h4("Available Signatures"),
                DT::dataTableOutput(ns("sig_table"))
              )
            )
          ),

          # Tab 2: Custom Metadata
          bslib::nav_panel(
            title = shiny::icon("file-upload"),
            "Custom Metadata",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Upload Custom Metadata"),
                  shiny::fileInput(
                    ns("meta_file"),
                    "Choose CSV File:",
                    accept = c("text/csv", ".csv")
                  ),
                  shiny::helpText("File should have 'sample' column and metadata columns"),
                  shiny::hr(),
                  shiny::h4("Manual Entry"),
                  shiny::textInput(
                    ns("meta_sample"),
                    "Sample ID:",
                    placeholder = "e.g., TCGA-XX-XXXX"
                  ),
                  shiny::textInput(
                    ns("meta_group"),
                    "Group:",
                    placeholder = "e.g., High, Low"
                  ),
                  shiny::actionButton(
                    ns("add_meta_btn"),
                    "Add Entry",
                    icon = shiny::icon("plus"),
                    class = "btn-info w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::h4("Metadata Preview"),
                DT::dataTableOutput(ns("meta_table")),
                shiny::hr(),
                shiny::downloadButton(
                  ns("download_meta_template"),
                  "Download Template",
                  class = "btn-outline-info"
                )
              )
            )
          ),

          # Tab 3: Sample Filtering
          bslib::nav_panel(
            title = shiny::icon("filter"),
            "Sample Filtering",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Filter Samples"),
                  shiny::selectInput(
                    ns("filter_cancer"),
                    "Cancer Type:",
                    choices = NULL,  # Populated server-side
                    multiple = TRUE
                  ),
                  shiny::selectInput(
                    ns("filter_sample_type"),
                    "Sample Type:",
                    choices = c("Tumor" = "tumor", "Normal" = "normal"),
                    selected = c("tumor", "normal"),
                    multiple = TRUE
                  ),
                  shiny::sliderInput(
                    ns("filter_expression"),
                    "Expression Range (percentile):",
                    min = 0,
                    max = 100,
                    value = c(0, 100)
                  ),
                  shiny::actionButton(
                    ns("apply_filter_btn"),
                    "Apply Filters",
                    icon = shiny::icon("filter"),
                    class = "btn-info w-100"
                  ),
                  shiny::hr(),
                  shiny::actionButton(
                    ns("reset_filter_btn"),
                    "Reset Filters",
                    icon = shiny::icon("undo"),
                    class = "btn-outline-secondary w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::h4("Filter Results"),
                shiny::verbatimTextOutput(ns("filter_stats")),
                shiny::hr(),
                shiny::h4("Filtered Samples"),
                DT::dataTableOutput(ns("filtered_samples_table"))
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
                  shiny::h4("Download Analysis Data"),
                  shiny::selectInput(
                    ns("download_data_type"),
                    "Data Type:",
                    choices = c(
                      "Expression Data" = "expression",
                      "Clinical Data" = "clinical",
                      "Mutation Data" = "mutation",
                      "Analysis Results" = "results"
                    )
                  ),
                  shiny::selectInput(
                    ns("download_format"),
                    "Format:",
                    choices = c("CSV" = "csv", "TSV" = "tsv", "RData" = "rdata"),
                    selected = "csv"
                  ),
                  shiny::checkboxInput(
                    ns("download_filtered_only"),
                    "Download filtered samples only",
                    value = TRUE
                  ),
                  shiny::hr(),
                  shiny::downloadButton(
                    ns("download_data_btn"),
                    "Download Data",
                    class = "btn-info w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::h4("Download Preview"),
                shiny::verbatimTextOutput(ns("download_preview")),
                shiny::hr(),
                shiny::h4("Data Summary"),
                shiny::tableOutput(ns("data_summary"))
              )
            )
          ),

          # Tab 5: Molecular Origin
          bslib::nav_panel(
            title = shiny::icon("database"),
            "Molecular Origin",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Set Data Source"),
                  shiny::selectInput(
                    ns("data_source"),
                    "Database:",
                    choices = c(
                      "TCGA (TOIL)" = "toil",
                      "PCAWG" = "pcawg",
                      "CCLE" = "ccle"
                    ),
                    selected = "toil"
                  ),
                  shiny::conditionalPanel(
                    condition = sprintf("input['%s'] == 'toil'", ns("data_source")),
                    shiny::selectInput(
                      ns("toil_data_type"),
                      "Data Type:",
                      choices = c(
                        "mRNA Expression" = "mRNA",
                        "Transcript" = "transcript",
                        "Protein" = "protein",
                        "Methylation" = "methylation",
                        "miRNA" = "miRNA",
                        "Mutation" = "mutation",
                        "CNV" = "cnv"
                      )
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = sprintf("input['%s'] == 'pcawg'", ns("data_source")),
                    shiny::selectInput(
                      ns("pcawg_data_type"),
                      "Data Type:",
                      choices = c(
                        "mRNA Expression" = "mRNA",
                        "Gene Fusion" = "fusion",
                        "Promoter Activity" = "promoter",
                        "miRNA" = "miRNA",
                        "APOBEC Mutagenesis" = "apobec"
                      )
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = sprintf("input['%s'] == 'ccle'", ns("data_source")),
                    shiny::selectInput(
                      ns("ccle_data_type"),
                      "Data Type:",
                      choices = c(
                        "mRNA Expression" = "mRNA",
                        "Protein" = "protein",
                        "CNV" = "cnv",
                        "Mutation" = "mutation"
                      )
                    )
                  ),
                  shiny::actionButton(
                    ns("set_origin_btn"),
                    "Set Origin",
                    icon = shiny::icon("check"),
                    class = "btn-info w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::h4("Current Settings"),
                shiny::verbatimTextOutput(ns("origin_settings")),
                shiny::hr(),
                shiny::h4("Available Data"),
                DT::dataTableOutput(ns("available_data_table"))
              )
            )
          ),

          # Tab 6: Sample Grouping
          bslib::nav_panel(
            title = shiny::icon("object-group"),
            "Sample Grouping",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::wellPanel(
                  shiny::h4("Define Groups"),
                  shiny::textInput(
                    ns("group_name"),
                    "Group Name:",
                    value = "Group_A"
                  ),
                  shiny::selectInput(
                    ns("group_by"),
                    "Group By:",
                    choices = c(
                      "Expression Level" = "expression",
                      "Clinical Feature" = "clinical",
                      "Custom" = "custom"
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = sprintf("input['%s'] == 'expression'", ns("group_by")),
                    shiny::textInput(
                      ns("group_gene"),
                      "Gene:",
                      value = "TP53"
                    ),
                    shiny::sliderInput(
                      ns("group_cutoff"),
                      "Cutoff (percentile):",
                      min = 0,
                      max = 100,
                      value = 50
                    )
                  ),
                  shiny::actionButton(
                    ns("create_group_btn"),
                    "Create Group",
                    icon = shiny::icon("plus"),
                    class = "btn-info w-100"
                  )
                )
              ),
              shiny::column(
                width = 8,
                shiny::h4("Group Summary"),
                shiny::verbatimTextOutput(ns("group_summary")),
                shiny::hr(),
                shiny::h4("Group Membership"),
                DT::dataTableOutput(ns("group_table"))
              )
            )
          )
        )
      )
    )
  )
}

#' TPC Functions Module Server
#'
#' @param id Module ID
#' @param app_state Shared reactive state
#' @param async_compute Async compute engine
#' @export
mod_tpc_functions_server <- function(id, app_state, async_compute) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive values for storing data
    rv <- shiny::reactiveValues(
      signatures = list(),
      metadata = data.frame(),
      filtered_samples = NULL,
      groups = list(),
      data_origin = list(database = "toil", data_type = "mRNA")
    )

    # Load cancer types for filtering
    shiny::observe({
      tcga_cancers <- load_data("tcga_cancers")
      shiny::updateSelectInput(session, "filter_cancer",
                               choices = tcga_cancers,
                               selected = tcga_cancers[1:3])
    })

    # ===== Custom Signature Tab =====

    shiny::observeEvent(input$create_sig_btn, {
      genes <- strsplit(input$sig_genes, "\n")[[1]]
      genes <- trimws(genes[genes != ""])

      if (length(genes) == 0) {
        shiny::showNotification("Please enter at least one gene", type = "error")
        return()
      }

      sig <- add_signature(
        name = input$sig_name,
        genes = genes,
        weights = rep(1, length(genes))
      )

      rv$signatures[[input$sig_name]] <- sig
      shiny::showNotification(paste("Signature", input$sig_name, "created!"), type = "message")
    })

    output$sig_preview <- shiny::renderPrint({
      if (length(rv$signatures) == 0) {
        cat("No signatures created yet.")
        return()
      }

      sig <- rv$signatures[[length(rv$signatures)]]
      if (!is.null(sig)) {
        cat("Signature Name:", sig$name, "\n")
        cat("Number of Genes:", sig$n_genes, "\n")
        cat("Genes:", paste(sig$genes, collapse = ", "), "\n")
        cat("Created:", format(sig$created), "\n")
      }
    })

    output$sig_table <- DT::renderDataTable({
      if (length(rv$signatures) == 0) {
        return(data.frame(Message = "No signatures available"))
      }

      sig_df <- purrr::map_dfr(rv$signatures, function(sig) {
        data.frame(
          Name = sig$name,
          Genes = sig$n_genes,
          Created = format(sig$created),
          stringsAsFactors = FALSE
        )
      })

      DT::datatable(sig_df, options = list(pageLength = 5))
    })

    # ===== Custom Metadata Tab =====

    shiny::observeEvent(input$meta_file, {
      req(input$meta_file)

      tryCatch({
        meta <- utils::read.csv(input$meta_file$datapath, stringsAsFactors = FALSE)

        if (!"sample" %in% colnames(meta)) {
          shiny::showNotification("Metadata must have 'sample' column", type = "error")
          return()
        }

        rv$metadata <- meta
        shiny::showNotification(paste("Loaded", nrow(meta), "metadata entries"), type = "message")
      }, error = function(e) {
        shiny::showNotification(paste("Error loading file:", conditionMessage(e)), type = "error")
      })
    })

    shiny::observeEvent(input$add_meta_btn, {
      if (input$meta_sample == "" || input$meta_group == "") {
        shiny::showNotification("Please enter both sample ID and group", type = "error")
        return()
      }

      new_entry <- data.frame(
        sample = input$meta_sample,
        group = input$meta_group,
        stringsAsFactors = FALSE
      )

      rv$metadata <- rbind(rv$metadata, new_entry)
      shiny::showNotification("Metadata entry added", type = "message")
    })

    output$meta_table <- DT::renderDataTable({
      if (nrow(rv$metadata) == 0) {
        return(data.frame(Message = "No metadata available"))
      }

      DT::datatable(rv$metadata, options = list(pageLength = 10, scrollX = TRUE))
    })

    output$download_meta_template <- shiny::downloadHandler(
      filename = function() "metadata_template.csv",
      content = function(file) {
        template <- data.frame(
          sample = c("TCGA-XX-XXXX-01", "TCGA-XX-XXXX-02"),
          group = c("High", "Low"),
          stringsAsFactors = FALSE
        )
        utils::write.csv(template, file, row.names = FALSE)
      }
    )

    # ===== Sample Filtering Tab =====

    shiny::observeEvent(input$apply_filter_btn, {
      # Get all samples
      tcga_gtex <- load_data("tcga_gtex")
      all_samples <- tcga_gtex$sample

      # Apply filters
      filters <- list(
        cancer_type = input$filter_cancer,
        sample_type = input$filter_sample_type
      )

      result <- filter_samples(all_samples, filters)
      rv$filtered_samples <- result$samples

      shiny::showNotification(
        paste("Filtered:", result$filtered_count, "of", result$original_count, "samples"),
        type = "message"
      )
    })

    shiny::observeEvent(input$reset_filter_btn, {
      rv$filtered_samples <- NULL
      shiny::showNotification("Filters reset", type = "message")
    })

    output$filter_stats <- shiny::renderPrint({
      if (is.null(rv$filtered_samples)) {
        cat("No filters applied.\n")
        cat("Use the filter controls to select samples.")
      } else {
        cat("Filtered Samples:", length(rv$filtered_samples), "\n")
        cat("Cancer Types:", paste(input$filter_cancer, collapse = ", "), "\n")
        cat("Sample Types:", paste(input$filter_sample_type, collapse = ", "))
      }
    })

    output$filtered_samples_table <- DT::renderDataTable({
      if (is.null(rv$filtered_samples)) {
        return(data.frame(Message = "No filtered samples"))
      }

      df <- data.frame(Sample = rv$filtered_samples, stringsAsFactors = FALSE)
      DT::datatable(df, options = list(pageLength = 10))
    })

    # ===== Data Download Tab =====

    output$download_preview <- shiny::renderPrint({
      cat("Download Settings:\n")
      cat("Data Type:", input$download_data_type, "\n")
      cat("Format:", input$download_format, "\n")
      cat("Filtered Only:", input$download_filtered_only, "\n")

      if (!is.null(rv$filtered_samples)) {
        cat("Sample Count:", length(rv$filtered_samples), "\n")
      }
    })

    output$data_summary <- shiny::renderTable({
      data.frame(
        Metric = c("Total Samples", "Filtered Samples", "Signatures", "Groups"),
        Count = c(
          nrow(load_data("tcga_gtex")),
          ifelse(is.null(rv$filtered_samples), 0, length(rv$filtered_samples)),
          length(rv$signatures),
          length(rv$groups)
        )
      )
    })

    output$download_data_btn <- shiny::downloadHandler(
      filename = function() {
        ext <- switch(input$download_format,
                      "csv" = ".csv",
                      "tsv" = ".tsv",
                      "rdata" = ".RData")
        paste0(input$download_data_type, "_data", ext)
      },
      content = function(file) {
        # Create sample data
        data <- data.frame(
          sample = if (!is.null(rv$filtered_samples)) rv$filtered_samples else c("sample1", "sample2"),
          value = c(1, 2),
          stringsAsFactors = FALSE
        )

        switch(input$download_format,
               "csv" = utils::write.csv(data, file, row.names = FALSE),
               "tsv" = utils::write.table(data, file, sep = "\t", row.names = FALSE),
               "rdata" = save(data, file = file)
        )
      }
    )

    # ===== Molecular Origin Tab =====

    shiny::observeEvent(input$set_origin_btn, {
      data_type <- switch(input$data_source,
                          "toil" = input$toil_data_type,
                          "pcawg" = input$pcawg_data_type,
                          "ccle" = input$ccle_data_type
      )

      rv$data_origin <- list(
        database = input$data_source,
        data_type = data_type
      )

      shiny::showNotification(
        paste("Data origin set to:", input$data_source, "-", data_type),
        type = "message"
      )
    })

    output$origin_settings <- shiny::renderPrint({
      cat("Current Data Origin:\n")
      cat("Database:", rv$data_origin$database, "\n")
      cat("Data Type:", rv$data_origin$data_type, "\n")
    })

    output$available_data_table <- DT::renderDataTable({
      data <- data.frame(
        Database = c("TCGA (TOIL)", "PCAWG", "CCLE"),
        DataTypes = c("mRNA, Protein, Methylation, miRNA, Mutation, CNV",
                      "mRNA, Fusion, Promoter, miRNA, APOBEC",
                      "mRNA, Protein, CNV, Mutation"),
        Samples = c("~11,000", "~2,800", "~1,000"),
        stringsAsFactors = FALSE
      )

      DT::datatable(data, options = list(pageLength = 5))
    })

    # ===== Sample Grouping Tab =====

    shiny::observeEvent(input$create_group_btn, {
      group_def <- list(
        name = input$group_name,
        by = input$group_by,
        samples = c("sample1", "sample2", "sample3")  # Placeholder
      )

      rv$groups[[input$group_name]] <- group_def
      shiny::showNotification(paste("Group", input$group_name, "created"), type = "message")
    })

    output$group_summary <- shiny::renderPrint({
      if (length(rv$groups) == 0) {
        cat("No groups defined yet.")
        return()
      }

      cat("Defined Groups:", length(rv$groups), "\n\n")
      for (name in names(rv$groups)) {
        group <- rv$groups[[name]]
        cat("-", name, "(", length(group$samples), "samples)\n")
      }
    })

    output$group_table <- DT::renderDataTable({
      if (length(rv$groups) == 0) {
        return(data.frame(Message = "No groups available"))
      }

      group_df <- purrr::map_dfr(names(rv$groups), function(name) {
        group <- rv$groups[[name]]
        data.frame(
          Group = name,
          DefinedBy = group$by,
          SampleCount = length(group$samples),
          stringsAsFactors = FALSE
        )
      })

      DT::datatable(group_df, options = list(pageLength = 5))
    })

    # Return reactive values for use by other modules
    return(rv)
  })
}
