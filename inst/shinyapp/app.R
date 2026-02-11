#' ZinaSuite Shiny Application
#'
#' @description
#' Main Shiny application for ZinaSuite - Interactive analysis of UCSC Xena data.
#' Features a modern bslib interface with comprehensive analysis modules.

# Source all module files
# Use tryCatch to handle different execution contexts
modules_dir <- tryCatch({
  file.path(dirname(sys.frame(1)$ofile), "modules")
}, error = function(e) {
  # Fallback for when sys.frame(1)$ofile is not available
  file.path(system.file("shinyapp", package = "ZinaSuite"), "modules")
})

if (!is.null(modules_dir) && dir.exists(modules_dir)) {
  module_files <- list.files(modules_dir, pattern = "\\.R$", full.names = TRUE)
  for (f in module_files) {
    source(f, local = TRUE)
  }
}

# UI Definition
ui <- bslib::page_navbar(
  title = htmltools::tags$span(
    htmltools::tags$img(src = "logo.png", height = "30px", style = "margin-right: 10px;"),
    "ZinaSuite"
  ),
  bg = "#2C3E50",
  fillable = TRUE,
  window_title = "ZinaSuite - UCSC Xena Data Analysis",
  navbar_options = bslib::navbar_options(
    bg = "#2C3E50",
    theme = "dark"
  ),

  # Home Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("home"),
      "Home"
    ),
    value = "home",
    mod_home_ui("home")
  ),

  # Quick Analysis Tab
  bslib::nav_menu(
    title = bslib::tooltip(
      shiny::icon("bolt"),
      "Quick Analysis"
    ),
    icon = shiny::icon("bolt"),
    bslib::nav_panel(
      title = "TCGA Quick",
      value = "quick_tcga",
      mod_quick_tcga_ui("quick_tcga")
    ),
    bslib::nav_panel(
      title = "PCAWG Quick",
      value = "quick_pcawg",
      mod_quick_pcawg_ui("quick_pcawg")
    ),
    bslib::nav_panel(
      title = "CCLE Quick",
      value = "quick_ccle",
      mod_quick_ccle_ui("quick_ccle")
    )
  ),

  # General Analysis Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("chart-line"),
      "General Analysis"
    ),
    value = "general",
    mod_general_analysis_ui("general")
  ),

  # Pan-Cancer Analysis Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("globe"),
      "Pan-Cancer"
    ),
    value = "pancan",
    mod_pancan_ui("pancan")
  ),

  # Immune Analysis Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("shield-virus"),
      "Immune"
    ),
    value = "immune",
    mod_immune_ui("immune")
  ),

  # Mutation Analysis Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("dna"),
      "Mutation"
    ),
    value = "mutation",
    mod_mutation_ui("mutation")
  ),

  # Dimension Reduction Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("project-diagram"),
      "Dimension"
    ),
    value = "dimension",
    mod_dimension_ui("dimension")
  ),

  # PharmacoGenomics Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("pills"),
      "PharmacoGenomics"
    ),
    value = "pharma",
    mod_pharmacogenomics_ui("pharma")
  ),

  # Batch Analysis Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("tasks"),
      "Batch"
    ),
    value = "batch",
    mod_batch_ui("batch")
  ),

  # Data Query Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("database"),
      "Data Query"
    ),
    value = "data_query",
    mod_data_query_ui("data_query")
  ),

  # About Tab
  bslib::nav_panel(
    title = bslib::tooltip(
      shiny::icon("info-circle"),
      "About"
    ),
    value = "about",
    mod_about_ui("about")
  ),

  # Sidebar with additional controls
  sidebar = bslib::sidebar(
    open = "closed",
    width = 300,
    shiny::h5("Analysis Settings"),
    shiny::hr(),
    shiny::selectInput(
      "global_source",
      "Default Data Source:",
      choices = c("TCGA" = "tcga", "PCAWG" = "pcawg", "CCLE" = "ccle"),
      selected = "tcga"
    ),
    shiny::selectInput(
      "global_cancer",
      "Default Cancer Type:",
      choices = c("All" = "all", "BRCA", "LUAD", "LUSC", "COAD", "READ"),
      selected = "all"
    ),
    shiny::hr(),
    shiny::h5("Cache Management"),
    shiny::actionButton(
      "clear_cache",
      "Clear Cache",
      icon = shiny::icon("trash"),
      class = "btn-warning btn-sm"
    ),
    shiny::hr(),
    shiny::h5("Help"),
    shiny::actionButton(
      "show_help",
      "Show Documentation",
      icon = shiny::icon("book"),
      class = "btn-info btn-sm"
    )
  ),

  # Footer
  footer = shiny::tags$footer(
    class = "bg-light py-3 mt-auto",
    shiny::tags$div(
      class = "container text-center text-muted",
      paste0("ZinaSuite v", as.character(packageVersion("ZinaSuite"))),
      " | Powered by ",
      shiny::tags$a(href = "https://bioconductor.org", "Bioconductor"),
      " | Data from ",
      shiny::tags$a(href = "https://xena.ucsc.edu", "UCSC Xena")
    )
  )
)

# Server Definition
server <- function(input, output, session) {
  # Initialize reactive values for shared state
  app_state <- shiny::reactiveValues(
    data = NULL,
    analysis_results = NULL,
    plots = NULL,
    async_tasks = list(),
    settings = list(
      source = "tcga",
      cancer = "all"
    )
  )

  # Initialize async compute engine
  async_compute <- NULL
  tryCatch({
    async_compute <- AsyncCompute$new(n_workers = 2)
    async_compute$start()
    shiny::showNotification(
      "Async compute engine initialized",
      type = "message",
      duration = 3
    )
  }, error = function(e) {
    shiny::showNotification(
      "Warning: Could not initialize async compute engine",
      type = "warning"
    )
  })

  # Cleanup on session end
  shiny::onSessionEnded(function() {
    if (!is.null(async_compute)) {
      async_compute$stop()
    }
  })

  # Observe global settings
  shiny::observeEvent(input$global_source, {
    app_state$settings$source <- input$global_source
  })

  shiny::observeEvent(input$global_cancer, {
    app_state$settings$cancer <- input$global_cancer
  })

  # Clear cache button
  shiny::observeEvent(input$clear_cache, {
    tryCatch({
      clear_zina_cache()
      shiny::showNotification("Cache cleared successfully", type = "message")
    }, error = function(e) {
      shiny::showNotification(paste("Error clearing cache:", e$message), type = "error")
    })
  })

  # Show help button
  shiny::observeEvent(input$show_help, {
    shiny::showModal(shiny::modalDialog(
      title = "ZinaSuite Documentation",
      shiny::includeMarkdown(system.file("doc", "introduction.Rmd", package = "ZinaSuite")),
      easyClose = TRUE,
      footer = shiny::modalButton("Close")
    ))
  })

  # Call module servers
  mod_home_server("home")
  mod_quick_tcga_server("quick_tcga", app_state, async_compute)
  mod_quick_pcawg_server("quick_pcawg", app_state, async_compute)
  mod_quick_ccle_server("quick_ccle", app_state, async_compute)
  mod_general_analysis_server("general", app_state, async_compute)
  mod_pancan_server("pancan", app_state, async_compute)
  mod_immune_server("immune", app_state, async_compute)
  mod_mutation_server("mutation", app_state, async_compute)
  mod_dimension_server("dimension", app_state, async_compute)
  mod_pharmacogenomics_server("pharma", app_state, async_compute)
  mod_batch_server("batch", app_state, async_compute)
  mod_data_query_server("data_query", app_state, async_compute)
  mod_about_server("about")
}

# Create and return shiny app object
shiny::shinyApp(ui = ui, server = server)
