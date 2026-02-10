#' ZinaSuite Shiny Application
#'
#' @description
#' Main Shiny application for ZinaSuite - Interactive analysis of UCSC Xena data.
#'
#' @export
#' @import shiny
#' @import shinyWidgets
#' @import shinydashboard
#'
#' @examples
#' \dontrun{
#' # Launch the Shiny app
#' run_zinasuite()
#' }
run_zinasuite <- function() {
  # Check required packages
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install with: install.packages('shiny')")
  }
  if (!requireNamespace("shinydashboard", quietly = TRUE)) {
    stop("Package 'shinydashboard' is required. Install with: install.packages('shinydashboard')")
  }

  # Source modules
  source_modules <- function() {
    modules_dir <- system.file("shinyapp", "modules", package = "ZinaSuite")
    if (dir.exists(modules_dir)) {
      module_files <- list.files(modules_dir, pattern = "\\.R$", full.names = TRUE)
      for (f in module_files) {
        source(f, local = TRUE)
      }
    }
  }

  # UI
  ui <- shinydashboard::dashboardPage(
    skin = "blue",

    # Header
    shinydashboard::dashboardHeader(
      title = "ZinaSuite",
      titleWidth = 250
    ),

    # Sidebar
    shinydashboard::dashboardSidebar(
      width = 250,
      shinydashboard::sidebarMenu(
        id = "sidebar_menu",
        shinydashboard::menuItem("Home", tabName = "home", icon = shiny::icon("home")),
        shinydashboard::menuItem("Data Query", tabName = "data_query", icon = shiny::icon("database")),
        shinydashboard::menuItem("Analysis", tabName = "analysis", icon = shiny::icon("chart-line")),
        shinydashboard::menuItem("Visualization", tabName = "visualization", icon = shiny::icon("chart-bar")),
        shinydashboard::menuItem("Batch Analysis", tabName = "batch", icon = shiny::icon("tasks")),
        shinydashboard::menuItem("About", tabName = "about", icon = shiny::icon("info-circle"))
      )
    ),

    # Body
    shinydashboard::dashboardBody(
      shiny::tags$head(
        shiny::tags$style(shiny::HTML("
          .content-wrapper { background-color: #f4f6f9; }
          .box { border-top: 3px solid #3c8dbc; }
        "))
      ),

      shinydashboard::tabItems(
        # Home tab
        shinydashboard::tabItem(
          tabName = "home",
          mod_home_ui("home")
        ),

        # Data Query tab
        shinydashboard::tabItem(
          tabName = "data_query",
          mod_data_query_ui("data_query")
        ),

        # Analysis tab
        shinydashboard::tabItem(
          tabName = "analysis",
          mod_analysis_ui("analysis")
        ),

        # Visualization tab
        shinydashboard::tabItem(
          tabName = "visualization",
          mod_visualization_ui("visualization")
        ),

        # Batch Analysis tab
        shinydashboard::tabItem(
          tabName = "batch",
          mod_batch_ui("batch")
        ),

        # About tab
        shinydashboard::tabItem(
          tabName = "about",
          mod_about_ui("about")
        )
      )
    )
  )

  # Server
  server <- function(input, output, session) {
    # Source modules in server context
    source_modules()

    # Initialize reactive values for shared state
    app_state <- shiny::reactiveValues(
      data = NULL,
      analysis_results = NULL,
      plots = NULL,
      async_tasks = list()
    )

    # Initialize async compute engine
    async_compute <- NULL
    tryCatch({
      async_compute <- AsyncCompute$new(n_workers = 2)
      async_compute$start()
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

    # Call module servers
    mod_home_server("home")
    mod_data_query_server("data_query", app_state)
    mod_analysis_server("analysis", app_state)
    mod_visualization_server("visualization", app_state)
    mod_batch_server("batch", app_state, async_compute)
    mod_about_server("about")
  }

  # Run app
  shiny::shinyApp(ui = ui, server = server)
}
