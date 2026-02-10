#' Run ZinaSuite Shiny Application
#'
#' @description
#' Launch the ZinaSuite Shiny application for interactive analysis of UCSC Xena data.
#'
#' @param port Port number for the Shiny application (default: random available port)
#' @param host Host IP address (default: "127.0.0.1")
#' @param launch.browser Whether to launch browser automatically (default: TRUE)
#' @param ... Additional arguments passed to shiny::runApp
#' @return Shiny application object (invisible)
#' @export
#'
#' @examples
#' \donttest{
#' # Run the application with default settings
#' run_zinasuite()
#'
#' # Run on a specific port
#' run_zinasuite(port = 8080)
#'
#' # Run without launching browser
#' run_zinasuite(launch.browser = FALSE)
#' }
run_zinasuite <- function(port = getOption("shiny.port"),
                          host = getOption("shiny.host", "127.0.0.1"),
                          launch.browser = TRUE,
                          ...) {
  # Check if shiny is installed
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the application. ",
         "Please install it with: install.packages('shiny')")
  }

  # Get app directory
  app_dir <- system.file("shinyapp", package = "ZinaSuite")

  if (app_dir == "") {
    stop("Shiny application not found. Please ensure ZinaSuite is properly installed.")
  }

  # Run the application
  shiny::runApp(
    appDir = app_dir,
    port = port,
    host = host,
    launch.browser = launch.browser,
    ...
  )
}

#' Get Shiny App Directory
#'
#' @description
#' Get the file path to the Shiny application directory.
#'
#' @return Character string with the app directory path
#' @export
#'
#' @examples
#' \donttest{
#' app_dir <- get_app_dir()
#' list.files(app_dir)
#' }
get_app_dir <- function() {
  app_dir <- system.file("shinyapp", package = "ZinaSuite")

  if (app_dir == "") {
    stop("Shiny application not found. Please ensure ZinaSuite is properly installed.")
  }

  app_dir
}
