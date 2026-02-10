#' Check if Suggested Package is Installed
#'
#' @description
#' Internal helper to check if a suggested package is installed.
#' Uses rlang::check_installed for better error messages.
#'
#' @param pkg Package name to check
#' @param reason Optional reason why the package is needed
#' @return TRUE if package is installed, otherwise throws error
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' check_suggested("survival")
#' }
check_suggested <- function(pkg, reason = NULL) {
  rlang::check_installed(pkg, reason = reason)
  TRUE
}

#' Check Multiple Suggested Packages
#'
#' @description
#' Check if multiple suggested packages are installed.
#'
#' @param pkgs Vector of package names
#' @param reason Optional reason why the packages are needed
#' @return TRUE if all packages are installed
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' check_suggested_pkgs(c("survival", "survminer"))
#' }
check_suggested_pkgs <- function(pkgs, reason = NULL) {
  for (pkg in pkgs) {
    check_suggested(pkg, reason)
  }
  TRUE
}

#' Check Shiny Dependencies
#'
#' @description
#' Check if all Shiny-related dependencies are installed.
#'
#' @return TRUE if all packages are installed
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' check_shiny_deps()
#' }
check_shiny_deps <- function() {
  shiny_pkgs <- c(
    "shiny", "shinydashboard", "shinyWidgets",
    "shinyjs", "shinyBS", "shinycssloaders", "waiter"
  )
  check_suggested_pkgs(
    shiny_pkgs,
    reason = "required to run the ZinaSuite Shiny application"
  )
}

#' Check Visualization Dependencies
#'
#' @description
#' Check if visualization-related dependencies are installed.
#'
#' @param type Type of visualization: "basic", "survival", "interactive", "heatmap"
#' @return TRUE if all packages are installed
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' check_vis_deps("survival")
#' }
check_vis_deps <- function(type = "basic") {
  vis_pkgs <- switch(type,
    "basic" = c("ggplot2"),
    "survival" = c("survival", "survminer", "ggpubr"),
    "interactive" = c("plotly"),
    "heatmap" = c("ComplexHeatmap", "RColorBrewer"),
    "stats" = c("ggstatsplot"),
    stop("Unknown visualization type: ", type)
  )

  check_suggested_pkgs(
    vis_pkgs,
    reason = paste("required for", type, "visualizations")
  )
}

#' Check Analysis Dependencies
#'
#' @description
#' Check if analysis-related dependencies are installed.
#'
#' @param type Type of analysis: "survival", "dimension", "table"
#' @return TRUE if all packages are installed
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' check_analysis_deps("survival")
#' }
check_analysis_deps <- function(type = "survival") {
  analysis_pkgs <- switch(type,
    "survival" = c("survival", "survminer"),
    "dimension" = c("Rtsne", "umap"),
    "table" = c("DT"),
    stop("Unknown analysis type: ", type)
  )

  check_suggested_pkgs(
    analysis_pkgs,
    reason = paste("required for", type, "analysis")
  )
}
