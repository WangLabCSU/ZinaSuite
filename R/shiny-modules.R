#' Shiny Module UI Functions
#'
#' @description
#' UI functions for all Shiny modules in ZinaSuite.
#' These are re-exported from the inst/shinyapp/modules/ directory.
#'
#' @name shiny-module-ui
#' @keywords internal
NULL

#' @rdname shiny-module-ui
#' @export
mod_home_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_home.R", package = "ZinaSuite"), local = TRUE)
  mod_home_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_data_query_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_data_query.R", package = "ZinaSuite"), local = TRUE)
  mod_data_query_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_analysis_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_analysis.R", package = "ZinaSuite"), local = TRUE)
  mod_analysis_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_visualization_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_visualization.R", package = "ZinaSuite"), local = TRUE)
  mod_visualization_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_pancan_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_pancan.R", package = "ZinaSuite"), local = TRUE)
  mod_pancan_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_mutation_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_mutation.R", package = "ZinaSuite"), local = TRUE)
  mod_mutation_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_dimension_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_dimension.R", package = "ZinaSuite"), local = TRUE)
  mod_dimension_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_immune_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_immune.R", package = "ZinaSuite"), local = TRUE)
  mod_immune_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_pharmacogenomics_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_pharmacogenomics.R", package = "ZinaSuite"), local = TRUE)
  mod_pharmacogenomics_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_batch_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_batch.R", package = "ZinaSuite"), local = TRUE)
  mod_batch_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_about_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_about.R", package = "ZinaSuite"), local = TRUE)
  mod_about_ui(id)
}

#' @rdname shiny-module-ui
#' @export
mod_quick_tcga_ui <- function(id) {
  source(system.file("shinyapp", "modules", "mod_quick_tcga.R", package = "ZinaSuite"), local = TRUE)
  mod_quick_tcga_ui(id)
}

#' Shiny Module Server Functions
#'
#' @description
#' Server functions for all Shiny modules in ZinaSuite.
#'
#' @name shiny-module-server
#' @keywords internal
NULL

#' @rdname shiny-module-server
#' @export
mod_home_server <- function(id) {
  source(system.file("shinyapp", "modules", "mod_home.R", package = "ZinaSuite"), local = TRUE)
  mod_home_server(id)
}

#' @rdname shiny-module-server
#' @export
mod_data_query_server <- function(id, app_state) {
  source(system.file("shinyapp", "modules", "mod_data_query.R", package = "ZinaSuite"), local = TRUE)
  mod_data_query_server(id, app_state)
}

#' @rdname shiny-module-server
#' @export
mod_analysis_server <- function(id, app_state) {
  source(system.file("shinyapp", "modules", "mod_analysis.R", package = "ZinaSuite"), local = TRUE)
  mod_analysis_server(id, app_state)
}

#' @rdname shiny-module-server
#' @export
mod_visualization_server <- function(id, app_state) {
  source(system.file("shinyapp", "modules", "mod_visualization.R", package = "ZinaSuite"), local = TRUE)
  mod_visualization_server(id, app_state)
}

#' @rdname shiny-module-server
#' @export
mod_pancan_server <- function(id, app_state) {
  source(system.file("shinyapp", "modules", "mod_pancan.R", package = "ZinaSuite"), local = TRUE)
  mod_pancan_server(id, app_state)
}

#' @rdname shiny-module-server
#' @export
mod_mutation_server <- function(id, app_state) {
  source(system.file("shinyapp", "modules", "mod_mutation.R", package = "ZinaSuite"), local = TRUE)
  mod_mutation_server(id, app_state)
}

#' @rdname shiny-module-server
#' @export
mod_dimension_server <- function(id, app_state) {
  source(system.file("shinyapp", "modules", "mod_dimension.R", package = "ZinaSuite"), local = TRUE)
  mod_dimension_server(id, app_state)
}

#' @rdname shiny-module-server
#' @export
mod_immune_server <- function(id, app_state) {
  source(system.file("shinyapp", "modules", "mod_immune.R", package = "ZinaSuite"), local = TRUE)
  mod_immune_server(id, app_state)
}

#' @rdname shiny-module-server
#' @export
mod_pharmacogenomics_server <- function(id, app_state) {
  source(system.file("shinyapp", "modules", "mod_pharmacogenomics.R", package = "ZinaSuite"), local = TRUE)
  mod_pharmacogenomics_server(id, app_state)
}

#' @rdname shiny-module-server
#' @export
mod_batch_server <- function(id, app_state, async_compute) {
  source(system.file("shinyapp", "modules", "mod_batch.R", package = "ZinaSuite"), local = TRUE)
  mod_batch_server(id, app_state, async_compute)
}

#' @rdname shiny-module-server
#' @export
mod_about_server <- function(id) {
  source(system.file("shinyapp", "modules", "mod_about.R", package = "ZinaSuite"), local = TRUE)
  mod_about_server(id)
}

#' @rdname shiny-module-server
#' @export
mod_quick_tcga_server <- function(id, app_state, async_compute) {
  source(system.file("shinyapp", "modules", "mod_quick_tcga.R", package = "ZinaSuite"), local = TRUE)
  mod_quick_tcga_server(id, app_state, async_compute)
}
