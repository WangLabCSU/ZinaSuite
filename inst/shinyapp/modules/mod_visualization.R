#' Visualization Module UI
mod_visualization_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shinydashboard::box(
      title = "Plot Settings", width = 4, status = "primary", solidHeader = TRUE,
      shiny::selectInput(ns("plot_type"), "Plot Type:",
        choices = c("Histogram" = "histogram", "Boxplot" = "boxplot", "Scatter" = "scatter")
      ),
      shiny::actionButton(ns("plot_btn"), "Generate Plot", icon = shiny::icon("chart-bar"), class = "btn-primary")
    ),
    shinydashboard::box(
      title = "Plot", width = 8, status = "info", solidHeader = TRUE,
      shiny::plotOutput(ns("plot"), height = "500px")
    )
  )
}

#' Visualization Module Server
mod_visualization_server <- function(id, app_state) {
  shiny::moduleServer(id, function(input, output, session) {
    output$plot <- shiny::renderPlot({
      shiny::req(input$plot_btn)
      shiny::req(app_state$data)
      
      values <- as.numeric(app_state$data$values)
      viz <- Visualization$new()
      
      switch(input$plot_type,
        "histogram" = viz$plot_distribution(values, type = "histogram"),
        "boxplot" = viz$plot_distribution(values, type = "boxplot"),
        "scatter" = {
          shiny::req(length(values) > 1)
          viz$plot_correlation(1:length(values), values, type = "scatter")
        }
      )
    })
  })
}
