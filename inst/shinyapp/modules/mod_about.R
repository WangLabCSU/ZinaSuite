# About Module UI
mod_about_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shinydashboard::box(
      title = "About ZinaSuite",
      width = 12,
      status = "primary",
      solidHeader = TRUE,

      shiny::h3("ZinaSuite: Interactive Analysis of UCSC Xena Data"),
      shiny::p("Version: 0.1.0"),
      shiny::p("A comprehensive R package for querying and analyzing cancer genomics data from UCSC Xena."),

      shiny::hr(),

      shiny::h4("Features:"),
      shiny::tags$ul(
        shiny::tags$li("Query gene expression, mutation, CNV data from TCGA, PCAWG, CCLE"),
        shiny::tags$li("Parallel processing using mirai for fast batch analysis"),
        shiny::tags$li("Correlation, survival, and group comparison analyses"),
        shiny::tags$li("Publication-ready visualizations")
      ),

      shiny::hr(),

      shiny::h4("Citation:"),
      shiny::p("If you use ZinaSuite in your research, please cite:"),
      shiny::tags$blockquote(
        "ZinaSuite: An R package for interactive analysis of UCSC Xena data."
      )
    )
  )
}

# About Module Server
mod_about_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # About module logic (if needed)
  })
}
