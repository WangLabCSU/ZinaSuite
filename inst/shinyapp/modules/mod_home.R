#' Home Module UI
#'
#' @param id Module ID
#' @return UI elements
mod_home_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::card(
      bslib::card_header(
        class = "bg-primary text-white",
        shiny::icon("home"), "Welcome to ZinaSuite"
      ),
      bslib::card_body(
        shiny::h3("Interactive Analysis of UCSC Xena Data"),
        shiny::p("ZinaSuite is a comprehensive R package for querying and analyzing cancer genomics data from UCSC Xena."),

        shiny::hr(),

        shiny::h4("Key Features:"),
        shiny::tags$ul(
          shiny::tags$li(shiny::strong("Data Query:"), " Access TCGA, PCAWG, and CCLE data"),
          shiny::tags$li(shiny::strong("Parallel Processing:"), " Fast analysis using mirai"),
          shiny::tags$li(shiny::strong("Analysis Tools:"), " Correlation, survival, and more"),
          shiny::tags$li(shiny::strong("Visualization:"), " Publication-ready plots")
        ),

        shiny::hr(),

        shiny::h4("Quick Start:"),
        shiny::p("Use the navigation menu to explore different modules:", shiny::tags$br(),
                 shiny::icon("bolt"), " ", shiny::strong("Quick Analysis"), " - Fast access to common analyses", shiny::tags$br(),
                 shiny::icon("chart-line"), " ", shiny::strong("General Analysis"), " - Statistical analysis tools", shiny::tags$br(),
                 shiny::icon("globe"), " ", shiny::strong("Pan-Cancer"), " - Cross-cancer analysis", shiny::tags$br(),
                 shiny::icon("tasks"), " ", shiny::strong("Batch Analysis"), " - Analyze multiple genes in parallel")
      )
    ),

    bslib::layout_column_wrap(
      width = 1/3,
      bslib::card(
        bslib::card_header("TCGA Data"),
        shiny::h5("The Cancer Genome Atlas"),
        shiny::p("Comprehensive cancer genomics dataset with over 11,000 patients across 33 cancer types."),
        shiny::tags$ul(
          shiny::tags$li("Gene expression (TPM/FPKM)"),
          shiny::tags$li("Mutation status"),
          shiny::tags$li("Copy number variation"),
          shiny::tags$li("Clinical data")
        )
      ),
      bslib::card(
        bslib::card_header("Analysis Capabilities"),
        shiny::h5("Available Analyses"),
        shiny::tags$ul(
          shiny::tags$li(shiny::strong("Correlation:"), " Gene-gene correlation analysis"),
          shiny::tags$li(shiny::strong("Survival:"), " Kaplan-Meier and Cox regression"),
          shiny::tags$li(shiny::strong("Comparison:"), " Group comparison tests"),
          shiny::tags$li(shiny::strong("Dimensionality Reduction:"), " PCA, t-SNE, UMAP")
        )
      ),
      bslib::card(
        bslib::card_header("Parallel Computing"),
        shiny::h5("Powered by mirai"),
        shiny::p("ZinaSuite uses the mirai package for asynchronous and parallel computing."),
        shiny::tags$ul(
          shiny::tags$li("Fast batch gene queries"),
          shiny::tags$li("Parallel correlation analysis"),
          shiny::tags$li("Async survival analysis"),
          shiny::tags$li("Progress tracking")
        )
      )
    )
  )
}

#' Home Module Server
#'
#' @param id Module ID
mod_home_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # Home module logic (if needed)
  })
}
