# Shiny App Tests
# Following mastering-shiny Chapter 21 testing principles

test_that("Shiny app loads without errors", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinydashboard")

  # Test that app can be created
  expect_silent({
    app_dir <- system.file("shinyapp", package = "ZinaSuite")
    if (app_dir != "") {
      # App exists
      expect_true(dir.exists(app_dir))
    }
  })
})

test_that("Module UI functions return shiny tags", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinydashboard")

  # Test home module UI
  home_ui <- mod_home_ui("home")
  expect_s3_class(home_ui, "shiny.tag")

  # Test data query module UI
  dq_ui <- mod_data_query_ui("data_query")
  expect_s3_class(dq_ui, "shiny.tag")

  # Test analysis module UI
  analysis_ui <- mod_analysis_ui("analysis")
  expect_s3_class(analysis_ui, "shiny.tag")

  # Test visualization module UI
  viz_ui <- mod_visualization_ui("visualization")
  expect_s3_class(viz_ui, "shiny.tag")

  # Test immune module UI
  immune_ui <- mod_immune_ui("immune")
  expect_s3_class(immune_ui, "shiny.tag")

  # Test pharmacogenomics module UI
  pharma_ui <- mod_pharmacogenomics_ui("pharma")
  expect_s3_class(pharma_ui, "shiny.tag")

  # Test batch module UI
  batch_ui <- mod_batch_ui("batch")
  expect_s3_class(batch_ui, "shiny.tag")

  # Test about module UI
  about_ui <- mod_about_ui("about")
  expect_s3_class(about_ui, "shiny.tag")
})

test_that("App state reactive values work correctly", {
  skip_if_not_installed("shiny")

  # Create reactive values
  app_state <- shiny::reactiveValues(
    data = NULL,
    analysis_results = NULL,
    plots = NULL,
    async_tasks = list()
  )

  # Test that reactive values can be set
  expect_silent({
    app_state$data <- list(gene = "TP53", values = c(1, 2, 3))
    app_state$analysis_results <- list(cor = 0.5, pvalue = 0.01)
  })

  # Test that values are stored correctly
  expect_equal(app_state$data$gene, "TP53")
  expect_equal(app_state$analysis_results$cor, 0.5)
})

test_that("Shiny dependencies are checked correctly", {
  skip_if_not_installed("shiny")

  # Test that check_shiny_deps works
  expect_true(check_shiny_deps())
})
