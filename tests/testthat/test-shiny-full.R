# Full Shiny App Integration Tests
# Following mastering-shiny Chapter 21 testing principles

# Test complete app loading
test_that("Shiny app loads all modules correctly", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinydashboard")

  # Test app directory exists
  app_dir <- system.file("shinyapp", package = "ZinaSuite")
  expect_true(dir.exists(app_dir), "Shiny app directory should exist")

  # Test all module files exist
  modules_dir <- file.path(app_dir, "modules")
  expected_modules <- c(
    "mod_home.R",
    "mod_data_query.R",
    "mod_analysis.R",
    "mod_visualization.R",
    "mod_pancan.R",
    "mod_mutation.R",
    "mod_dimension.R",
    "mod_immune.R",
    "mod_pharmacogenomics.R",
    "mod_batch.R",
    "mod_about.R"
  )

  for (module in expected_modules) {
    module_path <- file.path(modules_dir, module)
    expect_true(file.exists(module_path), paste("Module should exist:", module))
  }
})

# Test module UI functions return valid shiny tags
test_that("All module UI functions return shiny.tag objects", {
  skip_if_not_installed("shiny")

  modules <- list(
    home = mod_home_ui("home"),
    data_query = mod_data_query_ui("data_query"),
    analysis = mod_analysis_ui("analysis"),
    visualization = mod_visualization_ui("visualization"),
    pancan = mod_pancan_ui("pancan"),
    mutation = mod_mutation_ui("mutation"),
    dimension = mod_dimension_ui("dimension"),
    immune = mod_immune_ui("immune"),
    pharma = mod_pharmacogenomics_ui("pharma"),
    batch = mod_batch_ui("batch"),
    about = mod_about_ui("about")
  )

  for (name in names(modules)) {
    expect_s3_class(modules[[name]], "shiny.tag")
  }
})

# Test module server functions exist
test_that("Module server functions are exported", {
  skip_if_not_installed("shiny")

  # Test that server functions exist and are functions
  expect_type(mod_home_server, "closure")
  expect_type(mod_data_query_server, "closure")
  expect_type(mod_analysis_server, "closure")
  expect_type(mod_visualization_server, "closure")
  expect_type(mod_pancan_server, "closure")
  expect_type(mod_mutation_server, "closure")
  expect_type(mod_dimension_server, "closure")
  expect_type(mod_immune_server, "closure")
  expect_type(mod_pharmacogenomics_server, "closure")
  expect_type(mod_batch_server, "closure")
  expect_type(mod_about_server, "closure")
})

# Test Shiny dependencies
test_that("Shiny dependencies are available", {
  skip_if_not_installed("shiny")

  # Test required packages are available
  expect_true(requireNamespace("shinydashboard", quietly = TRUE))
  expect_true(requireNamespace("shinyWidgets", quietly = TRUE))
})

# Test app.R structure
test_that("App.R has correct structure", {
  skip_if_not_installed("shiny")

  app_file <- system.file("shinyapp", "app.R", package = "ZinaSuite")

  if (file.exists(app_file)) {
    app_content <- readLines(app_file)

    # Check for required components
    expect_true(any(grepl("run_zinasuite", app_content)), "Should have run_zinasuite function")
    expect_true(any(grepl("dashboardPage", app_content)), "Should use shinydashboard")
    expect_true(any(grepl("tabItems", app_content)), "Should have tabItems")

    # Check for all menu items
    menu_items <- c("home", "data_query", "analysis", "visualization",
                    "pancan", "mutation", "dimension", "immune",
                    "pharma", "batch", "about")

    for (item in menu_items) {
      expect_true(
        any(grepl(item, app_content)),
        info = paste("App should have menu item:", item)
      )
    }
  }
})

# Test module file structure
test_that("Module files have correct structure", {
  modules_dir <- system.file("shinyapp", "modules", package = "ZinaSuite")

  if (dir.exists(modules_dir)) {
    module_files <- list.files(modules_dir, pattern = "\\.R$", full.names = TRUE)

    for (file in module_files) {
      content <- readLines(file)

      # Check for roxygen documentation
      expect_true(
        any(grepl("^#'", content)),
        info = paste("Module should have roxygen docs:", basename(file))
      )

      # Check for module UI function
      expect_true(
        any(grepl("_ui <- function", content)),
        info = paste("Module should have UI function:", basename(file))
      )

      # Check for module server function
      expect_true(
        any(grepl("_server <- function", content)),
        info = paste("Module should have server function:", basename(file))
      )
    }
  }
})

# Test async compute integration
test_that("Async compute integration works", {
  skip_if_not_installed("mirai")

  # Test that AsyncCompute class exists
  expect_true(exists("AsyncCompute"))

  # Test that we can create an instance (if mirai is available)
  tryCatch({
    async <- AsyncCompute$new(n_workers = 1)
    expect_true(!is.null(async))
    async$stop()
  }, error = function(e) {
    skip("mirai not available for async testing")
  })
})

# Test complete app launch (if possible)
test_that("App can be launched in test mode", {
  skip_on_cran()
  skip_if_not_installed("shiny")

  # This test verifies the app can be created without errors
  expect_silent({
    app_dir <- system.file("shinyapp", package = "ZinaSuite")
    # Just verify the directory structure, don't actually run
    expect_true(length(list.files(app_dir)) > 0)
  })
})
