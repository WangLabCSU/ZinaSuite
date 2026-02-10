# Test Shiny Modules with shinytest2

# Skip if shinytest2 is not available
skip_if_not_installed("shinytest2")

# Module UI Tests -------------------------------------------------------------

test_that("mod_home_ui returns valid UI", {
  ui <- mod_home_ui("home")
  expect_type(ui, "list")
  expect_true("shiny.tag" %in% class(ui) || "shiny.tag.list" %in% class(ui))
})

test_that("mod_quick_tcga_ui returns valid UI", {
  ui <- mod_quick_tcga_ui("quick_tcga")
  expect_type(ui, "list")
})

test_that("mod_pancan_ui returns valid UI", {
  ui <- mod_pancan_ui("pancan")
  expect_type(ui, "list")
})

test_that("mod_pharmacogenomics_ui returns valid UI", {
  ui <- mod_pharmacogenomics_ui("pharma")
  expect_type(ui, "list")
})

test_that("mod_data_query_ui returns valid UI", {
  ui <- mod_data_query_ui("data_query")
  expect_type(ui, "list")
})

# Module Server Logic Tests ---------------------------------------------------

test_that("mod_home_server initializes correctly", {
  testServer(mod_home_server, args = list(id = "home"), {
    # Test that server initializes without error
    expect_true(TRUE)
  })
})

test_that("mod_about_server initializes correctly", {
  testServer(mod_about_server, args = list(id = "about"), {
    expect_true(TRUE)
  })
})

# Shiny Logic Function Tests --------------------------------------------------

test_that("shiny logic functions work independently", {
  # Test query parameter building
  params <- build_query_params("TP53", "mRNA", "tcga")
  expect_equal(params$gene, "TP53")
  expect_equal(params$data_type, "mRNA")

  # Test validation
  validation <- validate_query_input("TP53", "mRNA", "tcga")
  expect_true(validation$valid)

  validation_invalid <- validate_query_input("", "mRNA", "tcga")
  expect_false(validation_invalid$valid)
})

# App State Management Tests --------------------------------------------------

test_that("app state can be created and modified", {
  # Create mock app state
  app_state <- shiny::reactiveValues(
    data = NULL,
    analysis_results = NULL,
    settings = list(source = "tcga", cancer = "all")
  )

  expect_null(app_state$data)
  expect_equal(app_state$settings$source, "tcga")

  # Modify state
  app_state$data <- data.frame(test = 1:10)
  app_state$settings$source <- "pcawg"

  expect_equal(nrow(app_state$data), 10)
  expect_equal(app_state$settings$source, "pcawg")
})

# Async Compute Integration Tests ---------------------------------------------

test_that("AsyncCompute can be initialized", {
  skip_if_not_installed("mirai")

  # Create async compute instance
  async <- AsyncCompute$new(n_workers = 1)
  expect_equal(async$info()$n_workers, 1)
  expect_false(async$info()$is_running)

  # Start and stop
  async$start()
  expect_true(async$info()$is_running)

  async$stop()
  expect_false(async$info()$is_running)
})

# Shiny App Structure Tests ---------------------------------------------------

test_that("app.R has correct structure", {
  app_file <- system.file("shinyapp", "app.R", package = "ZinaSuite")

  if (file.exists(app_file)) {
    app_code <- readLines(app_file)

    # Check for required components
    expect_true(any(grepl("ui <-", app_code, fixed = TRUE)))
    expect_true(any(grepl("server <- function", app_code, fixed = TRUE)))
    expect_true(any(grepl("shinyApp", app_code, fixed = TRUE)))
  } else {
    skip("app.R not found in installed package")
  }
})

# Module Integration Tests ----------------------------------------------------

test_that("all modules can be loaded", {
  # List of all module files
  modules_dir <- system.file("shinyapp", "modules", package = "ZinaSuite")

  if (dir.exists(modules_dir)) {
    module_files <- list.files(modules_dir, pattern = "\\.R$", full.names = TRUE)
    expect_gt(length(module_files), 0)

    # Check each module file can be parsed
    for (file in module_files) {
      expect_error(parse(file), NA, info = paste("Failed to parse:", basename(file)))
    }
  } else {
    skip("modules directory not found")
  }
})
