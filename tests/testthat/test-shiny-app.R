# Shiny Application Tests - Updated for current UI structure
#'
#' Tests for basic Shiny app functionality

library(shinytest2)

# Test basic app loading
test_that("Shiny app loads successfully", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-basic",
    timeout = 30000
  )
  on.exit(app$stop())
  
  # Wait for app to be ready
  app$wait_for_idle()
  
  # Verify app is running
  expect_true(app$is_running())
})

# Test navigation to Quick Analysis - TCGA
test_that("Navigation to TCGA Quick Analysis works", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-nav-tcga",
    timeout = 30000
  )
  on.exit(app$stop())
  
  # Click on Quick Analysis menu first
  app$click("navbar-nav_menu_1")
  app$wait_for_idle()
  
  # Then click on TCGA Quick
  app$click("navbar-quick_tcga")
  app$wait_for_idle()
  
  # Verify we're on the right page by checking for expected elements
  expect_true(app$is_visible("quick_tcga-gene_input"))
})

# Test TCGA Quick Analysis - Tumor vs Normal
test_that("TCGA Quick - Tumor vs Normal analysis", {
  skip_if_offline()
  skip_on_cran()
  
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tcga-tvn",
    timeout = 60000
  )
  on.exit(app$stop())
  
  # Navigate to TCGA Quick
  app$click("navbar-nav_menu_1")
  app$wait_for_idle()
  app$click("navbar-quick_tcga")
  app$wait_for_idle()
  
  # Set analysis type to Tumor vs Normal
  app$set_inputs(`quick_tcga-analysis_type` = "tvn")
  app$set_inputs(`quick_tcga-gene` = "TP53")
  app$set_inputs(`quick_tcga-data_type` = "mRNA")
  
  # Run analysis
  app$click("quick_tcga-run_btn")
  app$wait_for_idle(timeout = 30000)
  
  # Verify results appear
  expect_true(app$is_visible("quick_tcga-plot_output"))
})

# Test TCGA Quick Analysis - Correlation
test_that("TCGA Quick - Correlation analysis", {
  skip_if_offline()
  skip_on_cran()
  
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tcga-cor",
    timeout = 60000
  )
  on.exit(app$stop())
  
  # Navigate to TCGA Quick
  app$click("navbar-nav_menu_1")
  app$wait_for_idle()
  app$click("navbar-quick_tcga")
  app$wait_for_idle()
  
  # Set analysis type to Correlation
  app$set_inputs(`quick_tcga-analysis_type` = "cor")
  app$set_inputs(`quick_tcga-gene1` = "TP53")
  app$set_inputs(`quick_tcga-gene2` = "BRCA1")
  app$set_inputs(`quick_tcga-data_type` = "mRNA")
  app$set_inputs(`quick_tcga-cor_method` = "spearman")
  
  # Run analysis
  app$click("quick_tcga-run_btn")
  app$wait_for_idle(timeout = 30000)
  
  # Verify results appear
  expect_true(app$is_visible("quick_tcga-plot_output"))
})

# Test navigation to other tabs
test_that("Navigation to all main tabs works", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-nav-all",
    timeout = 30000
  )
  on.exit(app$stop())
  
  # Test General Analysis
  app$click("navbar-general")
  app$wait_for_idle()
  expect_true(app$is_running())
  
  # Test Pan-Cancer
  app$click("navbar-pancan")
  app$wait_for_idle()
  expect_true(app$is_running())
  
  # Test Immune
  app$click("navbar-immune")
  app$wait_for_idle()
  expect_true(app$is_running())
  
  # Test Mutation
  app$click("navbar-mutation")
  app$wait_for_idle()
  expect_true(app$is_running())
  
  # Test About
  app$click("navbar-about")
  app$wait_for_idle()
  expect_true(app$is_running())
})

# Test Data Query functionality
test_that("Data Query tab works", {
  skip_if_offline()
  skip_on_cran()
  
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-data-query",
    timeout = 60000
  )
  on.exit(app$stop())
  
  # Navigate to Data Query
  app$click("navbar-data_query")
  app$wait_for_idle()
  
  # Set query parameters
  app$set_inputs(`data_query-source` = "tcga")
  app$set_inputs(`data_query-data_type` = "mRNA")
  app$set_inputs(`data_query-gene` = "TP53")
  
  # Run query
  app$click("data_query-query_btn")
  app$wait_for_idle(timeout = 30000)
  
  # Verify results table appears
  expect_true(app$is_visible("data_query-results_table"))
})
