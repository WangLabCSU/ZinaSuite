#' Shiny Test Suite for ZinaSuite
#'
#' Comprehensive shinytest2 tests covering all core functionality
#'

library(shinytest2)

# Test basic app initialization
test_that("App initializes correctly", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-init",
    timeout = 30000
  )
  on.exit(app$stop())

  # Check initial state
  app$expect_values()
})

# Test Home Module
test_that("Home module displays correctly", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-home",
    timeout = 30000
  )
  on.exit(app$stop())

  # Navigate to home
  app$click("nav-home")
  app$wait_for_idle()

  # Verify home content
  expect_true(app$get_html(".card-header") %>% grepl("Welcome to ZinaSuite", .))
})

# Test TCGA Quick Analysis Module
test_that("TCGA Quick Analysis - Tumor vs Normal", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tcga-tn",
    timeout = 60000
  )
  on.exit(app$stop())

  # Navigate to TCGA Quick tab
  app$click("nav-quick_tcga")
  app$wait_for_idle()

  # Set inputs
  app$set_inputs(`quick_tcga-gene` = "TP53")
  app$set_inputs(`quick_tcga-analysis_type` = "tumor_normal")
  app$set_inputs(`quick_tcga-data_type` = "mRNA")
  app$set_inputs(`quick_tcga-mode` = "pancan")

  # Run analysis
  app$click("quick_tcga-run_btn")
  app$wait_for_idle(timeout = 30000)

  # Verify results
  expect_true(app$is_visible("quick_tcga-plot_output"))
})

# Test TCGA Quick Analysis - Gene Correlation
test_that("TCGA Quick Analysis - Gene Correlation", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tcga-cor",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-quick_tcga")
  app$wait_for_idle()

  app$set_inputs(`quick_tcga-gene` = "TP53")
  app$set_inputs(`quick_tcga-gene2` = "BRCA1")
  app$set_inputs(`quick_tcga-analysis_type` = "correlation")
  app$set_inputs(`quick_tcga-data_type` = "mRNA")

  app$click("quick_tcga-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("quick_tcga-plot_output"))
})

# Test TCGA Quick Analysis - TIL Correlation
test_that("TCGA Quick Analysis - TIL Correlation", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tcga-til",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-quick_tcga")
  app$wait_for_idle()

  app$set_inputs(`quick_tcga-gene` = "TP53")
  app$set_inputs(`quick_tcga-analysis_type` = "til")
  app$set_inputs(`quick_tcga-data_type` = "mRNA")
  app$set_inputs(`quick_tcga-cor_method` = "spearman")

  app$click("quick_tcga-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("quick_tcga-plot_output"))
})

# Test TCGA Quick Analysis - Immune Correlation
test_that("TCGA Quick Analysis - Immune Correlation", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tcga-immune",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-quick_tcga")
  app$wait_for_idle()

  app$set_inputs(`quick_tcga-gene` = "TP53")
  app$set_inputs(`quick_tcga-analysis_type` = "immune")
  app$set_inputs(`quick_tcga-data_type` = "mRNA")
  app$set_inputs(`quick_tcga-cor_method` = "spearman")

  app$click("quick_tcga-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("quick_tcga-plot_output"))
})

# Test TCGA Quick Analysis - Survival KM
test_that("TCGA Quick Analysis - Survival KM", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tcga-km",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-quick_tcga")
  app$wait_for_idle()

  app$set_inputs(`quick_tcga-gene` = "TP53")
  app$set_inputs(`quick_tcga-analysis_type` = "survival_km")
  app$set_inputs(`quick_tcga-data_type` = "mRNA")
  app$set_inputs(`quick_tcga-surv_measure` = "OS")

  app$click("quick_tcga-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("quick_tcga-plot_output"))
})

# Test PCAWG Quick Analysis Module
test_that("PCAWG Quick Analysis - Tumor vs Normal", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pcawg-tn",
    timeout = 60000
  )
  on.exit(app$stop())

  # Navigate to PCAWG Quick tab
  app$click("nav-quick_pcawg")
  app$wait_for_idle()

  app$set_inputs(`quick_pcawg-gene` = "TP53")
  app$set_inputs(`quick_pcawg-analysis_type` = "tumor_normal")
  app$set_inputs(`quick_pcawg-data_type` = "mRNA")

  app$click("quick_pcawg-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("quick_pcawg-plot_output"))
})

# Test PCAWG Quick Analysis - Gene Correlation
test_that("PCAWG Quick Analysis - Gene Correlation", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pcawg-cor",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-quick_pcawg")
  app$wait_for_idle()

  app$set_inputs(`quick_pcawg-gene` = "TP53")
  app$set_inputs(`quick_pcawg-gene2` = "BRCA1")
  app$set_inputs(`quick_pcawg-analysis_type` = "correlation")
  app$set_inputs(`quick_pcawg-data_type` = "mRNA")

  app$click("quick_pcawg-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("quick_pcawg-plot_output"))
})

# Test CCLE Quick Analysis Module
test_that("CCLE Quick Analysis - Expression Distribution", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-ccle-dist",
    timeout = 60000
  )
  on.exit(app$stop())

  # Navigate to CCLE Quick tab
  app$click("nav-quick_ccle")
  app$wait_for_idle()

  app$set_inputs(`quick_ccle-gene` = "TP53")
  app$set_inputs(`quick_ccle-analysis_type` = "distribution")
  app$set_inputs(`quick_ccle-data_type` = "mRNA")

  app$click("quick_ccle-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("quick_ccle-plot_output"))
})

# Test CCLE Quick Analysis - Drug Sensitivity
test_that("CCLE Quick Analysis - Drug Sensitivity", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-ccle-drug",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-quick_ccle")
  app$wait_for_idle()

  app$set_inputs(`quick_ccle-gene` = "TP53")
  app$set_inputs(`quick_ccle-analysis_type` = "drug_sensitivity")
  app$set_inputs(`quick_ccle-data_type` = "mRNA")
  app$set_inputs(`quick_ccle-drug` = "cisplatin")

  app$click("quick_ccle-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("quick_ccle-plot_output"))
})

# Test Data Query Module
test_that("Data Query Module - TCGA", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-query-tcga",
    timeout = 60000
  )
  on.exit(app$stop())

  # Navigate to Data Query tab
  app$click("nav-data_query")
  app$wait_for_idle()

  app$set_inputs(`data_query-data_source` = "tcga")
  app$set_inputs(`data_query-data_type` = "mRNA")
  app$set_inputs(`data_query-gene_symbol` = "TP53")

  app$click("data_query-query_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("data_query-query_result"))
})

# Test Data Query Module - PCAWG
test_that("Data Query Module - PCAWG", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-query-pcawg",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-data_query")
  app$wait_for_idle()

  app$set_inputs(`data_query-data_source` = "pcawg")
  app$set_inputs(`data_query-data_type` = "mRNA")
  app$set_inputs(`data_query-gene_symbol` = "TP53")

  app$click("data_query-query_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("data_query-query_result"))
})

# Test Data Query Module - CCLE
test_that("Data Query Module - CCLE", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-query-ccle",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-data_query")
  app$wait_for_idle()

  app$set_inputs(`data_query-data_source` = "ccle")
  app$set_inputs(`data_query-data_type` = "mRNA")
  app$set_inputs(`data_query-gene_symbol` = "TP53")

  app$click("data_query-query_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("data_query-query_result"))
})

# Test Pan-Cancer Module
test_that("Pan-Cancer Module displays", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pancan",
    timeout = 30000
  )
  on.exit(app$stop())

  app$click("nav-pancan")
  app$wait_for_idle()

  # Verify module is visible
  expect_true(app$is_visible("pancan"))
})

# Test Immune Module
test_that("Immune Module displays", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-immune",
    timeout = 30000
  )
  on.exit(app$stop())

  app$click("nav-immune")
  app$wait_for_idle()

  expect_true(app$is_visible("immune"))
})

# Test Mutation Module
test_that("Mutation Module displays", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-mutation",
    timeout = 30000
  )
  on.exit(app$stop())

  app$click("nav-mutation")
  app$wait_for_idle()

  expect_true(app$is_visible("mutation"))
})

# Test Dimension Module
test_that("Dimension Module displays", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-dimension",
    timeout = 30000
  )
  on.exit(app$stop())

  app$click("nav-dimension")
  app$wait_for_idle()

  expect_true(app$is_visible("dimension"))
})

# Test PharmacoGenomics Module
test_that("PharmacoGenomics Module displays", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pharma",
    timeout = 30000
  )
  on.exit(app$stop())

  app$click("nav-pharma")
  app$wait_for_idle()

  expect_true(app$is_visible("pharma"))
})

# Test Batch Module
test_that("Batch Module displays", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-batch",
    timeout = 30000
  )
  on.exit(app$stop())

  app$click("nav-batch")
  app$wait_for_idle()

  expect_true(app$is_visible("batch"))
})

# Test Global Settings
test_that("Global settings update correctly", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-settings",
    timeout = 30000
  )
  on.exit(app$stop())

  # Open sidebar
  app$click("sidebar-toggle")
  app$wait_for_idle()

  # Change settings
  app$set_inputs(`global_source` = "pcawg")
  app$set_inputs(`global_cancer` = "BRCA")

  app$wait_for_idle()

  # Verify settings changed
  expect_equal(app$get_value("global_source"), "pcawg")
  expect_equal(app$get_value("global_cancer"), "BRCA")
})

# Test Download Functionality
test_that("Download buttons work", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-download",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-quick_tcga")
  app$wait_for_idle()

  app$set_inputs(`quick_tcga-gene` = "TP53")
  app$set_inputs(`quick_tcga-analysis_type` = "tumor_normal")
  app$click("quick_tcga-run_btn")
  app$wait_for_idle(timeout = 30000)

  # Test download buttons exist
  expect_true(app$is_visible("quick_tcga-download_plot"))
  expect_true(app$is_visible("quick_tcga-download_data"))
})

# Test Error Handling
test_that("Error handling works correctly", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-error",
    timeout = 30000
  )
  on.exit(app$stop())

  app$click("nav-quick_tcga")
  app$wait_for_idle()

  # Try to run with empty gene
  app$set_inputs(`quick_tcga-gene` = "")
  app$click("quick_tcga-run_btn")
  app$wait_for_idle()

  # Should show error notification
  expect_true(app$is_visible(".shiny-notification"))
})
