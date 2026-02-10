#' Cross-Omics Module Tests
#'
#' Tests for cross-omics analysis functionality

library(shinytest2)

# Test Cross-Omics - Gene Analysis
test_that("Cross-Omics - Gene Analysis", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-cross-gene",
    timeout = 60000
  )
  on.exit(app$stop())

  # Navigate to Cross-Omics tab
  app$click("nav-cross_omics")
  app$wait_for_idle()

  # Set gene analysis inputs
  app$set_inputs(`cross_omics-gene_id` = "TP53")
  app$set_inputs(`cross_omics-gene_omics_types` = c("mrna", "mutation", "cnv"))

  # Run analysis
  app$click("cross_omics-run_gene_analysis")
  app$wait_for_idle(timeout = 30000)

  # Verify results
  expect_true(app$is_visible("cross_omics-gene_heatmap"))
})

# Test Cross-Omics - Gene Correlation Tab
test_that("Cross-Omics - Gene Correlation", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-cross-gene-cor",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-cross_omics")
  app$wait_for_idle()

  app$set_inputs(`cross_omics-gene_id` = "TP53")
  app$click("cross_omics-run_gene_analysis")
  app$wait_for_idle(timeout = 30000)

  # Click on Correlation tab
  app$click("cross_omics-gene_cor_plot")
  app$wait_for_idle()

  expect_true(app$is_visible("cross_omics-gene_cor_plot"))
})

# Test Cross-Omics - Pathway Analysis
test_that("Cross-Omics - Pathway Analysis", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-cross-pathway",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-cross_omics")
  app$wait_for_idle()

  # Navigate to Pathway tab
  app$click("cross_omics-Pathway Cross-Omics")
  app$wait_for_idle()

  # Set pathway analysis inputs
  app$set_inputs(`cross_omics-pathway_id` = "HALLMARK_APOPTOSIS")
  app$set_inputs(`cross_omics-pathway_omics_types` = c("pathway", "mrna"))

  # Run analysis
  app$click("cross_omics-run_pathway_analysis")
  app$wait_for_idle(timeout = 30000)

  # Verify results
  expect_true(app$is_visible("cross_omics-pathway_heatmap"))
})

# Test Cross-Omics - Multi-Omics Integration
test_that("Cross-Omics - Multi-Omics Integration", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-cross-integration",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-cross_omics")
  app$wait_for_idle()

  # Navigate to Integration tab
  app$click("cross_omics-Multi-Omics Integration")
  app$wait_for_idle()

  # Set integration inputs
  app$set_inputs(`cross_omics-integ_gene` = "TP53")
  app$set_inputs(`cross_omics-integ_plot_type` = "circos")

  # Run integration
  app$click("cross_omics-run_integration")
  app$wait_for_idle(timeout = 30000)

  # Verify results
  expect_true(app$is_visible("cross_omics-integration_plot"))
})

# Test Cross-Omics - Data Table View
test_that("Cross-Omics - Data Table View", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-cross-table",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-cross_omics")
  app$wait_for_idle()

  app$set_inputs(`cross_omics-gene_id` = "TP53")
  app$click("cross_omics-run_gene_analysis")
  app$wait_for_idle(timeout = 30000)

  # Click on Data Table tab
  app$click("cross_omics-gene_data_table")
  app$wait_for_idle()

  expect_true(app$is_visible("cross_omics-gene_data_table"))
})
