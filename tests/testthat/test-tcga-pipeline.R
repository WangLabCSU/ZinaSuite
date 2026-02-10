#' TCGA Pipeline Module Tests
#'
#' Tests for TCGA deep analysis pipeline functionality

library(shinytest2)

# Test TCGA Pipeline - Correlation One-to-One
test_that("TCGA Pipeline - Correlation O2O", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pipeline-cor-o2o",
    timeout = 60000
  )
  on.exit(app$stop())

  # Navigate to TCGA Pipeline tab
  app$click("nav-tcga_pipeline")
  app$wait_for_idle()

  # Set inputs
  app$set_inputs(`tcga_pipeline-analysis_type` = "cor")
  app$set_inputs(`tcga_pipeline-mode` = "o2o")
  app$set_inputs(`tcga_pipeline-molecule_x` = "TP53")
  app$set_inputs(`tcga_pipeline-molecule_y` = "BRCA1")
  app$set_inputs(`tcga_pipeline-data_type_x` = "mRNA")
  app$set_inputs(`tcga_pipeline-data_type_y` = "mRNA")
  app$set_inputs(`tcga_pipeline-cor_method` = "spearman")

  # Run analysis
  app$click("tcga_pipeline-run_btn")
  app$wait_for_idle(timeout = 30000)

  # Verify results
  expect_true(app$is_visible("tcga_pipeline-plot_output"))
})

# Test TCGA Pipeline - Correlation One-to-Many
test_that("TCGA Pipeline - Correlation O2M", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pipeline-cor-o2m",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tcga_pipeline")
  app$wait_for_idle()

  app$set_inputs(`tcga_pipeline-analysis_type` = "cor")
  app$set_inputs(`tcga_pipeline-mode` = "o2m")
  app$set_inputs(`tcga_pipeline-molecule_x` = "TP53")
  app$set_inputs(`tcga_pipeline-molecule_y` = "BRCA1")
  app$set_inputs(`tcga_pipeline-data_type_x` = "mRNA")
  app$set_inputs(`tcga_pipeline-cor_method` = "pearson")

  app$click("tcga_pipeline-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("tcga_pipeline-plot_output"))
})

# Test TCGA Pipeline - Correlation Many-to-One
test_that("TCGA Pipeline - Correlation M2O", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pipeline-cor-m2o",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tcga_pipeline")
  app$wait_for_idle()

  app$set_inputs(`tcga_pipeline-analysis_type` = "cor")
  app$set_inputs(`tcga_pipeline-mode` = "m2o")
  app$set_inputs(`tcga_pipeline-molecules` = "TP53\nKRAS\nEGFR\nPTEN")
  app$set_inputs(`tcga_pipeline-data_type_x` = "mRNA")
  app$set_inputs(`tcga_pipeline-cor_method` = "spearman")

  app$click("tcga_pipeline-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("tcga_pipeline-plot_output"))
})

# Test TCGA Pipeline - Comparison One-to-One
test_that("TCGA Pipeline - Comparison O2O", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pipeline-comp-o2o",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tcga_pipeline")
  app$wait_for_idle()

  app$set_inputs(`tcga_pipeline-analysis_type` = "comp")
  app$set_inputs(`tcga_pipeline-mode` = "o2o")
  app$set_inputs(`tcga_pipeline-molecule_x` = "TP53")
  app$set_inputs(`tcga_pipeline-data_type_x` = "mRNA")

  app$click("tcga_pipeline-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("tcga_pipeline-plot_output"))
})

# Test TCGA Pipeline - Comparison One-to-Many
test_that("TCGA Pipeline - Comparison O2M", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pipeline-comp-o2m",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tcga_pipeline")
  app$wait_for_idle()

  app$set_inputs(`tcga_pipeline-analysis_type` = "comp")
  app$set_inputs(`tcga_pipeline-mode` = "o2m")
  app$set_inputs(`tcga_pipeline-molecule_x` = "TP53")
  app$set_inputs(`tcga_pipeline-data_type_x` = "mRNA")

  app$click("tcga_pipeline-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("tcga_pipeline-plot_output"))
})

# Test TCGA Pipeline - Survival One-to-One
test_that("TCGA Pipeline - Survival O2O", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pipeline-sur-o2o",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tcga_pipeline")
  app$wait_for_idle()

  app$set_inputs(`tcga_pipeline-analysis_type` = "sur")
  app$set_inputs(`tcga_pipeline-mode` = "o2o")
  app$set_inputs(`tcga_pipeline-molecule_x` = "TP53")
  app$set_inputs(`tcga_pipeline-data_type_x` = "mRNA")
  app$set_inputs(`tcga_pipeline-surv_measure` = "OS")
  app$set_inputs(`tcga_pipeline-cutoff_mode` = "median")

  app$click("tcga_pipeline-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("tcga_pipeline-plot_output"))
})

# Test TCGA Pipeline - Cross-Omics
test_that("TCGA Pipeline - Cross-Omics", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pipeline-cross",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tcga_pipeline")
  app$wait_for_idle()

  app$set_inputs(`tcga_pipeline-analysis_type` = "cross")
  app$set_inputs(`tcga_pipeline-molecule_x` = "TP53")
  app$set_inputs(`tcga_pipeline-molecule_y` = "BRCA1")
  app$set_inputs(`tcga_pipeline-data_type_x` = "mRNA")
  app$set_inputs(`tcga_pipeline-data_type_y` = "protein")

  app$click("tcga_pipeline-run_btn")
  app$wait_for_idle(timeout = 30000)

  expect_true(app$is_visible("tcga_pipeline-plot_output"))
})

# Test TCGA Pipeline - Data Table Tab
test_that("TCGA Pipeline - Data Table View", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pipeline-table",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tcga_pipeline")
  app$wait_for_idle()

  app$set_inputs(`tcga_pipeline-analysis_type` = "cor")
  app$set_inputs(`tcga_pipeline-mode` = "o2o")
  app$set_inputs(`tcga_pipeline-molecule_x` = "TP53")
  app$set_inputs(`tcga_pipeline-molecule_y` = "BRCA1")

  app$click("tcga_pipeline-run_btn")
  app$wait_for_idle(timeout = 30000)

  # Click on Data Table tab
  app$click("tcga_pipeline-data_table")
  app$wait_for_idle()

  expect_true(app$is_visible("tcga_pipeline-data_table"))
})

# Test TCGA Pipeline - Statistics Tab
test_that("TCGA Pipeline - Statistics View", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-pipeline-stats",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tcga_pipeline")
  app$wait_for_idle()

  app$set_inputs(`tcga_pipeline-analysis_type` = "cor")
  app$set_inputs(`tcga_pipeline-mode` = "o2o")
  app$set_inputs(`tcga_pipeline-molecule_x` = "TP53")
  app$set_inputs(`tcga_pipeline-molecule_y` = "BRCA1")

  app$click("tcga_pipeline-run_btn")
  app$wait_for_idle(timeout = 30000)

  # Click on Statistics tab
  app$click("tcga_pipeline-stats_output")
  app$wait_for_idle()

  expect_true(app$is_visible("tcga_pipeline-stats_output"))
})
