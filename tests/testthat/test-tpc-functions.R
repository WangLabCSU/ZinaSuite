#'为TPC功能模块添加测试

library(shinytest2)

# Test TPC Functions - Custom Signature
test_that("TPC Functions - Custom Signature", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tpc-signature",
    timeout = 60000
  )
  on.exit(app$stop())

  # Navigate to TPC Functions tab
  app$click("nav-tpc_functions")
  app$wait_for_idle()

  # Set signature inputs
  app$set_inputs(`tpc_functions-sig_name` = "Test_Signature")
  app$set_inputs(`tpc_functions-sig_data_type` = "mRNA")
  app$set_inputs(`tpc_functions-sig_genes` = "TP53\nKRAS\nEGFR")

  # Create signature
  app$click("tpc_functions-create_sig_btn")
  app$wait_for_idle()

  # Verify signature was created
  expect_true(app$is_visible("tpc_functions-sig_table"))
})

# Test TPC Functions - Custom Metadata Upload
test_that("TPC Functions - Custom Metadata", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tpc-metadata",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tpc_functions")
  app$wait_for_idle()

  # Navigate to Custom Metadata tab
  app$click("tpc_functions-Custom Metadata")
  app$wait_for_idle()

  # Add manual entry
  app$set_inputs(`tpc_functions-meta_sample` = "TCGA-TEST-01")
  app$set_inputs(`tpc_functions-meta_group` = "High")
  app$click("tpc_functions-add_meta_btn")
  app$wait_for_idle()

  # Verify metadata was added
  expect_true(app$is_visible("tpc_functions-meta_table"))
})

# Test TPC Functions - Sample Filtering
test_that("TPC Functions - Sample Filtering", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tpc-filter",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tpc_functions")
  app$wait_for_idle()

  # Navigate to Sample Filtering tab
  app$click("tpc_functions-Sample Filtering")
  app$wait_for_idle()

  # Set filter options
  app$set_inputs(`tpc_functions-filter_sample_type` = c("tumor", "normal"))
  app$set_inputs(`tpc_functions-filter_expression` = c(0, 100))

  # Apply filters
  app$click("tpc_functions-apply_filter_btn")
  app$wait_for_idle()

  # Verify filter results
  expect_true(app$is_visible("tpc_functions-filter_stats"))
})

# Test TPC Functions - Data Download
test_that("TPC Functions - Data Download", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tpc-download",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tpc_functions")
  app$wait_for_idle()

  # Navigate to Data Download tab
  app$click("tpc_functions-Data Download")
  app$wait_for_idle()

  # Set download options
  app$set_inputs(`tpc_functions-download_data_type` = "expression")
  app$set_inputs(`tpc_functions-download_format` = "csv")
  app$set_inputs(`tpc_functions-download_filtered_only` = TRUE)

  # Verify download button exists
  expect_true(app$is_visible("tpc_functions-download_data_btn"))
})

# Test TPC Functions - Molecular Origin
test_that("TPC Functions - Molecular Origin", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tpc-origin",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tpc_functions")
  app$wait_for_idle()

  # Navigate to Molecular Origin tab
  app$click("tpc_functions-Molecular Origin")
  app$wait_for_idle()

  # Set data source
  app$set_inputs(`tpc_functions-data_source` = "tcga")
  app$set_inputs(`tpc_functions-toil_data_type` = "mRNA")

  # Set origin
  app$click("tpc_functions-set_origin_btn")
  app$wait_for_idle()

  # Verify settings updated
  expect_true(app$is_visible("tpc_functions-origin_settings"))
})

# Test TPC Functions - Sample Grouping
test_that("TPC Functions - Sample Grouping", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-tpc-group",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-tpc_functions")
  app$wait_for_idle()

  # Navigate to Sample Grouping tab
  app$click("tpc_functions-Sample Grouping")
  app$wait_for_idle()

  # Set group options
  app$set_inputs(`tpc_functions-group_name` = "Test_Group")
  app$set_inputs(`tpc_functions-group_by` = "expression")
  app$set_inputs(`tpc_functions-group_gene` = "TP53")
  app$set_inputs(`tpc_functions-group_cutoff` = 50)

  # Create group
  app$click("tpc_functions-create_group_btn")
  app$wait_for_idle()

  # Verify group was created
  expect_true(app$is_visible("tpc_functions-group_table"))
})
