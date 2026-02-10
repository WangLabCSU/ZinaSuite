#' Other Pages Module Tests
#'
#' Tests for additional tools and features

library(shinytest2)

# Test Other Pages - Daily Gene
test_that("Other Pages - Daily Gene", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-other-daily",
    timeout = 60000
  )
  on.exit(app$stop())

  # Navigate to Other Pages tab
  app$click("nav-other_pages")
  app$wait_for_idle()

  # Verify daily gene is displayed
  expect_true(app$is_visible("other_pages-daily_gene_name"))
  expect_true(app$is_visible("other_pages-daily_gene_expr"))
  expect_true(app$is_visible("other_pages-daily_gene_survival"))
})

# Test Other Pages - Pan-Cancer Search
test_that("Other Pages - Pan-Cancer Search", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-other-search",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-other_pages")
  app$wait_for_idle()

  # Navigate to Pan-Cancer Search tab
  app$click("other_pages-Pan-Cancer Search")
  app$wait_for_idle()

  # Set search parameters
  app$set_inputs(`other_pages-search_query` = "TP53")
  app$set_inputs(`other_pages-search_type` = "gene")

  # Run search
  app$click("other_pages-run_search")
  app$wait_for_idle(timeout = 30000)

  # Verify results
  expect_true(app$is_visible("other_pages-search_results_table"))
})

# Test Other Pages - File Upload
test_that("Other Pages - File Upload", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-other-upload",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-other_pages")
  app$wait_for_idle()

  # Navigate to File Upload tab
  app$click("other_pages-File Upload")
  app$wait_for_idle()

  # Verify upload interface
  expect_true(app$is_visible("other_pages-upload_file"))
  expect_true(app$is_visible("other_pages-download_template"))
})

# Test Other Pages - Data Download
test_that("Other Pages - Data Download", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-other-download",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-other_pages")
  app$wait_for_idle()

  # Navigate to Data Download tab
  app$click("other_pages-Data Download")
  app$wait_for_idle()

  # Set download options
  app$set_inputs(`other_pages-download_dataset` = "tcga")
  app$set_inputs(`other_pages-download_data_type` = "expression")
  app$set_inputs(`other_pages-download_format` = "csv")

  # Verify download button
  expect_true(app$is_visible("other_pages-download_data_btn"))
})

# Test Other Pages - ID Query Help
test_that("Other Pages - ID Query Help", {
  app <- AppDriver$new(
    system.file("shinyapp", package = "ZinaSuite"),
    name = "zinasuite-other-id",
    timeout = 60000
  )
  on.exit(app$stop())

  app$click("nav-other_pages")
  app$wait_for_idle()

  # Navigate to ID Query Help tab
  app$click("other_pages-ID Query Help")
  app$wait_for_idle()

  # Set ID conversion parameters
  app$set_inputs(`other_pages-convert_id_input` = "TP53")
  app$set_inputs(`other_pages-convert_id_from` = "Gene Symbol")
  app$set_inputs(`other_pages-convert_id_to` = "Ensembl ID")

  # Run conversion
  app$click("other_pages-convert_id_btn")
  app$wait_for_idle()

  # Verify result
  expect_true(app$is_visible("other_pages-convert_result"))
})
