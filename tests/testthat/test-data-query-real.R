# Test real data queries from UCSC Xena
# These tests require internet connection

test_that("get_pancan_gene_value works with real data", {
  skip_if_offline()
  skip_if_not_installed("UCSCXenaTools")

  # Query TP53 expression
  expr <- get_pancan_gene_value("TP53")

  expect_type(expr, "double")
  expect_gt(length(expr), 10000)  # Should have many samples
  expect_true(all(names(expr) != ""))  # All samples should have names
  expect_false(all(is.na(expr)))  # Not all values should be NA
})

test_that("get_pancan_mutation_status works with real data", {
  skip_if_offline()
  skip_if_not_installed("UCSCXenaTools")

  # Query TP53 mutation status
  mut <- get_pancan_mutation_status("TP53")

  expect_type(mut, "double")
  expect_gt(length(mut), 1000)
  expect_true(all(mut %in% c(0, 1, NA)))  # Should be binary
})

test_that("get_pancan_cn_value works with real data", {
  skip_if_offline()
  skip_if_not_installed("UCSCXenaTools")

  # Query MYC CNV
  cnv <- get_pancan_cn_value("MYC")

  expect_type(cnv, "double")
  expect_gt(length(cnv), 1000)
  expect_true(all(cnv %in% c(-2, -1, 0, 1, 2, NA)))  # GISTIC2 values
})

test_that("query_pancan_value unified interface works", {
  skip_if_offline()
  skip_if_not_installed("UCSCXenaTools")

  # Test mRNA
  expr <- query_pancan_value("TP53", data_type = "mRNA")
  expect_type(expr, "double")
  expect_gt(length(expr), 10000)

  # Test mutation
  mut <- query_pancan_value("TP53", data_type = "mutation")
  expect_type(mut, "double")

  # Test CNV
  cnv <- query_pancan_value("MYC", data_type = "cnv")
  expect_type(cnv, "double")
})
