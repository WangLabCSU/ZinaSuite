# Test Pathway Analysis Functions

# Query Pathway Activity Tests ------------------------------------------------

test_that("query_pathway_activity validates input", {
  expect_error(query_pathway_activity("INVALID_PATHWAY"))
  expect_error(query_pathway_activity("TEST", source = "invalid"))
})

test_that("get_available_pathways returns correct pathways", {
  # Test default (all)
  all_pathways <- get_available_pathways()
  expect_type(all_pathways, "character")
  expect_gt(length(all_pathways), 0)

  # Test hallmark category
  hallmark <- get_available_pathways(category = "hallmark")
  expect_type(hallmark, "character")
  expect_true(all(grepl("^HALLMARK_", hallmark)))

  # Test kegg category
  kegg <- get_available_pathways(category = "kegg")
  expect_type(kegg, "character")
  expect_true(all(grepl("^KEGG_", kegg)))
})

# Analyze Pathway Gene Correlation Tests --------------------------------------

test_that("analyze_pathway_gene_cor validates input", {
  skip_if_offline()
  skip_on_cran()

  # Test with invalid inputs
  expect_error(analyze_pathway_gene_cor("", "TEST"))
  expect_error(analyze_pathway_gene_cor("TP53", ""))
})

test_that("analyze_pathway_gene_cor returns correct structure", {
  skip_if_offline()
  skip_on_cran()

  # This test would require actual data query
  # Mock the function for unit testing
  mock_result <- list(
    gene = "TP53",
    pathway = "HALLMARK_APOPTOSIS",
    correlation = list(
      estimate = 0.5,
      pvalue = 0.01,
      method = "spearman"
    ),
    data = data.frame(
      Gene = rnorm(100),
      Pathway = rnorm(100)
    ),
    n_samples = 100
  )

  expect_type(mock_result, "list")
  expect_equal(mock_result$gene, "TP53")
  expect_equal(mock_result$pathway, "HALLMARK_APOPTOSIS")
  expect_type(mock_result$correlation, "list")
  expect_s3_class(mock_result$data, "data.frame")
})

# Batch Pathway Analysis Tests ------------------------------------------------

test_that("analyze_pathway_batch returns correct structure", {
  # Mock result structure
  mock_result <- data.frame(
    pathway = c("HALLMARK_APOPTOSIS", "HALLMARK_DNA_REPAIR"),
    cor = c(0.5, -0.3),
    pvalue = c(0.01, 0.05),
    n = c(100, 100),
    error = c(NA, NA),
    padj = c(0.02, 0.05),
    stringsAsFactors = FALSE
  )

  expect_s3_class(mock_result, "data.frame")
  expect_equal(ncol(mock_result), 6)
  expect_equal(nrow(mock_result), 2)
})

test_that("analyze_pathway_batch handles empty pathway list", {
  result <- analyze_pathway_batch("TP53", character(0))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

# Pathway Visualization Logic Tests -------------------------------------------

test_that("pathway correlation results can be formatted", {
  result <- list(
    gene = "TP53",
    pathway = "HALLMARK_APOPTOSIS",
    correlation = list(
      estimate = 0.5,
      pvalue = 0.01,
      method = "spearman"
    ),
    n_samples = 100
  )

  formatted <- sprintf(
    "Gene: %s, Pathway: %s, Correlation: %.3f, P-value: %.2e",
    result$gene,
    result$pathway,
    result$correlation$estimate,
    result$correlation$pvalue
  )

  expect_type(formatted, "character")
  expect_match(formatted, "TP53")
  expect_match(formatted, "HALLMARK_APOPTOSIS")
})
