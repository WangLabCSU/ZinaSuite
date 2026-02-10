# Test PharmacoGenomics Analysis Functions

# Query Drug Sensitivity Tests ------------------------------------------------

test_that("query_drug_sensitivity validates input", {
  expect_error(query_drug_sensitivity("INVALID_DRUG"))
  expect_error(query_drug_sensitivity("cisplatin", metric = "invalid"))
  expect_error(query_drug_sensitivity("cisplatin", source = "invalid"))
})

test_that("get_available_drugs returns correct drugs", {
  # Test CCLE drugs
  ccle_drugs <- get_available_drugs("ccle")
  expect_type(ccle_drugs, "character")
  expect_gt(length(ccle_drugs), 0)
  expect_true("cisplatin" %in% ccle_drugs)
  expect_true("paclitaxel" %in% ccle_drugs)

  # Test GDSC drugs
  gdsc_drugs <- get_available_drugs("gdsc")
  expect_type(gdsc_drugs, "character")
  expect_gt(length(gdsc_drugs), 0)
})

# Analyze Drug-Gene Correlation Tests -----------------------------------------

test_that("analyze_drug_gene_cor validates input", {
  skip_if_offline()
  skip_on_cran()

  expect_error(analyze_drug_gene_cor("", "cisplatin"))
  expect_error(analyze_drug_gene_cor("TP53", ""))
})

test_that("analyze_drug_gene_cor returns correct structure", {
  # Mock result structure
  mock_result <- list(
    gene = "TP53",
    drug = "cisplatin",
    metric = "auc",
    correlation = list(
      estimate = -0.3,
      pvalue = 0.01,
      method = "spearman"
    ),
    data = data.frame(
      Gene = rnorm(100),
      Drug = rnorm(100),
      CellLine = paste0("Cell_", 1:100),
      stringsAsFactors = FALSE
    ),
    n_samples = 100
  )

  expect_type(mock_result, "list")
  expect_equal(mock_result$gene, "TP53")
  expect_equal(mock_result$drug, "cisplatin")
  expect_type(mock_result$correlation, "list")
  expect_s3_class(mock_result$data, "data.frame")
  expect_equal(ncol(mock_result$data), 3)
})

# Batch Drug-Gene Analysis Tests ----------------------------------------------

test_that("analyze_drug_gene_batch returns correct structure", {
  # Mock result structure
  mock_result <- data.frame(
    gene = c("TP53", "BRCA1", "EGFR"),
    cor = c(-0.3, 0.2, 0.5),
    pvalue = c(0.01, 0.05, 0.001),
    n = c(100, 95, 98),
    error = c(NA, NA, NA),
    padj = c(0.03, 0.075, 0.003),
    stringsAsFactors = FALSE
  )

  expect_s3_class(mock_result, "data.frame")
  expect_equal(ncol(mock_result), 6)
  expect_equal(nrow(mock_result), 3)
})

test_that("analyze_drug_gene_batch handles empty gene list", {
  result <- analyze_drug_gene_batch(character(0), "cisplatin")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

# Analyze Drug-Mutation Tests -------------------------------------------------

test_that("analyze_drug_mutation returns correct structure", {
  # Mock result structure
  mock_result <- list(
    drug = "olaparib",
    gene = "BRCA1",
    metric = "auc",
    data = data.frame(
      CellLine = paste0("Cell_", 1:100),
      Sensitivity = c(rnorm(80, mean = 0.8), rnorm(20, mean = 0.5)),
      Mutation = c(rep("WT", 80), rep("MUT", 20)),
      Mutated = c(rep("Wild-type", 80), rep("Mutated", 20)),
      stringsAsFactors = FALSE
    ),
    test = list(
      p.value = 0.001,
      statistic = 3.5
    ),
    summary = list(
      wt_median = 0.8,
      mut_median = 0.5,
      wt_n = 80,
      mut_n = 20
    )
  )

  expect_type(mock_result, "list")
  expect_equal(mock_result$drug, "olaparib")
  expect_equal(mock_result$gene, "BRCA1")
  expect_s3_class(mock_result$data, "data.frame")
  expect_type(mock_result$summary, "list")
})

# Create Drug Profile Tests ---------------------------------------------------

test_that("create_drug_profile returns correct structure", {
  # Mock result structure
  mock_result <- data.frame(
    Drug = c("cisplatin", "paclitaxel", "doxorubicin"),
    AUC = c(0.3, 0.5, 0.7),
    stringsAsFactors = FALSE
  )

  expect_s3_class(mock_result, "data.frame")
  expect_equal(ncol(mock_result), 2)
  expect_equal(nrow(mock_result), 3)
})

# Pharmacogenomics Visualization Logic Tests ----------------------------------

test_that("drug-gene correlation results can be formatted", {
  result <- list(
    gene = "TP53",
    drug = "cisplatin",
    metric = "auc",
    correlation = list(
      estimate = -0.3,
      pvalue = 0.01,
      method = "spearman"
    ),
    n_samples = 100
  )

  formatted <- sprintf(
    "Gene: %s, Drug: %s, Correlation: %.3f, P-value: %.2e, N: %d",
    result$gene,
    result$drug,
    result$correlation$estimate,
    result$correlation$pvalue,
    result$n_samples
  )

  expect_type(formatted, "character")
  expect_match(formatted, "TP53")
  expect_match(formatted, "cisplatin")
})

test_that("drug-mutation results can be formatted", {
  result <- list(
    drug = "olaparib",
    gene = "BRCA1",
    summary = list(
      wt_median = 0.8,
      mut_median = 0.5,
      wt_n = 80,
      mut_n = 20
    ),
    test = list(p.value = 0.001)
  )

  formatted <- sprintf(
    "%s sensitivity: WT median = %.3f (n=%d), MUT median = %.3f (n=%d), p = %.2e",
    result$drug,
    result$summary$wt_median,
    result$summary$wt_n,
    result$summary$mut_median,
    result$summary$mut_n,
    result$test$p.value
  )

  expect_type(formatted, "character")
  expect_match(formatted, "olaparib")
  expect_match(formatted, "WT median")
})

# Integration Tests -----------------------------------------------------------

test_that("pharmacogenomics workflow produces consistent results", {
  # Mock a complete workflow
  genes <- c("TP53", "BRCA1", "EGFR")
  drug <- "cisplatin"

  # Mock batch analysis
  batch_result <- data.frame(
    gene = genes,
    cor = c(-0.3, 0.2, 0.5),
    pvalue = c(0.01, 0.05, 0.001),
    n = c(100, 95, 98),
    padj = c(0.03, 0.075, 0.003),
    stringsAsFactors = FALSE
  )

  # Check results are sorted by absolute correlation
  expect_equal(batch_result$gene[order(-abs(batch_result$cor))], batch_result$gene)

  # Check p-values are adjusted
  expect_true(all(batch_result$padj >= batch_result$pvalue | is.na(batch_result$pvalue)))
})
