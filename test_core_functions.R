#!/usr/bin/env Rscript
# Test core ZinaSuite functions with actual data

cat("=== Testing ZinaSuite Core Functions ===\n\n")

# Load package
suppressPackageStartupMessages(library(ZinaSuite))

# Test 1: Data loading
cat("Test 1: Loading sample data...\n")
tryCatch({
  data("sample_data", package = "ZinaSuite")
  cat("  ✓ Sample data loaded successfully\n")
  cat("  - Expression matrix dimensions:", dim(sample_data$expression), "\n")
  cat("  - Number of samples:", ncol(sample_data$expression), "\n")
  cat("  - Number of genes:", nrow(sample_data$expression), "\n")
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

# Test 2: Query molecule (mock)
cat("\nTest 2: Query molecule function...\n")
tryCatch({
  # This will use mock data since we don't have actual Xena connection
  result <- query_molecule("TP53", data_type = "mRNA", source = "tcga")
  if (!is.null(result)) {
    cat("  ✓ Query function works\n")
    cat("  - Returned", length(result), "values\n")
  } else {
    cat("  ! Query returned NULL (expected without data)\n")
  }
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

# Test 3: TCGA correlation analysis
cat("\nTest 3: TCGA correlation analysis...\n")
tryCatch({
  result <- run_tcga_cor_o2o("TP53", "BRCA1")
  cat("  ✓ TCGA correlation analysis completed\n")
  cat("  - Correlation:", round(result$correlation, 3), "\n")
  cat("  - P-value:", format(result$p_value, digits = 3), "\n")
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

# Test 4: TCGA group comparison
cat("\nTest 4: TCGA group comparison...\n")
tryCatch({
  result <- run_tcga_comp_o2o("TP53", "BRCA")
  cat("  ✓ TCGA group comparison completed\n")
  if (!is.null(result$plot)) {
    cat("  - Plot generated successfully\n")
  }
  if (!is.null(result$data)) {
    cat("  - Data frame with", nrow(result$data), "rows\n")
  }
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

# Test 5: Survival analysis
cat("\nTest 5: TCGA survival analysis...\n")
tryCatch({
  result <- run_tcga_sur_o2o("TP53", "BRCA")
  cat("  ✓ TCGA survival analysis completed\n")
  if (!is.null(result$plot)) {
    cat("  - KM plot generated\n")
  }
  if (!is.null(result$stats)) {
    cat("  - Log-rank p-value:", format(result$stats$p_value, digits = 3), "\n")
  }
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

# Test 6: PCAWG analysis
cat("\nTest 6: PCAWG correlation analysis...\n")
tryCatch({
  result <- run_pcawg_cor_o2o("TP53", "BRCA1")
  cat("  ✓ PCAWG correlation analysis completed\n")
  cat("  - Correlation:", round(result$correlation, 3), "\n")
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

# Test 7: CCLE analysis
cat("\nTest 7: CCLE correlation analysis...\n")
tryCatch({
  result <- run_ccle_cor_o2o("TP53", "BRCA1")
  cat("  ✓ CCLE correlation analysis completed\n")
  cat("  - Correlation:", round(result$correlation, 3), "\n")
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

# Test 8: General analysis functions
cat("\nTest 8: General analysis - scatter correlation...\n")
tryCatch({
  x <- rnorm(100)
  y <- rnorm(100)
  result <- ga_scatter_correlation(x, y, "Gene1", "Gene2")
  cat("  ✓ Scatter correlation completed\n")
  cat("  - Correlation:", round(result$correlation, 3), "\n")
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

cat("\nTest 9: General analysis - matrix correlation...\n")
tryCatch({
  # Create sample matrix
  mat <- matrix(rnorm(200), nrow = 10)
  rownames(mat) <- paste0("Gene", 1:10)
  colnames(mat) <- paste0("Sample", 1:20)
  result <- ga_matrix_correlation(mat)
  cat("  ✓ Matrix correlation completed\n")
  cat("  - Correlation matrix dimensions:", dim(result$correlation_matrix), "\n")
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

cat("\nTest 10: Dimension reduction...\n")
tryCatch({
  mat <- matrix(rnorm(200), nrow = 10)
  rownames(mat) <- paste0("Gene", 1:10)
  colnames(mat) <- paste0("Sample", 1:20)
  result <- ga_dimension_reduction(mat, method = "pca")
  cat("  ✓ Dimension reduction completed\n")
  cat("  - PCA dimensions:", dim(result$coordinates), "\n")
}, error = function(e) {
  cat("  ✗ Error:", conditionMessage(e), "\n")
})

cat("\n=== Core Functions Test Complete ===\n")
