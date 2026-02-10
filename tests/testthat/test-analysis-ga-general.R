# Test General Analysis Functions for Shiny Modules

# Test data for dimension reduction
create_test_data <- function() {
  set.seed(123)
  data <- matrix(rnorm(200), nrow = 20)
  colnames(data) <- paste0("S", 1:10)
  rownames(data) <- paste0("G", 1:20)
  data
}

# Test ga_dimension_distribution ----
test_that("ga_dimension_distribution works with PCA", {
  data <- create_test_data()

  result <- ga_dimension_distribution(data, method = "pca")

  expect_type(result, "list")
  expect_named(result, c("coordinates", "model", "variance_explained", "plot"))
  expect_s3_class(result$plot, "ggplot")
  expect_equal(ncol(result$coordinates), 3)  # PC1, PC2, sample
  expect_equal(nrow(result$coordinates), 10)  # 10 samples
})

test_that("ga_dimension_distribution validates input", {
  expect_error(
    ga_dimension_distribution("not a matrix"),
    "data must be a matrix or data frame"
  )
})

test_that("ga_dimension_distribution handles color_by", {
  data <- create_test_data()
  groups <- rep(c("A", "B"), each = 5)

  result <- ga_dimension_distribution(
    data,
    method = "pca",
    color_by = groups
  )

  expect_true("color_group" %in% colnames(result$coordinates))
  expect_equal(result$coordinates$color_group, groups)
})

# Test ga_group_comparison ----
test_that("ga_group_comparison validates input", {
  expect_error(
    ga_group_comparison(grp_df = data.frame(x = 1)),
    "ncol\\(grp_df\\) > 1 is not TRUE"
  )
})

# Test ga_survival_analysis ----
test_that("ga_survival_analysis validates input columns", {
  # Wrong number of columns
  surv_df <- data.frame(
    sample = c("S1", "S2"),
    time = c(365, 200)
  )

  expect_error(
    ga_survival_analysis(surv_df = surv_df),
    "surv_df must have 4 columns"
  )
})

test_that("ga_survival_analysis checks sample size", {
  surv_df <- data.frame(
    sample = c("S1", "S2", "S3"),
    value = c(1.2, 2.3, 3.1),
    time = c(365, 200, 500),
    status = c(1, 0, 1),
    stringsAsFactors = FALSE
  )

  expect_error(
    ga_survival_analysis(surv_df = surv_df),
    "Insufficient samples for survival analysis"
  )
})

# Test ga_scatter_correlation ----
test_that("ga_scatter_correlation validates inputs", {
  expect_error(
    ga_scatter_correlation(
      dataset1 = "test",
      id1 = c("A", "B"),  # Should be length 1
      dataset2 = "test",
      id2 = "C"
    ),
    "length\\(id1\\) == 1 is not TRUE"
  )
})

# Test ga_matrix_correlation ----
test_that("ga_matrix_correlation validates input", {
  expect_error(
    ga_matrix_correlation(
      dataset = "test",
      ids = "TP53"  # Need at least 2
    ),
    "length\\(ids\\) >= 2 is not TRUE"
  )
})

# Test error handling ----
test_that("functions handle missing data gracefully", {
  # Test with empty data frame
  empty_df <- data.frame(
    sample = character(),
    group = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )

  expect_error(
    ga_group_comparison(grp_df = empty_df),
    "Insufficient samples"
  )
})

# Test dimension reduction with different methods ----
test_that("ga_dimension_distribution handles t-SNE", {
  skip_if_not_installed("Rtsne")

  data <- create_test_data()

  result <- ga_dimension_distribution(data, method = "tsne")

  expect_type(result, "list")
  expect_s3_class(result$plot, "ggplot")
})

test_that("ga_dimension_distribution handles UMAP", {
  skip_if_not_installed("umap")

  data <- create_test_data()

  result <- ga_dimension_distribution(data, method = "umap")

  expect_type(result, "list")
  expect_s3_class(result$plot, "ggplot")
})

# Test data imputation ----
test_that("ga_dimension_distribution handles NA values", {
  data <- create_test_data()
  data[1, 1] <- NA
  data[2, 2] <- NA

  result <- ga_dimension_distribution(data, method = "pca")

  expect_type(result, "list")
  expect_equal(nrow(result$coordinates), 10)
})
