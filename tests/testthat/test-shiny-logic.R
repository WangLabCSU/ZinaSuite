# Test Shiny Module Business Logic
# These tests verify UI/server-independent business logic functions

# Data Query Logic Tests --------------------------------------------------

test_that("build_query_params creates correct structure", {
  params <- build_query_params("TP53", "mRNA", "tcga")

  expect_type(params, "list")
  expect_equal(params$gene, "TP53")
  expect_equal(params$data_type, "mRNA")
  expect_equal(params$source, "tcga")
  expect_true(!is.null(params$timestamp))
})

test_that("validate_query_input validates correctly", {
  # Valid input
  result <- validate_query_input("TP53", "mRNA", "tcga")
  expect_true(result$valid)
  expect_equal(result$message, "Valid input")

  # Empty gene
  result <- validate_query_input("", "mRNA", "tcga")
  expect_false(result$valid)
  expect_match(result$message, "Gene symbol is required")

  # Invalid characters in gene
  result <- validate_query_input("TP53!", "mRNA", "tcga")
  expect_false(result$valid)
  expect_match(result$message, "invalid characters")

  # Invalid data type
  result <- validate_query_input("TP53", "invalid_type", "tcga")
  expect_false(result$valid)
  expect_match(result$message, "Invalid data type")

  # Invalid source
  result <- validate_query_input("TP53", "mRNA", "invalid_source")
  expect_false(result$valid)
  expect_match(result$message, "Invalid source")
})

test_that("format_query_results formats correctly", {
  # Successful query
  query_result <- list(
    success = TRUE,
    data = c(1, 2, 3, 4, 5),
    gene = "TP53",
    data_type = "mRNA",
    source = "tcga",
    sample_count = 5
  )

  formatted <- format_query_results(query_result)
  expect_type(formatted, "character")
  expect_match(formatted, "TP53")
  expect_match(formatted, "mRNA")
  expect_match(formatted, "TCGA")

  # Failed query
  query_result <- list(
    success = FALSE,
    error = "No data found"
  )

  formatted <- format_query_results(query_result)
  expect_match(formatted, "Error:")
  expect_match(formatted, "No data found")
})

test_that("create_query_dataframe creates correct structure", {
  # Successful query
  query_result <- list(
    success = TRUE,
    data = c(sample1 = 1.5, sample2 = 2.5, sample3 = 3.5),
    gene = "TP53",
    data_type = "mRNA",
    source = "tcga"
  )

  df <- create_query_dataframe(query_result)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
  expect_equal(colnames(df), c("Sample", "Value", "Gene", "DataType", "Source"))
  expect_equal(df$Gene[1], "TP53")

  # Failed query
  query_result <- list(
    success = FALSE,
    error = "Error message"
  )

  df <- create_query_dataframe(query_result)
  expect_s3_class(df, "data.frame")
  expect_equal(colnames(df), "Error")
})

# Analysis Logic Tests ----------------------------------------------------

test_that("build_analysis_params creates correct structure", {
  params <- build_analysis_params("correlation", "TP53", "BRCA1", "BRCA")

  expect_type(params, "list")
  expect_equal(params$analysis_type, "correlation")
  expect_equal(params$gene1, "TP53")
  expect_equal(params$gene2, "BRCA1")
  expect_equal(params$cancer, "BRCA")
  expect_true(!is.null(params$timestamp))
})

test_that("validate_analysis_input validates correctly", {
  # Valid input
  result <- validate_analysis_input("correlation", "TP53", "BRCA1")
  expect_true(result$valid)

  # Invalid analysis type
  result <- validate_analysis_input("invalid_type", "TP53", "BRCA1")
  expect_false(result$valid)
  expect_match(result$message, "'correlation' or 'survival'")

  # Empty gene1
  result <- validate_analysis_input("correlation", "", "BRCA1")
  expect_false(result$valid)
  expect_match(result$message, "Gene 1 is required")

  # Empty gene2
  result <- validate_analysis_input("correlation", "TP53", "")
  expect_false(result$valid)
  expect_match(result$message, "Gene 2 is required")
})

test_that("format_analysis_results formats correctly", {
  # Successful analysis
  analysis_result <- list(
    success = TRUE,
    result = list(
      estimate = 0.75,
      p.value = 0.001,
      method = "pearson"
    ),
    sample_count = 100,
    gene1 = list(mean = 5.5, sd = 1.2),
    gene2 = list(mean = 6.2, sd = 1.5)
  )

  formatted <- format_analysis_results(analysis_result)
  expect_type(formatted, "character")
  expect_match(formatted, "Sample Count: 100")
  expect_match(formatted, "0.75")
  expect_match(formatted, "pearson")

  # Failed analysis
  analysis_result <- list(
    success = FALSE,
    error = "Insufficient samples"
  )

  formatted <- format_analysis_results(analysis_result)
  expect_match(formatted, "Error:")
})

# Visualization Logic Tests -----------------------------------------------

test_that("build_plot_params creates correct structure", {
  params <- build_plot_params("histogram", "My Title", "X Label", "Y Label", "blue")

  expect_type(params, "list")
  expect_equal(params$plot_type, "histogram")
  expect_equal(params$title, "My Title")
  expect_equal(params$x_label, "X Label")
  expect_equal(params$y_label, "Y Label")
  expect_equal(params$color, "blue")
})

test_that("validate_plot_input validates correctly", {
  # Valid data
  result <- validate_plot_input(c(1, 2, 3, 4, 5), "histogram")
  expect_true(result$valid)

  # Empty data
  result <- validate_plot_input(NULL, "histogram")
  expect_false(result$valid)
  expect_match(result$message, "Data is empty")

  # All NA
  result <- validate_plot_input(c(NA, NA, NA), "histogram")
  expect_false(result$valid)
  expect_match(result$message, "only NA")

  # Invalid plot type
  result <- validate_plot_input(c(1, 2, 3), "invalid_type")
  expect_false(result$valid)
  expect_match(result$message, "Invalid plot type")
})

test_that("generate_distribution_plot returns ggplot object", {
  skip_if_not_installed("ggplot2")

  data <- rnorm(100)
  plot <- generate_distribution_plot(data, "Test Plot", "histogram")

  expect_s3_class(plot, "ggplot")
})

test_that("generate_correlation_plot returns ggplot object", {
  skip_if_not_installed("ggplot2")

  x <- rnorm(50)
  y <- rnorm(50)
  plot <- generate_correlation_plot(x, y, "Correlation Plot")

  expect_s3_class(plot, "ggplot")
})

# Batch Processing Logic Tests --------------------------------------------

test_that("create_batch_job creates correct structure", {
  params <- list(gene1 = "TP53", gene2 = "BRCA1")
  job <- create_batch_job("correlation", params, priority = 2)

  expect_type(job, "list")
  expect_match(job$id, "^job_")
  expect_equal(job$type, "correlation")
  expect_equal(job$params, params)
  expect_equal(job$priority, 2)
  expect_equal(job$status, "pending")
  expect_true(!is.null(job$created))
})

test_that("format_batch_results creates correct data frame", {
  jobs <- list(
    create_batch_job("query", list(gene = "TP53")),
    create_batch_job("correlation", list(gene1 = "TP53", gene2 = "BRCA1"))
  )

  # Simulate completion
  jobs[[1]]$status <- "completed"
  jobs[[1]]$started <- Sys.time()
  jobs[[1]]$completed <- Sys.time() + 1

  jobs[[2]]$status <- "failed"
  jobs[[2]]$started <- Sys.time()
  jobs[[2]]$completed <- Sys.time() + 2

  result <- format_batch_results(jobs)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(colnames(result), c("JobID", "Type", "Status", "Created", "Duration"))
  expect_equal(result$Status, c("completed", "failed"))
})

# Status and Notification Logic Tests -------------------------------------

test_that("create_status creates correct structure", {
  status <- create_status("running", "Processing data", list(step = 1))

  expect_type(status, "list")
  expect_equal(status$status, "running")
  expect_equal(status$message, "Processing data")
  expect_equal(status$details$step, 1)
  expect_true(!is.null(status$timestamp))
})

test_that("format_status formats correctly", {
  # With status
  status <- create_status("running", "Processing")
  formatted <- format_status(status)
  expect_type(formatted, "character")
  expect_match(formatted, "Processing")

  # NULL status
  formatted <- format_status(NULL)
  expect_equal(formatted, "Ready")
})

test_that("create_notification creates correct structure", {
  notif <- create_notification("success", "Done", "Analysis complete")

  expect_type(notif, "list")
  expect_equal(notif$type, "success")
  expect_equal(notif$title, "Done")
  expect_equal(notif$message, "Analysis complete")
  expect_true(!is.null(notif$timestamp))
})

# Progress Tracking Logic Tests -------------------------------------------

test_that("create_progress_tracker creates correct structure", {
  tracker <- create_progress_tracker(5, "Analysis")

  expect_type(tracker, "list")
  expect_equal(tracker$total, 5)
  expect_equal(tracker$current, 0)
  expect_equal(tracker$description, "Analysis")
  expect_equal(length(tracker$steps), 5)
  expect_true(!is.null(tracker$started))
})

test_that("update_progress updates correctly", {
  tracker <- create_progress_tracker(3)
  tracker <- update_progress(tracker, 1, "Step 1")

  expect_equal(tracker$current, 1)
  expect_equal(tracker$steps[1], "Step 1")
})

test_that("calculate_progress_pct calculates correctly", {
  tracker <- create_progress_tracker(10)
  expect_equal(calculate_progress_pct(tracker), 0)

  tracker$current <- 5
  expect_equal(calculate_progress_pct(tracker), 50)

  tracker$current <- 10
  expect_equal(calculate_progress_pct(tracker), 100)
})

test_that("format_progress formats correctly", {
  tracker <- create_progress_tracker(5, "Analysis")
  tracker <- update_progress(tracker, 2, "Processing data")

  formatted <- format_progress(tracker)
  expect_type(formatted, "character")
  expect_match(formatted, "Analysis:")
  expect_match(formatted, "40%")
  expect_match(formatted, "Processing data")
})
