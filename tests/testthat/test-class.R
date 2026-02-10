test_that("CacheManager works correctly", {
  skip_on_cran()

  # Create temporary cache
  cm <- CacheManager$new(cache_dir = tempfile())

  # Test set and get
  test_data <- list(a = 1, b = 2, c = "test")
  cm$set("key1", test_data)
  expect_equal(cm$get("key1"), test_data)

  # Test has
  expect_true(cm$has("key1"))
  expect_false(cm$has("nonexistent"))

  # Test memory cache
  cm$set("key2", "value2")
  expect_equal(cm$get("key2"), "value2")

  # Test remove
  cm$remove("key1")
  expect_false(cm$has("key1"))

  # Test clear
  cm$clear()
  expect_false(cm$has("key2"))

  # Test info
  info <- cm$info()
  expect_type(info, "list")
  expect_true("cache_dir" %in% names(info))
})

test_that("AsyncCompute can submit and collect single tasks", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  engine <- AsyncCompute$new(n_workers = 2)
  engine$start()
  on.exit(engine$stop(), add = TRUE)

  expect_true(engine$info()$is_running)

  # Submit a simple task - just return a constant
  task_id <- engine$submit(42)

  expect_type(task_id, "character")

  # Wait for completion and collect
  result <- engine$collect(task_id, wait = TRUE)
  expect_equal(result, 42)
})

test_that("AsyncCompute batch processing works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  engine <- AsyncCompute$new(n_workers = 2)
  engine$start()
  on.exit(engine$stop(), add = TRUE)

  # Submit batch task
  items <- 1:5
  task_id <- engine$submit_batch(
    items = items,
    fn = function(x) x^2,
    .progress = FALSE
  )

  expect_type(task_id, "character")

  # Collect results
  results <- engine$collect_batch(task_id, wait = TRUE)
  expect_equal(results, list(1, 4, 9, 16, 25))
})

test_that("AsyncCompute progress tracking works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  engine <- AsyncCompute$new(n_workers = 2)
  engine$start()
  on.exit(engine$stop(), add = TRUE)

  # Submit batch task
  task_id <- engine$submit_batch(
    items = 1:10,
    fn = function(x) {
      Sys.sleep(0.05)
      x * 2
    },
    .progress = FALSE
  )

  # Check progress
  progress <- engine$get_progress(task_id)
  expect_type(progress, "list")
  expect_true("ready" %in% names(progress))
  expect_true("total" %in% names(progress))
  expect_equal(progress$total, 10)

  # Wait for completion
  results <- engine$collect_batch(task_id, wait = TRUE)
  expect_equal(length(results), 10)
})

test_that("AsyncCompute cleanup works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  engine <- AsyncCompute$new(n_workers = 2)
  engine$start()

  # Submit a task
  task_id <- engine$submit({
    Sys.sleep(0.1)
    "test"
  })

  # Stop should clean up
  engine$stop()
  expect_false(engine$info()$is_running)
})

test_that("XenaData initialization works", {
  skip_if_not_installed("UCSCXenaTools")

  # Should work with valid host
  xd <- XenaData$new(host = "toilHub")
  expect_equal(xd$get_info()$name, "Xena")

  # Should fail with invalid host
  expect_error(XenaData$new(host = "invalidHost"))
})

test_that("XenaData can query gene expression", {
  skip_on_cran()
  skip_if_not_installed("UCSCXenaTools")

  xd <- XenaData$new(host = "toilHub")

  # Query TP53 expression
  expr <- xd$get_gene_expression("TP53")

  expect_type(expr, "double")
  expect_gt(length(expr), 10000)
  expect_true(all(names(expr) != ""))
})

test_that("XenaData can query mutation status", {
  skip_on_cran()
  skip_if_not_installed("UCSCXenaTools")

  xd <- XenaData$new(host = "toilHub")

  # Query TP53 mutation
  mut <- xd$get_mutation_status("TP53")

  expect_type(mut, "double")
  expect_true(all(mut %in% c(0, 1, NA), na.rm = TRUE))
})

test_that("XenaData can query CNV", {
  skip_on_cran()
  skip_if_not_installed("UCSCXenaTools")

  xd <- XenaData$new(host = "toilHub")

  # Query MYC CNV
  cnv <- xd$get_cnv("MYC")

  expect_type(cnv, "double")
  expect_true(all(cnv %in% c(-2, -1, 0, 1, 2, NA), na.rm = TRUE))
})

test_that("XenaData parallel batch query works", {
  skip_on_cran()
  skip_if_not_installed("mirai")
  skip_if_not_installed("UCSCXenaTools")

  xd <- XenaData$new(host = "toilHub")

  # Query multiple genes in parallel
  genes <- c("TP53", "BRCA1", "EGFR")
  results <- xd$query_batch_parallel(genes, data_type = "mRNA", n_workers = 2, .progress = FALSE)

  expect_type(results, "list")
  expect_equal(names(results), genes)
  expect_true(all(vapply(results, is.numeric, logical(1))))
})
