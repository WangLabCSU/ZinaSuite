test_that("mirai basic functionality works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  # Start daemons
  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  # Create a simple mirai task
  m <- mirai::mirai({
    x <- 10
    y <- 20
    x + y
  })

  # Check it's unresolved initially
  expect_true(mirai::unresolved(m))

  # Collect result
  result <- mirai::call_mirai(m)$data
  expect_equal(result, 30)

  # Check it's resolved
  expect_false(mirai::unresolved(m))
})

test_that("mirai with arguments works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  # Create mirai with arguments
  m <- mirai::mirai({
    x * y
  }, x = 5, y = 6)

  result <- mirai::call_mirai(m)$data
  expect_equal(result, 30)
})

test_that("mirai_map parallel processing works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  # Parallel map
  results <- mirai::mirai_map(
    1:5,
    function(x) x^2,
    .progress = FALSE
  )

  # Collect results
  collected <- results[]
  expect_type(collected, "list")
  expect_length(collected, 5)
  expect_equal(collected, list(1, 4, 9, 16, 25))
})

test_that("mirai_map with complex function works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  # More complex computation
  results <- mirai::mirai_map(
    1:10,
    function(x) {
      # Simulate some work
      sum <- 0
      for (i in 1:x) {
        sum <- sum + i
      }
      sum
    },
    .progress = FALSE
  )

  collected <- results[]
  expect_type(collected, "list")
  expect_length(collected, 10)

  # Verify results (sum of 1:n = n*(n+1)/2)
  expected <- as.list(1:10 * 2:11 / 2)
  expect_equal(collected, expected)
})

test_that("mirai error handling works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  # Create a task that will error
  m <- mirai::mirai({
    stop("Test error")
  })

  # Collect - the error is captured in the result
  result <- mirai::call_mirai(m)$data

  # Result should be an errorValue object
  expect_true(inherits(result, "errorValue") || inherits(result, "try-error"))
})

test_that("mirai daemons configuration works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  # Start with specific number of daemons
  info <- mirai::daemons(4)
  on.exit(mirai::daemons(0), add = TRUE)

  # Should be active
  expect_true(length(info) > 0)

  # Reset daemons
  mirai::daemons(0)

  # Should be inactive
  info <- mirai::daemons(0)
})

test_that("mirai non-blocking check works", {
  skip_on_cran()
  skip_if_not_installed("mirai")

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  # Create a task that takes some time
  m <- mirai::mirai({
    Sys.sleep(0.5)
    42
  })

  # Initially should be unresolved
  expect_true(mirai::unresolved(m))

  # Wait a bit and check again
  Sys.sleep(0.1)
  # Still unresolved
  expect_true(mirai::unresolved(m))

  # Wait for completion
  result <- mirai::call_mirai(m)$data
  expect_equal(result, 42)
  expect_false(mirai::unresolved(m))
})
