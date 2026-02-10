#' AsyncCompute R6 Class
#'
#' @description
#' Asynchronous computation engine based on mirai. Provides parallel processing
#' capabilities for batch operations.
#'
#' @examples
#' \dontrun{
#' # Create async compute instance
#' ac <- AsyncCompute$new(n_workers = 4)
#'
#' # Start workers
#' ac$start()
#'
#' # Submit single task
#' task_id <- ac$submit({
#'   Sys.sleep(1)
#'   rnorm(100)
#' })
#'
#' # Check if ready
#' ac$is_ready(task_id)
#'
#' # Collect result
#' result <- ac$collect(task_id)
#'
#' # Submit batch tasks
#' batch_id <- ac$submit_batch(
#'   items = c("TP53", "BRCA1", "EGFR"),
#'   fn = function(gene) {
#'     list(gene = gene, length = nchar(gene))
#'   }
#' )
#'
#' # Collect batch results
#' results <- ac$collect_batch(batch_id)
#'
#' # Stop workers
#' ac$stop()
#' }
#'
#' @export
AsyncCompute <- R6::R6Class(
  "AsyncCompute",

  public = list(
    #' @description
    #' Initialize a new AsyncCompute instance
    #' @param n_workers Number of parallel workers (default: detectCores - 1)
    initialize = function(n_workers = parallel::detectCores() - 1) {
      private$n_workers <- max(1, n_workers)
      private$task_list <- list()
      private$is_running <- FALSE
    },

    #' @description
    #' Start the async compute workers
    #' @return Invisible self
    start = function() {
      if (!private$is_running) {
        mirai::daemons(n = private$n_workers)
        private$is_running <- TRUE
      }
      invisible(self)
    },

    #' @description
    #' Stop the async compute workers
    #' @return Invisible self
    stop = function() {
      if (private$is_running) {
        mirai::daemons(0)
        private$is_running <- FALSE
        private$task_list <- list()
      }
      invisible(self)
    },

    #' @description
    #' Submit a single async task
    #' @param expr Expression to evaluate (quoted)
    #' @return Task ID
    submit = function(expr) {
      private$ensure_running()
      task_id <- private$generate_id()

      # Use do.call to properly pass the expression
      m <- do.call(mirai::mirai, list(expr))

      private$task_list[[task_id]] <- list(
        type = "single",
        mirai = m,
        submitted = Sys.time()
      )

      task_id
    },

    #' @description
    #' Submit batch tasks using mirai_map
    #' @param items List/vector of items to process
    #' @param fn Function to apply to each item
    #' @param ... Additional arguments passed to fn
    #' @param .progress Show progress bar (default: TRUE)
    #' @return Task ID
    submit_batch = function(items, fn, ..., .progress = TRUE) {
      private$ensure_running()
      task_id <- private$generate_id()

      mp <- mirai::mirai_map(
        items,
        fn,
        ...,
        .progress = .progress
      )

      private$task_list[[task_id]] <- list(
        type = "batch",
        mirai = mp,
        submitted = Sys.time(),
        total = length(items)
      )

      task_id
    },

    #' @description
    #' Collect result from a single task
    #' @param task_id Task ID
    #' @param wait Whether to wait for completion (default: TRUE)
    #' @param timeout Timeout in seconds (default: Inf)
    #' @return Task result or NULL if not ready and wait=FALSE
    collect = function(task_id, wait = TRUE, timeout = Inf) {
      task <- private$task_list[[task_id]]
      if (is.null(task)) {
        stop("Task not found: ", task_id)
      }

      if (task$type != "single") {
        stop("Use collect_batch() for batch tasks")
      }

      m <- task$mirai

      if (wait) {
        if (is.finite(timeout)) {
          start_time <- Sys.time()
          while (mirai::unresolved(m) && as.numeric(Sys.time() - start_time) < timeout) {
            Sys.sleep(0.1)
          }
          if (mirai::unresolved(m)) {
            stop("Timeout waiting for task: ", task_id)
          }
        }
        result <- mirai::call_mirai(m)$data
        private$task_list[[task_id]] <- NULL  # Remove completed task
        return(result)
      } else {
        if (mirai::unresolved(m)) {
          return(NULL)
        } else {
          result <- m$data
          private$task_list[[task_id]] <- NULL
          return(result)
        }
      }
    },

    #' @description
    #' Collect results from a batch task
    #' @param task_id Task ID
    #' @param wait Whether to wait for completion (default: TRUE)
    #' @return List of results or NULL if not ready and wait=FALSE
    collect_batch = function(task_id, wait = TRUE) {
      task <- private$task_list[[task_id]]
      if (is.null(task)) {
        stop("Task not found: ", task_id)
      }

      if (task$type != "batch") {
        stop("Use collect() for single tasks")
      }

      mp <- task$mirai

      if (wait) {
        results <- mp[]
        private$task_list[[task_id]] <- NULL
        return(results)
      } else {
        if (all(!mirai::unresolved(mp))) {
          results <- mp[]
          private$task_list[[task_id]] <- NULL
          return(results)
        } else {
          return(NULL)
        }
      }
    },

    #' @description
    #' Check if a task is ready
    #' @param task_id Task ID
    #' @return TRUE if ready, FALSE otherwise
    is_ready = function(task_id) {
      task <- private$task_list[[task_id]]
      if (is.null(task)) return(NA)

      m <- task$mirai
      if (task$type == "batch") {
        all(!mirai::unresolved(m))
      } else {
        !mirai::unresolved(m)
      }
    },

    #' @description
    #' Get progress of a batch task
    #' @param task_id Task ID
    #' @return List with progress information
    get_progress = function(task_id) {
      task <- private$task_list[[task_id]]
      if (is.null(task)) {
        return(NULL)
      }

      m <- task$mirai
      if (task$type == "batch") {
        ready <- sum(!mirai::unresolved(m))
        total <- length(m)
        list(
          ready = ready,
          total = total,
          percent = ready / total,
          status = if (ready == total) "completed" else "running"
        )
      } else {
        ready <- as.integer(!mirai::unresolved(m))
        list(
          ready = ready,
          total = 1,
          percent = as.numeric(ready),
          status = if (ready) "completed" else "running"
        )
      }
    },

    #' @description
    #' Get information about the async compute instance
    #' @return List with information
    info = function() {
      list(
        n_workers = private$n_workers,
        is_running = private$is_running,
        active_tasks = length(private$task_list)
      )
    }
  ),

  private = list(
    # Finalizer - stop workers when object is garbage collected
    finalize = function() {
      if (private$is_running) {
        mirai::daemons(0)
      }
    },

    n_workers = NULL,
    is_running = FALSE,
    task_list = NULL,
    task_counter = 0,

    ensure_running = function() {
      if (!private$is_running) {
        self$start()
      }
    },

    generate_id = function() {
      private$task_counter <- private$task_counter + 1
      paste0("task_", private$task_counter, "_", format(Sys.time(), "%H%M%S"), "_", sample(1000, 1))
    }
  )
)
