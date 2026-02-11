#' ZinaSuite Utility Functions
#'
#' @description
#' Internal utility functions for the ZinaSuite package.
#'
#' @name zinasuite-utils
#' @keywords internal
NULL

#' Null-default operator
#'
#' @param lhs First value
#' @param rhs Default value if lhs is NULL
#' @return lhs if not NULL, otherwise rhs
#' @name null-coalescing
#' @keywords internal
`%||%` <- function(lhs, rhs) if (is.null(lhs)) rhs else lhs

#' Check if required packages are installed
#'
#' @param pkgs Character vector of package names
#' @param stop_on_missing Whether to stop if packages are missing
#' @return Logical indicating if all packages are available
#' @keywords internal
check_packages <- function(pkgs, stop_on_missing = TRUE) {
  missing <- !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)

  if (any(missing)) {
    msg <- paste("Required package(s) not installed:",
                 paste(pkgs[missing], collapse = ", "))
    if (stop_on_missing) {
      stop(msg)
    } else {
      warning(msg)
      return(FALSE)
    }
  }
  TRUE
}

#' Get data source instance
#'
#' @param source Data source name: "tcga", "pcawg", "ccle"
#' @return DataSource instance
#' @keywords internal
get_data_source <- function(source = c("tcga", "pcawg", "ccle")) {
  source <- match.arg(source)

  switch(source,
    tcga = XenaData$new(host = "toilHub"),
    pcawg = PCAWGData$new(),
    ccle = CCLEData$new()
  )
}

#' Get or create global async compute engine
#'
#' @param n_workers Number of workers
#' @return AsyncCompute instance
#' @keywords internal
get_async_engine <- function(n_workers = parallel::detectCores() - 1) {
  if (!exists(".zina_async_engine", envir = .ZinaSuiteEnv)) {
    assign(".zina_async_engine",
           AsyncCompute$new(n_workers = n_workers),
           envir = .ZinaSuiteEnv)
  }
  get(".zina_async_engine", envir = .ZinaSuiteEnv)
}

#' Clean up global async compute engine
#'
#' @keywords internal
cleanup_async_engine <- function() {
  if (exists(".zina_async_engine", envir = .ZinaSuiteEnv)) {
    engine <- get(".zina_async_engine", envir = .ZinaSuiteEnv)
    engine$stop()
    rm(".zina_async_engine", envir = .ZinaSuiteEnv)
  }
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{pipe}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Null-coalescing operator
#'
#' Returns the left-hand side if not NULL and has length > 0, otherwise returns the right-hand side.
#'
#' @param lhs Left-hand side value
#' @param rhs Right-hand side value
#' @return lhs if not NULL and has length > 0, otherwise rhs
#' @name null_coalesce
#' @rdname null-coalescing
#' @keywords internal
#' @export
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs) && length(lhs) > 0) lhs else rhs
}

#' Format bytes to human-readable string
#'
#' @param bytes Number of bytes
#' @return Formatted string
#' @keywords internal
format_bytes <- function(bytes) {
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- floor(log(bytes, 1024)) + 1
  i <- min(i, length(units))
  sprintf("%.2f %s", bytes / (1024^(i-1)), units[i])
}

#' Retry a function call with exponential backoff
#'
#' @param expr Expression to evaluate
#' @param max_tries Maximum number of attempts
#' @param initial_delay Initial delay in seconds
#' @return Result of expression evaluation
#' @keywords internal
retry <- function(expr, max_tries = 3, initial_delay = 1) {
  for (i in seq_len(max_tries)) {
    result <- tryCatch(
      eval(expr),
      error = function(e) {
        if (i == max_tries) stop(e)
        delay <- initial_delay * (2^(i-1))
        message(sprintf("Attempt %d failed, retrying in %.1f seconds...", i, delay))
        Sys.sleep(delay)
        NULL
      }
    )
    if (!is.null(result)) return(result)
  }
}

#' Validate gene symbol
#'
#' @param gene Gene symbol to validate
#' @return TRUE if valid, FALSE otherwise
#' @keywords internal
is_valid_gene <- function(gene) {
  if (!is.character(gene) || length(gene) != 1) return(FALSE)
  if (nchar(gene) == 0) return(FALSE)
  # Basic check for valid gene symbol format
  grepl("^[A-Za-z0-9_-]+$", gene)
}

#' Validate sample IDs
#'
#' @param sample_ids Character vector of sample IDs
#' @return TRUE if valid, FALSE otherwise
#' @keywords internal
is_valid_samples <- function(sample_ids) {
  is.character(sample_ids) && length(sample_ids) > 0 && all(nchar(sample_ids) > 0)
}

#' Convert TCGA barcode to sample type
#'
#' @param barcode TCGA barcode
#' @return Sample type: "tumor" or "normal"
#' @keywords internal
tcga_barcode_to_type <- function(barcode) {
  # Extract sample type code (positions 14-15)
  type_code <- substr(barcode, 14, 15)

  # 01-09 are tumor samples, 10-19 are normal samples
  code_num <- as.numeric(type_code)
  ifelse(code_num <= 9, "tumor", "normal")
}

#' Extract cancer type from TCGA barcode
#'
#' @param barcode TCGA barcode
#' @return Cancer type abbreviation
#' @keywords internal
tcga_barcode_to_cancer <- function(barcode) {
  # Extract cancer type code (positions 1-2)
  toupper(substr(barcode, 1, 2))
}

#' Create a progress callback function
#'
#' @param total Total number of items
#' @param message Progress message
#' @return Callback function
#' @keywords internal
make_progress_callback <- function(total, message = "Processing") {
  count <- 0
  function(...) {
    count <<- count + 1
    if (count %% max(1, floor(total / 10)) == 0) {
      pct <- round(100 * count / total)
      message(sprintf("%s: %d%% (%d/%d)", message, pct, count, total))
    }
  }
}
