#' Data Validation and Cleaning Tools
#'
#' @description
#' Comprehensive data validation and cleaning utilities for genomic data analysis.
#'
#' @name data-validation
NULL

#' Validate Gene Symbol
#'
#' @param gene Character vector of gene symbols
#' @param allow_multiple Whether to allow multiple genes
#' @return TRUE if valid, throws error otherwise
#' @export
#'
#' @examples
#' validate_gene_symbol("TP53")
#' validate_gene_symbol(c("TP53", "BRCA1", "EGFR"), allow_multiple = TRUE)
validate_gene_symbol <- function(gene, allow_multiple = FALSE) {
  if (is.null(gene) || length(gene) == 0) {
    stop("Gene symbol cannot be NULL or empty")
  }

  if (!allow_multiple && length(gene) > 1) {
    stop("Only single gene symbol is allowed. Set allow_multiple = TRUE for multiple genes.")
  }

  # Check for valid gene symbol format
  invalid_genes <- gene[!grepl("^[A-Za-z0-9_-]+$", gene)]
  if (length(invalid_genes) > 0) {
    stop("Invalid gene symbol(s): ", paste(invalid_genes, collapse = ", "))
  }

  # Check length (most gene symbols are < 20 characters)
  long_genes <- gene[nchar(gene) > 30]
  if (length(long_genes) > 0) {
    warning("Gene symbol(s) unusually long: ", paste(long_genes, collapse = ", "))
  }

  invisible(TRUE)
}

#' Validate Cancer Type
#'
#' @param cancer Character vector of cancer types
#' @param valid_cancers Optional vector of valid cancer types
#' @return TRUE if valid, throws error otherwise
#' @export
#'
#' @examples
#' validate_cancer_type("BRCA")
#' validate_cancer_type(c("BRCA", "LUAD", "LUSC"))
validate_cancer_type <- function(cancer, valid_cancers = NULL) {
  if (is.null(cancer) || length(cancer) == 0) {
    stop("Cancer type cannot be NULL or empty")
  }

  # Default TCGA cancer types
  if (is.null(valid_cancers)) {
    valid_cancers <- c(
      "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
      "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
      "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
      "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
    )
  }

  invalid_cancers <- cancer[!cancer %in% valid_cancers]
  if (length(invalid_cancers) > 0) {
    stop("Invalid cancer type(s): ", paste(invalid_cancers, collapse = ", "),
         "\nValid types: ", paste(valid_cancers, collapse = ", "))
  }

  invisible(TRUE)
}

#' Validate Sample IDs
#'
#' @param sample_ids Character vector of sample IDs
#' @param pattern Regex pattern for validation
#' @return TRUE if valid, throws error otherwise
#' @export
#'
#' @examples
#' validate_sample_ids(c("TCGA-01-1234-01", "TCGA-02-5678-01"))
validate_sample_ids <- function(sample_ids, pattern = NULL) {
  if (is.null(sample_ids) || length(sample_ids) == 0) {
    stop("Sample IDs cannot be NULL or empty")
  }

  # Default TCGA sample ID pattern
  if (is.null(pattern)) {
    pattern <- "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]$"
  }

  invalid_samples <- sample_ids[!grepl(pattern, sample_ids)]
  if (length(invalid_samples) > 0) {
    stop("Invalid sample ID(s): ", paste(head(invalid_samples, 5), collapse = ", "),
         ifelse(length(invalid_samples) > 5, "...", ""))
  }

  invisible(TRUE)
}

#' Clean Expression Data
#'
#' @param data Expression data matrix or data frame
#' @param remove_na Whether to remove rows with NA values
#' @param remove_inf Whether to remove infinite values
#' @param log_transform Whether to log2 transform (if data is not already)
#' @param min_value Minimum value threshold
#' @return Cleaned data
#' @export
#'
#' @examples
#' \dontrun{
#' data <- query_gene_expression("TP53")
#' clean_data <- clean_expression_data(data, remove_na = TRUE)
#' }
clean_expression_data <- function(data, remove_na = TRUE, remove_inf = TRUE,
                                   log_transform = FALSE, min_value = 0) {
  if (is.null(data)) {
    stop("Data cannot be NULL")
  }

  # Convert to numeric matrix if data frame
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  # Remove NA rows
  if (remove_na) {
    na_rows <- apply(data, 1, function(x) all(is.na(x)))
    data <- data[!na_rows, , drop = FALSE]

    # Remove rows with any NA
    na_rows <- apply(data, 1, function(x) any(is.na(x)))
    if (any(na_rows)) {
      warning("Removing ", sum(na_rows), " rows with NA values")
      data <- data[!na_rows, , drop = FALSE]
    }
  }

  # Remove infinite values
  if (remove_inf) {
    inf_rows <- apply(data, 1, function(x) any(is.infinite(x)))
    if (any(inf_rows)) {
      warning("Removing ", sum(inf_rows), " rows with infinite values")
      data <- data[!inf_rows, , drop = FALSE]
    }
  }

  # Log transform if requested
  if (log_transform) {
    # Check if data might already be log transformed
    if (max(data, na.rm = TRUE) > 100) {
      data <- log2(data + 1)
    }
  }

  # Apply minimum value threshold
  if (min_value > 0) {
    data[data < min_value] <- min_value
  }

  data
}

#' Detect Outliers
#'
#' @param x Numeric vector
#' @param method Outlier detection method: "iqr", "zscore", "mad"
#' @param threshold Threshold for outlier detection
#' @return Logical vector indicating outliers
#' @export
#'
#' @examples
#' data <- rnorm(100)
#' data[1:5] <- data[1:5] + 10  # Add outliers
#' outliers <- detect_outliers(data)
#' sum(outliers)
detect_outliers <- function(x, method = c("iqr", "zscore", "mad"), threshold = NULL) {
  method <- match.arg(method)

  if (is.null(threshold)) {
    threshold <- switch(method,
                        "iqr" = 1.5,
                        "zscore" = 3,
                        "mad" = 3)
  }

  switch(method,
         "iqr" = {
           q1 <- quantile(x, 0.25, na.rm = TRUE)
           q3 <- quantile(x, 0.75, na.rm = TRUE)
           iqr <- q3 - q1
           lower <- q1 - threshold * iqr
           upper <- q3 + threshold * iqr
           x < lower | x > upper
         },
         "zscore" = {
           z <- scale(x)
           abs(z) > threshold
         },
         "mad" = {
           med <- median(x, na.rm = TRUE)
           mad_val <- mad(x, na.rm = TRUE)
           if (mad_val == 0) mad_val <- 0.0001  # Avoid division by zero
           abs(x - med) / mad_val > threshold
         }
  )
}

#' Validate Correlation Input
#'
#' @param x First numeric vector
#' @param y Second numeric vector
#' @param min_observations Minimum number of valid observations
#' @return TRUE if valid, throws error otherwise
#' @export
#'
#' @examples
#' x <- rnorm(50)
#' y <- rnorm(50)
#' validate_correlation_input(x, y)
validate_correlation_input <- function(x, y, min_observations = 10) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length")
  }

  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors")
  }

  # Count valid observations
  valid_obs <- sum(!is.na(x) & !is.na(y))
  if (valid_obs < min_observations) {
    stop("Insufficient valid observations: ", valid_obs,
         " (minimum required: ", min_observations, ")")
  }

  # Check for constant values
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    stop("Cannot compute correlation with constant values")
  }

  invisible(TRUE)
}

#' Data Quality Report
#'
#' @param data Data to check
#' @return List with quality metrics
#' @export
#'
#' @examples
#' data <- data.frame(
#'   gene1 = c(1, 2, NA, 4, 5),
#'   gene2 = c(2, NA, 4, 5, 6)
#' )
#' report <- data_quality_report(data)
#' print(report)
data_quality_report <- function(data) {
  if (is.null(data)) {
    return(list(error = "Data is NULL"))
  }

  # Convert to data frame if matrix
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  report <- list(
    dimensions = dim(data),
    total_cells = prod(dim(data)),
    missing_values = sum(is.na(data)),
    missing_percentage = round(sum(is.na(data)) / prod(dim(data)) * 100, 2),
    infinite_values = sum(is.infinite(as.matrix(data))),
    zero_values = sum(as.matrix(data) == 0, na.rm = TRUE),
    negative_values = sum(as.matrix(data) < 0, na.rm = TRUE),
    duplicate_rows = sum(duplicated(data)),
    timestamp = Sys.time()
  )

  # Add column-wise statistics if data frame
  if (is.data.frame(data)) {
    report$column_stats <- lapply(data, function(col) {
      if (is.numeric(col)) {
        list(
          type = "numeric",
          mean = mean(col, na.rm = TRUE),
          sd = sd(col, na.rm = TRUE),
          min = min(col, na.rm = TRUE),
          max = max(col, na.rm = TRUE),
          na_count = sum(is.na(col))
        )
      } else {
        list(
          type = class(col)[1],
          unique_values = length(unique(col)),
          na_count = sum(is.na(col))
        )
      }
    })
  }

  class(report) <- c("DataQualityReport", "list")
  report
}

#' Print Data Quality Report
#'
#' @param x DataQualityReport object
#' @param ... Additional arguments
#' @export
print.DataQualityReport <- function(x, ...) {
  cat("=== Data Quality Report ===\n\n")
  cat("Dimensions: ", paste(x$dimensions, collapse = " x "), "\n")
  cat("Total cells: ", x$total_cells, "\n")
  cat("Missing values: ", x$missing_values, " (", x$missing_percentage, "%)\n", sep = "")
  cat("Infinite values: ", x$infinite_values, "\n", sep = "")
  cat("Zero values: ", x$zero_values, "\n", sep = "")
  cat("Negative values: ", x$negative_values, "\n", sep = "")
  cat("Duplicate rows: ", x$duplicate_rows, "\n", sep = "")
  cat("\nGenerated: ", format(x$timestamp), "\n", sep = "")
}

#' Batch Validate Data
#'
#' @param data_list Named list of datasets
#' @return Data frame with validation results
#' @export
#'
#' @examples
#' data_list <- list(
#'   dataset1 = data.frame(a = 1:10, b = rnorm(10)),
#'   dataset2 = data.frame(c = letters[1:5], d = rnorm(5))
#' )
#' results <- batch_validate_data(data_list)
batch_validate_data <- function(data_list) {
  if (!is.list(data_list) || is.null(names(data_list))) {
    stop("data_list must be a named list")
  }

  results <- lapply(names(data_list), function(name) {
    data <- data_list[[name]]

    result <- tryCatch({
      report <- data_quality_report(data)
      list(
        name = name,
        status = "OK",
        rows = report$dimensions[1],
        cols = report$dimensions[2],
        missing_pct = report$missing_percentage,
        message = "Valid"
      )
    }, error = function(e) {
      list(
        name = name,
        status = "ERROR",
        rows = NA,
        cols = NA,
        missing_pct = NA,
        message = conditionMessage(e)
      )
    })

    result
  })

  do.call(rbind, lapply(results, as.data.frame))
}
