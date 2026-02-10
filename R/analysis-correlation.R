#' Correlation Analysis
#'
#' @description
#' Perform correlation analysis between molecular features.
#'
#' @param x First variable (numeric vector)
#' @param y Second variable (numeric vector)
#' @param method Correlation method: "pearson", "spearman", "kendall"
#' @param use_purity Whether to adjust for tumor purity (default: FALSE)
#' @return List with correlation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Query TP53 and BRCA1 expression
#' tp53_expr <- query_gene_expression("TP53")
#' brca1_expr <- query_gene_expression("BRCA1")
#'
#' # Calculate correlation
#' result <- analyze_correlation(tp53_expr, brca1_expr, method = "pearson")
#' print(result$estimate)  # Correlation coefficient
#' print(result$pvalue)    # P-value
#' }
analyze_correlation <- function(x, y,
                                method = c("pearson", "spearman", "kendall"),
                                use_purity = FALSE) {
  method <- match.arg(method)

  # Remove NA values
  complete_cases <- stats::complete.cases(x, y)
  x <- x[complete_cases]
  y <- y[complete_cases]

  if (length(x) < 3) {
    stop("Insufficient complete cases for correlation analysis (need at least 3)")
  }

  # Calculate correlation
  cor_result <- stats::cor.test(x, y, method = method)

  list(
    estimate = unname(cor_result$estimate),
    pvalue = cor_result$p.value,
    statistic = unname(cor_result$statistic),
    method = method,
    n = length(x),
    conf_int = if (!is.null(cor_result$conf.int)) {
      c(lower = cor_result$conf.int[1], upper = cor_result$conf.int[2])
    } else {
      NULL
    }
  )
}

#' Batch Correlation Analysis with mirai Parallel Processing
#'
#' @description
#' Perform correlation analysis between a target gene and multiple candidate genes
#' using parallel processing via mirai.
#'
#' @param target_gene Target gene symbol
#' @param candidate_genes Vector of candidate gene symbols to correlate with
#' @param data_type Type of molecular data: "mRNA", "protein", etc.
#' @param source Data source: "tcga", "pcawg", "ccle"
#' @param method Correlation method: "pearson", "spearman", "kendall"
#' @param n_workers Number of parallel workers (default: 4)
#' @param adjust_method Multiple testing correction method (default: "fdr")
#' @param .progress Show progress bar (default: TRUE)
#' @return Data frame with correlation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Parallel correlation analysis of TP53 with multiple genes
#' results <- analyze_correlation_batch(
#'   target_gene = "TP53",
#'   candidate_genes = c("BRCA1", "EGFR", "KRAS", "MYC", "CDKN2A"),
#'   data_type = "mRNA",
#'   source = "tcga",
#'   method = "pearson",
#'   n_workers = 4,
#'   .progress = TRUE
#' )
#'
#' # View top correlated genes
#' head(results[order(-abs(results$cor)), ])
#' }
analyze_correlation_batch <- function(target_gene,
                                      candidate_genes,
                                      data_type = "mRNA",
                                      source = "tcga",
                                      method = "pearson",
                                      n_workers = 4,
                                      adjust_method = "fdr",
                                      .progress = TRUE) {
  # Get target gene data
  target_data <- query_gene_expression(target_gene, source = source)

  if (is.null(target_data)) {
    stop("Failed to retrieve data for target gene: ", target_gene)
  }

  # Remove NA values from target data
  valid_samples <- !is.na(target_data)
  target_data_clean <- target_data[valid_samples]

  # Start mirai daemons
  mirai::daemons(n_workers)
  on.exit(mirai::daemons(0), add = TRUE)

  # Prepare argument list for each gene
  arg_list <- lapply(candidate_genes, function(gene) {
    list(
      gene = gene,
      target_data = target_data_clean,
      valid_samples = names(target_data_clean),
      data_type = data_type,
      source = source,
      method = method
    )
  })

  # Parallel correlation analysis
  results <- mirai::mirai_map(
    arg_list,
    function(args) {
      # Load required packages in worker
      if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
        stop("UCSCXenaTools not available")
      }

      # Get gene data
      xe <- UCSCXenaTools::XenaData
      host_url <- unique(xe$XenaHosts[xe$XenaHostNames == "toilHub"])
      dataset <- "TcgaTargetGtex_rsem_gene_tpm"

      gene_result <- UCSCXenaTools::fetch_dense_values(
        host = host_url[1],
        dataset = dataset,
        identifiers = args$gene,
        check = FALSE,
        use_probeMap = TRUE
      )

      if (is.null(gene_result) || nrow(gene_result) == 0) {
        return(list(
          gene = args$gene,
          cor = NA,
          pvalue = NA,
          n = 0,
          error = "Failed to retrieve gene data"
        ))
      }

      # Convert to named vector
      gene_data <- as.numeric(gene_result[1, ])
      names(gene_data) <- colnames(gene_result)

      # Match samples with target data
      common_samples <- intersect(args$valid_samples, names(gene_data))
      if (length(common_samples) < 3) {
        return(list(
          gene = args$gene,
          cor = NA,
          pvalue = NA,
          n = length(common_samples),
          error = "Insufficient common samples"
        ))
      }

      x <- args$target_data[common_samples]
      y <- gene_data[common_samples]

      # Calculate correlation
      cor_result <- stats::cor.test(x, y, method = args$method)

      list(
        gene = args$gene,
        cor = unname(cor_result$estimate),
        pvalue = cor_result$p.value,
        n = length(common_samples),
        error = NULL
      )
    },
    .progress = .progress
  )

  # Collect results
  results_list <- results[]

  # Convert to data frame
  results_df <- do.call(rbind, lapply(results_list, function(r) {
    data.frame(
      gene = r$gene,
      cor = r$cor,
      pvalue = r$pvalue,
      n = r$n,
      stringsAsFactors = FALSE
    )
  }))

  # Adjust p-values
  valid_pvalues <- !is.na(results_df$pvalue)
  if (any(valid_pvalues)) {
    results_df$padj <- NA
    results_df$padj[valid_pvalues] <- stats::p.adjust(
      results_df$pvalue[valid_pvalues],
      method = adjust_method
    )
  }

  # Sort by absolute correlation
  results_df <- results_df[order(-abs(results_df$cor)), ]
  rownames(results_df) <- NULL

  results_df
}

#' Correlation Matrix Analysis
#'
#' @description
#' Compute correlation matrix for multiple genes.
#'
#' @param genes Vector of gene symbols
#' @param data_type Type of molecular data
#' @param source Data source
#' @param method Correlation method
#' @return List containing correlation matrix and p-value matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Correlation matrix for DNA repair genes
#' genes <- c("BRCA1", "BRCA2", "TP53", "ATM", "CHEK2")
#' result <- analyze_correlation_matrix(genes, data_type = "mRNA")
#'
#' # Access correlation matrix
#' print(result$cor_matrix)
#'
#' # Visualize
#' corrplot::corrplot(result$cor_matrix)
#' }
analyze_correlation_matrix <- function(genes,
                                       data_type = "mRNA",
                                       source = "tcga",
                                       method = "pearson") {
  # Query all genes
  message("Querying data for ", length(genes), " genes...")
  gene_data_list <- query_molecules(
    identifiers = genes,
    data_type = data_type,
    source = source,
    n_workers = 4,
    .progress = TRUE
  )

  # Filter out failed queries
  valid_genes <- names(gene_data_list)[!vapply(gene_data_list, is.null, logical(1))]
  if (length(valid_genes) < 2) {
    stop("At least 2 genes with valid data are required")
  }

  if (length(valid_genes) < length(genes)) {
    warning("Failed to retrieve data for ", length(genes) - length(valid_genes), " genes")
  }

  # Find common samples
  sample_lists <- lapply(gene_data_list[valid_genes], names)
  common_samples <- Reduce(intersect, sample_lists)

  if (length(common_samples) < 10) {
    stop("Insufficient common samples across genes")
  }

  # Create data matrix
  data_matrix <- do.call(cbind, lapply(gene_data_list[valid_genes], function(x) {
    x[common_samples]
  }))
  colnames(data_matrix) <- valid_genes

  # Remove rows with NA
  complete_rows <- stats::complete.cases(data_matrix)
  data_matrix <- data_matrix[complete_rows, ]

  # Compute correlation matrix
  cor_matrix <- stats::cor(data_matrix, method = method, use = "pairwise.complete.obs")

  # Compute p-value matrix
  n <- nrow(data_matrix)
  p_matrix <- matrix(NA, nrow = ncol(data_matrix), ncol = ncol(data_matrix))
  colnames(p_matrix) <- rownames(p_matrix) <- valid_genes

  for (i in seq_along(valid_genes)) {
    for (j in seq_along(valid_genes)) {
      if (i != j) {
        test_result <- stats::cor.test(data_matrix[, i], data_matrix[, j], method = method)
        p_matrix[i, j] <- test_result$p.value
      } else {
        p_matrix[i, j] <- 0
      }
    }
  }

  list(
    cor_matrix = cor_matrix,
    p_matrix = p_matrix,
    n_samples = n,
    genes = valid_genes
  )
}

#' Partial Correlation Analysis
#'
#' @description
#' Perform partial correlation analysis adjusting for confounding variables.
#'
#' @param x First variable
#' @param y Second variable
#' @param z Confounding variable(s) to adjust for
#' @param method Correlation method
#' @return Partial correlation results
#' @export
#'
#' @examples
#' \donttest{
#' # Partial correlation adjusting for tumor purity
#' # Note: x, y, z must have the same length
#' set.seed(123)
#' n <- 100
#' x <- rnorm(n)
#' y <- rnorm(n) + 0.5 * x
#' z <- rnorm(n) + 0.3 * x + 0.3 * y
#'
#' result <- analyze_partial_correlation(
#'   x, y, z,
#'   method = "pearson"
#' )
#' print(result$estimate)
#' }
analyze_partial_correlation <- function(x, y, z,
                                        method = c("pearson", "spearman")) {
  method <- match.arg(method)

  # Ensure z is a matrix/data.frame
  if (is.vector(z)) {
    z <- as.matrix(z)
  }

  # Find complete cases
  complete_cases <- stats::complete.cases(x, y, z)
  x <- x[complete_cases]
  y <- y[complete_cases]
  z <- z[complete_cases, , drop = FALSE]

  if (method == "spearman") {
    x <- rank(x)
    y <- rank(y)
    z <- apply(z, 2, rank)
  }

  # Residualize x and y on z
  x_resid <- stats::residuals(stats::lm(x ~ z))
  y_resid <- stats::residuals(stats::lm(y ~ z))

  # Correlate residuals
  cor_result <- stats::cor.test(x_resid, y_resid, method = "pearson")

  list(
    estimate = unname(cor_result$estimate),
    pvalue = cor_result$p.value,
    statistic = unname(cor_result$statistic),
    n = length(x),
    method = method,
    adjusted_for = colnames(z) %||% "confounder"
  )
}
