#' AnalysisEngine R6 Class
#'
#' @description
#' R6 class for performing various bioinformatics analyses including
#' correlation, survival, comparison, and dimensionality reduction.
#'
#' @export
#' @examples
#' \dontrun{
#' # Create analysis engine
#' engine <- AnalysisEngine$new()
#'
#' # Perform correlation analysis
#' tp53_expr <- query_gene_expression("TP53")
#' brca1_expr <- query_gene_expression("BRCA1")
#' result <- engine$correlation(tp53_expr, brca1_expr)
#'
#' # Perform survival analysis
#' surv_data <- load_data("tcga_surv")
#' result <- engine$survival(
#'   time = surv_data$OS.time,
#'   status = surv_data$OS,
#'   group = surv_data$Type
#' )
#' }
AnalysisEngine <- R6::R6Class(
  "AnalysisEngine",

  public = list(
    #' @description
    #' Initialize a new AnalysisEngine instance
    #' @param n_workers Number of parallel workers for async operations
    initialize = function(n_workers = 4) {
      private$n_workers <- n_workers
      private$async_compute <- AsyncCompute$new(n_workers = n_workers)
    },

    #' @description
    #' Perform correlation analysis between two variables
    #' @param x First variable (numeric vector)
    #' @param y Second variable (numeric vector)
    #' @param method Correlation method: "pearson", "spearman", "kendall"
    #' @return List with correlation results
    correlation = function(x, y, method = c("pearson", "spearman", "kendall")) {
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
    },

    #' @description
    #' Perform batch correlation analysis in parallel
    #' @param target_gene Target gene symbol
    #' @param candidate_genes Vector of candidate gene symbols
    #' @param data_type Type of molecular data
    #' @param source Data source
    #' @param method Correlation method
    #' @param adjust_method Multiple testing correction method
    #' @return Data frame with correlation results
    correlation_batch = function(target_gene,
                                  candidate_genes,
                                  data_type = "mRNA",
                                  source = "tcga",
                                  method = "pearson",
                                  adjust_method = "fdr") {
      analyze_correlation_batch(
        target_gene = target_gene,
        candidate_genes = candidate_genes,
        data_type = data_type,
        source = source,
        method = method,
        n_workers = private$n_workers,
        adjust_method = adjust_method,
        .progress = TRUE
      )
    },

    #' @description
    #' Perform survival analysis
    #' @param time Survival time (in days)
    #' @param status Event status (0 = censored, 1 = event/death)
    #' @param group Optional grouping variable
    #' @param analysis_type Type of analysis: "km", "cox", "both"
    #' @return List with survival analysis results
    survival = function(time, status, group = NULL,
                        analysis_type = c("km", "cox", "both")) {
      analyze_survival(time, status, group, analysis_type)
    },

    #' @description
    #' Perform survival analysis by gene expression
    #' @param gene Gene symbol
    #' @param cutoff_method Method to define groups
    #' @param cutoff_value Custom cutoff value
    #' @param data_type Type of molecular data
    #' @param source Data source
    #' @param analysis_type Type of survival analysis
    #' @return List with survival analysis results
    survival_by_expression = function(gene,
                                       cutoff_method = "median",
                                       cutoff_value = NULL,
                                       data_type = "mRNA",
                                       source = "tcga",
                                       analysis_type = "both") {
      analyze_survival_by_expression(
        gene = gene,
        cutoff_method = cutoff_method,
        cutoff_value = cutoff_value,
        data_type = data_type,
        source = source,
        analysis_type = analysis_type
      )
    },

    #' @description
    #' Perform univariate Cox regression for multiple genes
    #' @param genes Vector of gene symbols
    #' @param data_type Type of molecular data
    #' @param source Data source
    #' @param adjust_method Multiple testing correction method
    #' @return Data frame with Cox regression results
    unicox_batch = function(genes,
                            data_type = "mRNA",
                            source = "tcga",
                            adjust_method = "fdr") {
      analyze_unicox_batch(
        genes = genes,
        data_type = data_type,
        source = source,
        n_workers = private$n_workers,
        adjust_method = adjust_method,
        .progress = TRUE
      )
    },

    #' @description
    #' Perform group comparison analysis
    #' @param data Numeric vector of data to compare
    #' @param group Grouping variable
    #' @param method Test method: "t.test", "wilcox", "anova"
    #' @return List with comparison results
    comparison = function(data, group, method = c("t.test", "wilcox", "anova")) {
      method <- match.arg(method)

      # Remove NA values
      valid <- stats::complete.cases(data, group)
      data <- data[valid]
      group <- group[valid]

      if (length(unique(group)) < 2) {
        stop("Need at least 2 groups for comparison")
      }

      result <- switch(method,
        "t.test" = {
          if (length(unique(group)) != 2) {
            stop("t.test requires exactly 2 groups")
          }
          stats::t.test(data ~ group)
        },
        "wilcox" = {
          if (length(unique(group)) != 2) {
            stop("Wilcoxon test requires exactly 2 groups")
          }
          stats::wilcox.test(data ~ group)
        },
        "anova" = {
          fit <- stats::aov(data ~ group)
          summary(fit)
        }
      )

      list(
        test = method,
        result = result,
        n_groups = length(unique(group)),
        n_total = length(data)
      )
    },

    #' @description
    #' Perform dimensionality reduction
    #' @param data Matrix or data frame with samples as rows
    #' @param method Reduction method: "pca", "tsne", "umap"
    #' @param n_components Number of components/dimensions
    #' @param ... Additional parameters passed to reduction method
    #' @return List with reduced coordinates and model
    dimension_reduction = function(data,
                                    method = c("pca", "tsne", "umap"),
                                    n_components = 2,
                                    ...) {
      method <- match.arg(method)

      # Remove rows with NA
      data <- data[stats::complete.cases(data), ]

      # Check for required packages
      if (method == "tsne" && !requireNamespace("Rtsne", quietly = TRUE)) {
        stop("Package 'Rtsne' is required for t-SNE")
      }
      if (method == "umap" && !requireNamespace("umap", quietly = TRUE)) {
        stop("Package 'umap' is required for UMAP")
      }

      result <- switch(method,
        "pca" = {
          # PCA
          pca_result <- stats::prcomp(data, scale. = TRUE, ...)
          list(
            coordinates = pca_result$x[, 1:min(n_components, ncol(pca_result$x))],
            model = pca_result,
            variance_explained = summary(pca_result)$importance[2, 1:min(n_components, ncol(pca_result$x))]
          )
        },
        "tsne" = {
          # t-SNE
          tsne_result <- Rtsne::Rtsne(data, dims = n_components, ...)
          list(
            coordinates = tsne_result$Y,
            model = tsne_result
          )
        },
        "umap" = {
          # UMAP
          umap_result <- umap::umap(data, n_components = n_components, ...)
          list(
            coordinates = umap_result$layout,
            model = umap_result
          )
        }
      )

      result$method <- method
      result$n_components <- n_components
      result$n_samples <- nrow(data)

      result
    },

    #' @description
    #' Get analysis engine info
    #' @return List with engine information
    info = function() {
      list(
        n_workers = private$n_workers,
        async_compute_running = private$async_compute$info()$is_running
      )
    }
  ),

  private = list(
    n_workers = NULL,
    async_compute = NULL
  )
)
