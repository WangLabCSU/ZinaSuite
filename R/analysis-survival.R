#' Survival Analysis
#'
#' @description
#' Perform survival analysis using Kaplan-Meier and Cox regression.
#'
#' @param time Survival time (in days)
#' @param status Event status (0 = censored, 1 = event/death)
#' @param group Optional grouping variable for stratified analysis
#' @param analysis_type Type of analysis: "km" (Kaplan-Meier), "cox" (Cox regression), "both"
#' @return List with survival analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' # Load survival data
#' surv_data <- load_data("tcga_surv")
#' clinical <- load_data("tcga_clinical")
#'
#' # Get TP53 expression and create groups
#' tp53_expr <- query_gene_expression("TP53")
#' common_samples <- intersect(names(tp53_expr), surv_data$sample)
#'
#' # Create high/low expression groups
#' tp53_values <- tp53_expr[common_samples]
#' group <- ifelse(tp53_values > median(tp53_values), "High", "Low")
#'
#' # Perform survival analysis
#' result <- analyze_survival(
#'   time = surv_data$OS.time[match(common_samples, surv_data$sample)],
#'   status = surv_data$OS[match(common_samples, surv_data$sample)],
#'   group = group,
#'   analysis_type = "both"
#' )
#'
#' # Access results
#' print(result$km$pvalue)  # Log-rank test p-value
#' print(result$cox$hr)     # Hazard ratio
#' }
analyze_survival <- function(time,
                              status,
                              group = NULL,
                              analysis_type = c("km", "cox", "both")) {
  analysis_type <- match.arg(analysis_type)

  # Check for survival package
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for survival analysis")
  }

  # Create survival object
  surv_obj <- survival::Surv(time, status)

  results <- list()

  # Kaplan-Meier analysis
  if (analysis_type %in% c("km", "both")) {
    if (is.null(group)) {
      # Overall survival
      fit <- survival::survfit(surv_obj ~ 1)
      km_result <- list(
        fit = fit,
        pvalue = NULL,
        comparison = "Overall survival"
      )
    } else {
      # Stratified analysis
      fit <- survival::survfit(surv_obj ~ group)

      # Log-rank test
      logrank <- survival::survdiff(surv_obj ~ group)
      pvalue <- 1 - stats::pchisq(logrank$chisq, df = length(logrank$n) - 1)

      km_result <- list(
        fit = fit,
        pvalue = pvalue,
        comparison = logrank,
        group_levels = unique(group)
      )
    }
    results$km <- km_result
  }

  # Cox regression analysis
  if (analysis_type %in% c("cox", "both")) {
    if (is.null(group)) {
      # Null model (just intercept)
      cox_fit <- survival::coxph(surv_obj ~ 1)
      cox_result <- list(
        fit = cox_fit,
        hr = NULL,
        ci = NULL,
        pvalue = NULL,
        summary = summary(cox_fit)
      )
    } else {
      # Cox model with group
      cox_fit <- survival::coxph(surv_obj ~ group)
      cox_summary <- summary(cox_fit)

      cox_result <- list(
        fit = cox_fit,
        hr = cox_summary$conf.int[, "exp(coef)"],
        ci = cox_summary$conf.int[, c("lower .95", "upper .95")],
        pvalue = cox_summary$coefficients[, "Pr(>|z|)"],
        summary = cox_summary
      )
    }
    results$cox <- cox_result
  }

  results$n <- length(time)
  results$events <- sum(status, na.rm = TRUE)

  results
}

#' Survival Analysis by Gene Expression
#'
#' @description
#' Perform survival analysis comparing groups based on gene expression levels.
#'
#' @param gene Gene symbol
#' @param cutoff_method Method to define groups: "median", "tertile", "quartile", "custom"
#' @param cutoff_value Custom cutoff value (if cutoff_method = "custom")
#' @param data_type Type of molecular data
#' @param source Data source
#' @param analysis_type Type of survival analysis
#' @return List with survival analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze survival by TP53 expression (median split)
#' result <- analyze_survival_by_expression(
#'   gene = "TP53",
#'   cutoff_method = "median",
#'   analysis_type = "both"
#' )
#'
#' # View hazard ratio
#' print(result$survival$cox$hr)
#' }
analyze_survival_by_expression <- function(gene,
                                           cutoff_method = c("median", "tertile", "quartile", "custom"),
                                           cutoff_value = NULL,
                                           data_type = "mRNA",
                                           source = "tcga",
                                           analysis_type = "both") {
  cutoff_method <- match.arg(cutoff_method)

  # Query gene expression
  gene_expr <- query_gene_expression(gene, source = source)

  if (is.null(gene_expr)) {
    stop("Failed to retrieve expression data for gene: ", gene)
  }

  # Load survival data
  surv_data <- load_data("tcga_surv")

  # Find common samples using barcode matching
  # (expression and survival data may have slightly different ID formats)
  match_result <- match_samples(names(gene_expr), surv_data$Sample, "tcga", "tcga", match_by = "barcode")

  if (match_result$n_matched < 10) {
    stop("Insufficient samples with both expression and survival data")
  }

  # Get expression values for common samples
  expr_values <- gene_expr[match_result$idx1]
  names(expr_values) <- match_result$common_ids

  # Create groups based on cutoff method
  group <- switch(cutoff_method,
    "median" = {
      cutoff <- stats::median(expr_values, na.rm = TRUE)
      ifelse(expr_values > cutoff, "High", "Low")
    },
    "tertile" = {
      tertiles <- stats::quantile(expr_values, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
      cut(expr_values, breaks = tertiles, labels = c("Low", "Medium", "High"), include.lowest = TRUE)
    },
    "quartile" = {
      quartiles <- stats::quantile(expr_values, probs = c(0, 0.25, 0.75, 1), na.rm = TRUE)
      ifelse(expr_values <= quartiles[2], "Low",
        ifelse(expr_values >= quartiles[3], "High", "Medium")
      )
    },
    "custom" = {
      if (is.null(cutoff_value)) {
        stop("cutoff_value must be provided when cutoff_method = 'custom'")
      }
      ifelse(expr_values > cutoff_value, "High", "Low")
    }
  )

  # Get survival data for common samples
  time <- surv_data$OS.time[match_result$idx2]
  status <- surv_data$OS[match_result$idx2]

  # Remove NA values
  valid <- stats::complete.cases(time, status, group)
  time <- time[valid]
  status <- status[valid]
  group <- group[valid]

  # Perform survival analysis
  surv_result <- analyze_survival(time, status, group, analysis_type)

  list(
    gene = gene,
    cutoff_method = cutoff_method,
    cutoff_value = if (cutoff_method == "custom") cutoff_value else stats::median(expr_values, na.rm = TRUE),
    n_samples = match_result$n_matched,
    n_events = sum(status, na.rm = TRUE),
    group_sizes = table(group),
    survival = surv_result,
    expression = expr_values
  )
}

#' Univariate Cox Analysis for Multiple Genes
#'
#' @description
#' Perform univariate Cox regression analysis for multiple genes in parallel.
#'
#' @param genes Vector of gene symbols
#' @param data_type Type of molecular data
#' @param source Data source
#' @param n_workers Number of parallel workers
#' @param adjust_method Multiple testing correction method
#' @param .progress Show progress bar
#' @return Data frame with Cox regression results
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze multiple genes
#' genes <- c("TP53", "BRCA1", "EGFR", "MYC")
#' results <- analyze_unicox_batch(
#'   genes = genes,
#'   n_workers = 4,
#'   .progress = TRUE
#' )
#'
#' # View significant genes
#' subset(results, padj < 0.05)
#' }
analyze_unicox_batch <- function(genes,
                                  data_type = "mRNA",
                                  source = "tcga",
                                  n_workers = 4,
                                  adjust_method = "fdr",
                                  .progress = TRUE) {
  # Check for survival package
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required")
  }

  # Load survival data
  surv_data <- load_data("tcga_surv")

  # Start mirai daemons
  mirai::daemons(n_workers)
  on.exit(mirai::daemons(0), add = TRUE)

  # Prepare argument list
  arg_list <- lapply(genes, function(gene) {
    list(
      gene = gene,
      surv_data = surv_data,
      data_type = data_type,
      source = source
    )
  })

  # Parallel Cox analysis
  results <- mirai::mirai_map(
    arg_list,
    function(args) {
      # Load required packages
      if (!requireNamespace("UCSCXenaTools", quietly = TRUE) ||
          !requireNamespace("survival", quietly = TRUE)) {
        return(list(
          gene = args$gene,
          hr = NA,
          lower = NA,
          upper = NA,
          pvalue = NA,
          n = 0,
          error = "Required packages not available"
        ))
      }

      # Get gene expression
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
          hr = NA,
          lower = NA,
          upper = NA,
          pvalue = NA,
          n = 0,
          error = "Failed to retrieve gene data"
        ))
      }

      # Convert to named vector
      gene_data <- as.numeric(gene_result[1, ])
      names(gene_data) <- colnames(gene_result)

      # Match with survival data
      common_samples <- intersect(names(gene_data), args$surv_data$sample)
      if (length(common_samples) < 10) {
        return(list(
          gene = args$gene,
          hr = NA,
          lower = NA,
          upper = NA,
          pvalue = NA,
          n = length(common_samples),
          error = "Insufficient common samples"
        ))
      }

      expr <- gene_data[common_samples]
      time <- args$surv_data$OS.time[match(common_samples, args$surv_data$sample)]
      status <- args$surv_data$OS[match(common_samples, args$surv_data$sample)]

      # Remove NA
      valid <- stats::complete.cases(expr, time, status)
      if (sum(valid) < 10) {
        return(list(
          gene = args$gene,
          hr = NA,
          lower = NA,
          upper = NA,
          pvalue = NA,
          n = sum(valid),
          error = "Insufficient valid samples after removing NA"
        ))
      }

      expr <- expr[valid]
      time <- time[valid]
      status <- status[valid]

      # Cox regression
      surv_obj <- survival::Surv(time, status)
      cox_fit <- survival::coxph(surv_obj ~ expr)
      cox_summary <- summary(cox_fit)

      list(
        gene = args$gene,
        hr = cox_summary$conf.int[1, "exp(coef)"],
        lower = cox_summary$conf.int[1, "lower .95"],
        upper = cox_summary$conf.int[1, "upper .95"],
        pvalue = cox_summary$coefficients[1, "Pr(>|z|)"],
        n = sum(valid),
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
      hr = r$hr,
      lower = r$lower,
      upper = r$upper,
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

  # Sort by p-value
  results_df <- results_df[order(results_df$pvalue), ]
  rownames(results_df) <- NULL

  results_df
}
