#' Visualize Dimensionality Reduction
#'
#' @description
#' Create scatter plot showing dimensionality reduction results (PCA, t-SNE, UMAP).
#'
#' @param data Matrix or data frame with samples as rows
#' @param method Reduction method: "pca", "tsne", "umap"
#' @param n_components Number of components to compute
#' @param color_by Vector for coloring points (default: NULL)
#' @param shape_by Vector for shaping points (default: NULL)
#' @param title Plot title
#' @param ... Additional parameters passed to reduction method
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Query multiple genes
#' genes <- c("TP53", "BRCA1", "EGFR", "MYC")
#' expr_list <- query_molecules(genes, n_workers = 4)
#' 
#' # Create expression matrix
#' expr_matrix <- do.call(rbind, lapply(expr_list, function(x) {
#'   as.numeric(x)
#' }))
#' 
#' # PCA plot
#' p <- vis_identifier_dim_dist(expr_matrix, method = "pca")
#' print(p)
#' }
vis_identifier_dim_dist <- function(data,
                                    method = c("pca", "tsne", "umap"),
                                    n_components = 2,
                                    color_by = NULL,
                                    shape_by = NULL,
                                    title = NULL,
                                    ...) {
  method <- match.arg(method)
  
  # Check for required packages
  if (method == "tsne" && !requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Package 'Rtsne' is required for t-SNE. Install with: install.packages('Rtsne')")
  }
  if (method == "umap" && !requireNamespace("umap", quietly = TRUE)) {
    stop("Package 'umap' is required for UMAP. Install with: install.packages('umap')")
  }
  
  # Remove rows with NA
  complete_rows <- complete.cases(data)
  data <- data[complete_rows, ]
  
  if (nrow(data) < 3) {
    stop("Insufficient samples after removing NA values")
  }
  
  # Perform dimensionality reduction
  coords <- switch(method,
    "pca" = {
      pca_result <- stats::prcomp(data, scale. = TRUE, ...)
      pca_result$x[, 1:min(n_components, ncol(pca_result$x))]
    },
    "tsne" = {
      check_analysis_deps("dimension")
      tsne_result <- Rtsne::Rtsne(data, dims = n_components, ...)
      tsne_result$Y
    },
    "umap" = {
      check_analysis_deps("dimension")
      umap_result <- umap::umap(data, n_components = n_components, ...)
      umap_result$layout
    }
  )
  
  # Prepare data frame
  df <- as.data.frame(coords)
  colnames(df) <- paste0("Dim", 1:ncol(df))
  
  # Add color and shape variables
  if (!is.null(color_by)) {
    df$Color <- color_by[complete_rows]
  }
  if (!is.null(shape_by)) {
    df$Shape <- shape_by[complete_rows]
  }
  
  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2))
  
  if (!is.null(color_by) && !is.null(shape_by)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Color, shape = .data$Shape), alpha = 0.7, size = 3)
  } else if (!is.null(color_by)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Color), alpha = 0.7, size = 3)
  } else if (!is.null(shape_by)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(shape = .data$Shape), alpha = 0.7, size = 3)
  } else {
    p <- p + ggplot2::geom_point(alpha = 0.7, size = 3, color = "steelblue")
  }
  
  # Add variance explained for PCA
  if (method == "pca") {
    pca_result <- prcomp(data, scale. = TRUE)
    var_exp <- summary(pca_result)$importance[2, 1:2] * 100
    x_label <- sprintf("PC1 (%.1f%%)", var_exp[1])
    y_label <- sprintf("PC2 (%.1f%%)", var_exp[2])
  } else {
    x_label <- "Dimension 1"
    y_label <- "Dimension 2"
  }
  
  p <- p +
    ggplot2::labs(
      title = title %||% paste(toupper(method), "Plot"),
      x = x_label,
      y = y_label
    ) +
    theme_zinasuite()
  
  p
}


#' Visualize Group Comparison
#'
#' @description
#' Create boxplot or violin plot comparing groups.
#'
#' @param data Numeric vector of values
#' @param group Grouping variable
#' @param plot_type Type of plot: "boxplot", "violin", "bar"
#' @param test_method Statistical test method: "t.test", "wilcox", "anova"
#' @param title Plot title
#' @param palette Color palette name
#' @param ... Additional parameters
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Query gene expression
#' tp53_expr <- query_gene_expression("TP53")
#' 
#' # Load sample information
#' sample_info <- load_data("tcga_gtex")
#' 
#' # Match samples
#' common <- intersect(names(tp53_expr), sample_info$Sample)
#' 
#' # Create comparison plot
#' p <- vis_identifier_grp_comparison(
#'   data = as.numeric(tp53_expr[common]),
#'   group = sample_info$Type[match(common, sample_info$Sample)],
#'   plot_type = "boxplot",
#'   test_method = "wilcox"
#' )
#' print(p)
#' }
vis_identifier_grp_comparison <- function(data,
                                          group,
                                          plot_type = c("boxplot", "violin", "bar"),
                                          test_method = c("t.test", "wilcox", "anova"),
                                          title = NULL,
                                          palette = "default",
                                          ...) {
  plot_type <- match.arg(plot_type)
  test_method <- match.arg(test_method)
  
  # Prepare data
  df <- data.frame(
    Value = data,
    Group = group,
    stringsAsFactors = FALSE
  )
  
  # Remove NA
  df <- df[complete.cases(df), ]
  
  # Perform statistical test
  test_result <- switch(test_method,
    "t.test" = {
      if (length(unique(df$Group)) != 2) {
        warning("t.test requires exactly 2 groups, using ANOVA instead")
        aov(Value ~ Group, data = df)
      } else {
        t.test(Value ~ Group, data = df)
      }
    },
    "wilcox" = {
      if (length(unique(df$Group)) != 2) {
        warning("Wilcoxon test requires exactly 2 groups, using Kruskal-Wallis instead")
        kruskal.test(Value ~ Group, data = df)
      } else {
        wilcox.test(Value ~ Group, data = df)
      }
    },
    "anova" = aov(Value ~ Group, data = df)
  )
  
  # Extract p-value
  p_value <- if (inherits(test_result, "htest")) {
    test_result$p.value
  } else {
    summary(test_result)[[1]][["Pr(>F)"]][1]
  }
  
  # Create subtitle with test results
  subtitle <- sprintf("%s: p = %.2e", test_method, p_value)
  
  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Group, y = .data$Value, fill = .data$Group))
  
  p <- switch(plot_type,
    "boxplot" = p + ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5),
    "violin" = p + ggplot2::geom_violin(alpha = 0.7) + 
                 ggplot2::geom_boxplot(width = 0.1, alpha = 0.5),
    "bar" = {
      # Calculate mean and SE
      summary_df <- aggregate(Value ~ Group, data = df, 
                             FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
      summary_df <- do.call(data.frame, summary_df)
      colnames(summary_df) <- c("Group", "Mean", "SE")
      
      ggplot2::ggplot(summary_df, ggplot2::aes(x = .data$Group, y = .data$Mean, fill = .data$Group)) +
        ggplot2::geom_bar(stat = "identity", alpha = 0.7) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$Mean - .data$SE, ymax = .data$Mean + .data$SE), width = 0.2)
    }
  )
  
  if (plot_type != "bar") {
    p <- p +
      ggplot2::scale_fill_manual(values = get_palette(length(unique(df$Group)), palette)) +
      ggplot2::labs(
        title = title %||% "Group Comparison",
        subtitle = subtitle,
        x = "Group",
        y = "Value"
      ) +
      theme_zinasuite()
  }
  
  p
}


#' Visualize Group Survival
#'
#' @description
#' Create Kaplan-Meier survival curve by group.
#'
#' @param time Survival time
#' @param status Event status (0 = censored, 1 = event)
#' @param group Grouping variable
#' @param title Plot title
#' @param conf.int Show confidence intervals
#' @param risk.table Show risk table
#' @param ... Additional parameters
#' @return ggsurvplot object or ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Load survival data
#' surv_data <- load_data("tcga_surv")
#' 
#' # Create survival plot by cancer type
#' p <- vis_identifier_grp_surv(
#'   time = surv_data$OS.time,
#'   status = surv_data$OS,
#'   group = surv_data$Type,
#'   conf.int = TRUE,
#'   risk.table = TRUE
#' )
#' print(p)
#' }
vis_identifier_grp_surv <- function(time,
                                    status,
                                    group,
                                    title = NULL,
                                    conf.int = TRUE,
                                    risk.table = TRUE,
                                    ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required")
  }
  
  # Create survival object
  surv_obj <- survival::Surv(time, status)
  
  # Fit Kaplan-Meier model
  km_fit <- survival::survfit(surv_obj ~ group)
  
  # Log-rank test
  logrank <- survival::survdiff(surv_obj ~ group)
  p_value <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  
  # Create title
  plot_title <- title %||% sprintf("Survival by Group (p = %.3f)", p_value)
  
  # Check if survminer is available
  if (requireNamespace("survminer", quietly = TRUE)) {
    survminer::ggsurvplot(
      km_fit,
      conf.int = conf.int,
      risk.table = risk.table,
      title = plot_title,
      xlab = "Time",
      ylab = "Survival Probability",
      pval = TRUE,
      pval.method = TRUE,
      ...
    )
  } else {
    # Basic ggplot2 version
    df <- data.frame(
      time = km_fit$time,
      surv = km_fit$surv,
      lower = km_fit$lower,
      upper = km_fit$upper,
      group = rep(names(km_fit$strata), km_fit$strata)
    )
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$surv, color = .data$group)) +
      ggplot2::geom_step(size = 1) +
      ggplot2::labs(
        title = plot_title,
        x = "Time",
        y = "Survival Probability"
      ) +
      theme_zinasuite()
    
    if (conf.int) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$group), 
                                    alpha = 0.2, color = NA)
    }
    
    p
  }
}
