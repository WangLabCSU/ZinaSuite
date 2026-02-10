#' Visualize Gene Expression in Tumor vs Normal
#'
#' @description
#' Create boxplot or violin plot comparing gene expression between tumor and normal samples.
#'
#' @param gene Gene symbol
#' @param data_type Type of molecular data (default: "mRNA")
#' @param plot_type Type of plot: "boxplot" or "violin"
#' @param source Data source (default: "tcga")
#' @param palette Color palette name
#' @param title Plot title (auto-generated if NULL)
#' @param ... Additional parameters
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create tumor vs normal boxplot for TP53
#' p <- vis_toil_TvsN("TP53", plot_type = "boxplot")
#' print(p)
#'
#' # Create violin plot
#' p <- vis_toil_TvsN("TP53", plot_type = "violin")
#' print(p)
#' }
vis_toil_TvsN <- function(gene,
                          plot_type = c("boxplot", "violin"),
                          source = "tcga",
                          palette = "cancer",
                          title = NULL,
                          ...) {
  plot_type <- match.arg(plot_type)

  # Query gene expression
  gene_expr <- query_gene_expression(gene, source = source)

  # Load sample information
  sample_info <- load_data("tcga_gtex")

  # Match samples
  common_samples <- intersect(names(gene_expr), sample_info$Sample)

  if (length(common_samples) == 0) {
    stop("No matching samples found between expression data and sample info")
  }

  # Prepare data
  df <- data.frame(
    Sample = common_samples,
    Expression = as.numeric(gene_expr[common_samples]),
    stringsAsFactors = FALSE
  )

  # Add type information
  df$Type <- sample_info$Type[match(df$Sample, sample_info$Sample)]
  df$Tissue <- sample_info$Tissue[match(df$Sample, sample_info$Sample)]

  # Filter for tumor and normal only
  df <- df[df$Type %in% c("Tumor", "Normal"), ]

  if (nrow(df) == 0) {
    stop("No tumor or normal samples found")
  }

  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Tissue, y = .data$Expression, fill = .data$Type))

  if (plot_type == "boxplot") {
    p <- p + ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5)
  } else {
    p <- p + ggplot2::geom_violin(alpha = 0.7) +
      ggplot2::geom_boxplot(width = 0.1, alpha = 0.5)
  }

  p <- p +
    ggplot2::scale_fill_manual(values = get_palette(2, palette)) +
    ggplot2::labs(
      title = title %||% paste(gene, "Expression: Tumor vs Normal"),
      x = "Tissue",
      y = "Expression (log2 TPM)"
    ) +
    theme_zinasuite() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}


#' Visualize Univariate Cox Regression Results as Forest Plot
#'
#' @description
#' Create forest plot showing hazard ratios from univariate Cox regression.
#'
#' @param data Data frame with columns: gene, hr, lower, upper, pvalue
#' @param title Plot title
#' @param palette Color palette name
#' @param sort_by Sort results by: "hr" or "pvalue"
#' @param ... Additional parameters
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Run batch Cox analysis
#' results <- analyze_unicox_batch(c("TP53", "BRCA1", "EGFR"), n_workers = 2)
#'
#' # Create forest plot
#' p <- vis_unicox_tree(results)
#' print(p)
#' }
vis_unicox_tree <- function(data,
                            title = "Univariate Cox Regression",
                            palette = "default",
                            sort_by = c("hr", "pvalue"),
                            ...) {
  sort_by <- match.arg(sort_by)

  # Ensure required columns exist
  required_cols <- c("gene", "hr", "lower", "upper")
  if (!all(required_cols %in% colnames(data))) {
    stop("Data must contain columns: gene, hr, lower, upper")
  }

  # Sort data
  if (sort_by == "hr") {
    data <- data[order(data$hr), ]
  } else {
    data <- data[order(data$pvalue), ]
  }

  # Add significance indicator
  data$significant <- data$pvalue < 0.05

  # Create plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(gene, hr), y = hr))

  p <- p +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper, color = significant), width = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      x = "Gene",
      y = "Hazard Ratio (95% CI)",
      color = "Significant (p < 0.05)"
    ) +
    theme_zinasuite()

  p
}


#' Visualize Gene Expression Across Cancer Types (Anatomy Plot)
#'
#' @description
#' Create heatmap showing gene expression across different cancer types.
#'
#' @param gene Gene symbol
#' @param data_type Type of molecular data (default: "mRNA")
#' @param source Data source (default: "tcga")
#' @param title Plot title (auto-generated if NULL)
#' @param ... Additional parameters
#' @return ggplot object or pheatmap object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create anatomy heatmap for TP53
#' p <- vis_pancan_anatomy("TP53")
#' print(p)
#' }
vis_pancan_anatomy <- function(gene,
                               source = "tcga",
                               title = NULL,
                               ...) {
  # Query gene expression
  gene_expr <- query_gene_expression(gene, source = source)

  # Load sample information
  sample_info <- load_data("tcga_gtex")

  # Match samples
  common_samples <- intersect(names(gene_expr), sample_info$Sample)

  # Prepare data
  df <- data.frame(
    Sample = common_samples,
    Expression = as.numeric(gene_expr[common_samples]),
    Cancer = sample_info$Cancer[match(common_samples, sample_info$Sample)],
    Type = sample_info$Type[match(common_samples, sample_info$Sample)],
    stringsAsFactors = FALSE
  )

  # Calculate median expression by cancer type
  expr_summary <- stats::aggregate(
    Expression ~ Cancer + Type,
    data = df,
    FUN = median,
    na.rm = TRUE
  )

  # Reshape for heatmap
  expr_wide <- tidyr::pivot_wider(
    expr_summary,
    names_from = Type,
    values_from = Expression
  )

  # Create heatmap using ggplot2
  p <- ggplot2::ggplot(expr_summary, ggplot2::aes(x = Type, y = Cancer, fill = Expression)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = get_palette(100, "gradient")) +
    ggplot2::labs(
      title = title %||% paste(gene, "Expression Across Cancer Types"),
      x = "Sample Type",
      y = "Cancer Type"
    ) +
    theme_zinasuite()

  p
}


#' Visualize Gene-Immune Correlation
#'
#' @description
#' Create heatmap showing correlation between gene expression and immune cell infiltration.
#'
#' @param gene Gene symbol
#' @param immune_features Vector of immune feature names (default: common immune cell types)
#' @param data_type Type of molecular data (default: "mRNA")
#' @param source Data source (default: "tcga")
#' @param title Plot title (auto-generated if NULL)
#' @param ... Additional parameters
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create immune correlation heatmap for TP53
#' p <- vis_gene_immune_cor("TP53")
#' print(p)
#' }
vis_gene_immune_cor <- function(gene,
                                immune_features = NULL,
                                source = "tcga",
                                title = NULL,
                                ...) {
  # Query gene expression
  gene_expr <- query_gene_expression(gene, source = source)

  # Load immune data (TIL)
  immune_data <- load_data("tcga_TIL")

  # Match samples
  common_samples <- intersect(names(gene_expr), immune_data$Sample)

  if (length(common_samples) < 10) {
    stop("Insufficient samples with immune data")
  }

  # Get immune features if not specified
  if (is.null(immune_features)) {
    immune_features <- setdiff(colnames(immune_data), "Sample")
    immune_features <- immune_features[1:min(10, length(immune_features))]
  }

  # Calculate correlations
  cor_results <- data.frame(
    Feature = character(),
    Correlation = numeric(),
    Pvalue = numeric(),
    stringsAsFactors = FALSE
  )

  for (feature in immune_features) {
    if (feature %in% colnames(immune_data)) {
      feature_values <- immune_data[match(common_samples, immune_data$Sample), feature]
      gene_values <- as.numeric(gene_expr[common_samples])

      valid <- stats::complete.cases(feature_values, gene_values)
      if (sum(valid) > 10) {
        cor_test <- stats::cor.test(gene_values[valid], feature_values[valid])
        cor_results <- rbind(cor_results, data.frame(
          Feature = feature,
          Correlation = cor_test$estimate,
          Pvalue = cor_test$p.value
        ))
      }
    }
  }

  if (nrow(cor_results) == 0) {
    stop("No valid correlations calculated")
  }

  # Create heatmap
  p <- ggplot2::ggplot(cor_results, ggplot2::aes(x = "Gene", y = stats::reorder(Feature, Correlation), fill = Correlation)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggplot2::labs(
      title = title %||% paste(gene, "vs Immune Infiltration"),
      x = "",
      y = "Immune Cell Type"
    ) +
    theme_zinasuite()

  p
}


#' Visualize Gene Correlation Scatter Plot
#'
#' @description
#' Create scatter plot showing correlation between two genes.
#'
#' @param gene1 First gene symbol
#' @param gene2 Second gene symbol
#' @param data_type Type of molecular data (default: "mRNA")
#' @param source Data source (default: "tcga")
#' @param color_by Variable to color points by (default: NULL)
#' @param title Plot title (auto-generated if NULL)
#' @param ... Additional parameters
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create correlation scatter plot
#' p <- vis_gene_cor("TP53", "BRCA1")
#' print(p)
#' }
vis_gene_cor <- function(gene1,
                         gene2,
                         source = "tcga",
                         color_by = NULL,
                         title = NULL,
                         ...) {
  # Query gene expressions
  expr1 <- query_gene_expression(gene1, source = source)
  expr2 <- query_gene_expression(gene2, source = source)

  # Find common samples
  common_samples <- intersect(names(expr1), names(expr2))

  if (length(common_samples) < 10) {
    stop("Insufficient common samples between genes")
  }

  # Prepare data
  df <- data.frame(
    Gene1 = as.numeric(expr1[common_samples]),
    Gene2 = as.numeric(expr2[common_samples]),
    Sample = common_samples,
    stringsAsFactors = FALSE
  )

  # Add color variable if specified
  if (!is.null(color_by)) {
    if (color_by == "Type") {
      sample_info <- load_data("tcga_gtex")
      df$Color <- sample_info$Type[match(df$Sample, sample_info$Sample)]
    }
  }

  # Calculate correlation
  cor_test <- stats::cor.test(df$Gene1, df$Gene2)
  cor_label <- sprintf("r = %.3f, p = %.2e", cor_test$estimate, cor_test$p.value)

  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Gene1, y = Gene2))

  if (!is.null(color_by) && "Color" %in% colnames(df)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = Color), alpha = 0.6, size = 2)
  } else {
    p <- p + ggplot2::geom_point(alpha = 0.6, size = 2, color = "steelblue")
  }

  p <- p +
    ggplot2::geom_smooth(method = "lm", color = "red", se = TRUE) +
    ggplot2::labs(
      title = title %||% paste(gene1, "vs", gene2),
      subtitle = cor_label,
      x = paste(gene1, "Expression"),
      y = paste(gene2, "Expression")
    ) +
    theme_zinasuite()

  p
}


#' Visualize Gene-TMB Correlation
#'
#' @description
#' Create scatter plot showing correlation between gene expression and tumor mutation burden (TMB).
#'
#' @param gene Gene symbol
#' @param data_type Type of molecular data (default: "mRNA")
#' @param source Data source (default: "tcga")
#' @param title Plot title (auto-generated if NULL)
#' @param ... Additional parameters
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create TMB correlation plot for TP53
#' p <- vis_gene_tmb_cor("TP53")
#' print(p)
#' }
vis_gene_tmb_cor <- function(gene,
                             source = "tcga",
                             title = NULL,
                             ...) {
  # Query gene expression
  gene_expr <- query_gene_expression(gene, source = source)

  # Load TMB data
  tmb_data <- load_data("tcga_tmb")

  # Match samples
  common_samples <- intersect(names(gene_expr), tmb_data$Sample)

  if (length(common_samples) < 10) {
    stop("Insufficient samples with TMB data")
  }

  # Prepare data
  df <- data.frame(
    Gene = as.numeric(gene_expr[common_samples]),
    TMB = tmb_data$TMB[match(common_samples, tmb_data$Sample)],
    stringsAsFactors = FALSE
  )

  # Remove NA
  df <- df[stats::complete.cases(df), ]

  # Calculate correlation
  cor_test <- stats::cor.test(df$Gene, df$TMB)
  cor_label <- sprintf("r = %.3f, p = %.2e", cor_test$estimate, cor_test$p.value)

  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Gene, y = TMB)) +
    ggplot2::geom_point(alpha = 0.6, size = 2, color = "steelblue") +
    ggplot2::geom_smooth(method = "lm", color = "red", se = TRUE) +
    ggplot2::labs(
      title = title %||% paste(gene, "vs TMB"),
      subtitle = cor_label,
      x = paste(gene, "Expression"),
      y = "Tumor Mutation Burden"
    ) +
    theme_zinasuite()

  p
}


#' Visualize Gene-MSI Correlation
#'
#' @description
#' Create boxplot comparing gene expression between MSI-high and MSI-low samples.
#'
#' @param gene Gene symbol
#' @param data_type Type of molecular data (default: "mRNA")
#' @param source Data source (default: "tcga")
#' @param title Plot title (auto-generated if NULL)
#' @param ... Additional parameters
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create MSI comparison plot for TP53
#' p <- vis_gene_msi_cor("TP53")
#' print(p)
#' }
vis_gene_msi_cor <- function(gene,
                             source = "tcga",
                             title = NULL,
                             ...) {
  # Query gene expression
  gene_expr <- query_gene_expression(gene, source = source)

  # Load MSI data
  msi_data <- load_data("tcga_MSI")

  # Match samples
  common_samples <- intersect(names(gene_expr), msi_data$Sample)

  if (length(common_samples) < 10) {
    stop("Insufficient samples with MSI data")
  }

  # Prepare data
  df <- data.frame(
    Gene = as.numeric(gene_expr[common_samples]),
    MSI = msi_data$MSI[match(common_samples, msi_data$Sample)],
    stringsAsFactors = FALSE
  )

  # Remove NA
  df <- df[stats::complete.cases(df), ]

  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = MSI, y = Gene, fill = MSI)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::scale_fill_manual(values = get_palette(2, "default")) +
    ggplot2::labs(
      title = title %||% paste(gene, "Expression by MSI Status"),
      x = "MSI Status",
      y = paste(gene, "Expression")
    ) +
    theme_zinasuite()

  p
}
