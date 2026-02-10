#' Visualize Gene-Immune Cell Correlation
#'
#' @description
#' Create a heatmap showing correlation between gene expression and immune cell
#' infiltration levels across cancer types.
#'
#' @param gene Gene symbol to analyze
#' @param data_type Molecular data type (default: "mRNA")
#' @param immune_features Vector of immune features to include. If NULL, uses all available.
#' @param cancers Vector of cancer types to include. If NULL, uses all available.
#' @param method Correlation method: "pearson", "spearman", or "kendall"
#' @param adjust_method Multiple testing correction method
#' @param plot_type Type of plot: "heatmap", "dotplot", or "barplot"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' # Create immune correlation heatmap for TP53
#' p <- vis_gene_immune_cor("TP53")
#' print(p)
#'
#' # Analyze specific immune features
#' p <- vis_gene_immune_cor(
#'   gene = "BRCA1",
#'   immune_features = c("CD8_T_cells", "Macrophages", "Tregs"),
#'   cancers = c("BRCA", "LUAD", "LUSC")
#' )
#' }
vis_gene_immune_cor <- function(gene,
                                 data_type = "mRNA",
                                 immune_features = NULL,
                                 cancers = NULL,
                                 method = c("spearman", "pearson", "kendall"),
                                 adjust_method = "fdr",
                                 plot_type = c("heatmap", "dotplot", "barplot")) {
  method <- match.arg(method)
  plot_type <- match.arg(plot_type)

  # Get gene expression data
  gene_expr <- query_gene_expression(gene, source = "tcga")

  if (is.null(gene_expr)) {
    stop("Failed to retrieve expression data for gene: ", gene)
  }

  # Get immune infiltration data
  immune_data <- load_data("tcga_TIL")

  if (is.null(immune_features)) {
    immune_features <- setdiff(colnames(immune_data), c("Sample", "Cancer"))
  }

  # Get sample cancer types
  sample_info <- load_data("tcga_gtex")
  sample_cancer <- setNames(
    sample_info$tissue[match(names(gene_expr), sample_info$sample)],
    names(gene_expr)
  )

  # Filter by cancers if specified
  if (!is.null(cancers)) {
    keep <- sample_cancer %in% cancers
    gene_expr <- gene_expr[keep]
    sample_cancer <- sample_cancer[keep]
  }

  # Calculate correlations by cancer type
  results <- list()

  for (cancer in unique(sample_cancer)) {
    if (is.na(cancer)) next

    # Get samples for this cancer
    cancer_samples <- names(sample_cancer)[sample_cancer == cancer]

    # Get gene expression for these samples
    cancer_expr <- gene_expr[names(gene_expr) %in% cancer_samples]

    if (length(cancer_expr) < 10) next

    # Get immune data for these samples
    cancer_immune <- immune_data[immune_data$Sample %in% cancer_samples, ]

    # Calculate correlations for each immune feature
    for (feature in immune_features) {
      if (!feature %in% colnames(cancer_immune)) next

      # Match samples
      common_samples <- intersect(names(cancer_expr), cancer_immune$Sample)
      if (length(common_samples) < 10) next

      x <- cancer_expr[common_samples]
      y <- setNames(cancer_immune[[feature]], cancer_immune$Sample)[common_samples]

      # Remove NA
      valid <- complete.cases(x, y)
      if (sum(valid) < 10) next

      cor_result <- cor.test(x[valid], y[valid], method = method)

      results[[length(results) + 1]] <- data.frame(
        Cancer = cancer,
        Feature = feature,
        Correlation = cor_result$estimate,
        Pvalue = cor_result$p.value,
        N = sum(valid),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) == 0) {
    stop("No valid correlations could be calculated")
  }

  result_df <- do.call(rbind, results)
  result_df$Padj <- p.adjust(result_df$Pvalue, method = adjust_method)
  result_df$Significance <- cut(
    result_df$Padj,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
  )

  # Create plot
  switch(plot_type,
    heatmap = plot_immune_heatmap(result_df, gene),
    dotplot = plot_immune_dotplot(result_df, gene),
    barplot = plot_immune_barplot(result_df, gene)
  )
}

#' Plot Immune Correlation Heatmap
#'
#' @param data Correlation results data frame
#' @param gene Gene name for title
#' @return ggplot object
#' @keywords internal
plot_immune_heatmap <- function(data, gene) {
  ggplot2::ggplot(data, ggplot2::aes(x = Cancer, y = Feature, fill = Correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = c(-1, 1),
      name = "Correlation"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Significance),
      color = "black", size = 3
    ) +
    ggplot2::labs(
      title = paste("Gene-Immune Correlation:", gene),
      x = "Cancer Type",
      y = "Immune Feature"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

#' Plot Immune Correlation Dotplot
#'
#' @param data Correlation results data frame
#' @param gene Gene name for title
#' @return ggplot object
#' @keywords internal
plot_immune_dotplot <- function(data, gene) {
  ggplot2::ggplot(data, ggplot2::aes(x = Cancer, y = Feature)) +
    ggplot2::geom_point(
      ggplot2::aes(size = abs(Correlation), color = Correlation)
    ) +
    ggplot2::scale_color_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = c(-1, 1)
    ) +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::labs(
      title = paste("Gene-Immune Correlation:", gene),
      x = "Cancer Type",
      y = "Immune Feature",
      size = "|Correlation|"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Plot Immune Correlation Barplot
#'
#' @param data Correlation results data frame
#' @param gene Gene name for title
#' @return ggplot object
#' @keywords internal
plot_immune_barplot <- function(data, gene) {
  # Select top correlations per cancer
  data_top <- data %>%
    dplyr::group_by(Cancer) %>%
    dplyr::slice_max(order_by = abs(Correlation), n = 5)

  ggplot2::ggplot(data_top, ggplot2::aes(x = reorder(Feature, Correlation), y = Correlation)) +
    ggplot2::geom_col(ggplot2::aes(fill = Correlation > 0)) +
    ggplot2::facet_wrap(~Cancer, scales = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
      labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
      name = "Direction"
    ) +
    ggplot2::labs(
      title = paste("Top Gene-Immune Correlations:", gene),
      x = "Immune Feature",
      y = "Correlation Coefficient"
    ) +
    ggplot2::theme_minimal()
}

#' Visualize Gene-TIL Correlation
#'
#' @description
#' Visualize correlation between gene expression and tumor-infiltrating
#' lymphocyte (TIL) fractions.
#'
#' @param gene Gene symbol to analyze
#' @param cell_type Specific TIL cell type or NULL for all types
#' @param cancers Vector of cancer types to include
#' @param method Correlation method
#' @param plot_type Type of visualization
#'
#' @return A ggplot object or list of plots
#' @export
#'
#' @examples
#' \donttest{
#' # Analyze TP53 correlation with all TIL types
#' p <- vis_gene_TIL_cor("TP53")
#'
#' # Focus on CD8 T cells
#' p <- vis_gene_TIL_cor("TP53", cell_type = "CD8_T_cells")
#' }
vis_gene_TIL_cor <- function(gene,
                              cell_type = NULL,
                              cancers = NULL,
                              method = c("spearman", "pearson"),
                              plot_type = c("heatmap", "scatter", "boxplot")) {
  method <- match.arg(method)
  plot_type <- match.arg(plot_type)

  # Get gene expression
  gene_expr <- query_gene_expression(gene, source = "tcga")

  # Get TIL data
  til_data <- load_data("tcga_TIL")

  # Determine cell types to analyze
  cell_types <- if (is.null(cell_type)) {
    setdiff(colnames(til_data), c("Sample", "Cancer"))
  } else {
    cell_type
  }

  # Get sample cancer info
  sample_info <- load_data("tcga_gtex")
  sample_cancer <- setNames(
    sample_info$tissue[match(names(gene_expr), sample_info$sample)],
    names(gene_expr)
  )

  # Filter by cancers
  if (!is.null(cancers)) {
    keep <- sample_cancer %in% cancers
    gene_expr <- gene_expr[keep]
    sample_cancer <- sample_cancer[keep]
    til_data <- til_data[til_data$Cancer %in% cancers, ]
  }

  switch(plot_type,
    heatmap = {
      # Use the immune correlation function
      vis_gene_immune_cor(
        gene = gene,
        immune_features = cell_types,
        cancers = cancers,
        method = method,
        plot_type = "heatmap"
      )
    },
    scatter = {
      # Create scatter plot for specific cell type
      if (is.null(cell_type) || length(cell_type) > 1) {
        stop("For scatter plot, specify a single cell_type")
      }

      # Merge data
      plot_data <- data.frame(
        Sample = til_data$Sample,
        TIL = til_data[[cell_type]],
        Cancer = til_data$Cancer,
        stringsAsFactors = FALSE
      )

      plot_data$Expression <- gene_expr[plot_data$Sample]
      plot_data <- plot_data[complete.cases(plot_data), ]

      ggplot2::ggplot(plot_data, ggplot2::aes(x = Expression, y = TIL)) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::geom_smooth(method = "lm", color = "red") +
        ggplot2::facet_wrap(~Cancer, scales = "free") +
        ggplot2::labs(
          title = paste(gene, "vs", cell_type),
          x = paste(gene, "Expression"),
          y = paste(cell_type, "Fraction")
        ) +
        ggplot2::theme_minimal()
    },
    boxplot = {
      # Create boxplot of TIL by gene expression quartiles
      if (is.null(cell_type) || length(cell_type) > 1) {
        stop("For boxplot, specify a single cell_type")
      }

      plot_data <- data.frame(
        Sample = til_data$Sample,
        TIL = til_data[[cell_type]],
        Cancer = til_data$Cancer,
        stringsAsFactors = FALSE
      )

      plot_data$Expression <- gene_expr[plot_data$Sample]
      plot_data <- plot_data[complete.cases(plot_data), ]

      # Create expression quartiles by cancer
      plot_data <- plot_data %>%
        dplyr::group_by(Cancer) %>%
        dplyr::mutate(
          Expr_Quartile = cut(
            Expression,
            breaks = quantile(Expression, probs = c(0, 0.25, 0.5, 0.75, 1)),
            labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
            include.lowest = TRUE
          )
        )

      ggplot2::ggplot(plot_data, ggplot2::aes(x = Expr_Quartile, y = TIL)) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = Expr_Quartile)) +
        ggplot2::facet_wrap(~Cancer, scales = "free_y") +
        ggplot2::labs(
          title = paste(cell_type, "by", gene, "Expression Quartiles"),
          x = paste(gene, "Expression Quartile"),
          y = paste(cell_type, "Fraction")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )
    }
  )
}

#' Visualize Gene-TMB Correlation
#'
#' @description
#' Visualize correlation between gene expression and tumor mutation burden (TMB).
#'
#' @param gene Gene symbol to analyze
#' @param cancers Vector of cancer types to include
#' @param method Correlation method
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' # Analyze TP53 correlation with TMB
#' p <- vis_gene_tmb_cor("TP53")
#'
#' # Focus on specific cancers
#' p <- vis_gene_tmb_cor("TP53", cancers = c("LUAD", "LUSC", "BRCA"))
#' }
vis_gene_tmb_cor <- function(gene,
                              cancers = NULL,
                              method = c("spearman", "pearson")) {
  method <- match.arg(method)

  # Get gene expression
  gene_expr <- query_gene_expression(gene, source = "tcga")

  # Get TMB data
  tmb_data <- load_data("tcga_tmb")

  # Merge data
  plot_data <- data.frame(
    Sample = tmb_data$Sample,
    TMB = tmb_data$tmb,
    Cancer = tmb_data$Cancer,
    stringsAsFactors = FALSE
  )

  plot_data$Expression <- gene_expr[plot_data$Sample]
  plot_data <- plot_data[complete.cases(plot_data), ]

  # Filter by cancers
  if (!is.null(cancers)) {
    plot_data <- plot_data[plot_data$Cancer %in% cancers, ]
  }

  # Calculate correlation by cancer
  cor_by_cancer <- plot_data %>%
    dplyr::group_by(Cancer) %>%
    dplyr::summarise(
      Correlation = cor(Expression, log1p(TMB), method = method),
      Pvalue = cor.test(Expression, log1p(TMB), method = method)$p.value,
      N = dplyr::n(),
      .groups = "drop"
    )

  # Create combined plot
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Expression, y = log1p(TMB))) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_smooth(method = "lm", color = "red") +
    ggplot2::facet_wrap(~Cancer, scales = "free") +
    ggplot2::labs(
      title = paste(gene, "Expression vs TMB"),
      x = paste(gene, "Expression"),
      y = "log(TMB + 1)"
    ) +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(
    cor_by_cancer,
    ggplot2::aes(x = reorder(Cancer, Correlation), y = Correlation, fill = Correlation > 0)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac")
    ) +
    ggplot2::labs(
      title = "Correlation by Cancer Type",
      x = "Cancer",
      y = "Correlation Coefficient"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # Return both plots (user can arrange with patchwork or cowplot)
  list(scatter = p1, bar = p2, correlation = cor_by_cancer)
}

#' Visualize Gene-MSI Correlation
#'
#' @description
#' Visualize correlation between gene expression and microsatellite instability (MSI).
#'
#' @param gene Gene symbol to analyze
#' @param cancers Vector of cancer types to include
#' @param plot_type Type of plot: "boxplot" or "violin"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' # Analyze TP53 expression by MSI status
#' p <- vis_gene_msi_cor("TP53")
#' }
vis_gene_msi_cor <- function(gene,
                              cancers = NULL,
                              plot_type = c("boxplot", "violin")) {
  plot_type <- match.arg(plot_type)

  # Get gene expression
  gene_expr <- query_gene_expression(gene, source = "tcga")

  # Get MSI data
  msi_data <- load_data("tcga_MSI")

  # Merge data
  plot_data <- data.frame(
    Sample = msi_data$Sample,
    MSI_Score = msi_data$msi_score,
    MSI_Status = msi_data$msi_status,
    Cancer = msi_data$Cancer,
    stringsAsFactors = FALSE
  )

  plot_data$Expression <- gene_expr[plot_data$Sample]
  plot_data <- plot_data[complete.cases(plot_data), ]

  # Filter by cancers
  if (!is.null(cancers)) {
    plot_data <- plot_data[plot_data$Cancer %in% cancers, ]
  }

  # Remove NA MSI status
  plot_data <- plot_data[!is.na(plot_data$MSI_Status), ]

  if (plot_type == "boxplot") {
    ggplot2::ggplot(plot_data, ggplot2::aes(x = MSI_Status, y = Expression, fill = MSI_Status)) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~Cancer, scales = "free_y") +
      ggplot2::labs(
        title = paste(gene, "Expression by MSI Status"),
        x = "MSI Status",
        y = paste(gene, "Expression")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  } else {
    ggplot2::ggplot(plot_data, ggplot2::aes(x = MSI_Status, y = Expression, fill = MSI_Status)) +
      ggplot2::geom_violin() +
      ggplot2::geom_boxplot(width = 0.2) +
      ggplot2::facet_wrap(~Cancer, scales = "free_y") +
      ggplot2::labs(
        title = paste(gene, "Expression by MSI Status"),
        x = "MSI Status",
        y = paste(gene, "Expression")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  }
}

#' Visualize Gene-Stemness Correlation
#'
#' @description
#' Visualize correlation between gene expression and tumor stemness scores.
#'
#' @param gene Gene symbol to analyze
#' @param cancers Vector of cancer types to include
#' @param stemness_type Type of stemness score: "mRNA", "methylation", or "combined"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' # Analyze TP53 correlation with stemness
#' p <- vis_gene_stemness_cor("TP53")
#'
#' # Use methylation-based stemness
#' p <- vis_gene_stemness_cor("TP53", stemness_type = "methylation")
#' }
vis_gene_stemness_cor <- function(gene,
                                   cancers = NULL,
                                   stemness_type = c("combined", "mRNA", "methylation")) {
  stemness_type <- match.arg(stemness_type)

  # Get gene expression
  gene_expr <- query_gene_expression(gene, source = "tcga")

  # Get stemness data
  stem_data <- load_data("tcga_stemness")

  # Select stemness column
  stem_col <- switch(stemness_type,
    combined = "stemness",
    mRNA = "stemness_mRNA",
    methylation = "stemness_methylation"
  )

  # Merge data
  plot_data <- data.frame(
    Sample = stem_data$Sample,
    Stemness = stem_data[[stem_col]],
    Cancer = stem_data$Cancer,
    stringsAsFactors = FALSE
  )

  plot_data$Expression <- gene_expr[plot_data$Sample]
  plot_data <- plot_data[complete.cases(plot_data), ]

  # Filter by cancers
  if (!is.null(cancers)) {
    plot_data <- plot_data[plot_data$Cancer %in% cancers, ]
  }

  # Calculate correlation
  cor_result <- cor.test(plot_data$Expression, plot_data$Stemness, method = "spearman")

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Expression, y = Stemness)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_smooth(method = "lm", color = "red") +
    ggplot2::facet_wrap(~Cancer, scales = "free") +
    ggplot2::labs(
      title = paste(gene, "vs", stemness_type, "Stemness"),
      subtitle = paste(
        "Overall correlation:",
        sprintf("r = %.3f, p = %.3e", cor_result$estimate, cor_result$p.value)
      ),
      x = paste(gene, "Expression"),
      y = paste(stemness_type, "Stemness Score")
    ) +
    ggplot2::theme_minimal()
}
