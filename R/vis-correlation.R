#' Correlation Visualization
#'
#' @description
#' Create scatter plots and heatmaps for correlation visualization.
#'
#' @param x First variable (numeric vector) or data frame for heatmap
#' @param y Second variable (numeric vector, optional for heatmap)
#' @param type Plot type: "scatter" or "heatmap"
#' @param method Correlation method for annotation
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param ... Additional arguments passed to ggplot
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Query gene expression data
#' tp53_expr <- query_gene_expression("TP53", source = "tcga")
#' brca1_expr <- query_gene_expression("BRCA1", source = "tcga")
#'
#' # Create scatter plot
#' p <- vis_correlation(tp53_expr, brca1_expr, type = "scatter")
#' print(p)
#'
#' # Create correlation heatmap for multiple genes
#' genes <- c("TP53", "BRCA1", "EGFR", "KRAS")
#' expr_data <- sapply(genes, query_gene_expression, source = "tcga")
#' p <- vis_correlation(expr_data, type = "heatmap")
#' print(p)
#' }
vis_correlation <- function(x, y = NULL,
                            type = c("scatter", "heatmap"),
                            method = c("pearson", "spearman"),
                            title = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            ...) {
  type <- match.arg(type)
  method <- match.arg(method)

  if (type == "scatter") {
    # Scatter plot
    if (is.null(y)) {
      stop("y is required for scatter plot")
    }

    # Remove NA values
    complete_cases <- stats::complete.cases(x, y)
    x <- x[complete_cases]
    y <- y[complete_cases]

    # Calculate correlation
    cor_result <- stats::cor.test(x, y, method = method)
    cor_value <- cor_result$estimate
    p_value <- cor_result$p.value

    # Create data frame
    data <- data.frame(x = x, y = y)

    # Build plot
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_point(alpha = 0.5, ...) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red") +
      ggplot2::labs(
        title = title %||% paste("Correlation Plot (", method, ")"),
        subtitle = sprintf("r = %.3f, p = %.2e", cor_value, p_value),
        x = xlab %||% "X",
        y = ylab %||% "Y"
      ) +
      ggplot2::theme_minimal()

  } else {
    # Heatmap
    if (is.null(y)) {
      # x is a data frame or matrix
      cor_matrix <- stats::cor(x, method = method, use = "pairwise.complete.obs")
    } else {
      # Combine x and y into matrix
      data <- cbind(x, y)
      cor_matrix <- stats::cor(data, method = method, use = "pairwise.complete.obs")
    }

    # Convert to long format
    cor_df <- as.data.frame(as.table(cor_matrix))
    colnames(cor_df) <- c("Var1", "Var2", "Correlation")

    # Build plot
    p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$Correlation)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$Correlation)), size = 3) +
      ggplot2::scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, limits = c(-1, 1)
      ) +
      ggplot2::labs(
        title = title %||% paste("Correlation Matrix (", method, ")"),
        x = NULL, y = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text.y = ggplot2::element_text()
      )
  }

  p
}

#' Gene Correlation Plot
#'
#' @description
#' Create a scatter plot showing correlation between two genes.
#'
#' @param gene1 First gene symbol
#' @param gene2 Second gene symbol
#' @param data_type Molecular data type (default: "mRNA")
#' @param source Data source (default: "tcga")
#' @param color_by Optional variable to color points by (e.g., cancer type)
#' @param method Correlation method (default: "pearson")
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot correlation between TP53 and BRCA1
#' p <- vis_gene_correlation("TP53", "BRCA1", method = "pearson")
#' print(p)
#'
#' # Color by cancer type
#' p <- vis_gene_correlation("TP53", "BRCA1", color_by = "cancer")
#' print(p)
#' }
vis_gene_correlation <- function(gene1, gene2,
                                 data_type = "mRNA",
                                 source = "tcga",
                                 color_by = NULL,
                                 method = c("pearson", "spearman")) {
  method <- match.arg(method)

  # Query gene expression
  expr1 <- query_gene_expression(gene1, source = source)
  expr2 <- query_gene_expression(gene2, source = source)

  if (is.null(expr1) || is.null(expr2)) {
    stop("Failed to retrieve data for one or both genes")
  }

  # Match samples using barcode matching
  match_result <- match_samples(names(expr1), names(expr2), "tcga", "tcga", match_by = "barcode")

  if (match_result$n_matched == 0) {
    stop("No common samples found between the two genes")
  }

  # Create data frame
  data <- data.frame(
    gene1 = as.numeric(expr1[match_result$idx1]),
    gene2 = as.numeric(expr2[match_result$idx2]),
    sample = match_result$common_ids
  )

  # Add color variable if specified
  if (!is.null(color_by)) {
    if (color_by == "cancer") {
      sample_info <- load_data("tcga_gtex")
      # Match using barcode
      sample_bc <- substr(sample_info$Sample, 1, 12)
      data_bc <- substr(data$sample, 1, 12)
      data$cancer <- sample_info$Tissue[match(data_bc, sample_bc)]
    }
  }

  # Calculate correlation
  cor_result <- stats::cor.test(data$gene1, data$gene2, method = method)

  # Build plot
  if (!is.null(color_by)) {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$gene1, y = .data$gene2, color = .data[[color_by]])) +
      ggplot2::geom_point(alpha = 0.6)
  } else {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$gene1, y = .data$gene2)) +
      ggplot2::geom_point(alpha = 0.5, color = "steelblue")
  }

  p <- p +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red") +
    ggplot2::labs(
      title = paste(gene1, "vs", gene2),
      subtitle = sprintf("%s correlation: r = %.3f, p = %.2e", method, cor_result$estimate, cor_result$p.value),
      x = paste(gene1, "expression"),
      y = paste(gene2, "expression")
    ) +
    ggplot2::theme_minimal()

  p
}

#' Pan-Cancer Correlation Plot
#'
#' @description
#' Create a plot showing correlation between two genes across cancer types.
#'
#' @param gene1 First gene symbol
#' @param gene2 Second gene symbol
#' @param data_type Molecular data type (default: "mRNA")
#' @param source Data source (default: "tcga")
#' @param method Correlation method (default: "pearson")
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot pan-cancer correlation
#' p <- vis_pancan_correlation("TP53", "BRCA1")
#' print(p)
#' }
vis_pancan_correlation <- function(gene1, gene2,
                                   data_type = "mRNA",
                                   source = "tcga",
                                   method = c("pearson", "spearman")) {
  method <- match.arg(method)

  # Query gene expression
  expr1 <- query_gene_expression(gene1, source = source)
  expr2 <- query_gene_expression(gene2, source = source)

  if (is.null(expr1) || is.null(expr2)) {
    stop("Failed to retrieve data for one or both genes")
  }

  # Load sample information
  sample_info <- load_data("tcga_gtex")

  # Match samples
  common_samples <- intersect(names(expr1), names(expr2))
  common_samples <- intersect(common_samples, sample_info$sample)

  expr1 <- expr1[common_samples]
  expr2 <- expr2[common_samples]

  # Get cancer types
  cancer_types <- sample_info$tissue[match(common_samples, sample_info$sample)]

  # Create data frame
  data <- data.frame(
    gene1 = as.numeric(expr1),
    gene2 = as.numeric(expr2),
    cancer = cancer_types,
    sample = common_samples
  )

  # Calculate correlation by cancer type
  cor_by_cancer <- data %>%
    dplyr::group_by(.data$cancer) %>%
    dplyr::summarise(
      cor = cor(.data$gene1, .data$gene2, method = method, use = "complete.obs"),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n >= 10) %>%
    dplyr::arrange(.data$cor)

  # Create plot
  p <- ggplot2::ggplot(cor_by_cancer, ggplot2::aes(x = stats::reorder(.data$cancer, .data$cor), y = .data$cor)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Pan-Cancer Correlation:", gene1, "vs", gene2),
      subtitle = paste("Method:", method),
      x = "Cancer Type",
      y = "Correlation Coefficient"
    ) +
    ggplot2::theme_minimal()

  p
}
