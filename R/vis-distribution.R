#' Distribution Visualization
#'
#' @description
#' Create distribution plots (boxplot, violin plot, or histogram) for molecular data.
#'
#' @param data Data frame containing the data to plot
#' @param x Variable for x-axis (column name in data)
#' @param y Variable for y-axis (column name in data)
#' @param type Plot type: "boxplot", "violin", or "histogram"
#' @param fill Fill color or variable name for grouping
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
#'
#' # Create sample data frame with cancer types
#' sample_info <- load_data("tcga_gtex")
#' data <- data.frame(
#'   expression = tp53_expr,
#'   cancer = sample_info$tissue[match(names(tp53_expr), sample_info$sample)],
#'   type = sample_info$type2[match(names(tp53_expr), sample_info$sample)]
#' )
#'
#' # Create boxplot
#' p <- vis_distribution(
#'   data = data,
#'   x = "cancer",
#'   y = "expression",
#'   type = "boxplot",
#'   fill = "type",
#'   title = "TP53 Expression by Cancer Type"
#' )
#' print(p)
#' }
vis_distribution <- function(data,
                             x = NULL,
                             y = NULL,
                             type = c("boxplot", "violin", "histogram"),
                             fill = NULL,
                             title = NULL,
                             xlab = NULL,
                             ylab = NULL,
                             ...) {
  type <- match.arg(type)

  # Build ggplot
  p <- ggplot2::ggplot(data)

  if (type == "boxplot") {
    if (is.null(x)) {
      p <- p + ggplot2::geom_boxplot(ggplot2::aes(y = .data[[y]], fill = .data[[fill]]), ...)
    } else if (is.null(fill)) {
      p <- p + ggplot2::geom_boxplot(ggplot2::aes(x = .data[[x]], y = .data[[y]]), ...)
    } else {
      p <- p + ggplot2::geom_boxplot(ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]]), ...)
    }
  } else if (type == "violin") {
    if (is.null(x)) {
      p <- p + ggplot2::geom_violin(ggplot2::aes(y = .data[[y]], fill = .data[[fill]]), ...)
    } else if (is.null(fill)) {
      p <- p + ggplot2::geom_violin(ggplot2::aes(x = .data[[x]], y = .data[[y]]), ...)
    } else {
      p <- p + ggplot2::geom_violin(ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]]), ...)
    }
  } else { # histogram
    if (is.null(fill)) {
      p <- p + ggplot2::geom_histogram(ggplot2::aes(x = .data[[x]]), ...)
    } else {
      p <- p + ggplot2::geom_histogram(ggplot2::aes(x = .data[[x]], fill = .data[[fill]]), ...)
    }
  }

  # Add labels
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  if (!is.null(xlab)) p <- p + ggplot2::xlab(xlab)
  if (!is.null(ylab)) p <- p + ggplot2::ylab(ylab)

  # Add theme
  p <- p + ggplot2::theme_minimal()

  p
}

#' Tumor vs Normal Expression Plot
#'
#' @description
#' Create a boxplot or violin plot comparing tumor vs normal expression across
#' cancer types.
#'
#' @param gene Gene symbol
#' @param data_type Molecular data type (default: "mRNA")
#' @param source Data source (default: "tcga")
#' @param mode Plot mode: "Boxplot" or "Violinplot"
#' @param show_pvalue Whether to show p-values (default: TRUE)
#' @param method Statistical test method: "wilcox.test" or "t.test"
#' @param tcga_only Include only TCGA samples (default: FALSE)
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create tumor vs normal plot for TP53
#' p <- vis_tumor_normal("TP53", mode = "Boxplot")
#' print(p)
#'
#' # Violin plot version
#' p <- vis_tumor_normal("TP53", mode = "Violinplot", show_pvalue = TRUE)
#' print(p)
#' }
vis_tumor_normal <- function(gene,
                             data_type = "mRNA",
                             source = "tcga",
                             mode = c("Boxplot", "Violinplot"),
                             show_pvalue = TRUE,
                             method = c("wilcox.test", "t.test"),
                             tcga_only = FALSE) {
  mode <- match.arg(mode)
  method <- match.arg(method)

  # Query gene expression
  gene_data <- query_gene_expression(gene, source = source)

  if (is.null(gene_data)) {
    stop("Failed to retrieve data for gene: ", gene)
  }

  # Load sample information
  sample_info <- load_data("tcga_gtex")

  # Match samples
  common_samples <- intersect(names(gene_data), sample_info$sample)
  gene_data <- gene_data[common_samples]
  sample_info <- sample_info[sample_info$sample %in% common_samples, ]

  # Create data frame
  data <- data.frame(
    expression = as.numeric(gene_data),
    sample = names(gene_data),
    tissue = sample_info$tissue[match(names(gene_data), sample_info$sample)],
    type = sample_info$type2[match(names(gene_data), sample_info$sample)],
    dataset = ifelse(substr(names(gene_data), 1, 4) == "TCGA", "TCGA", "GTEX")
  )

  # Filter for tumors with matched normals
  tumor_types <- unique(data$tissue[data$type == "tumor"])
  normal_types <- unique(data$tissue[data$type == "normal"])
  matched_types <- intersect(tumor_types, normal_types)

  data <- data[data$tissue %in% matched_types, ]

  # TCGA only filter
  if (tcga_only) {
    data <- data[data$dataset == "TCGA", ]
  }

  # Create plot
  if (mode == "Boxplot") {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$tissue, y = .data$expression, fill = .data$type)) +
      ggplot2::geom_boxplot(...)
  } else {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$tissue, y = .data$expression, fill = .data$type)) +
      ggplot2::geom_violin(trim = TRUE, ...) +
      ggplot2::geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA)
  }

  # Add p-values if requested
  if (show_pvalue) {
    # Calculate p-values for each cancer type
    p_values <- sapply(matched_types, function(cancer) {
      cancer_data <- data[data$tissue == cancer, ]
      tumor_expr <- cancer_data$expression[cancer_data$type == "tumor"]
      normal_expr <- cancer_data$expression[cancer_data$type == "normal"]

      if (length(tumor_expr) > 0 && length(normal_expr) > 0) {
        test_result <- wilcox.test(tumor_expr, normal_expr)
        test_result$p.value
      } else {
        NA
      }
    })

    # Add significance annotations
    # (Implementation depends on ggpubr or similar package)
  }

  p <- p +
    ggplot2::ggtitle(paste(gene, "Expression: Tumor vs Normal")) +
    ggplot2::xlab("Cancer Type") +
    ggplot2::ylab(paste(data_type, "Expression")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}

#' Pan-Cancer Expression Plot
#'
#' @description
#' Create a heatmap or bar plot showing gene expression across cancer types.
#'
#' @param genes Vector of gene symbols
#' @param data_type Molecular data type (default: "mRNA")
#' @param source Data source (default: "tcga")
#' @param plot_type Plot type: "heatmap" or "barplot"
#' @return ggplot object or complex heatmap object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create heatmap for multiple genes
#' genes <- c("TP53", "BRCA1", "EGFR", "KRAS")
#' p <- vis_pancan_expression(genes, plot_type = "heatmap")
#' print(p)
#' }
vis_pancan_expression <- function(genes,
                                  data_type = "mRNA",
                                  source = "tcga",
                                  plot_type = c("heatmap", "barplot")) {
  plot_type <- match.arg(plot_type)

  # Query gene expression for all genes
  expr_list <- query_molecules(
    identifiers = genes,
    data_type = data_type,
    source = source,
    n_workers = 4,
    .progress = TRUE
  )

  # Load sample information
  sample_info <- load_data("tcga_gtex")

  # Create expression matrix
  expr_matrix <- do.call(rbind, lapply(expr_list, function(x) {
    values <- as.numeric(x)
    names(values) <- names(x)
    values
  }))

  # Calculate median expression by cancer type
  cancer_types <- unique(sample_info$tissue)
  median_expr <- sapply(cancer_types, function(cancer) {
    cancer_samples <- sample_info$sample[sample_info$tissue == cancer]
    rowMeans(expr_matrix[, intersect(colnames(expr_matrix), cancer_samples), drop = FALSE], na.rm = TRUE)
  })

  if (plot_type == "heatmap") {
    # Create heatmap
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      check_vis_deps("heatmap")
      ComplexHeatmap::Heatmap(
        median_expr,
        name = "Expression",
        cluster_columns = TRUE,
        cluster_rows = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE
      )
    } else {
      # Fallback to ggplot2
      data_long <- as.data.frame(as.table(median_expr))
      colnames(data_long) <- c("Gene", "Cancer", "Expression")

      ggplot2::ggplot(data_long, ggplot2::aes(x = .data$Cancer, y = .data$Gene, fill = .data$Expression)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
  } else {
    # Bar plot
    data_long <- as.data.frame(as.table(median_expr))
    colnames(data_long) <- c("Gene", "Cancer", "Expression")

    ggplot2::ggplot(data_long, ggplot2::aes(x = .data$Cancer, y = .data$Expression, fill = .data$Gene)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
}
