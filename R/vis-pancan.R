#' Visualize Gene Expression in Tumor vs Normal
#'
#' @description
#' Create boxplot or violin plot comparing gene expression between tumor and normal samples.
#'
#' @param gene Gene symbol
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
  p <- ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(.data$gene, .data$hr), y = .data$hr))

  p <- p +
    ggplot2::geom_point(ggplot2::aes(color = .data$significant), size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper, color = .data$significant), width = 0.2) +
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


#' Visualize Gene Correlation Scatter Plot
#'
#' @description
#' Create scatter plot showing correlation between two genes.
#'
#' @param gene1 First gene symbol
#' @param gene2 Second gene symbol
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
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Gene1, y = .data$Gene2))

  if (!is.null(color_by) && "Color" %in% colnames(df)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Color), alpha = 0.6, size = 2)
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

#' Visualize Gene and Pathway Correlation
#'
#' @description
#' Create scatter plot showing correlation between gene expression and pathway activity score.
#'
#' @param gene Gene symbol (e.g., "TP53")
#' @param data_type Data type: "mRNA", "transcript", "protein", "methylation", "miRNA"
#' @param pw_name Pathway name (e.g., "HALLMARK_APOPTOSIS")
#' @param cancer_choose Cancer type(s) to include
#' @param cor_method Correlation method: "spearman" or "pearson"
#' @param use_regline Whether to add regression line (default: TRUE)
#' @param alpha Point transparency (default: 0.5)
#' @param color Point color (default: "#000000")
#' @param filter_tumor Whether to filter to tumor samples only (default: TRUE)
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Visualize TP53 correlation with apoptosis pathway
#' p <- vis_gene_pw_cor(
#'   gene = "TP53",
#'   data_type = "mRNA",
#'   pw_name = "HALLMARK_APOPTOSIS",
#'   cancer_choose = "BRCA"
#' )
#' print(p)
#' }
vis_gene_pw_cor <- function(gene = "TP53",
                            data_type = "mRNA",
                            pw_name = "HALLMARK_APOPTOSIS",
                            cancer_choose = "BRCA",
                            cor_method = c("spearman", "pearson"),
                            use_regline = TRUE,
                            alpha = 0.5,
                            color = "#000000",
                            filter_tumor = TRUE) {
  cor_method <- match.arg(cor_method)

  # Load pathway data
  tcga_pw <- load_data("tcga_PW")
  tcga_pw_meta <- load_data("tcga_PW_meta")

  # Validate pathway name
  if (!is.null(pw_name)) {
    if (!(pw_name %in% tcga_pw_meta$ID)) {
      stop("Invalid pathway name. See load_data('tcga_PW_meta') for available pathways.")
    }
  }

  # Query gene expression
  gene_expr <- query_molecule_value(gene, data_type = data_type, source = "tcga")

  if (is.null(gene_expr) || all(is.na(gene_expr))) {
    warning("No gene expression data available for ", gene)
    return(NULL)
  }

  message("Retrieved expression data for ", gene)

  # Load sample information
  tcga_gtex <- load_data("tcga_gtex")

  # Filter samples
  if (filter_tumor) {
    sample_filter <- tcga_gtex$type2 == "tumor"
  } else {
    sample_filter <- TRUE
  }

  # Filter by cancer type
  cancer_filter <- tcga_gtex$tissue %in% cancer_choose

  # Get filtered samples
  filtered_samples <- tcga_gtex$sample[sample_filter & cancer_filter]

  # Prepare gene expression data
  df <- data.frame(
    Sample = names(gene_expr),
    Gene = as.numeric(gene_expr),
    stringsAsFactors = FALSE
  )

  # Add pathway scores
  pw_samples <- rownames(tcga_pw)
  common_samples <- intersect(df$Sample, pw_samples)
  common_samples <- intersect(common_samples, filtered_samples)

  if (length(common_samples) < 10) {
    stop("Insufficient samples for analysis (need at least 10)")
  }

  # Extract pathway scores
  pw_scores <- tcga_pw[common_samples, pw_name]

  # Prepare final data
  plot_data <- data.frame(
    Sample = common_samples,
    Gene = df$Gene[match(common_samples, df$Sample)],
    Pathway = as.numeric(pw_scores),
    Cancer = tcga_gtex$tissue[match(common_samples, tcga_gtex$sample)],
    stringsAsFactors = FALSE
  )

  # Remove NAs
  plot_data <- plot_data[stats::complete.cases(plot_data), ]

  if (nrow(plot_data) < 10) {
    stop("Insufficient valid samples after removing NAs")
  }

  # Calculate correlation
  cor_res <- stats::cor.test(plot_data$Gene, plot_data$Pathway, method = cor_method)
  cor_label <- sprintf("%s: r=%.3f, p=%.2e", cor_method, cor_res$estimate, cor_res$p.value)

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Gene, y = .data$Pathway)) +
    ggplot2::geom_point(alpha = alpha, color = color, size = 3) +
    ggplot2::labs(
      x = paste(gene, data_type),
      y = pw_name,
      title = paste(gene, "vs", pw_name),
      subtitle = cor_label
    ) +
    theme_zinasuite()

  if (use_regline) {
    p <- p + ggplot2::geom_smooth(method = "lm", color = "red", se = TRUE)
  }

  # Store data as attribute for download
  attr(p, "data") <- plot_data

  p
}

