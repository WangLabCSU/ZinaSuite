#' Visualize Gene Expression on Anatomical Map
#'
#' @description
#' Create a visualization of gene expression levels across different anatomical
#' locations (organs) based on TCGA cancer types. This provides an intuitive
#' overview of where a gene is highly expressed in the human body.
#'
#' @param gene Gene symbol to visualize
#' @param data_type Type of molecular data: "mRNA", "protein", "cnv"
#' @param plot_type Type of plot: "heatmap", "dotplot", or "barplot"
#' @param cancers Vector of cancer types to include. If NULL, uses all available.
#' @param show_normal Whether to show normal tissue data from GTEx
#' @param color_palette Color palette for the plot
#' @param text_size Size of text labels
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create anatomical expression map for TP53
#' p <- vis_pancan_anatomy("TP53")
#' print(p)
#'
#' # Show protein expression
#' p <- vis_pancan_anatomy("TP53", data_type = "protein")
#'
#' # Include normal tissue comparison
#' p <- vis_pancan_anatomy("BRCA1", show_normal = TRUE)
#' }
vis_pancan_anatomy <- function(gene,
                                data_type = c("mRNA", "protein", "cnv"),
                                plot_type = c("heatmap", "dotplot", "barplot"),
                                cancers = NULL,
                                show_normal = FALSE,
                                color_palette = c("#2166ac", "#f7f7f7", "#b2182b"),
                                text_size = 10) {
  data_type <- match.arg(data_type)
  plot_type <- match.arg(plot_type)

  # Get gene expression data
  gene_expr <- switch(data_type,
    "mRNA" = query_gene_expression(gene, source = "tcga"),
    "protein" = query_protein(gene, source = "tcga"),
    "cnv" = query_cnv(gene, source = "tcga"),
    stop("Unsupported data_type: ", data_type)
  )

  if (is.null(gene_expr)) {
    stop("Failed to retrieve expression data for gene: ", gene)
  }

  # Get sample cancer types using barcode matching
  sample_info <- load_data("tcga_gtex")
  match_result <- match_samples(names(gene_expr), sample_info$Sample, "tcga", "tcga", match_by = "barcode")

  sample_cancer <- stats::setNames(
    sample_info$Tissue[match_result$idx2],
    match_result$common_ids
  )

  # Get organ mapping
  organ_map <- load_data("TCGA.organ")

  # Calculate median expression by cancer type
  cancer_types <- unique(sample_cancer)
  cancer_types <- cancer_types[!is.na(cancer_types)]

  if (!is.null(cancers)) {
    cancer_types <- intersect(cancer_types, cancers)
  }

  expr_by_cancer <- sapply(cancer_types, function(cancer) {
    samples <- names(sample_cancer)[sample_cancer == cancer]
    median(gene_expr[samples], na.rm = TRUE)
  })

  # Create plot data
  plot_data <- data.frame(
    Cancer = cancer_types,
    Expression = expr_by_cancer,
    stringsAsFactors = FALSE
  )

  # Merge with organ mapping
  plot_data <- merge(plot_data, organ_map, by = "Cancer", all.x = TRUE)

  # Remove rows without organ mapping
  plot_data <- plot_data[!is.na(plot_data$Organ), ]

  if (nrow(plot_data) == 0) {
    stop("No valid data after merging with organ mapping")
  }

  # Create the plot
  switch(plot_type,
    heatmap = plot_anatomy_heatmap(plot_data, gene, color_palette, text_size),
    dotplot = plot_anatomy_dotplot(plot_data, gene, color_palette, text_size),
    barplot = plot_anatomy_barplot(plot_data, gene, color_palette, text_size)
  )
}

#' Plot Anatomy Heatmap
#'
#' @param data Plot data frame
#' @param gene Gene name
#' @param color_palette Color palette
#' @param text_size Text size
#' @return ggplot object
#' @keywords internal
plot_anatomy_heatmap <- function(data, gene, color_palette, text_size) {
  # Order by organ system
  data$System <- factor(data$System, levels = unique(data$System))
  data$Organ <- factor(data$Organ, levels = unique(data$Organ))

  ggplot2::ggplot(data, ggplot2::aes(x = .data$Cancer, y = .data$Organ, fill = .data$Expression)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = color_palette[1],
      mid = color_palette[2],
      high = color_palette[3],
      midpoint = stats::median(data$Expression, na.rm = TRUE),
      name = paste(gene, "Expression")
    ) +
    ggplot2::facet_grid(.data$System ~ ., scales = "free_y", space = "free") +
    ggplot2::labs(
      title = paste("Anatomical Expression Map:", gene),
      x = "Cancer Type",
      y = "Anatomical Location"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = text_size),
      axis.text.y = ggplot2::element_text(size = text_size),
      strip.text = ggplot2::element_text(face = "bold", size = text_size + 2),
      panel.grid = ggplot2::element_blank()
    )
}

#' Plot Anatomy Dotplot
#'
#' @param data Plot data frame
#' @param gene Gene name
#' @param color_palette Color palette
#' @param text_size Text size
#' @return ggplot object
#' @keywords internal
plot_anatomy_dotplot <- function(data, gene, color_palette, text_size) {
  # Calculate expression quartile for sizing
  data$Size <- cut(
    data$Expression,
    breaks = stats::quantile(data$Expression, probs = c(0, 0.25, 0.5, 0.75, 1)),
    labels = c("Q1", "Q2", "Q3", "Q4"),
    include.lowest = TRUE
  )

  ggplot2::ggplot(data, ggplot2::aes(x = .data$Cancer, y = .data$Organ)) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$Size, color = .data$Expression)
    ) +
    ggplot2::scale_color_gradient2(
      low = color_palette[1],
      mid = color_palette[2],
      high = color_palette[3],
      midpoint = stats::median(data$Expression, na.rm = TRUE),
      name = paste(gene, "Expression")
    ) +
    ggplot2::scale_size_manual(
      values = c("Q1" = 2, "Q2" = 4, "Q3" = 6, "Q4" = 8),
      name = "Expression\nQuartile"
    ) +
    ggplot2::facet_grid(.data$System ~ ., scales = "free_y", space = "free") +
    ggplot2::labs(
      title = paste("Anatomical Expression Map:", gene),
      x = "Cancer Type",
      y = "Anatomical Location"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = text_size),
      axis.text.y = ggplot2::element_text(size = text_size),
      strip.text = ggplot2::element_text(face = "bold", size = text_size + 2)
    )
}

#' Plot Anatomy Barplot
#'
#' @param data Plot data frame
#' @param gene Gene name
#' @param color_palette Color palette
#' @param text_size Text size
#' @return ggplot object
#' @keywords internal
plot_anatomy_barplot <- function(data, gene, color_palette, text_size) {
  # Order by expression level
  data$Cancer <- stats::reorder(data$Cancer, data$Expression)

  ggplot2::ggplot(data, ggplot2::aes(x = .data$Cancer, y = .data$Expression, fill = .data$Expression)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient2(
      low = color_palette[1],
      mid = color_palette[2],
      high = color_palette[3],
      midpoint = stats::median(data$Expression, na.rm = TRUE),
      name = paste(gene, "Expression")
    ) +
    ggplot2::facet_grid(.data$System ~ ., scales = "free_y", space = "free") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Expression by Anatomical Location:", gene),
      x = "Cancer Type",
      y = "Median Expression"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = text_size),
      strip.text = ggplot2::element_text(face = "bold", size = text_size + 2),
      legend.position = "right"
    )
}

#' Get Organ Information for Cancer Types
#'
#' @description
#' Retrieve anatomical and system information for TCGA cancer types.
#'
#' @param cancers Vector of cancer types. If NULL, returns all.
#' @return Data frame with Cancer, Organ, and System columns
#' @export
#'
#' @examples
#' \dontrun{
#' # Get organ info for all cancers
#' organ_info <- get_organ_info()
#'
#' # Get specific cancer info
#' brca_info <- get_organ_info("BRCA")
#' }
get_organ_info <- function(cancers = NULL) {
  organ_map <- load_data("TCGA.organ")

  if (!is.null(cancers)) {
    organ_map <- organ_map[organ_map$Cancer %in% cancers, ]
  }

  organ_map
}
