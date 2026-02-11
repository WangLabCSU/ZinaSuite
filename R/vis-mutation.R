#' Visualize Mutation vs Wild-Type Comparison
#'
#' @description
#' Compare molecular features (expression, CNV, etc.) between mutant and wild-type
#' samples for a given gene. Supports visualization of tumor vs normal comparisons
#' as well.
#'
#' @param gene Gene symbol to analyze
#' @param target_gene Gene for which to compare expression/molecular features.
#'   If NULL, uses the same as `gene` (self-comparison).
#' @param data_type Type of molecular data to compare: "mRNA", "protein", "cnv", "methylation"
#' @param cancers Vector of cancer types to include. If NULL, uses all available.
#' @param plot_type Type of plot: "boxplot", "violin", "barplot", or "dotplot"
#' @param test_method Statistical test method: "wilcox.test" or "t.test"
#' @param show_pvalue Whether to show p-values on the plot
#' @param show_ns Whether to show "ns" for non-significant results
#' @param color_palette Color palette for the groups
#' @param facet_by Whether to facet by "cancer" or "none"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare TP53 expression in TP53 mutant vs wild-type samples
#' p <- vis_toil_Mut("TP53")
#'
#' # Compare BRCA1 expression in TP53 mutant vs wild-type
#' p <- vis_toil_Mut("TP53", target_gene = "BRCA1")
#'
#' # Use violin plot for specific cancers
#' p <- vis_toil_Mut(
#'   gene = "TP53",
#'   cancers = c("BRCA", "LUAD", "LUSC"),
#'   plot_type = "violin"
#' )
#' }
vis_toil_Mut <- function(gene,
                          target_gene = NULL,
                          data_type = c("mRNA", "protein", "cnv", "methylation"),
                          cancers = NULL,
                          plot_type = c("boxplot", "violin", "barplot", "dotplot"),
                          test_method = c("wilcox.test", "t.test"),
                          show_pvalue = TRUE,
                          show_ns = FALSE,
                          color_palette = c("#DF2020", "#2E5AAC"),
                          facet_by = c("cancer", "none")) {
  data_type <- match.arg(data_type)
  plot_type <- match.arg(plot_type)
  test_method <- match.arg(test_method)
  facet_by <- match.arg(facet_by)

  # Use gene as target if not specified
  if (is.null(target_gene)) {
    target_gene <- gene
  }

  # Get mutation status
  mut_status <- query_mutation(gene, source = "tcga")

  if (is.null(mut_status)) {
    stop("Failed to retrieve mutation data for gene: ", gene)
  }

  # Convert mutation status to binary
  mut_binary <- ifelse(mut_status == "0" | is.na(mut_status), "WT", "Mutant")

  # Get target gene expression/feature
  target_data <- switch(data_type,
    "mRNA" = query_gene_expression(target_gene, source = "tcga"),
    "protein" = query_protein(target_gene, source = "tcga"),
    "cnv" = query_cnv(target_gene, source = "tcga"),
    "methylation" = query_methylation(target_gene, source = "tcga"),
    stop("Unsupported data_type: ", data_type)
  )

  if (is.null(target_data)) {
    stop("Failed to retrieve ", data_type, " data for gene: ", target_gene)
  }

  # Get sample cancer types using barcode matching
  sample_info <- load_data("tcga_gtex")
  match_result <- match_samples(names(target_data), sample_info$Sample, "tcga", "tcga", match_by = "barcode")

  if (match_result$n_matched == 0) {
    stop("No matching samples found")
  }

  # Create plot data
  plot_data <- data.frame(
    Sample = match_result$common_ids,
    Value = as.numeric(target_data[match_result$idx1]),
    Mutation = mut_binary[match(match_result$common_ids, names(mut_binary))],
    Cancer = sample_info$Tissue[match_result$idx2],
    stringsAsFactors = FALSE
  )

  # Remove samples with missing mutation status
  plot_data <- plot_data[!is.na(plot_data$Mutation), ]

  # Filter by cancers
  if (!is.null(cancers)) {
    plot_data <- plot_data[plot_data$Cancer %in% cancers, ]
  }

  # Remove NA values
  plot_data <- plot_data[complete.cases(plot_data), ]

  if (nrow(plot_data) == 0) {
    stop("No valid data after filtering")
  }

  # Calculate statistics by cancer
  stats_by_cancer <- plot_data %>%
    dplyr::group_by(.data$Cancer) %>%
    dplyr::filter(length(unique(.data$Mutation)) == 2) %>%
    dplyr::summarise(
      test = list(broom::tidy(switch(test_method,
        "wilcox.test" = stats::wilcox.test(.data$Value ~ .data$Mutation),
        "t.test" = stats::t.test(.data$Value ~ .data$Mutation)
      ))),
      n_wt = sum(.data$Mutation == "WT"),
      n_mut = sum(.data$Mutation == "Mutant"),
      .groups = "drop"
    ) %>%
    tidyr::unnest(.data$test)

  # Determine y-position for p-value labels
  y_max <- max(plot_data$Value, na.rm = TRUE)
  y_min <- min(plot_data$Value, na.rm = TRUE)
  y_range <- y_max - y_min

  # Create base plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$Mutation, y = .data$Value, fill = .data$Mutation)
  )

  # Add geometry based on plot_type
  p <- switch(plot_type,
    boxplot = p + ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA),
    violin = p +
      ggplot2::geom_violin(alpha = 0.7, trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA),
    barplot = p +
      ggplot2::stat_summary(fun = mean, geom = "bar", alpha = 0.7) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.2),
    dotplot = p +
      ggplot2::geom_jitter(ggplot2::aes(color = Mutation), alpha = 0.5, width = 0.2) +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 3, shape = 18)
  )

  # Add raw points for boxplot/violin
  if (plot_type %in% c("boxplot", "violin")) {
    p <- p + ggplot2::geom_jitter(width = 0.2, alpha = 0.3, size = 0.5)
  }

  # Set colors
  p <- p +
    ggplot2::scale_fill_manual(values = setNames(color_palette, c("Mutant", "WT"))) +
    ggplot2::scale_color_manual(values = setNames(color_palette, c("Mutant", "WT")))

  # Add faceting
  if (facet_by == "cancer") {
    p <- p + ggplot2::facet_wrap(~Cancer, scales = "free_y")
  }

  # Add p-values
  if (show_pvalue && nrow(stats_by_cancer) > 0) {
    # Create p-value label
    stats_by_cancer$p_label <- ifelse(
      stats_by_cancer$p.value < 0.001,
      "***",
      ifelse(
        stats_by_cancer$p.value < 0.01,
        "**",
        ifelse(stats_by_cancer$p.value < 0.05, "*", ifelse(show_ns, "ns", ""))
      )
    )

    if (facet_by == "cancer") {
      # Add p-value as text on each facet
      p <- p + ggplot2::geom_text(
        data = stats_by_cancer[stats_by_cancer$p_label != "", ],
        ggplot2::aes(x = 1.5, y = y_max + y_range * 0.1, label = p_label),
        inherit.aes = FALSE,
        vjust = 0,
        size = 4
      )
    } else {
      # Overall comparison
      overall_test <- switch(test_method,
        "wilcox.test" = stats::wilcox.test(.data$Value ~ .data$Mutation, data = plot_data),
        "t.test" = stats::t.test(.data$Value ~ .data$Mutation, data = plot_data)
      )

      p_label <- ifelse(
        overall_test$p.value < 0.001,
        "***",
        ifelse(
          overall_test$p.value < 0.01,
          "**",
          ifelse(overall_test$p.value < 0.05, "*", ifelse(show_ns, "ns", ""))
        )
      )

      if (p_label != "") {
        p <- p + ggplot2::annotate(
          "text",
          x = 1.5,
          y = y_max + y_range * 0.1,
          label = p_label,
          vjust = 0,
          size = 5
        )
      }
    }
  }

  # Add labels
  y_label <- switch(data_type,
    "mRNA" = paste(target_gene, "Expression (TPM)"),
    "protein" = paste(target_gene, "Protein Level"),
    "cnv" = paste(target_gene, "Copy Number (GISTIC2)"),
    "methylation" = paste(target_gene, "Methylation (Beta)"),
    "Value"
  )

  p <- p + ggplot2::labs(
    title = paste(target_gene, "by", gene, "Mutation Status"),
    subtitle = ifelse(
      gene == target_gene,
      "Self-comparison (cis-effect)",
      "Trans-effect analysis"
    ),
    x = paste(gene, "Status"),
    y = y_label,
    fill = "Status"
  )

  # Apply theme
  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  # Return plot with statistics
  structure(
    p,
    stats = stats_by_cancer,
    class = c("vis_mut_plot", class(p))
  )
}

#' Visualize Mutation Frequency
#'
#' @description
#' Create a barplot showing mutation frequency across cancer types.
#'
#' @param genes Vector of gene symbols to analyze
#' @param cancers Vector of cancer types to include
#' @param min_freq Minimum mutation frequency to include (0-1)
#' @param plot_type Type of plot: "bar", "heatmap", or "oncoplot"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Show mutation frequency of key genes
#' p <- vis_mutation_frequency(c("TP53", "KRAS", "PIK3CA", "PTEN"))
#'
#' # Focus on specific cancers
#' p <- vis_mutation_frequency(
#'   genes = c("TP53", "BRCA1", "BRCA2"),
#'   cancers = c("BRCA", "OV", "UCEC")
#' )
#' }
vis_mutation_frequency <- function(genes,
                                    cancers = NULL,
                                    min_freq = 0,
                                    plot_type = c("bar", "heatmap", "oncoplot")) {
  plot_type <- match.arg(plot_type)

  # Get mutation data for all genes
  mut_list <- lapply(genes, function(g) {
    mut <- query_mutation(g, source = "tcga")
    if (is.null(mut)) return(NULL)
    data.frame(
      Sample = names(mut),
      Gene = g,
      Mutated = mut != "0" & !is.na(mut),
      stringsAsFactors = FALSE
    )
  })

  mut_data <- do.call(rbind, mut_list)

  if (is.null(mut_data) || nrow(mut_data) == 0) {
    stop("No mutation data retrieved")
  }

  # Get cancer types
  sample_info <- load_data("tcga_gtex")
  mut_data$Cancer <- sample_info$tissue[match(mut_data$Sample, sample_info$sample)]

  # Filter by cancers
  if (!is.null(cancers)) {
    mut_data <- mut_data[mut_data$Cancer %in% cancers, ]
  }

  # Calculate frequencies
  freq_data <- mut_data %>%
    dplyr::group_by(.data$Cancer, .data$Gene) %>%
    dplyr::summarise(
      Frequency = mean(.data$Mutated, na.rm = TRUE),
      N = sum(.data$Mutated, na.rm = TRUE),
      Total = dplyr::n(),
      .groups = "drop"
    )

  # Filter by minimum frequency
  if (min_freq > 0) {
    genes_keep <- freq_data %>%
      dplyr::group_by(.data$Gene) %>%
      dplyr::summarise(max_freq = max(.data$Frequency)) %>%
      dplyr::filter(.data$max_freq >= min_freq) %>%
      dplyr::pull(.data$Gene)

    freq_data <- freq_data[freq_data$Gene %in% genes_keep, ]
  }

  # Create plot
  switch(plot_type,
    bar = {
      ggplot2::ggplot(freq_data, ggplot2::aes(x = stats::reorder(.data$Gene, -.data$Frequency), y = .data$Frequency)) +
        ggplot2::geom_col(ggplot2::aes(fill = .data$Cancer)) +
        ggplot2::facet_wrap(~.data$Cancer, scales = "free_x") +
        ggplot2::labs(
          title = "Mutation Frequency by Cancer Type",
          x = "Gene",
          y = "Mutation Frequency"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    },
    heatmap = {
      ggplot2::ggplot(freq_data, ggplot2::aes(x = .data$Gene, y = .data$Cancer, fill = .data$Frequency)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, 1)) +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.2f", .data$Frequency)),
          color = "black", size = 3
        ) +
        ggplot2::labs(
          title = "Mutation Frequency Heatmap",
          x = "Gene",
          y = "Cancer Type"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          panel.grid = ggplot2::element_blank()
        )
    },
    oncoplot = {
      # Create binary mutation matrix
      mut_matrix <- mut_data %>%
        dplyr::select(Sample, Gene, Mutated) %>%
        tidyr::pivot_wider(
          names_from = Gene,
          values_from = Mutated,
          values_fill = FALSE
        )

      # This would need more complex implementation for a proper oncoplot
      # For now, return a simplified version
      message("Oncoplot implementation requires additional dependencies (ComplexHeatmap)")
      vis_mutation_frequency(genes, cancers, min_freq, "heatmap")
    }
  )
}

#' Visualize Co-occurrence of Mutations
#'
#' @description
#' Analyze and visualize co-occurrence or mutual exclusivity of mutations
#' in two genes.
#'
#' @param gene1 First gene symbol
#' @param gene2 Second gene symbol
#' @param cancers Vector of cancer types to include
#' @param plot_type Type of plot: "mosaic", "bar", or "table"
#'
#' @return A ggplot object or contingency table
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze TP53 and KRAS co-occurrence
#' result <- vis_mutation_cooccurrence("TP53", "KRAS")
#' print(result$plot)
#' print(result$table)
#' }
vis_mutation_cooccurrence <- function(gene1,
                                       gene2,
                                       cancers = NULL,
                                       plot_type = c("mosaic", "bar", "table")) {
  plot_type <- match.arg(plot_type)

  # Get mutation data
  mut1 <- query_mutation(gene1, source = "tcga")
  mut2 <- query_mutation(gene2, source = "tcga")

  # Create data frame
  samples <- intersect(names(mut1), names(mut2))

  cooc_data <- data.frame(
    Sample = samples,
    Gene1 = mut1[samples] != "0" & !is.na(mut1[samples]),
    Gene2 = mut2[samples] != "0" & !is.na(mut2[samples]),
    stringsAsFactors = FALSE
  )

  # Get cancer types
  sample_info <- load_data("tcga_gtex")
  cooc_data$Cancer <- sample_info$tissue[match(cooc_data$Sample, sample_info$sample)]

  # Filter by cancers
  if (!is.null(cancers)) {
    cooc_data <- cooc_data[cooc_data$Cancer %in% cancers, ]
  }

  # Create contingency table
  cont_table <- table(
    Mutation_Gene1 = cooc_data$Gene1,
    Mutation_Gene2 = cooc_data$Gene2
  )

  # Fisher's exact test
  fisher_test <- fisher.test(cont_table)

  # Create plot data
  plot_data <- as.data.frame(cont_table)
  colnames(plot_data) <- c("Gene1_Mut", "Gene2_Mut", "Count")

  # Create labels
  plot_data$Label <- paste0(
    ifelse(plot_data$Gene1_Mut == "TRUE", gene1, paste0(gene1, " WT")),
    " / ",
    ifelse(plot_data$Gene2_Mut == "TRUE", gene2, paste0(gene2, " WT"))
  )

  if (plot_type == "table") {
    return(list(
      table = cont_table,
      fisher_test = fisher_test,
      interpretation = ifelse(
        fisher_test$p.value < 0.05,
        ifelse(fisher_test$estimate > 1, "Co-occurring", "Mutually exclusive"),
        "No significant association"
      )
    ))
  }

  p <- switch(plot_type,
    mosaic = {
      # Simplified mosaic plot using bar
      ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Label, y = .data$Count, fill = .data$Label)) +
        ggplot2::geom_col() +
        ggplot2::labs(
          title = paste("Mutation Co-occurrence:", gene1, "vs", gene2),
          subtitle = sprintf(
            "Fisher's exact test: OR = %.2f, p = %.3e",
            fisher_test$estimate, fisher_test$p.value
          ),
          x = "Mutation Status",
          y = "Number of Samples"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )
    },
    bar = {
      ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Label, y = .data$Count, fill = .data$Label)) +
        ggplot2::geom_col() +
        ggplot2::geom_text(
          ggplot2::aes(label = .data$Count),
          vjust = -0.5
        ) +
        ggplot2::labs(
          title = paste("Mutation Co-occurrence:", gene1, "vs", gene2),
          subtitle = sprintf(
            "Fisher's exact test: OR = %.2f, p = %.3e",
            fisher_test$estimate, fisher_test$p.value
          ),
          x = "Mutation Status",
          y = "Number of Samples"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )
    }
  )

  structure(
    p,
    table = cont_table,
    fisher_test = fisher_test,
    class = c("vis_cooc_plot", class(p))
  )
}
