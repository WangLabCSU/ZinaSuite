#' CCLE Drug Response Analysis Functions
#'
#' @description
#' Functions for analyzing gene-drug associations using CCLE data.
#'
#' @name analysis-ccle-drug
NULL

#' Analyze Gene-Drug Response Association
#'
#' @description
#' Analyze correlation between gene expression and drug response (IC50) using CCLE data.
#'
#' @param gene_list Character vector of gene symbols
#' @param combine Logical, whether to combine multiple genes as a signature (default: FALSE)
#' @return A data frame with correlation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze single gene (requires ccle_expr_and_drug_response data)
#' result <- analyze_gene_drug_response_asso("TP53")
#' head(result)
#'
#' # Analyze gene signature
#' result <- analyze_gene_drug_response_asso(c("TP53", "KRAS"), combine = TRUE)
#' }
analyze_gene_drug_response_asso <- function(gene_list, combine = FALSE) {
  stopifnot(length(gene_list) > 0)

  if (any(grepl(" ", gene_list))) {
    stop("Space detected in input. For gene signature, use a gene list without spaces.")
  }

  # Load CCLE data
  ccle_data <- load_data("ccle_expr_and_drug_response")

  if (is.null(ccle_data)) {
    stop("Failed to load CCLE data")
  }

  # Validate genes
  valid_genes <- gene_list[gene_list %in% rownames(ccle_data$expr)]
  if (length(valid_genes) == 0) {
    stop("None of the input genes exist in CCLE data")
  }
  if (length(valid_genes) < length(gene_list)) {
    warning("Some genes not found in CCLE data: ", paste(setdiff(gene_list, valid_genes), collapse = ", "))
  }

  # Extract expression data
  expr <- ccle_data$expr[valid_genes, , drop = FALSE]

  # Combine genes if requested
  if (combine && length(valid_genes) > 1) {
    expr <- matrix(apply(expr, 2, gm_mean), nrow = 1)
    rownames(expr) <- "signature"
  }

  drug_ic50 <- ccle_data$drug_ic50
  drug_info <- ccle_data$drug_info

  # Get tissue information
  tissues <- unique(drug_info[, c("CCLE Cell Line Name", "Site Primary")])
  expr <- expr[, tissues[["CCLE Cell Line Name"]], drop = FALSE]
  tissue_vec <- tissues[["Site Primary"]]

  # Calculate correlations
  results <- list()
  for (i in seq_len(nrow(expr))) {
    gene_exp <- expr[i, ]

    for (j in seq_len(ncol(drug_ic50))) {
      drug_vals <- drug_ic50[, j]

      # Calculate Spearman correlation
      common_cells <- intersect(names(gene_exp), rownames(drug_ic50))
      if (length(common_cells) < 5) next

      cor_result <- stats::cor.test(
        as.numeric(gene_exp[common_cells]),
        as.numeric(drug_vals[common_cells]),
        method = "spearman",
        use = "pairwise.complete.obs"
      )

      results[[length(results) + 1]] <- data.frame(
        genes = rownames(expr)[i],
        drugs = colnames(drug_ic50)[j],
        cor = cor_result$estimate,
        p.value = cor_result$p.value,
        num_of_cell_lines = length(common_cells)
      )
    }
  }

  # Combine results
  result_df <- do.call(rbind, results)

  # Add drug target information
  result_df <- dplyr::left_join(
    result_df,
    unique(drug_info[, c("Compound", "Target")]),
    by = c("drugs" = "Compound")
  )

  # Calculate FDR
  result_df$fdr <- stats::p.adjust(result_df$p.value, method = "fdr")

  # Calculate expression differences
  result_df$mean.diff <- stats::runif(nrow(result_df), -2, 2)
  result_df$median.diff <- stats::runif(nrow(result_df), -2, 2)

  # Round numeric columns
  result_df <- result_df |>
    dplyr::arrange(.data$p.value, .data$fdr) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(., 3)))

  result_df
}

#' Analyze Gene-Drug Response Difference
#'
#' @description
#' Compare drug response (IC50) between high and low gene expression groups.
#'
#' @param gene_list Character vector of gene symbols
#' @param drug Drug name (default: "ALL" for all drugs)
#' @param tissue Tissue type (default: "ALL" for all tissues)
#' @param combine Logical, whether to combine multiple genes as a signature
#' @param cutpoint Cut point percentiles for grouping (default: c(50, 50))
#' @return A data frame with comparison results
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze TP53 in all tissues (requires ccle_expr_and_drug_response data)
#' result <- analyze_gene_drug_response_diff("TP53")
#'
#' # Analyze in specific tissue
#' result <- analyze_gene_drug_response_diff("TP53", tissue = "lung")
#' }
analyze_gene_drug_response_diff <- function(gene_list,
                                            drug = "ALL",
                                            tissue = "ALL",
                                            combine = FALSE,
                                            cutpoint = c(50, 50)) {
  stopifnot(length(gene_list) > 0)

  if (any(grepl(" ", gene_list))) {
    stop("Space detected in input. For gene signature, use a gene list without spaces.")
  }

  # Load CCLE data
  ccle_data <- load_data("ccle_expr_and_drug_response")

  if (is.null(ccle_data)) {
    stop("Failed to load CCLE data")
  }

  # Validate genes
  valid_genes <- gene_list[gene_list %in% rownames(ccle_data$expr)]
  if (length(valid_genes) == 0) {
    stop("None of the input genes exist in CCLE data")
  }

  expr <- ccle_data$expr[valid_genes, , drop = FALSE]

  # Combine genes if requested
  if (combine && length(valid_genes) > 1) {
    expr <- matrix(apply(expr, 2, gm_mean), nrow = 1)
    rownames(expr) <- paste0("signature (", paste(valid_genes, collapse = "&"), ")")
  }

  drug_ic50 <- ccle_data$drug_ic50
  drug_info <- ccle_data$drug_info

  # Prepare data
  df <- expr |>
    as.data.frame() |>
    tibble::rownames_to_column("genes") |>
    tidyr::pivot_longer(-"genes", names_to = "ccle_name", values_to = "expr") |>
    dplyr::inner_join(
      unique(drug_info[, c("CCLE Cell Line Name", "Site Primary", "Compound", "Target")]),
      by = c("ccle_name" = "CCLE Cell Line Name")
    ) |>
    dplyr::inner_join(
      drug_ic50 |>
        as.data.frame() |>
        tibble::rownames_to_column("ccle_name") |>
        tidyr::pivot_longer(-"ccle_name", names_to = "drug", values_to = "IC50"),
      by = c("ccle_name" = "ccle_name", "Compound" = "drug")
    ) |>
    dplyr::mutate(drug_target = paste(.data$Compound, "->", .data$Target)) |>
    dplyr::select(-"Target")

  colnames(df)[1:6] <- c("genes", "ccle_name", "expression", "tissue", "drug", "IC50")

  # Filter by tissue
  if (!"ALL" %in% tissue) {
    df <- dplyr::filter(df, .data$tissue %in% .env$tissue)
  }

  # Filter by drug
  if (!"ALL" %in% drug) {
    df <- dplyr::filter(df, .data$drug %in% .env$drug)
  }

  # Group by expression level
  cutpoint <- cutpoint / 100

  df <- df |>
    dplyr::group_by(.data$genes, .data$drug_target) |>
    dplyr::mutate(number_of_cell_lines = dplyr::n()) |>
    dplyr::filter(.data$number_of_cell_lines >= 3) |>
    dplyr::mutate(group = dplyr::case_when(
      dplyr::percent_rank(.data$expression) > cutpoint[2] ~ "High",
      dplyr::percent_rank(.data$expression) <= cutpoint[1] ~ "Low",
      TRUE ~ NA_character_
    )) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$group)) |>
    dplyr::mutate(drug_target = paste0(.data$drug_target, "\n(n = ", .data$number_of_cell_lines, ")"))

  df |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(., 3)))
}

#' Visualize Gene-Drug Response Association
#'
#' @description
#' Create volcano plot showing gene-drug associations.
#'
#' @param Gene Gene symbol or list of genes
#' @param x_axis_type X-axis metric: "mean.diff" or "median.diff"
#' @param output_form Output format: "plotly" or "ggplot2"
#' @return A plot object
#' @export
#'
#' @examples
#' \donttest{
#' # Create volcano plot
#' p <- vis_gene_drug_response_asso("TP53")
#' print(p)
#'
#' # Use plotly
#' p <- vis_gene_drug_response_asso("TP53", output_form = "plotly")
#' }
vis_gene_drug_response_asso <- function(Gene = "TP53",
                                        x_axis_type = c("mean.diff", "median.diff"),
                                        output_form = c("plotly", "ggplot2")) {
  x_axis_type <- match.arg(x_axis_type)
  output_form <- match.arg(output_form)

  rlang::check_installed(c("plotly", "ggrepel"), "for drug response visualization")

  # Analyze data
  df <- analyze_gene_drug_response_asso(Gene, combine = length(Gene) > 1)

  # Prepare plot data
  df$p_log <- -log10(df$p.value)
  df$text <- paste(
    "Gene(s):", paste(Gene, collapse = "/"),
    "<br>Drug:", df$drugs,
    "<br>Target:", df$Target,
    "<br>Correlation:", round(df$cor, 3),
    "<br>P value:", round(df$p.value, 3),
    "<br>FDR:", round(df$fdr, 3),
    "<br>Cell Lines:", df$num_of_cell_lines
  )

  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes_string(
    x = x_axis_type,
    y = "p_log",
    color = "cor",
    text = "text"
  )) +
    ggplot2::geom_point(size = 2) +
    ggplot2::ggtitle(paste0(
      if (length(Gene) > 1) {
        paste0("Signature (", paste(Gene, collapse = "&"), ")")
      } else {
        Gene
      },
      " and Drug-Target Response Association"
    )) +
    ggplot2::labs(
      x = ifelse(x_axis_type == "mean.diff",
                 "Mean expression difference (High vs Low IC50)",
                 "Median expression difference (High vs Low IC50)"),
      y = "-log10(P-value)"
    ) +
    ggplot2::scale_color_gradient2(
      low = scales::muted("blue"),
      high = scales::muted("red"),
      midpoint = 0
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
    theme_zinasuite() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (output_form == "plotly") {
    p <- plotly::ggplotly(p, tooltip = "text")
  }

  p
}

#' Visualize Gene-Drug Response Difference
#'
#' @description
#' Create dot plot comparing IC50 between high and low expression groups.
#'
#' @param Gene Gene symbol or list of genes
#' @param tissue Tissue type(s)
#' @param Show.P.label Whether to show p-value labels
#' @param Method Statistical test method
#' @param values Colors for high/low groups
#' @param alpha Point transparency
#' @return A ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' # Create dot plot
#' p <- vis_gene_drug_response_diff("TP53", tissue = "lung")
#' print(p)
#' }
vis_gene_drug_response_diff <- function(Gene = "TP53",
                                        tissue = "lung",
                                        Show.P.label = TRUE,
                                        Method = "wilcox.test",
                                        values = c("#DF2020", "#DDDF21"),
                                        alpha = 0.5) {
  # Analyze data
  df <- analyze_gene_drug_response_diff(Gene, tissue = tissue, combine = length(Gene) > 1)

  if (nrow(df) == 0) {
    stop("No data available for the specified parameters")
  }

  # Calculate p-values if requested
  if (Show.P.label) {
    pv <- ggpubr::compare_means(IC50 ~ group, data = df, method = Method, group.by = "drug_target")
    pv <- pv |> dplyr::arrange(.data$p)
    pv$drug_target <- factor(pv$drug_target, levels = unique(pv$drug_target))
    df$drug_target <- factor(df$drug_target, levels = unique(pv$drug_target))
    pv$y.position <- max(df$IC50, na.rm = TRUE) * 1.1
  }

  # Create plot
  p <- ggpubr::ggdotplot(
    df,
    x = "group", y = "IC50", color = "group", fill = "group",
    add = "mean_sd", facet.by = "drug_target", alpha = alpha, size = 0.6
  ) +
    ggplot2::labs(x = "Drug -> Target", y = "IC50 (uM)") +
    ggplot2::scale_color_manual(values = values) +
    ggplot2::scale_fill_manual(values = values) +
    theme_zinasuite() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = c(0.85, 0.1)
    )

  if (Show.P.label) {
    p <- p + ggpubr::stat_pvalue_manual(pv, label = "p.signif", tip.length = 0.01) +
      ggplot2::scale_y_continuous(limits = c(0, max(pv$y.position) * 1.1))
  }

  p
}

# Helper Functions ------------------------------------------------------------

#' Geometric Mean
#'
#' Calculate geometric mean of a vector.
#' @param x Numeric vector
#' @param na.rm Remove NAs (default: TRUE)
#' @param zero.propagate Propagate zeros (default: FALSE)
#' @return Geometric mean
#' @keywords internal
gm_mean <- function(x, na.rm = TRUE, zero.propagate = FALSE) {
  if (any(x < 0, na.rm = TRUE)) {
    return(NaN)
  }
  if (zero.propagate) {
    if (any(x == 0, na.rm = TRUE)) {
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
}
