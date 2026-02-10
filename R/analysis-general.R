#' General Analysis Functions
#'
#' @description
#' This file contains general analysis functions extracted from UCSCXenaShiny
#' for use in ZinaSuite Shiny applications.
#'
#' @name analysis-general
NULL

#' Scatter Correlation Analysis
#'
#' @description
#' Perform scatter correlation analysis between two molecular identifiers.
#' This function extracts the core logic from UCSCXenaShiny's general analysis module.
#'
#' @param dataset1 Character string specifying the first dataset name
#' @param id1 Character string specifying the first molecule identifier
#' @param dataset2 Character string specifying the second dataset name
#' @param id2 Character string specifying the second molecule identifier
#' @param samples Optional character vector of sample names to filter
#' @param use_ggstats Logical indicating whether to use ggstatsplot (default: FALSE)
#' @param use_simple_axis_label Logical indicating whether to use simple axis labels (default: TRUE)
#' @param line_color Character string for regression line color (default: "blue")
#' @param alpha Numeric value for point transparency (default: 0.5)
#' @param ... Additional parameters passed to plotting functions
#'
#' @return A ggplot object showing the scatter correlation plot
#' @export
#'
#' @examples
#' \donttest{
#' # Basic correlation plot
#' p <- analyze_scatter_correlation(
#'   dataset1 = "TcgaTargetGtex_rsem_gene_tpm",
#'   id1 = "TP53",
#'   dataset2 = "TcgaTargetGtex_rsem_gene_tpm",
#'   id2 = "KRAS"
#' )
#' print(p)
#'
#' # With specific samples
#' samples <- c("TCGA-D5-5538-01", "TCGA-VM-A8C8-01")
#' p <- analyze_scatter_correlation(
#'   dataset1 = "TcgaTargetGtex_rsem_gene_tpm",
#'   id1 = "TP53",
#'   dataset2 = "TcgaTargetGtex_rsem_gene_tpm",
#'   id2 = "KRAS",
#'   samples = samples
#' )
#' }
analyze_scatter_correlation <- function(dataset1, id1, dataset2, id2,
                                        samples = NULL,
                                        use_ggstats = FALSE,
                                        use_simple_axis_label = TRUE,
                                        line_color = "blue",
                                        alpha = 0.5, ...) {
  # Validate inputs
  stopifnot(length(id1) == 1, length(id2) == 1)
  stopifnot(is.character(dataset1), is.character(dataset2))

  # Query molecule values
  id1_value <- query_molecule_value(dataset1, id1)
  id2_value <- query_molecule_value(dataset2, id2)

  if (is.null(id1_value) || is.null(id2_value)) {
    stop("Failed to retrieve data for one or both identifiers")
  }

  # Create data frame
  df <- dplyr::inner_join(
    dplyr::tibble(
      sample = names(id1_value),
      X = as.numeric(id1_value)
    ),
    dplyr::tibble(
      sample = names(id2_value),
      Y = as.numeric(id2_value)
    ),
    by = "sample"
  )

  # Filter samples if specified
  if (!is.null(samples)) {
    df <- dplyr::filter(df, .data$sample %in% samples)
  }

  if (nrow(df) < 3) {
    stop("Insufficient samples for correlation analysis (need at least 3)")
  }

  # Build axis labels
  xlab <- if (use_simple_axis_label) {
    id1
  } else {
    paste0(id1, "(", attr(id1_value, "label"), ")")
  }

  ylab <- if (use_simple_axis_label) {
    id2
  } else {
    paste0(id2, "(", attr(id2_value, "label"), ")")
  }

  # Create plot
  if (!use_ggstats) {
    rlang::check_installed("ggpubr", "for scatter correlation plots")

    p <- ggpubr::ggscatter(
      data = df,
      x = "X",
      y = "Y",
      xlab = xlab,
      ylab = ylab,
      alpha = alpha,
      add = "reg.line",
      add.params = list(color = line_color, fill = "lightgray"),
      cor.coef = TRUE,
      ...
    )
  } else {
    rlang::check_installed("ggstatsplot", "for enhanced scatter plots")

    p <- ggstatsplot::ggscatterstats(
      data = df,
      x = "X",
      y = "Y",
      xlab = xlab,
      ylab = ylab,
      ...
    )
  }

  p
}

#' Matrix Correlation Analysis
#'
#' @description
#' Perform correlation analysis for multiple molecular identifiers.
#' Creates a correlation matrix visualization.
#'
#' @param dataset Character string specifying the dataset name
#' @param ids Character vector of molecule identifiers
#' @param samples Optional character vector of sample names to filter
#' @param matrix.type Character string: "full", "upper", or "lower" (default: "full")
#' @param type Character string: "parametric", "nonparametric", "robust", or "bayes" (default: "parametric")
#' @param partial Logical indicating whether to compute partial correlation (default: FALSE)
#' @param sig.level Numeric significance level (default: 0.05)
#' @param p.adjust.method Character string for p-value adjustment method (default: "fdr")
#' @param color_low Character string for low correlation color (default: "#E69F00")
#' @param color_high Character string for high correlation color (default: "#009E73")
#' @param ... Additional parameters passed to ggcorrmat
#'
#' @return A ggplot object showing the correlation matrix
#' @export
#'
#' @examples
#' \donttest{
#' # Correlation matrix for multiple genes
#' p <- analyze_matrix_correlation(
#'   dataset = "TcgaTargetGtex_rsem_gene_tpm",
#'   ids = c("TP53", "KRAS", "PTEN", "BRCA1")
#' )
#' print(p)
#'
#' # Upper triangle only
#' p <- analyze_matrix_correlation(
#'   dataset = "TcgaTargetGtex_rsem_gene_tpm",
#'   ids = c("TP53", "KRAS", "PTEN"),
#'   matrix.type = "upper"
#' )
#' }
analyze_matrix_correlation <- function(dataset, ids,
                                       samples = NULL,
                                       matrix.type = c("full", "upper", "lower"),
                                       type = c("parametric", "nonparametric", "robust", "bayes"),
                                       partial = FALSE,
                                       sig.level = 0.05,
                                       p.adjust.method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                       color_low = "#E69F00",
                                       color_high = "#009E73",
                                       ...) {
  # Validate inputs
  stopifnot(length(ids) >= 2)

  matrix.type <- match.arg(matrix.type)
  type <- match.arg(type)
  p.adjust.method <- match.arg(p.adjust.method)

  # Query data for all identifiers
  df <- purrr::map(ids, function(x) {
    message("Querying data of identifier ", x, " from dataset: ", dataset)
    data <- query_molecule_value(dataset, x)

    if (is.null(data)) {
      warning("Failed to retrieve data for identifier: ", x)
      return(NULL)
    }

    dplyr::tibble(
      sample = names(data),
      y = as.numeric(data)
    ) |>
      dplyr::rename(!!x := .data$y)
  }) |>
    purrr::compact() |>
    purrr::reduce(dplyr::full_join, by = "sample")

  if (ncol(df) < 3) {  # sample column + at least 2 genes
    stop("At least 2 valid identifiers are required for correlation matrix")
  }

  # Filter samples if specified
  if (!is.null(samples)) {
    df <- dplyr::filter(df, .data$sample %in% samples)
  }

  # Check for sufficient data
  if (nrow(df) < 5) {
    stop("Insufficient samples after filtering (need at least 5)")
  }

  # Create correlation plot
  rlang::check_installed("ggstatsplot", "for correlation matrix visualization")

  colors <- c(color_low, "white", color_high)

  p <- ggstatsplot::ggcorrmat(
    data = df,
    matrix.type = matrix.type,
    type = type,
    partial = partial,
    sig.level = sig.level,
    p.adjust.method = p.adjust.method,
    colors = colors,
    ...
  )

  p
}

#' Group Comparison Analysis
#'
#' @description
#' Perform group comparison analysis for a molecular identifier.
#' Supports both between-group and within-group comparisons.
#'
#' @param dataset Optional character string specifying the dataset name
#' @param id Optional character string specifying the molecule identifier
#' @param grp_df Data frame with grouping information:
#'   - If dataset and id provided: 2-3 columns (sample, group, optional facet)
#'   - If no dataset/id: 3-4 columns (sample, value, group, optional facet)
#' @param samples Optional character vector of sample names to filter
#' @param fun_type Character string: "betweenstats" or "withinstats" (default: "betweenstats")
#' @param type Character string: "parametric", "nonparametric", "robust", or "bayes" (default: "parametric")
#' @param pairwise.comparisons Logical indicating whether to show pairwise comparisons (default: TRUE)
#' @param p.adjust.method Character string for p-value adjustment method (default: "fdr")
#' @param ggtheme ggplot2 theme object (default: cowplot::theme_cowplot())
#' @param ... Additional parameters passed to ggbetweenstats or ggwithinstats
#'
#' @return A ggplot object showing the group comparison
#' @export
#'
#' @examples
#' \donttest{
#' # Create sample grouping data
#' grp_df <- data.frame(
#'   sample = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'   group = c("A", "A", "B", "B", "C", "C"),
#'   value = c(1.2, 1.5, 2.3, 2.1, 3.1, 2.9)
#' )
#'
#' # Direct comparison with value column
#' p <- analyze_group_comparison(grp_df = grp_df)
#' print(p)
#'
#' # Using dataset and identifier
#' # p <- analyze_group_comparison(
#' #   dataset = "TcgaTargetGtex_rsem_gene_tpm",
#' #   id = "TP53",
#' #   grp_df = grp_df[, c("sample", "group")]
#' # )
#' }
analyze_group_comparison <- function(dataset = NULL,
                                     id = NULL,
                                     grp_df,
                                     samples = NULL,
                                     fun_type = c("betweenstats", "withinstats"),
                                     type = c("parametric", "nonparametric", "robust", "bayes"),
                                     pairwise.comparisons = TRUE,
                                     p.adjust.method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                     ggtheme = ggplot2::theme_bw(),
                                     ...) {
  # Validate inputs
  stopifnot(ncol(grp_df) > 1)

  fun_type <- match.arg(fun_type)
  type <- match.arg(type)
  p.adjust.method <- match.arg(p.adjust.method)

  colnames(grp_df)[1] <- "sample"

  # Determine data source
  if (!is.null(dataset) && !is.null(id)) {
    message("Querying data of identifier ", id, " from dataset ", dataset)
    id_value <- query_molecule_value(dataset, id)

    if (is.null(id_value)) {
      stop("Failed to retrieve data for identifier: ", id)
    }

    df <- dplyr::tibble(
      sample = names(id_value),
      X = as.numeric(id_value)
    )
    colnames(df)[2] <- id

    df <- dplyr::inner_join(df, grp_df, by = "sample")
    do_grp <- ncol(grp_df) >= 3
  } else {
    message("Directly using 'grp_df' for comparison analysis")
    df <- grp_df
    do_grp <- ncol(grp_df) >= 4
  }

  # Filter samples if specified
  if (!is.null(samples)) {
    df <- dplyr::filter(df, .data$sample %in% samples)
  }

  if (nrow(df) < 3) {
    stop("Insufficient samples for group comparison (need at least 3)")
  }

  # Create comparison plot
  rlang::check_installed("ggstatsplot", "for group comparison visualization")

  if (do_grp) {
    # Grouped comparison
    fun <- if (fun_type == "betweenstats") {
      ggstatsplot::grouped_ggbetweenstats
    } else {
      ggstatsplot::grouped_ggwithinstats
    }

    p <- fun(
      data = df,
      x = !!rlang::sym(colnames(df)[3]),
      y = !!rlang::sym(colnames(df)[2]),
      grouping.var = !!rlang::sym(colnames(df)[4]),
      type = type,
      pairwise.comparisons = pairwise.comparisons,
      p.adjust.method = p.adjust.method,
      ggtheme = ggtheme,
      ...
    )
  } else {
    # Simple comparison
    fun <- if (fun_type == "betweenstats") {
      ggstatsplot::ggbetweenstats
    } else {
      ggstatsplot::ggwithinstats
    }

    p <- fun(
      data = df,
      x = !!rlang::sym(colnames(df)[3]),
      y = !!rlang::sym(colnames(df)[2]),
      type = type,
      pairwise.comparisons = pairwise.comparisons,
      p.adjust.method = p.adjust.method,
      ggtheme = ggtheme,
      ...
    )
  }

  p
}

#' Survival Analysis
#'
#' @description
#' Perform survival analysis for a molecular identifier.
#' Supports automatic or custom cutoff for grouping samples.
#'
#' @param dataset Optional character string specifying the dataset name
#' @param id Optional character string specifying the molecule identifier
#' @param surv_df Data frame with survival information:
#'   - If dataset and id provided: 3 columns (sample, time, status)
#'   - If no dataset/id: 4 columns (sample, value, time, status)
#' @param samples Optional character vector of sample names to filter
#' @param cutoff_mode Character string: "Auto", "Custom", or "None" (default: "Auto")
#' @param cutpoint Numeric vector of cutpoints for Custom mode (default: c(50, 50))
#' @param palette Character string for color palette (default: "aaas")
#' @param ... Additional parameters passed to survival plotting functions
#'
#' @return A ggplot object showing the survival analysis
#' @export
#'
#' @examples
#' \donttest{
#' # Create sample survival data
#' surv_df <- data.frame(
#'   sample = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'   value = c(1.2, 2.3, 3.1, 1.5, 2.1, 2.9),
#'   time = c(365, 200, 500, 180, 400, 300),
#'   status = c(1, 0, 1, 1, 0, 1)
#' )
#'
#' # Direct survival analysis
#' p <- analyze_survival(surv_df = surv_df)
#' print(p)
#'
#' # With custom cutpoints
#' p <- analyze_survival(
#'   surv_df = surv_df,
#'   cutoff_mode = "Custom",
#'   cutpoint = c(25, 75)
#' )
#' }
analyze_survival <- function(dataset = NULL,
                             id = NULL,
                             surv_df,
                             samples = NULL,
                             cutoff_mode = c("Auto", "Custom", "None"),
                             cutpoint = c(50, 50),
                             palette = "aaas",
                             ...) {
  cutoff_mode <- match.arg(cutoff_mode)

  # Determine data source
  if (!is.null(dataset) && !is.null(id)) {
    message("Querying data of identifier ", id, " from dataset ", dataset)
    id_value <- query_molecule_value(dataset, id)

    if (is.null(id_value)) {
      stop("Failed to retrieve data for identifier: ", id)
    }

    df <- dplyr::tibble(
      sample = names(id_value),
      value = as.numeric(id_value)
    )

    if (ncol(surv_df) == 3) {
      colnames(surv_df) <- c("sample", "time", "status")
    } else {
      stop("When using dataset and id, surv_df must have 3 columns: sample, time, status")
    }

    df <- dplyr::inner_join(df, surv_df, by = "sample")
  } else {
    message("Directly using 'surv_df' for survival analysis")
    df <- surv_df

    if (ncol(df) == 4) {
      colnames(df) <- c("sample", "value", "time", "status")
    } else {
      stop("When not using dataset/id, surv_df must have 4 columns: sample, value, time, status")
    }
  }

  # Filter samples if specified
  if (!is.null(samples)) {
    df <- dplyr::filter(df, .data$sample %in% samples)
  }

  if (nrow(df) < 10) {
    stop("Insufficient samples for survival analysis (need at least 10)")
  }

  # Create survival plot
  rlang::check_installed(c("survival", "survminer"), "for survival analysis")

  if (cutoff_mode != "None") {
    p <- sur_plot(df, cutoff_mode, cutpoint, palette = palette, ...)
  } else {
    colnames(df)[2] <- "group"
    p <- p_survplot(df, palette = palette, ...)
  }

  p
}

#' Dimension Distribution Analysis
#'
#' @description
#' Perform dimension reduction and distribution analysis.
#' Supports PCA, t-SNE, and UMAP methods.
#'
#' @param data Matrix or data frame with features as rows and samples as columns
#' @param method Character string: "pca", "tsne", or "umap" (default: "pca")
#' @param n_components Integer number of components to compute (default: 2)
#' @param color_by Optional vector for coloring points
#' @param shape_by Optional vector for shaping points
#' @param ... Additional parameters passed to dimension reduction functions
#'
#' @return A list containing:
#'   - plot: ggplot object
#'   - coordinates: data frame with dimension coordinates
#'   - model: fitted model object
#' @export
#'
#' @examples
#' \donttest{
#' # Create sample data
#' set.seed(123)
#' data <- matrix(rnorm(1000), nrow = 100)
#' colnames(data) <- paste0("S", 1:10)
#'
#' # PCA analysis
#' result <- analyze_dimension_distribution(data, method = "pca")
#' print(result$plot)
#'
#' # With coloring
#' groups <- rep(c("A", "B"), each = 5)
#' result <- analyze_dimension_distribution(data, method = "pca", color_by = groups)
#' }
analyze_dimension_distribution <- function(data,
                                           method = c("pca", "tsne", "umap"),
                                           n_components = 2,
                                           color_by = NULL,
                                           shape_by = NULL,
                                           ...) {
  method <- match.arg(method)

  # Validate data
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("data must be a matrix or data frame")
  }

  # Transpose if needed (samples should be rows)
  if (ncol(data) < nrow(data)) {
    data <- t(data)
  }

  # Remove rows/columns with all NA
  data <- data[apply(data, 1, function(x) !all(is.na(x))), ]
  data <- data[, apply(data, 2, function(x) !all(is.na(x)))]

  # Impute remaining NAs
  if (any(is.na(data))) {
    data <- apply(data, 2, function(x) {
      x[is.na(x)] <- median(x, na.rm = TRUE)
      x
    })
  }

  result <- switch(method,
    pca = {
      # PCA analysis
      pca_result <- stats::prcomp(data, scale. = TRUE, ...)

      coords <- as.data.frame(pca_result$x[, 1:n_components])
      colnames(coords) <- paste0("PC", 1:n_components)

      list(
        coordinates = coords,
        model = pca_result,
        variance_explained = summary(pca_result)$importance[2, 1:n_components]
      )
    },

    tsne = {
      # t-SNE analysis
      rlang::check_installed("Rtsne", "for t-SNE dimension reduction")

      tsne_result <- Rtsne::Rtsne(
        data,
        dims = n_components,
        perplexity = min(30, nrow(data) - 1),
        verbose = FALSE,
        ...
      )

      coords <- as.data.frame(tsne_result$Y)
      colnames(coords) <- paste0("tSNE", 1:n_components)

      list(
        coordinates = coords,
        model = tsne_result
      )
    },

    umap = {
      # UMAP analysis
      rlang::check_installed("umap", "for UMAP dimension reduction")

      umap_result <- umap::umap(data, n_components = n_components, ...)

      coords <- as.data.frame(umap_result$layout)
      colnames(coords) <- paste0("UMAP", 1:n_components)

      list(
        coordinates = coords,
        model = umap_result
      )
    }
  )

  # Add metadata for plotting
  result$coordinates$sample <- rownames(data)

  if (!is.null(color_by)) {
    result$coordinates$color_group <- color_by
  }

  if (!is.null(shape_by)) {
    result$coordinates$shape_group <- shape_by
  }

  # Create plot
  x_col <- colnames(result$coordinates)[1]
  y_col <- colnames(result$coordinates)[2]

  p <- ggplot2::ggplot(
    result$coordinates,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = if (!is.null(color_by)) .data$color_group else NULL,
      shape = if (!is.null(shape_by)) .data$shape_group else NULL
    )
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::labs(
      title = paste(toupper(method), "Dimension Reduction"),
      x = x_col,
      y = y_col
    ) +
    ggplot2::theme_minimal()

  # Add variance explained to PCA plot
  if (method == "pca") {
    ve <- result$variance_explained
    p <- p + ggplot2::labs(
      x = sprintf("%s (%.1f%%)", x_col, ve[1] * 100),
      y = sprintf("%s (%.1f%%)", y_col, ve[2] * 100)
    )
  }

  result$plot <- p

  result
}

#' Helper function for survival plotting with cutoff
#' @keywords internal
sur_plot <- function(df, cutoff_mode, cutpoint, palette, ...) {
  # Calculate groups based on cutoff mode
  if (cutoff_mode == "Auto") {
    # Use median as cutoff
    median_val <- median(df$value, na.rm = TRUE)
    df$group <- ifelse(df$value >= median_val, "High", "Low")
  } else if (cutoff_mode == "Custom") {
    # Use custom percentiles
    q <- quantile(df$value, probs = cutpoint / 100, na.rm = TRUE)
    df$group <- dplyr::case_when(
      df$value <= q[1] ~ "Low",
      df$value >= q[2] ~ "High",
      TRUE ~ "Medium"
    )
    # Remove medium group for binary comparison
    df <- df[df$group != "Medium", ]
  }

  p_survplot(df, palette, ...)
}

#' Helper function for survival plotting
#' @keywords internal
p_survplot <- function(df, palette, ...) {
  # Create survival object
  surv_obj <- survival::Surv(df$time, df$status)

  # Fit survival curve
  fit <- survival::survfit(surv_obj ~ group, data = df)

  # Create plot
  p <- survminer::ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    palette = palette,
    ...
  )

  p
}
