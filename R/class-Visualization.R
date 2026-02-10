#' Visualization R6 Class
#'
#' @description
#' R6 class for creating various bioinformatics visualizations using ggplot2.
#'
#' @export
#' @examples
#' \dontrun{
#' # Create visualization engine
#' viz <- Visualization$new()
#'
#' # Create distribution plot
#' tp53_expr <- query_gene_expression("TP53")
#' p <- viz$plot_distribution(tp53_expr, type = "boxplot")
#' print(p)
#'
#' # Create correlation plot
#' brca1_expr <- query_gene_expression("BRCA1")
#' p <- viz$plot_correlation(tp53_expr, brca1_expr, type = "scatter")
#' print(p)
#' }
Visualization <- R6::R6Class(
  "Visualization",

  public = list(
    #' @description
    #' Initialize a new Visualization instance
    #' @param theme Default theme to use
    initialize = function(theme = "default") {
      private$theme <- theme
    },

    #' @description
    #' Create distribution plot
    #' @param data Numeric vector or data frame
    #' @param type Plot type: "boxplot", "violin", "histogram", "density"
    #' @param group Optional grouping variable
    #' @param title Plot title
    #' @param xlab X-axis label
    #' @param ylab Y-axis label
    #' @return ggplot object
    plot_distribution = function(data,
                                  type = c("boxplot", "violin", "histogram", "density"),
                                  group = NULL,
                                  title = NULL,
                                  xlab = "Value",
                                  ylab = NULL) {
      type <- match.arg(type)

      # Prepare data
      if (is.vector(data)) {
        df <- data.frame(value = data)
        if (!is.null(group)) {
          df$group <- group
        }
      } else {
        df <- as.data.frame(data)
      }

      # Create base plot
      p <- ggplot2::ggplot(df)

      # Add geometry based on type
      p <- switch(type,
        "boxplot" = {
          if (is.null(group)) {
            p + ggplot2::geom_boxplot(ggplot2::aes(y = value), fill = "steelblue", alpha = 0.7)
          } else {
            p + ggplot2::geom_boxplot(ggplot2::aes(x = group, y = value, fill = group), alpha = 0.7)
          }
        },
        "violin" = {
          if (is.null(group)) {
            p + ggplot2::geom_violin(ggplot2::aes(x = "", y = value), fill = "steelblue", alpha = 0.7)
          } else {
            p + ggplot2::geom_violin(ggplot2::aes(x = group, y = value, fill = group), alpha = 0.7)
          }
        },
        "histogram" = {
          p + ggplot2::geom_histogram(ggplot2::aes(x = value), bins = 30, fill = "steelblue", alpha = 0.7)
        },
        "density" = {
          if (is.null(group)) {
            p + ggplot2::geom_density(ggplot2::aes(x = value), fill = "steelblue", alpha = 0.7)
          } else {
            p + ggplot2::geom_density(ggplot2::aes(x = value, color = group, fill = group), alpha = 0.3)
          }
        }
      )

      # Add labels and theme
      p <- p +
        ggplot2::labs(title = title, x = xlab, y = ylab %||% type) +
        theme_zinasuite()

      p
    },

    #' @description
    #' Create correlation plot
    #' @param x First variable
    #' @param y Second variable
    #' @param type Plot type: "scatter", "hex", "density2d"
    #' @param add_regression Whether to add regression line
    #' @param title Plot title
    #' @return ggplot object
    plot_correlation = function(x, y,
                                 type = c("scatter", "hex", "density2d"),
                                 add_regression = TRUE,
                                 title = NULL) {
      type <- match.arg(type)

      # Prepare data
      df <- data.frame(x = x, y = y)
      df <- df[stats::complete.cases(df), ]

      # Calculate correlation
      cor_result <- stats::cor.test(df$x, df$y)
      cor_label <- sprintf("r = %.3f, p = %.3e", cor_result$estimate, cor_result$p.value)

      # Create base plot
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y))

      # Add geometry
      p <- switch(type,
        "scatter" = p + ggplot2::geom_point(alpha = 0.5, color = "steelblue"),
        "hex" = {
          if (requireNamespace("hexbin", quietly = TRUE)) {
            p + ggplot2::geom_hex(bins = 30)
          } else {
            p + ggplot2::geom_point(alpha = 0.5, color = "steelblue")
          }
        },
        "density2d" = p + ggplot2::geom_density_2d(color = "steelblue")
      )

      # Add regression line
      if (add_regression) {
        p <- p + ggplot2::geom_smooth(method = "lm", color = "red", se = TRUE)
      }

      # Add labels
      p <- p +
        ggplot2::labs(
          title = title %||% "Correlation Plot",
          subtitle = cor_label,
          x = "X",
          y = "Y"
        ) +
        theme_zinasuite()

      p
    },

    #' @description
    #' Create survival plot (Kaplan-Meier)
    #' @param surv_result Survival analysis result from analyze_survival()
    #' @param conf.int Show confidence intervals
    #' @param risk.table Show risk table
    #' @param title Plot title
    #' @return ggplot object (or ggsurvplot if survminer available)
    plot_survival = function(surv_result,
                              conf.int = TRUE,
                              risk.table = TRUE,
                              title = NULL) {
      if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Package 'survival' is required for survival plots")
      }

      km_fit <- surv_result$km$fit

      # Check if survminer is available for better plots
      if (requireNamespace("survminer", quietly = TRUE)) {
        survminer::ggsurvplot(
          km_fit,
          conf.int = conf.int,
          risk.table = risk.table,
          title = title %||% "Kaplan-Meier Survival Curve",
          xlab = "Time",
          ylab = "Survival Probability",
          pval = !is.null(surv_result$km$pvalue),
          pval.method = TRUE
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

        p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = surv, color = group)) +
          ggplot2::geom_step() +
          ggplot2::labs(
            title = title %||% "Kaplan-Meier Survival Curve",
            x = "Time",
            y = "Survival Probability"
          ) +
          theme_zinasuite()

        if (conf.int) {
          p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2)
        }

        p
      }
    },

    #' @description
    #' Create heatmap
    #' @param data Matrix or data frame
    #' @param scale Scale data: "row", "column", "none"
    #' @param clustering Clustering: "both", "row", "column", "none"
    #' @param title Plot title
    #' @param colors Color palette
    #' @return ggplot object (or pheatmap if available)
    plot_heatmap = function(data,
                            scale = c("none", "row", "column"),
                            clustering = c("both", "row", "column", "none"),
                            title = NULL,
                            colors = NULL) {
      scale <- match.arg(scale)
      clustering <- match.arg(clustering)

      # Use pheatmap if available
      if (requireNamespace("pheatmap", quietly = TRUE)) {
        cluster_rows <- clustering %in% c("both", "row")
        cluster_cols <- clustering %in% c("both", "column")

        pheatmap::pheatmap(
          data,
          scale = scale,
          cluster_rows = cluster_rows,
          cluster_cols = cluster_cols,
          main = title,
          color = colors %||% get_palette(100, "gradient")
        )
      } else {
        # Basic ggplot2 version
        df <- as.data.frame(as.table(as.matrix(data)))
        colnames(df) <- c("row", "col", "value")

        ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = value)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradientn(colors = colors %||% get_palette(100, "gradient")) +
          ggplot2::labs(title = title) +
          theme_zinasuite() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      }
    },

    #' @description
    #' Create volcano plot
    #' @param data Data frame with logFC and pvalue columns
    #' @param logfc_col Column name for log fold change
    #' @param pval_col Column name for p-value
    #' @param title Plot title
    #' @return ggplot object
    plot_volcano = function(data,
                            logfc_col = "logFC",
                            pval_col = "pvalue",
                            title = NULL) {
      # Calculate -log10(p-value)
      data$neg_log10_pval <- -log10(data[[pval_col]])

      # Define significance
      data$significance <- "Not significant"
      data$significance[data[[pval_col]] < 0.05 & abs(data[[logfc_col]]) > 1] <- "Significant"
      data$significance[data[[pval_col]] < 0.05 & data[[logfc_col]] > 1] <- "Up-regulated"
      data$significance[data[[pval_col]] < 0.05 & data[[logfc_col]] < -1] <- "Down-regulated"

      ggplot2::ggplot(data, ggplot2::aes_string(x = logfc_col, y = "neg_log10_pval", color = "significance")) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_color_manual(values = c(
          "Not significant" = "grey",
          "Significant" = "blue",
          "Up-regulated" = "red",
          "Down-regulated" = "green"
        )) +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
        ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
        ggplot2::labs(
          title = title %||% "Volcano Plot",
          x = "Log Fold Change",
          y = "-log10(p-value)"
        ) +
        theme_zinasuite()
    },

    #' @description
    #' Create forest plot for Cox regression results
    #' @param data Data frame with hr, lower, upper, pvalue columns
    #' @param title Plot title
    #' @return ggplot object
    plot_forest = function(data, title = NULL) {
      # Ensure required columns exist
      required_cols <- c("hr", "lower", "upper")
      if (!all(required_cols %in% colnames(data))) {
        stop("Data must contain columns: hr, lower, upper")
      }

      # Add gene/variable column if not present
      if (!"gene" %in% colnames(data) && !"variable" %in% colnames(data)) {
        data$variable <- rownames(data)
      }

      var_col <- if ("gene" %in% colnames(data)) "gene" else "variable"

      ggplot2::ggplot(data, ggplot2::aes_string(x = var_col, y = "hr")) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = title %||% "Forest Plot",
          x = "Variable",
          y = "Hazard Ratio (95% CI)"
        ) +
        theme_zinasuite()
    },

    #' @description
    #' Set theme for plots
    #' @param theme Theme name: "default", "publication", "minimal"
    set_theme = function(theme = c("default", "publication", "minimal")) {
      theme <- match.arg(theme)
      private$theme <- theme
    },

    #' @description
    #' Save plot to file
    #' @param plot ggplot object
    #' @param filename Output filename
    #' @param width Width in inches
    #' @param height Height in inches
    #' @param dpi Resolution
    save_plot = function(plot, filename, width = 8, height = 6, dpi = 300) {
      save_plot(plot, filename, width, height, dpi)
    }
  ),

  private = list(
    theme = "default"
  )
)
