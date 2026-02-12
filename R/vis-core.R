#' Visualization Core Functions
#'
#' @description
#' Core visualization functions for ZinaSuite based on ggplot2.
#'
#' @name vis-core
#' @keywords internal
NULL

#' Set ZinaSuite Theme
#'
#' @description
#' Apply a consistent theme for all ZinaSuite visualizations.
#'
#' @param base_size Base font size
#' @param base_family Base font family
#' @return ggplot2 theme object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = mpg, y = wt)) +
#'   geom_point() +
#'   theme_zinasuite()
#' }
theme_zinasuite <- function(base_size = 12, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.position = "right",
      legend.box = "vertical"
    )
}

#' Get Color Palette
#'
#' @description
#' Get a color palette for visualizations.
#'
#' @param n Number of colors
#' @param palette Palette name: "default", "cancer", "gradient"
#' @return Vector of colors
#' @export
#'
#' @examples
#' # Get 5 colors from default palette
#' colors <- get_palette(5)
#'
#' # Get gradient palette
#' grad_colors <- get_palette(10, palette = "gradient")
get_palette <- function(n, palette = c("default", "cancer", "gradient", "diverging", "nature")) {
  palette <- match.arg(palette)

  switch(palette,
    "default" = {
      if (n <= 8) {
        RColorBrewer::brewer.pal(max(n, 3), "Set2")[1:n]
      } else {
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
      }
    },
    "cancer" = , "nature" = {
      # TCGA cancer type colors (Nature style)
      tcga_colors <- c(
        "#ED2891", "#F7941D", "#8DC63F", "#00AEEF", "#F9ED32",
        "#C4A4CC", "#00A651", "#F58220", "#6E298D", "#D2232A",
        "#00A99D", "#F15A29", "#8E76B5", "#FDB913", "#00AE4D",
        "#0072BC", "#F15A22", "#8CC63F", "#ED1C24", "#2E3192",
        "#00A651", "#F7941D", "#92278F", "#00AEEF", "#8DC63F"
      )
      if (n <= length(tcga_colors)) {
        tcga_colors[1:n]
      } else {
        grDevices::colorRampPalette(tcga_colors)(n)
      }
    },
    "gradient" = {
      grDevices::colorRampPalette(c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"))(n)
    },
    "diverging" = {
      grDevices::colorRampPalette(c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020"))(n)
    }
  )
}

#' Validate Plot Data
#'
#' @description
#' Internal function to validate data before plotting.
#'
#' @param data Data frame to validate
#' @param required_cols Required column names
#' @keywords internal
validate_plot_data <- function(data, required_cols) {
  if (!is.data.frame(data)) {
    stop("Data must be a data frame")
  }

  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  invisible(TRUE)
}

#' Add Statistical Annotation
#'
#' @description
#' Add statistical test results to plots.
#'
#' @param p ggplot object
#' @param test_result Statistical test result
#' @param x X position for annotation
#' @param y Y position for annotation
#' @return ggplot object with annotation
#' @keywords internal
add_stat_annotation <- function(p, test_result, x = NULL, y = NULL) {
  if (is.null(test_result) || is.null(test_result$pvalue)) {
    return(p)
  }

  pval <- test_result$pvalue
  if (pval < 0.001) {
    label <- "p < 0.001***"
  } else if (pval < 0.01) {
    label <- sprintf("p = %.3f**", pval)
  } else if (pval < 0.05) {
    label <- sprintf("p = %.3f*", pval)
  } else {
    label <- sprintf("p = %.3f (ns)", pval)
  }

  if (is.null(x) || is.null(y)) {
    # Auto-position based on plot limits
    plot_build <- ggplot2::ggplot_build(p)
    y_range <- plot_build$layout$panel_params[[1]]$y.range
    x_range <- plot_build$layout$panel_params[[1]]$x.range
    x <- mean(x_range)
    y <- y_range[2] * 0.95
  }

  p + ggplot2::annotate(
    "text",
    x = x,
    y = y,
    label = label,
    size = 4,
    fontface = "bold",
    hjust = 0.5
  )
}

#' Save Plot
#'
#' @description
#' Save a ggplot with consistent settings.
#'
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution
#' @param ... Additional arguments to ggsave
#' @export
#'
#' @examples
#' \dontrun{
#' p <- ggplot(mtcars, aes(x = mpg, y = wt)) + geom_point()
#' save_plot(p, "myplot.pdf", width = 8, height = 6)
#' }
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300, ...) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )
}
