#' Survival Visualization
#'
#' @description
#' Create Kaplan-Meier survival curves and forest plots for survival analysis.
#'
#' @param surv_result Result from analyze_survival function
#' @param type Plot type: "km" (Kaplan-Meier) or "forest"
#' @param title Plot title
#' @param ... Additional arguments
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Load survival data
#' surv_data <- load_data("tcga_surv")
#'
#' # Perform survival analysis
#' result <- analyze_survival(
#'   time = surv_data$OS.time,
#'   status = surv_data$OS,
#'   group = surv_data$Type,
#'   analysis_type = "both"
#' )
#'
#' # Plot KM curve
#' p <- vis_survival(result, type = "km", title = "Overall Survival")
#' print(p)
#' }
vis_survival <- function(surv_result,
                         type = c("km", "forest"),
                         title = NULL,
                         ...) {
  type <- match.arg(type)

  if (type == "km") {
    vis_kaplan_meier(surv_result$km, title, ...)
  } else {
    vis_forest_plot(surv_result$cox, title, ...)
  }
}

#' Kaplan-Meier Survival Curve
#'
#' @description
#' Create a Kaplan-Meier survival curve plot.
#'
#' @param km_result KM analysis result from analyze_survival
#' @param title Plot title
#' @param conf.int Whether to show confidence intervals
#' @param risk.table Whether to show risk table
#' @return ggplot object or ggsurvplot object
#' @keywords internal
vis_kaplan_meier <- function(km_result, title = NULL, conf.int = TRUE, risk.table = FALSE) {
  fit <- km_result$fit
  
  # Always use basic ggplot2 version for reliability
  vis_kaplan_meier_basic(fit, title, conf.int, km_result$pvalue)
}

#' Basic Kaplan-Meier plot using ggplot2
#' @keywords internal
vis_kaplan_meier_basic <- function(fit, title = NULL, conf.int = TRUE, pvalue = NULL) {
  # Extract survival data from fit using summary
  surv_summary <- summary(fit, times = fit$time)
  
  # Create data frame for plotting
  if (is.null(fit$strata)) {
    # Single curve
    plot_data <- data.frame(
      time = surv_summary$time,
      surv = surv_summary$surv,
      lower = surv_summary$lower,
      upper = surv_summary$upper,
      group = "All"
    )
  } else {
    # Multiple strata - extract strata information from summary
    strata_vec <- as.character(surv_summary$strata)
    # Clean up strata names (remove "group=" prefix if present)
    strata_vec <- gsub("^group=", "", strata_vec)
    plot_data <- data.frame(
      time = surv_summary$time,
      surv = surv_summary$surv,
      lower = surv_summary$lower,
      upper = surv_summary$upper,
      group = strata_vec
    )
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$time, y = .data$surv, color = .data$group)) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::labs(
      title = title %||% "Kaplan-Meier Survival Curve",
      x = "Time",
      y = "Survival Probability",
      color = "Group"
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()
  
  # Add confidence intervals if requested
  if (conf.int && !is.null(plot_data$lower)) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$group), 
                                   alpha = 0.2, color = NA)
  }
  
  # Add p-value if available
  if (!is.null(pvalue)) {
    p <- p + ggplot2::annotate("text", x = max(plot_data$time) * 0.1, y = 0.1, 
                               label = paste("p =", format(pvalue, digits = 3)),
                               hjust = 0, size = 4)
  }
  
  p
}

#' Forest Plot for Cox Regression
#'
#' @description
#' Create a forest plot for Cox regression results.
#'
#' @param cox_result Cox regression result from analyze_cox_regression
#' @param title Plot title
#' @return ggplot object
#' @keywords internal
vis_forest_plot <- function(cox_result, title = NULL) {
  if (is.null(cox_result) || is.null(cox_result$coefficients)) {
    stop("No Cox regression results to plot")
  }

  # Extract coefficients and confidence intervals
  coef_data <- as.data.frame(cox_result$coefficients)
  coef_data$variable <- rownames(coef_data)
  colnames(coef_data) <- c("coef", "hr", "se", "z", "pvalue", "variable")

  # Get confidence intervals
  ci_data <- as.data.frame(cox_result$conf.int)
  ci_data$variable <- rownames(ci_data)

  # Merge data
  plot_data <- merge(coef_data, ci_data, by = "variable")
  plot_data$variable <- gsub("group", "", plot_data$variable)

  # Create forest plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$hr, y = stats::reorder(.data$variable, .data$hr))) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$`lower .95`, xmax = .data$`upper .95`), height = 0.2) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = title %||% "Cox Regression Forest Plot",
      x = "Hazard Ratio",
      y = NULL
    ) +
    ggplot2::theme_minimal()

  p
}

#' Univariate Cox Forest Plot
#'
#' @description
#' Create a forest plot for univariate Cox regression results across multiple genes.
#'
#' @param unicox_result Result from analyze_unicox_batch function
#' @param title Plot title
#' @param top_n Show top N significant genes (default: all)
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform univariate Cox analysis
#' genes <- c("TP53", "BRCA1", "EGFR", "KRAS")
#' surv_data <- load_data("tcga_surv")
#' result <- analyze_unicox_batch(genes, surv_data)
#'
#' # Create forest plot
#' p <- vis_unicox_forest(result, title = "Univariate Cox Analysis")
#' print(p)
#' }
vis_unicox_forest <- function(unicox_result, title = NULL, top_n = NULL) {
  # Check if input is a gene list (character vector) or analysis result (data frame)
  if (is.character(unicox_result)) {
    # Input is a gene list, need to run analysis first
    message("Running univariate Cox analysis for ", length(unicox_result), " genes...")
    data <- analyze_unicox_batch(unicox_result)
  } else {
    # Input is already analysis result
    data <- unicox_result
  }

  # Filter valid results
  data <- data[!is.na(data$hr), ]

  # Select top N if specified
  if (!is.null(top_n)) {
    data <- head(data, top_n)
  }

  # Create forest plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$hr, y = stats::reorder(.data$gene, .data$hr))) +
    ggplot2::geom_point(ggplot2::aes(color = .data$pvalue < 0.05), size = 3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lower, xmax = .data$upper), height = 0.2) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
    ggplot2::labs(
      title = title %||% "Univariate Cox Regression",
      subtitle = "Red points indicate significant genes (p < 0.05)",
      x = "Hazard Ratio",
      y = "Gene"
    ) +
    ggplot2::theme_minimal()

  p
}

#' Survival Analysis by Gene Expression
#'
#' @description
#' Create Kaplan-Meier curves comparing survival between high and low expression groups.
#'
#' @param gene Gene symbol
#' @param source Data source (default: "tcga")
#' @param cancer TCGA cancer type code (e.g., "BRCA", "LUAD"). If NULL, uses all TCGA samples.
#' @param cutoff_method Cutoff method: "median", "tertile", "quartile"
#' @param title Plot title
#' @return ggsurvplot object or ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create survival plot by TP53 expression in BRCA
#' p <- vis_survival_by_gene(
#'   gene = "TP53",
#'   cancer = "BRCA",
#'   cutoff_method = "median",
#'   title = "Survival by TP53 Expression"
#' )
#' print(p)
#' }
vis_survival_by_gene <- function(gene,
                                 cancer = NULL,
                                 source = "tcga",
                                 cutoff_method = c("median", "tertile", "quartile"),
                                 title = NULL) {
  cutoff_method <- match.arg(cutoff_method)

  # Perform analysis
  result <- analyze_survival_by_expression(
    gene = gene,
    cancer = cancer,
    cutoff_method = cutoff_method,
    source = source
  )

  # Create plot
  vis_kaplan_meier(
    result$survival$km,
    title = title %||% paste("Survival by", gene, "Expression (", cutoff_method, "cutoff)")
  )
}
