#' TCGA Pipeline Analysis Functions
#'
#' @description
#' Deep analysis pipeline functions for TCGA data, supporting correlation,
#' comparison, survival, and cross-omics analysis with multiple modes.
#'
#' @name analysis-tcga-pipeline
NULL

# Helper function for gene pair correlation visualization
#' @keywords internal
vis_gene_pair_cor <- function(gene1, gene2, data_type1 = "mRNA", data_type2 = "mRNA",
                               cancer_choose = NULL, cor_method = "pearson",
                               use_regline = TRUE, alpha = 0.6) {
  # Use the existing vis_gene_correlation function
  vis_gene_correlation(
    gene1 = gene1,
    gene2 = gene2,
    data_type1 = data_type1,
    data_type2 = data_type2,
    cancer = cancer_choose,
    method = cor_method,
    use_regline = use_regline,
    alpha = alpha
  )
}

# TCGA Correlation Analysis Functions -----------------------------------------

#' Run TCGA Correlation Analysis - One-to-One
#'
#' @param molecule1 First molecule identifier
#' @param molecule2 Second molecule identifier
#' @param cancer Cancer type code
#' @param data_type1 Data type for first molecule
#' @param data_type2 Data type for second molecule
#' @param cor_method Correlation method (pearson, spearman, kendall)
#' @param show_regline Whether to show regression line
#' @param point_alpha Point transparency
#' @return List containing plot, data, and stats
#' @export
run_tcga_cor_o2o <- function(molecule1, molecule2, cancer, data_type1, data_type2, cor_method, show_regline = TRUE, point_alpha = 0.6) {
  # Use existing vis functions for consistency
  plot <- vis_gene_pair_cor(
    gene1 = molecule1,
    gene2 = molecule2,
    data_type1 = data_type1,
    data_type2 = data_type2,
    cancer_choose = cancer,
    cor_method = cor_method,
    use_regline = show_regline,
    alpha = point_alpha
  )

  # Extract data from plot attribute
  data <- attr(plot, "data")

  # Calculate statistics
  cor_test <- stats::cor.test(data$Gene1, data$Gene2, method = cor_method)

  stats <- list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    method = cor_method,
    n_samples = nrow(data)
  )

  list(plot = plot, data = data, stats = stats)
}

#' Run TCGA Correlation Analysis - One-to-Many
#'
#' @param molecule1 First molecule identifier
#' @param molecule2 Second molecule identifier
#' @param data_type1 Data type for first molecule
#' @param data_type2 Data type for second molecule
#' @param cor_method Correlation method
#' @return List containing plot, data, and stats
#' @export
run_tcga_cor_o2m <- function(molecule1, molecule2, data_type1, data_type2, cor_method) {
  # Get all TCGA cancers
  tcga_cancers <- load_data("tcga_cancers")

  # Calculate correlation for each cancer
  cor_results <- purrr::map_dfr(tcga_cancers, function(cancer) {
    tryCatch({
      plot <- vis_gene_pair_cor(
        gene1 = molecule1,
        gene2 = molecule2,
        data_type1 = data_type1,
        data_type2 = data_type2,
        cancer_choose = cancer,
        cor_method = cor_method,
        use_regline = FALSE,
        alpha = 0.5
      )

      data <- attr(plot, "data")
      if (is.null(data) || nrow(data) < 10) {
        return(data.frame(Cancer = cancer, Correlation = NA, P_value = NA, N = 0))
      }

      cor_test <- stats::cor.test(data$Gene1, data$Gene2, method = cor_method)

      data.frame(
        Cancer = cancer,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value,
        N = nrow(data)
      )
    }, error = function(e) {
      data.frame(Cancer = cancer, Correlation = NA, P_value = NA, N = 0)
    })
  })

  # Remove NAs
  cor_results <- cor_results[!is.na(cor_results$Correlation), ]

  if (nrow(cor_results) == 0) {
    stop("No valid correlation results across cancers")
  }

  # Create forest plot
  plot <- ggplot2::ggplot(cor_results, ggplot2::aes(x = stats::reorder(.data$Cancer, .data$Correlation), y = .data$Correlation)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$N, color = .data$P_value < 0.05)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(molecule1, "vs", molecule2, "across TCGA cancers"),
      x = "Cancer Type",
      y = "Correlation Coefficient",
      size = "Sample Size",
      color = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = cor_results, stats = summary(cor_results$Correlation))
}

#' Run TCGA Correlation Analysis - Many-to-One
#'
#' @param molecules_text Text with multiple molecules (one per line)
#' @param cancer Cancer type code
#' @param data_type Data type
#' @param cor_method Correlation method
#' @return List containing plot, data, and stats
#' @export
run_tcga_cor_m2o <- function(molecules_text, cancer, data_type, cor_method) {
  # Parse molecules
  molecules <- strsplit(molecules_text, "\n")[[1]]
  molecules <- trimws(molecules[molecules != ""])

  if (length(molecules) < 2) {
    stop("At least 2 molecules required for many-to-one analysis")
  }

  # Use first molecule as reference
  ref_molecule <- molecules[1]

  # Calculate correlation for each molecule
  cor_results <- purrr::map_dfr(molecules, function(mol) {
    tryCatch({
      plot <- vis_gene_pair_cor(
        gene1 = ref_molecule,
        gene2 = mol,
        data_type1 = data_type,
        data_type2 = data_type,
        cancer_choose = cancer,
        cor_method = cor_method,
        use_regline = FALSE,
        alpha = 0.5
      )

      data <- attr(plot, "data")
      if (is.null(data) || nrow(data) < 5) {
        return(data.frame(Molecule = mol, Correlation = NA, P_value = NA, N = 0))
      }

      cor_test <- stats::cor.test(data$Gene1, data$Gene2, method = cor_method)

      data.frame(
        Molecule = mol,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value,
        N = nrow(data)
      )
    }, error = function(e) {
      data.frame(Molecule = mol, Correlation = NA, P_value = NA, N = 0)
    })
  })

  # Remove NAs
  cor_results <- cor_results[!is.na(cor_results$Correlation), ]

  # Create bar plot
  plot <- ggplot2::ggplot(cor_results, ggplot2::aes(x = stats::reorder(.data$Molecule, .data$Correlation), y = .data$Correlation)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = .data$P_value < 0.05)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "steelblue")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Correlation with", ref_molecule, "in", cancer),
      x = "Molecule",
      y = "Correlation Coefficient",
      fill = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = cor_results, stats = summary(cor_results$Correlation))
}

# TCGA Comparison Analysis Functions -----------------------------------------

#' Run TCGA Comparison Analysis - One-to-One
#'
#' @param molecule Molecule identifier
#' @param cancer Cancer type code
#' @param data_type Data type
#' @return List containing plot, data, and stats
#' @export
run_tcga_comp_o2o <- function(molecule, cancer, data_type) {
  # Use existing vis function
  plot <- vis_toil_TvsN(
    gene = molecule,
    data_type = data_type,
    cancer_choose = cancer
  )

  # Extract data from plot attribute
  data <- attr(plot, "data")

  if (is.null(data)) {
    stop("No data available for comparison")
  }

  # Calculate statistics
  tumor_vals <- data$Expression[data$Type == "tumor"]
  normal_vals <- data$Expression[data$Type == "normal"]

  if (length(tumor_vals) > 0 && length(normal_vals) > 0) {
    test_result <- stats::wilcox.test(tumor_vals, normal_vals)
    p_value <- test_result$p.value
  } else {
    p_value <- NA
  }

  stats <- list(
    p_value = p_value,
    tumor_median = stats::median(tumor_vals, na.rm = TRUE),
    normal_median = stats::median(normal_vals, na.rm = TRUE),
    fold_change = stats::median(tumor_vals, na.rm = TRUE) / stats::median(normal_vals, na.rm = TRUE)
  )

  list(plot = plot, data = data, stats = stats)
}

#' Run TCGA Comparison Analysis - One-to-Many
#'
#' @param molecule Molecule identifier
#' @param data_type Data type
#' @return List containing plot, data, and stats
#' @export
run_tcga_comp_o2m <- function(molecule, data_type) {
  # Get all TCGA cancers
  tcga_cancers <- load_data("tcga_cancers")

  # Calculate for each cancer
  comp_results <- purrr::map_dfr(tcga_cancers, function(cancer) {
    tryCatch({
      plot <- vis_toil_TvsN(
        gene = molecule,
        data_type = data_type,
        cancer_choose = cancer
      )

      data <- attr(plot, "data")
      if (is.null(data)) {
        return(data.frame(Cancer = cancer, Tumor_Median = NA, Normal_Median = NA, P_value = NA, Fold_Change = NA))
      }

      tumor_vals <- data$Expression[data$Type == "tumor"]
      normal_vals <- data$Expression[data$Type == "normal"]

      if (length(tumor_vals) > 0 && length(normal_vals) > 0) {
        test_result <- stats::wilcox.test(tumor_vals, normal_vals)
        data.frame(
          Cancer = cancer,
          Tumor_Median = stats::median(tumor_vals, na.rm = TRUE),
          Normal_Median = stats::median(normal_vals, na.rm = TRUE),
          P_value = test_result$p.value,
          Fold_Change = stats::median(tumor_vals, na.rm = TRUE) / stats::median(normal_vals, na.rm = TRUE)
        )
      } else {
        data.frame(Cancer = cancer, Tumor_Median = NA, Normal_Median = NA, P_value = NA, Fold_Change = NA)
      }
    }, error = function(e) {
      data.frame(Cancer = cancer, Tumor_Median = NA, Normal_Median = NA, P_value = NA, Fold_Change = NA)
    })
  })

  # Remove NAs
  comp_results <- comp_results[!is.na(comp_results$Tumor_Median), ]

  if (nrow(comp_results) == 0) {
    stop("No valid comparison results across cancers")
  }

  # Create plot
  plot <- ggplot2::ggplot(comp_results, ggplot2::aes(x = stats::reorder(.data$Cancer, .data$Fold_Change), y = .data$Fold_Change)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = .data$P_value < 0.05)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "steelblue")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(molecule, "Tumor vs Normal across TCGA cancers"),
      x = "Cancer Type",
      y = "Fold Change (Tumor/Normal)",
      fill = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = comp_results, stats = summary(comp_results$Fold_Change))
}

#' Run TCGA Comparison Analysis - Many-to-One
#'
#' @param molecules_text Text with multiple molecules
#' @param cancer Cancer type code
#' @param data_type Data type
#' @return List containing plot, data, and stats
#' @export
run_tcga_comp_m2o <- function(molecules_text, cancer, data_type) {
  # Parse molecules
  molecules <- strsplit(molecules_text, "\n")[[1]]
  molecules <- trimws(molecules[molecules != ""])

  if (length(molecules) < 2) {
    stop("At least 2 molecules required for many-to-one analysis")
  }

  # Calculate for each molecule
  comp_results <- purrr::map_dfr(molecules, function(mol) {
    tryCatch({
      plot <- vis_toil_TvsN(
        gene = mol,
        data_type = data_type,
        cancer_choose = cancer
      )

      data <- attr(plot, "data")
      if (is.null(data)) {
        return(data.frame(Molecule = mol, Tumor_Median = NA, Normal_Median = NA, P_value = NA, Fold_Change = NA))
      }

      tumor_vals <- data$Expression[data$Type == "tumor"]
      normal_vals <- data$Expression[data$Type == "normal"]

      if (length(tumor_vals) > 0 && length(normal_vals) > 0) {
        test_result <- stats::wilcox.test(tumor_vals, normal_vals)
        data.frame(
          Molecule = mol,
          Tumor_Median = stats::median(tumor_vals, na.rm = TRUE),
          Normal_Median = stats::median(normal_vals, na.rm = TRUE),
          P_value = test_result$p.value,
          Fold_Change = stats::median(tumor_vals, na.rm = TRUE) / stats::median(normal_vals, na.rm = TRUE)
        )
      } else {
        data.frame(Molecule = mol, Tumor_Median = NA, Normal_Median = NA, P_value = NA, Fold_Change = NA)
      }
    }, error = function(e) {
      data.frame(Molecule = mol, Tumor_Median = NA, Normal_Median = NA, P_value = NA, Fold_Change = NA)
    })
  })

  # Remove NAs
  comp_results <- comp_results[!is.na(comp_results$Tumor_Median), ]

  # Create plot
  plot <- ggplot2::ggplot(comp_results, ggplot2::aes(x = stats::reorder(.data$Molecule, .data$Fold_Change), y = .data$Fold_Change)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = .data$P_value < 0.05)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "steelblue")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Tumor vs Normal in", cancer),
      x = "Molecule",
      y = "Fold Change (Tumor/Normal)",
      fill = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = comp_results, stats = summary(comp_results$Fold_Change))
}

# TCGA Survival Analysis Functions -------------------------------------------

#' Run TCGA Survival Analysis - One-to-One
#'
#' @param molecule Molecule identifier
#' @param cancer Cancer type code
#' @param data_type Data type
#' @param surv_measure Survival measure (OS or PFI)
#' @param cutoff_mode Cutoff mode for grouping
#' @return List containing plot, data, and stats
#' @export
run_tcga_sur_o2o <- function(molecule, cancer, data_type, surv_measure = "OS", cutoff_mode = "Auto") {
  # Use existing vis function
  plot <- vis_survival(
    gene = molecule,
    data_type = data_type,
    cancer = cancer,
    measure = surv_measure,
    cutoff_mode = cutoff_mode
  )

  # Extract data from plot attribute
  data <- attr(plot, "data")

  if (is.null(data)) {
    stop("No survival data available")
  }

  # Calculate statistics
  fit <- survival::survdiff(survival::Surv(time, status) ~ group, data = data)
  p_value <- 1 - stats::pchisq(fit$chisq, length(fit$n) - 1)

  stats <- list(
    p_value = p_value,
    log_rank_chisq = fit$chisq,
    n_high = sum(data$group == "high"),
    n_low = sum(data$group == "low")
  )

  list(plot = plot, data = data, stats = stats)
}

#' Run TCGA Survival Analysis - One-to-Many
#'
#' @param molecule Molecule identifier
#' @param data_type Data type
#' @param surv_measure Survival measure
#' @param cutoff_mode Cutoff mode
#' @return List containing plot, data, and stats
#' @export
run_tcga_sur_o2m <- function(molecule, data_type, surv_measure = "OS", cutoff_mode = "Auto") {
  # Get all TCGA cancers
  tcga_cancers <- load_data("tcga_cancers")

  # Calculate for each cancer
  surv_results <- purrr::map_dfr(tcga_cancers, function(cancer) {
    tryCatch({
      plot <- vis_survival(
        gene = molecule,
        data_type = data_type,
        cancer = cancer,
        measure = surv_measure,
        cutoff_mode = cutoff_mode
      )

      data <- attr(plot, "data")
      if (is.null(data) || length(unique(data$group)) < 2) {
        return(data.frame(Cancer = cancer, HR = NA, P_value = NA, N = 0))
      }

      # Fit Cox model
      fit <- survival::coxph(survival::Surv(time, status) ~ group, data = data)
      summary_fit <- summary(fit)

      data.frame(
        Cancer = cancer,
        HR = summary_fit$conf.int[1, "exp(coef)"],
        CI_lower = summary_fit$conf.int[1, "lower .95"],
        CI_upper = summary_fit$conf.int[1, "upper .95"],
        P_value = summary_fit$coefficients[1, "Pr(>|z|)"],
        N = nrow(data)
      )
    }, error = function(e) {
      data.frame(Cancer = cancer, HR = NA, CI_lower = NA, CI_upper = NA, P_value = NA, N = 0)
    })
  })

  # Remove NAs
  surv_results <- surv_results[!is.na(surv_results$HR), ]

  if (nrow(surv_results) == 0) {
    stop("No valid survival results across cancers")
  }

  # Create forest plot
  plot <- ggplot2::ggplot(surv_results, ggplot2::aes(x = stats::reorder(.data$Cancer, .data$HR), y = .data$HR)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$N, color = .data$P_value < 0.05)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$CI_lower, ymax = .data$CI_upper), width = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(molecule, "Survival across TCGA cancers"),
      x = "Cancer Type",
      y = "Hazard Ratio (95% CI)",
      size = "Sample Size",
      color = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = surv_results, stats = summary(surv_results$HR))
}

#' Run TCGA Survival Analysis - Many-to-One
#'
#' @param molecules_text Text with multiple molecules
#' @param cancer Cancer type code
#' @param data_type Data type
#' @param surv_measure Survival measure
#' @param cutoff_mode Cutoff mode
#' @return List containing plot, data, and stats
#' @export
run_tcga_sur_m2o <- function(molecules_text, cancer, data_type, surv_measure = "OS", cutoff_mode = "Auto") {
  # Parse molecules
  molecules <- strsplit(molecules_text, "\n")[[1]]
  molecules <- trimws(molecules[molecules != ""])

  if (length(molecules) < 2) {
    stop("At least 2 molecules required for many-to-one analysis")
  }

  # Calculate for each molecule
  surv_results <- purrr::map_dfr(molecules, function(mol) {
    tryCatch({
      plot <- vis_survival(
        gene = mol,
        data_type = data_type,
        cancer = cancer,
        measure = surv_measure,
        cutoff_mode = cutoff_mode
      )

      data <- attr(plot, "data")
      if (is.null(data) || length(unique(data$group)) < 2) {
        return(data.frame(Molecule = mol, HR = NA, P_value = NA, N = 0))
      }

      # Fit Cox model
      fit <- survival::coxph(survival::Surv(time, status) ~ group, data = data)
      summary_fit <- summary(fit)

      data.frame(
        Molecule = mol,
        HR = summary_fit$conf.int[1, "exp(coef)"],
        CI_lower = summary_fit$conf.int[1, "lower .95"],
        CI_upper = summary_fit$conf.int[1, "upper .95"],
        P_value = summary_fit$coefficients[1, "Pr(>|z|)"],
        N = nrow(data)
      )
    }, error = function(e) {
      data.frame(Molecule = mol, HR = NA, P_value = NA, N = 0)
    })
  })

  # Remove NAs
  surv_results <- surv_results[!is.na(surv_results$HR), ]

  # Create forest plot
  plot <- ggplot2::ggplot(surv_results, ggplot2::aes(x = stats::reorder(.data$Molecule, .data$HR), y = .data$HR)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$N, color = .data$P_value < 0.05)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$CI_lower, ymax = .data$CI_upper), width = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Survival in", cancer),
      x = "Molecule",
      y = "Hazard Ratio (95% CI)",
      size = "Sample Size",
      color = "Significant (p<0.05)"
    ) +
    theme_zinasuite()

  list(plot = plot, data = surv_results, stats = summary(surv_results$HR))
}

# TCGA Cross-Omics Analysis Functions ----------------------------------------

#' Run TCGA Cross-Omics Analysis - One-to-Many
#'
#' @param molecule1 First molecule identifier
#' @param molecule2 Second molecule identifier
#' @param cancer Cancer type code
#' @param data_type1 Data type for first molecule
#' @param data_type2 Data type for second molecule
#' @return List containing plot, data, and stats
#' @export
run_tcga_cross_o2m <- function(molecule1, molecule2, cancer, data_type1, data_type2) {
  # Use existing cross-omics function
  plot <- analyze_cross_gene(
    gene1 = molecule1,
    gene2 = molecule2,
    cancer = cancer,
    data_type1 = data_type1,
    data_type2 = data_type2
  )

  # Extract data from plot attribute
  data <- attr(plot, "data")

  if (is.null(data)) {
    stop("No cross-omics data available")
  }

  # Calculate statistics
  cor_test <- stats::cor.test(data$Value1, data$Value2, method = "spearman")

  stats <- list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    method = "spearman",
    n_samples = nrow(data)
  )

  list(plot = plot, data = data, stats = stats)
}
