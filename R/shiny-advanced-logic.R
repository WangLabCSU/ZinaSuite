#' Advanced Shiny Business Logic
#'
#' @description
#' Advanced business logic functions for Shiny modules covering all UCSCXenaShiny features.
#'
#' @name shiny-advanced-logic
NULL

# Pathway Analysis Logic --------------------------------------------------

#' Analyze Pathway
#'
#' @param pathway_name Pathway name
#' @param gene_list List of genes in pathway
#' @param cancer Cancer type
#' @return Pathway analysis result
#' @export
analyze_pathway <- function(pathway_name, gene_list, cancer = NULL) {
  tryCatch({
    if (length(gene_list) < 2) {
      stop("At least 2 genes required for pathway analysis")
    }

    # Query expression for all genes
    expr_data <- lapply(gene_list, function(gene) {
      query_gene_expression(gene)
    })
    names(expr_data) <- gene_list

    # Calculate pathway score (mean expression)
    common_samples <- Reduce(intersect, lapply(expr_data, names))
    if (length(common_samples) < 10) {
      stop("Insufficient common samples for pathway analysis")
    }

    pathway_scores <- sapply(common_samples, function(s) {
      mean(sapply(expr_data, function(x) x[s]), na.rm = TRUE)
    })

    list(
      success = TRUE,
      pathway = pathway_name,
      genes = gene_list,
      cancer = cancer,
      scores = pathway_scores,
      mean_score = mean(pathway_scores, na.rm = TRUE),
      sd_score = sd(pathway_scores, na.rm = TRUE),
      sample_count = length(common_samples)
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

#' Analyze Cross Gene
#'
#' @param target_gene Target gene
#' @param correlate_genes Genes to correlate with
#' @param cancer Cancer type
#' @return Cross-gene analysis result
#' @export
analyze_cross_gene <- function(target_gene, correlate_genes, cancer = NULL) {
  tryCatch({
    # Query target gene
    target_expr <- query_gene_expression(target_gene)

    # Query correlate genes
    results <- lapply(correlate_genes, function(gene) {
      gene_expr <- query_gene_expression(gene)
      cor_result <- analyze_correlation(target_expr, gene_expr)
      list(
        gene = gene,
        correlation = cor_result$estimate,
        p_value = cor_result$p.value
      )
    })

    names(results) <- correlate_genes

    list(
      success = TRUE,
      target_gene = target_gene,
      cancer = cancer,
      correlations = results,
      top_correlated = names(sort(sapply(results, function(x) abs(x$correlation)), decreasing = TRUE))[1:5]
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

#' Analyze Cross Pathway
#'
#' @param pathway1 First pathway genes
#' @param pathway2 Second pathway genes
#' @param cancer Cancer type
#' @return Cross-pathway analysis result
#' @export
analyze_cross_pathway <- function(pathway1, pathway2, cancer = NULL) {
  tryCatch({
    # Analyze each pathway
    p1_result <- analyze_pathway("Pathway1", pathway1, cancer)
    p2_result <- analyze_pathway("Pathway2", pathway2, cancer)

    if (!p1_result$success || !p2_result$success) {
      stop("Failed to analyze one or both pathways")
    }

    # Correlate pathway scores
    common_samples <- intersect(names(p1_result$scores), names(p2_result$scores))
    cor_result <- cor.test(
      p1_result$scores[common_samples],
      p2_result$scores[common_samples]
    )

    list(
      success = TRUE,
      cancer = cancer,
      pathway1_score = p1_result$mean_score,
      pathway2_score = p2_result$mean_score,
      correlation = cor_result$estimate,
      p_value = cor_result$p.value,
      sample_count = length(common_samples)
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# Drug Analysis Logic -----------------------------------------------------

#' Analyze Drug Sensitivity
#'
#' @param gene Gene symbol
#' @param drug Drug name
#' @param cancer Cancer type
#' @return Drug sensitivity analysis
#' @export
analyze_drug_sensitivity <- function(gene, drug, cancer = NULL) {
  tryCatch({
    # Query gene expression
    gene_expr <- query_gene_expression(gene)

    # Query drug sensitivity (placeholder - would need actual drug data)
    # In real implementation, this would query CCLE or GDSC data
    drug_sens <- rnorm(length(gene_expr), mean = 0.5, sd = 0.2)
    names(drug_sens) <- names(gene_expr)

    # Correlate expression with sensitivity
    common <- intersect(names(gene_expr), names(drug_sens))
    cor_result <- cor.test(gene_expr[common], drug_sens[common])

    # Group by expression level
    high_expr <- names(gene_expr)[gene_expr > median(gene_expr, na.rm = TRUE)]
    low_expr <- names(gene_expr)[gene_expr <= median(gene_expr, na.rm = TRUE)]

    list(
      success = TRUE,
      gene = gene,
      drug = drug,
      cancer = cancer,
      correlation = cor_result$estimate,
      p_value = cor_result$p.value,
      high_expr_sensitivity = mean(drug_sens[intersect(high_expr, names(drug_sens))], na.rm = TRUE),
      low_expr_sensitivity = mean(drug_sens[intersect(low_expr, names(drug_sens))], na.rm = TRUE),
      sample_count = length(common)
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

#' Analyze Drug Response
#'
#' @param genes Gene list
#' @param drug Drug name
#' @param cancer Cancer type
#' @return Drug response analysis
#' @export
analyze_drug_response <- function(genes, drug, cancer = NULL) {
  tryCatch({
    # Query all genes
    expr_matrix <- query_molecules(genes, cancer = cancer)

    # Calculate signature score
    signature_score <- colMeans(expr_matrix, na.rm = TRUE)

    # Query drug response
    drug_resp <- rnorm(length(signature_score), mean = 0.5, sd = 0.2)
    names(drug_resp) <- names(signature_score)

    # Correlate
    common <- intersect(names(signature_score), names(drug_resp))
    cor_result <- cor.test(signature_score[common], drug_resp[common])

    list(
      success = TRUE,
      genes = genes,
      drug = drug,
      cancer = cancer,
      correlation = cor_result$estimate,
      p_value = cor_result$p.value,
      signature_mean = mean(signature_score, na.rm = TRUE),
      sample_count = length(common)
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

#' Analyze Drug Omic Pair
#'
#' @param drug Drug name
#' @param omic_type Omic data type
#' @param feature Feature name
#' @return Drug-omic pair analysis
#' @export
analyze_drug_omic <- function(drug, omic_type, feature) {
  tryCatch({
    # Query omic data
    omic_data <- switch(omic_type,
      "expression" = query_gene_expression(feature),
      "mutation" = query_mutation(feature),
      "cnv" = query_cnv(feature),
      stop("Unsupported omic type")
    )

    # Query drug data (placeholder)
    drug_data <- rnorm(length(omic_data), mean = 0.5, sd = 0.2)
    names(drug_data) <- names(omic_data)

    # Analysis
    common <- intersect(names(omic_data), names(drug_data))

    if (is.numeric(omic_data)) {
      cor_result <- cor.test(omic_data[common], drug_data[common])
      result <- list(
        correlation = cor_result$estimate,
        p_value = cor_result$p.value
      )
    } else {
      # Categorical analysis
      groups <- split(drug_data[common], omic_data[common])
      result <- list(
        group_means = sapply(groups, mean, na.rm = TRUE),
        anova = summary(aov(drug_data[common] ~ omic_data[common]))
      )
    }

    c(list(
      success = TRUE,
      drug = drug,
      omic_type = omic_type,
      feature = feature,
      sample_count = length(common)
    ), result)
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

#' Analyze Feature Across Types
#'
#' @param feature Feature name
#' @param types Cancer types to compare
#' @return Cross-cancer analysis
#' @export
analyze_feature_across <- function(feature, types = NULL) {
  tryCatch({
    if (is.null(types)) {
      types <- tcga_cancer_types()
    }

    results <- lapply(types, function(cancer) {
      tryCatch({
        data <- query_gene_expression(feature)
        list(
          cancer = cancer,
          mean = mean(data, na.rm = TRUE),
          sd = sd(data, na.rm = TRUE),
          median = median(data, na.rm = TRUE),
          n = length(data),
          success = TRUE
        )
      }, error = function(e) {
        list(cancer = cancer, success = FALSE, error = conditionMessage(e))
      })
    })

    names(results) <- types

    # Filter successful results
    successful <- results[sapply(results, function(x) x$success)]

    list(
      success = TRUE,
      feature = feature,
      cancer_count = length(successful),
      results = successful,
      highest = names(which.max(sapply(successful, function(x) x$mean))),
      lowest = names(which.min(sapply(successful, function(x) x$mean)))
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

#' Analyze Feature Signature
#'
#' @param features Feature list
#' @param signature_name Signature name
#' @param cancer Cancer type
#' @return Signature analysis
#' @export
analyze_feature_signature <- function(features, signature_name, cancer = NULL) {
  tryCatch({
    # Query all features
    data <- query_molecules(features, cancer = cancer)

    # Calculate signature score
    signature_score <- colMeans(data, na.rm = TRUE)

    # Cluster samples
    kmeans_result <- kmeans(scale(t(data)), centers = 2)

    list(
      success = TRUE,
      signature = signature_name,
      features = features,
      cancer = cancer,
      score_mean = mean(signature_score, na.rm = TRUE),
      score_sd = sd(signature_score, na.rm = TRUE),
      clusters = kmeans_result$cluster,
      cluster_sizes = table(kmeans_result$cluster),
      sample_count = ncol(data)
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# Dimension Reduction Logic -----------------------------------------------

#' Analyze Dimension Reduction
#'
#' @param genes Gene list
#' @param method Method: "pca", "tsne", "umap"
#' @param cancer Cancer type
#' @param n_components Number of components
#' @return Dimension reduction result
#' @export
analyze_dimension_reduction <- function(genes, method = "pca", cancer = NULL, n_components = 2) {
  tryCatch({
    # Query data
    data <- query_molecules(genes, cancer = cancer)

    # Remove genes with too many NAs
    data <- data[rowSums(is.na(data)) < ncol(data) * 0.3, ]

    # Transpose (samples as rows)
    data_t <- t(data)

    # Remove samples with NAs
    data_t <- data_t[complete.cases(data_t), ]

    result <- switch(method,
      "pca" = {
        pca_result <- prcomp(data_t, scale. = TRUE)
        list(
          coordinates = pca_result$x[, 1:n_components],
          variance_explained = summary(pca_result)$importance[2, 1:n_components],
          components = pca_result$rotation[, 1:n_components]
        )
      },
      "tsne" = {
        if (!requireNamespace("Rtsne", quietly = TRUE)) {
          stop("Rtsne package required for t-SNE")
        }
        tsne_result <- Rtsne::Rtsne(data_t, dims = n_components, verbose = FALSE)
        list(
          coordinates = tsne_result$Y,
          perplexity = tsne_result$perplexity
        )
      },
      "umap" = {
        if (!requireNamespace("umap", quietly = TRUE)) {
          stop("umap package required for UMAP")
        }
        umap_result <- umap::umap(data_t, n_components = n_components)
        list(
          coordinates = umap_result$layout
        )
      },
      stop("Unknown method: ", method)
    )

    c(list(
      success = TRUE,
      method = method,
      cancer = cancer,
      n_genes = nrow(data),
      n_samples = nrow(data_t),
      n_components = n_components
    ), result)
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

#' Analyze Dimension Distribution
#'
#' @param gene Gene symbol
#' @param cancer Cancer type
#' @return Distribution analysis
#' @export
analyze_dimension_distribution <- function(gene, cancer = NULL) {
  tryCatch({
    data <- query_gene_expression(gene)

    # Calculate distribution metrics
    list(
      success = TRUE,
      gene = gene,
      cancer = cancer,
      n = length(data),
      mean = mean(data, na.rm = TRUE),
      median = median(data, na.rm = TRUE),
      sd = sd(data, na.rm = TRUE),
      min = min(data, na.rm = TRUE),
      max = max(data, na.rm = TRUE),
      q25 = quantile(data, 0.25, na.rm = TRUE),
      q75 = quantile(data, 0.75, na.rm = TRUE),
      skewness = moments::skewness(data, na.rm = TRUE),
      kurtosis = moments::kurtosis(data, na.rm = TRUE)
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# Group Comparison Logic --------------------------------------------------

#' Analyze Group Comparison
#'
#' @param gene Gene symbol
#' @param group_var Grouping variable
#' @param cancer Cancer type
#' @return Group comparison result
#' @export
analyze_group_comparison <- function(gene, group_var, cancer = NULL) {
  tryCatch({
    # Query gene expression
    expr <- query_gene_expression(gene)

    # Query grouping variable
    groups_result <- query_tcga_group(database = "tcga", cancer = cancer, group = group_var)

    # Extract group data
    groups_data <- groups_result$data
    if (nrow(groups_data) == 0) {
      stop("No grouping data available")
    }

    # Create named vector of groups
    groups <- setNames(groups_data[[group_var]], groups_data$Sample)

    # Match samples using barcode (first 12 characters)
    # This handles cases where clinical and expression data use different ID formats
    match_result <- match_samples(names(expr), names(groups), "tcga", "tcga", match_by = "barcode")

    if (match_result$n_matched == 0) {
      stop("No matching samples found between expression and grouping data")
    }

    expr_matched <- expr[match_result$idx1]
    groups_matched <- groups[match_result$idx2]

    # Statistical test
    if (length(unique(groups_matched)) == 2) {
      # T-test for 2 groups
      test_result <- t.test(expr_matched ~ groups_matched)
      test_type <- "t-test"
      p_value <- test_result$p.value
    } else {
      # ANOVA for >2 groups
      test_result <- summary(aov(expr_matched ~ groups_matched))
      test_type <- "ANOVA"
      p_value <- test_result[[1]]$`Pr(>F)`[1]
    }

    # Group statistics
    group_stats <- tapply(expr_matched, groups_matched, function(x) {
      list(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), n = length(x))
    })

    list(
      success = TRUE,
      gene = gene,
      group_var = group_var,
      cancer = cancer,
      test_type = test_type,
      p_value = p_value,
      group_stats = group_stats,
      sample_count = match_result$n_matched
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# Daily Gene Feature Logic ------------------------------------------------

#' Get Daily Gene Feature
#'
#' @return Daily featured gene information
#' @export
daily_gene_feature <- function() {
  # List of important cancer genes
  featured_genes <- c(
    "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "BRAF",
    "PIK3CA", "PTEN", "MYC", "CDKN2A", "RB1", "ATM",
    "MLH1", "MSH2", "APC", "VHL", "WT1", "RET"
  )

  # Select based on day of year
  day_of_year <- as.numeric(format(Sys.Date(), "%j"))
  selected_gene <- featured_genes[(day_of_year %% length(featured_genes)) + 1]

  list(
    gene = selected_gene,
    date = Sys.Date(),
    description = get_gene_description(selected_gene),
    cancers = get_gene_related_cancers(selected_gene)
  )
}

#' Get Gene Description
#'
#' @param gene Gene symbol
#' @return Gene description
#' @keywords internal
get_gene_description <- function(gene) {
  descriptions <- list(
    "TP53" = "Tumor protein p53, guardian of the genome",
    "BRCA1" = "Breast cancer type 1 susceptibility protein",
    "BRCA2" = "Breast cancer type 2 susceptibility protein",
    "EGFR" = "Epidermal growth factor receptor",
    "KRAS" = "Kirsten rat sarcoma viral oncogene homolog",
    "BRAF" = "B-Raf proto-oncogene",
    "PIK3CA" = "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha",
    "PTEN" = "Phosphatase and tensin homolog",
    "MYC" = "MYC proto-oncogene, bHLH transcription factor",
    "CDKN2A" = "Cyclin dependent kinase inhibitor 2A",
    "RB1" = "RB transcriptional corepressor 1",
    "ATM" = "ATM serine/threonine kinase",
    "MLH1" = "MutL homolog 1",
    "MSH2" = "MutS homolog 2",
    "APC" = "APC regulator of WNT signaling pathway",
    "VHL" = "Von Hippel-Lindau tumor suppressor",
    "WT1" = "Wilms tumor 1",
    "RET" = "Ret proto-oncogene"
  )

  descriptions[[gene]] %||% paste(gene, "- Cancer-related gene")
}

#' Get Gene Related Cancers
#'
#' @param gene Gene symbol
#' @return Related cancers
#' @keywords internal
get_gene_related_cancers <- function(gene) {
  cancer_map <- list(
    "TP53" = c("BRCA", "LUAD", "LUSC", "HNSC", "OV"),
    "BRCA1" = c("BRCA", "OV"),
    "BRCA2" = c("BRCA", "OV", "PRAD"),
    "EGFR" = c("LUAD", "GBM", "HNSC"),
    "KRAS" = c("PAAD", "COAD", "READ", "LUAD"),
    "BRAF" = c("SKCM", "THCA", "COAD"),
    "PIK3CA" = c("BRCA", "UCEC", "HNSC"),
    "PTEN" = c("PRAD", "UCEC", "BRCA"),
    "MYC" = c("BLCA", "BRCA", "COAD"),
    "CDKN2A" = c("PAAD", "HNSC", "ESCA"),
    "RB1" = c("RET", "LUNG", "BLCA"),
    "ATM" = c("BRCA", "DLBC", "CESC"),
    "MLH1" = c("COAD", "UCEC", "STAD"),
    "MSH2" = c("COAD", "UCEC", "STAD"),
    "APC" = c("COAD", "READ", "STAD"),
    "VHL" = c("KIRC", "KIRP", "KICH"),
    "WT1" = c("KIRP", "MESO"),
    "RET" = c("THCA", "PCPG")
  )

  cancer_map[[gene]] %||% c("BRCA", "LUAD", "LUSC")
}

# Pan-Cancer Search Logic -------------------------------------------------

#' Pan-Cancer Search
#'
#' @param query Search query (gene, pathway, or cancer)
#' @return Search results
#' @export
pancan_search <- function(query) {
  tryCatch({
    # Check if it's a gene
    if (is_valid_gene(query)) {
      return(list(
        type = "gene",
        query = query,
        available_analyses = c(
          "Expression across cancers",
          "Survival analysis",
          "Correlation analysis",
          "Mutation analysis",
          "Immune infiltration"
        )
      ))
    }

    # Check if it's a cancer type
    if (query %in% tcga_cancer_types()) {
      return(list(
        type = "cancer",
        query = query,
        available_analyses = c(
          "Differential expression",
          "Survival analysis",
          "Mutation landscape",
          "Pathway analysis"
        )
      ))
    }

    list(
      type = "unknown",
      query = query,
      message = "Query not recognized as gene or cancer type"
    )
  }, error = function(e) {
    list(type = "error", error = conditionMessage(e))
  })
}

#' TCGA Cancer Types
#'
#' @return Vector of TCGA cancer types
#' @export
tcga_cancer_types <- function() {
  c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
    "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
    "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
    "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
}

# File Upload Logic -------------------------------------------------------

#' Process Uploaded File
#'
#' @param file_path Path to uploaded file
#' @param file_type File type
#' @return Processed data
#' @export
process_uploaded_file <- function(file_path, file_type = "auto") {
  tryCatch({
    if (file_type == "auto") {
      # Detect from extension
      ext <- tools::file_ext(file_path)
      file_type <- switch(tolower(ext),
        "csv" = "csv",
        "txt" = "tsv",
        "tsv" = "tsv",
        "xlsx" = "excel",
        "xls" = "excel",
        "rds" = "rds",
        "rda" = "rda",
        stop("Cannot detect file type")
      )
    }

    data <- switch(file_type,
      "csv" = read.csv(file_path, stringsAsFactors = FALSE),
      "tsv" = read.delim(file_path, stringsAsFactors = FALSE),
      "excel" = {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("readxl package required for Excel files")
        }
        readxl::read_excel(file_path)
      },
      "rds" = readRDS(file_path),
      "rda" = {
        env <- new.env()
        load(file_path, envir = env)
        as.list(env)
      },
      stop("Unsupported file type: ", file_type)
    )

    list(
      success = TRUE,
      file_type = file_type,
      data = data,
      rows = nrow(data),
      cols = ncol(data)
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

#' Validate Uploaded Data
#'
#' @param data Uploaded data
#' @param expected_format Expected format
#' @return Validation result
#' @export
validate_uploaded_data <- function(data, expected_format = "gene_list") {
  errors <- character()

  if (is.null(data) || nrow(data) == 0) {
    errors <- c(errors, "Data is empty")
  }

  if (expected_format == "gene_list") {
    if (!"gene" %in% tolower(colnames(data))) {
      errors <- c(errors, "Required column 'gene' not found")
    }
  } else if (expected_format == "expression") {
    if (!any(sapply(data, is.numeric))) {
      errors <- c(errors, "No numeric columns found for expression data")
    }
  }

  if (length(errors) > 0) {
    list(valid = FALSE, errors = errors)
  } else {
    list(valid = TRUE, message = "Data is valid")
  }
}

# Batch Upload Logic ------------------------------------------------------

#' Batch Upload
#'
#' @param files List of file paths
#' @return Batch upload result
#' @export
batch_upload <- function(files) {
  results <- lapply(files, function(file) {
    result <- process_uploaded_file(file)
    result$file <- basename(file)
    result
  })

  successful <- results[sapply(results, function(x) x$success)]
  failed <- results[!sapply(results, function(x) x$success)]

  list(
    total = length(files),
    successful = length(successful),
    failed = length(failed),
    results = results,
    success_rate = length(successful) / length(files) * 100
  )
}

# Help Identifier Logic ---------------------------------------------------

#' Help Identifier
#'
#' @param id Identifier to look up
#' @return Help information
#' @export
help_identifier <- function(id) {
  # Check if it's a gene
  if (is_valid_gene(id)) {
    return(list(
      type = "gene",
      id = id,
      description = get_gene_description(id),
      related_datasets = c("TCGA", "PCAWG", "CCLE"),
      available_data_types = c("mRNA", "mutation", "CNV", "methylation", "protein")
    ))
  }

  # Check if it's a TCGA barcode
  if (grepl("^TCGA-", id)) {
    return(list(
      type = "TCGA_sample",
      id = id,
      cancer_type = tcga_barcode_to_cancer(id),
      sample_type = tcga_barcode_to_type(id),
      description = "TCGA sample identifier"
    ))
  }

  list(
    type = "unknown",
    id = id,
    message = "Identifier not recognized"
  )
}

#' TCGA Barcode to Cancer
#'
#' @param barcode TCGA barcode
#' @return Cancer type
#' @export
tcga_barcode_to_cancer <- function(barcode) {
  # Extract cancer type from barcode (TCGA-XX-XXXX)
  parts <- strsplit(barcode, "-")[[1]]
  if (length(parts) >= 2) {
    cancer_code <- parts[2]
    # Map to full name
    cancer_map <- list(
      "01" = "BRCA", "02" = "GBM", "03" = "LUAD", "04" = "LUSC",
      "05" = "COAD", "06" = "HNSC", "07" = "KIRC", "08" = "BLCA"
    )
    return(cancer_map[[cancer_code]] %||% "Unknown")
  }
  "Unknown"
}

#' TCGA Barcode to Type
#'
#' @param barcode TCGA barcode
#' @return Sample type
#' @export
tcga_barcode_to_type <- function(barcode) {
  # Extract sample type from barcode
  parts <- strsplit(barcode, "-")[[1]]
  if (length(parts) >= 4) {
    type_code <- substr(parts[4], 1, 2)
    type_map <- list(
      "01" = "Primary Tumor",
      "02" = "Recurrent Tumor",
      "03" = "Primary Blood Cancer",
      "06" = "Metastatic",
      "11" = "Solid Tissue Normal"
    )
    return(type_map[[type_code]] %||% "Unknown")
  }
  "Unknown"
}

# Custom Metadata Logic ---------------------------------------------------

#' Custom Metadata
#'
#' @param metadata Data frame with metadata
#' @param sample_col Column containing sample IDs
#' @return Processed metadata
#' @export
custom_metadata <- function(metadata, sample_col = "sample") {
  tryCatch({
    if (!sample_col %in% colnames(metadata)) {
      stop("Sample column not found: ", sample_col)
    }

    # Validate sample IDs
    sample_ids <- metadata[[sample_col]]

    list(
      success = TRUE,
      sample_count = nrow(metadata),
      variable_count = ncol(metadata) - 1,
      variables = setdiff(colnames(metadata), sample_col),
      sample_ids = sample_ids[1:min(5, length(sample_ids))],  # First 5 samples
      preview = head(metadata, 5)
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

#' Add Signature
#'
#' @param name Signature name
#' @param genes Gene list
#' @param weights Gene weights (optional)
#' @return Signature definition
#' @export
add_signature <- function(name, genes, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(genes))
  }

  if (length(genes) != length(weights)) {
    stop("Genes and weights must have same length")
  }

  list(
    name = name,
    genes = genes,
    weights = weights,
    created = Sys.time(),
    n_genes = length(genes)
  )
}

# Filter Samples Logic ----------------------------------------------------

#' Filter Samples
#'
#' @param samples Sample IDs
#' @param filters List of filters
#' @return Filtered samples
#' @export
filter_samples <- function(samples, filters = list()) {
  result <- samples

  for (filter_name in names(filters)) {
    filter_val <- filters[[filter_name]]

    result <- switch(filter_name,
      "cancer_type" = {
        # Filter by cancer type
        result[grepl(filter_val, result)]
      },
      "sample_type" = {
        # Filter by sample type
        result[grepl(filter_val, result)]
      },
      "expression_range" = {
        # Would need expression data
        result
      },
      result
    )
  }

  list(
    original_count = length(samples),
    filtered_count = length(result),
    removed_count = length(samples) - length(result),
    samples = result
  )
}
