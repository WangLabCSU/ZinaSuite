#' Pathway Analysis Functions
#'
#' @description
#' Functions for pathway enrichment and activity analysis.
#'
#' @name pathway-analysis
NULL

#' Query Pathway Activity
#'
#' @description
#' Query pathway activity scores from UCSC Xena databases.
#'
#' @param pathway_name Pathway name or identifier
#' @param source Data source: "tcga" (default), "pcawg"
#' @param dataset Specific dataset name (optional)
#' @return Named vector of pathway activity scores
#' @export
#'
#' @examples
#' \dontrun{
#' # Query PI3K-Akt pathway activity
#' pi3k_scores <- query_pathway_activity("PI3K-Akt signaling pathway")
#' summary(pi3k_scores)
#' }
query_pathway_activity <- function(pathway_name,
                                   source = c("tcga", "pcawg"),
                                   dataset = NULL) {
  source <- match.arg(source)

  # Map source to host and default dataset
  host <- switch(source,
    "tcga" = "toilHub",
    "pcawg" = "pcawgHub"
  )

  default_dataset <- switch(source,
    "tcga" = "TcgaTargetGtex_rsem_gene_tpm",
    "pcawg" = "PCAWG_specimen"
  )

  dataset <- dataset %||% default_dataset

  # Query pathway activity from Xena
  xe <- UCSCXenaTools::XenaData
  host_url <- unique(xe$XenaHosts[xe$XenaHostNames == host])

  result <- UCSCXenaTools::fetch_dense_values(
    host = host_url[1],
    dataset = dataset,
    identifiers = pathway_name,
    check = FALSE
  )

  if (is.null(result) || nrow(result) == 0) {
    stop("Pathway not found: ", pathway_name)
  }

  # Convert to named vector
  activity <- as.numeric(result[1, ])
  names(activity) <- colnames(result)

  activity
}

#' Analyze Pathway Correlation with Gene
#'
#' @description
#' Analyze correlation between pathway activity and gene expression.
#'
#' @param gene Gene symbol
#' @param pathway Pathway name
#' @param source Data source
#' @param method Correlation method
#' @return List with correlation results and plot data
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze correlation between TP53 and apoptosis pathway
#' result <- analyze_pathway_gene_cor("TP53", "Apoptosis")
#' print(result$correlation)
#' }
analyze_pathway_gene_cor <- function(gene,
                                     pathway,
                                     source = "tcga",
                                     method = c("spearman", "pearson", "kendall")) {
  method <- match.arg(method)

  # Query gene expression and pathway activity
  gene_expr <- query_gene_expression(gene, source = source)
  pathway_act <- query_pathway_activity(pathway, source = source)

  # Find common samples
  common_samples <- intersect(names(gene_expr), names(pathway_act))

  if (length(common_samples) < 10) {
    stop("Insufficient common samples for correlation analysis")
  }

  # Prepare data
  df <- data.frame(
    Gene = as.numeric(gene_expr[common_samples]),
    Pathway = as.numeric(pathway_act[common_samples]),
    stringsAsFactors = FALSE
  )

  # Calculate correlation
  cor_result <- analyze_correlation(df$Gene, df$Pathway, method = method)

  list(
    gene = gene,
    pathway = pathway,
    correlation = cor_result,
    data = df,
    n_samples = length(common_samples)
  )
}

#' Get Available Pathways
#'
#' @description
#' Get list of available pathways for analysis.
#'
#' @param source Data source
#' @param category Pathway category (optional)
#' @return Character vector of pathway names
#' @export
#'
#' @examples
#' \dontrun{
#' pathways <- get_available_pathways()
#' head(pathways)
#' }
get_available_pathways <- function(source = "tcga",
                                   category = c("all", "hallmark", "kegg", "reactome", "go")) {
  category <- match.arg(category)

  # Common cancer-related pathways
  hallmark_pathways <- c(
    "HALLMARK_APOPTOSIS",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_HYPOXIA",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_ANDROGEN_RESPONSE",
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_UV_RESPONSE_UP",
    "HALLMARK_UV_RESPONSE_DN",
    "HALLMARK_XENOBIOTIC_METABOLISM",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_PROTEIN_SECRETION",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_APICAL_SURFACE",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_HEME_METABOLISM",
    "HALLMARK_COAGULATION",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "HALLMARK_HEDGEHOG_SIGNALING",
    "HALLMARK_NOTCH_SIGNALING",
    "HALLMARK_ANGIOGENESIS"
  )

  kegg_pathways <- c(
    "KEGG_CELL_CYCLE",
    "KEGG_APOPTOSIS",
    "KEGG_P53_SIGNALING_PATHWAY",
    "KEGG_PI3K_AKT_SIGNALING_PATHWAY",
    "KEGG_MAPK_SIGNALING_PATHWAY",
    "KEGG_JAK_STAT_SIGNALING_PATHWAY",
    "KEGG_TGF_BETA_SIGNALING_PATHWAY",
    "KEGG_WNT_SIGNALING_PATHWAY",
    "KEGG_NOTCH_SIGNALING_PATHWAY",
    "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
    "KEGG_VEGF_SIGNALING_PATHWAY",
    "KEGG_ERBB_SIGNALING_PATHWAY",
    "KEGG_INSULIN_SIGNALING_PATHWAY",
    "KEGG_MTOR_SIGNALING_PATHWAY",
    "KEGG_FOCAL_ADHESION",
    "KEGG_ECM_RECEPTOR_INTERACTION",
    "KEGG_ADHERENS_JUNCTION",
    "KEGG_TIGHT_JUNCTION",
    "KEGG_GAP_JUNCTION",
    "KEGG_CELL_ADHESION_MOLECULES",
    "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
    "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
    "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
    "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
    "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
    "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY",
    "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
    "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
    "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
    "KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY",
    "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
    "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION",
    "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION",
    "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"
  )

  switch(category,
    "all" = c(hallmark_pathways, kegg_pathways),
    "hallmark" = hallmark_pathways,
    "kegg" = kegg_pathways,
    hallmark_pathways
  )
}

#' Batch Pathway Analysis
#'
#' @description
#' Analyze correlation between a gene and multiple pathways.
#'
#' @param gene Gene symbol
#' @param pathways Vector of pathway names
#' @param source Data source
#' @param method Correlation method
#' @param n_workers Number of parallel workers
#' @param .progress Show progress bar
#' @return Data frame with correlation results
#' @export
#'
#' @examples
#' \dontrun{
#' pathways <- get_available_pathways(category = "hallmark")
#' results <- analyze_pathway_batch("TP53", pathways[1:10])
#' head(results)
#' }
analyze_pathway_batch <- function(gene,
                                  pathways,
                                  source = "tcga",
                                  method = "spearman",
                                  n_workers = 4,
                                  .progress = TRUE) {
  # Query gene expression once
  gene_expr <- query_gene_expression(gene, source = source)

  # Start mirai daemons
  mirai::daemons(n_workers)
  on.exit(mirai::daemons(0), add = TRUE)

  # Prepare argument list
  arg_list <- lapply(pathways, function(pathway) {
    list(
      pathway = pathway,
      gene_expr = gene_expr,
      source = source,
      method = method
    )
  })

  # Parallel analysis
  results <- mirai::mirai_map(
    arg_list,
    function(args) {
      tryCatch({
        # Query pathway activity
        xe <- UCSCXenaTools::XenaData
        host_url <- unique(xe$XenaHosts[xe$XenaHostNames == "toilHub"])
        dataset <- "TcgaTargetGtex_rsem_gene_tpm"

        pathway_result <- UCSCXenaTools::fetch_dense_values(
          host = host_url[1],
          dataset = dataset,
          identifiers = args$pathway,
          check = FALSE
        )

        if (is.null(pathway_result) || nrow(pathway_result) == 0) {
          return(list(
            pathway = args$pathway,
            cor = NA,
            pvalue = NA,
            n = 0,
            error = "Pathway not found"
          ))
        }

        pathway_act <- as.numeric(pathway_result[1, ])
        names(pathway_act) <- colnames(pathway_result)

        # Find common samples
        common_samples <- intersect(names(args$gene_expr), names(pathway_act))

        if (length(common_samples) < 10) {
          return(list(
            pathway = args$pathway,
            cor = NA,
            pvalue = NA,
            n = length(common_samples),
            error = "Insufficient common samples"
          ))
        }

        # Calculate correlation
        x <- as.numeric(args$gene_expr[common_samples])
        y <- as.numeric(pathway_act[common_samples])

        cor_result <- stats::cor.test(x, y, method = args$method)

        list(
          pathway = args$pathway,
          cor = unname(cor_result$estimate),
          pvalue = cor_result$p.value,
          n = length(common_samples),
          error = NULL
        )
      }, error = function(e) {
        list(
          pathway = args$pathway,
          cor = NA,
          pvalue = NA,
          n = 0,
          error = conditionMessage(e)
        )
      })
    },
    .progress = .progress
  )

  # Collect results
  results_list <- results[]

  # Convert to data frame
  results_df <- do.call(rbind, lapply(results_list, function(r) {
    data.frame(
      pathway = r$pathway,
      cor = r$cor,
      pvalue = r$pvalue,
      n = r$n,
      error = r$error,
      stringsAsFactors = FALSE
    )
  }))

  # Adjust p-values
  valid_pvalues <- !is.na(results_df$pvalue)
  if (any(valid_pvalues)) {
    results_df$padj <- NA
    results_df$padj[valid_pvalues] <- stats::p.adjust(
      results_df$pvalue[valid_pvalues],
      method = "fdr"
    )
  }

  # Sort by absolute correlation
  results_df <- results_df[order(-abs(results_df$cor)), ]
  rownames(results_df) <- NULL

  results_df
}
