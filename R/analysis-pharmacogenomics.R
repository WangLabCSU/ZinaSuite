#' PharmacoGenomics Analysis Functions
#'
#' @description
#' Functions for drug sensitivity and pharmacogenomic analysis using CCLE and GDSC data.
#'
#' @name pharmacogenomics-analysis
NULL

#' Query Drug Sensitivity Data
#'
#' @description
#' Query drug sensitivity data (IC50, AUC, etc.) from CCLE or GDSC databases.
#'
#' @param drug_name Drug name or identifier
#' @param metric Sensitivity metric: "ic50", "auc", "ec50"
#' @param source Data source: "ccle" (default), "gdsc"
#' @return Named vector of drug sensitivity values
#' @export
#'
#' @examples
#' \dontrun{
#' # Query cisplatin IC50 values from CCLE
#' cisplatin_ic50 <- query_drug_sensitivity("cisplatin", metric = "ic50")
#' summary(cisplatin_ic50)
#' }
query_drug_sensitivity <- function(drug_name,
                                   metric = c("ic50", "auc", "ec50"),
                                   source = c("ccle", "gdsc")) {
  metric <- match.arg(metric)
  source <- match.arg(source)

  # Map source to dataset
  dataset <- switch(source,
    "ccle" = "ccle\xC3\x82\xC2\xB1Drug_AUC",
    "gdsc" = "GDSC_v7_AUC"
  )

  # Query drug sensitivity from Xena
  xe <- UCSCXenaTools::XenaData
  host_url <- unique(xe$XenaHosts[xe$XenaHostNames == "publicHub"])

  result <- UCSCXenaTools::fetch_dense_values(
    host = host_url[1],
    dataset = dataset,
    identifiers = drug_name,
    check = FALSE
  )

  if (is.null(result) || nrow(result) == 0) {
    stop("Drug not found: ", drug_name)
  }

  # Convert to named vector
  sensitivity <- as.numeric(result[1, ])
  names(sensitivity) <- colnames(result)

  sensitivity
}

#' Analyze Drug-Gene Correlation
#'
#' @description
#' Analyze correlation between gene expression and drug sensitivity.
#'
#' @param gene Gene symbol
#' @param drug Drug name
#' @param metric Sensitivity metric
#' @param cor_method Correlation method
#' @return List with correlation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze correlation between EGFR expression and erlotinib sensitivity
#' result <- analyze_drug_gene_cor("EGFR", "erlotinib")
#' print(result$correlation)
#' }
analyze_drug_gene_cor <- function(gene,
                                  drug,
                                  metric = c("ic50", "auc"),
                                  cor_method = c("spearman", "pearson", "kendall")) {
  metric <- match.arg(metric)
  cor_method <- match.arg(cor_method)

  # Query gene expression and drug sensitivity
  gene_expr <- query_gene_expression(gene, source = "ccle")
  drug_sens <- query_drug_sensitivity(drug, metric = metric, source = "ccle")

  # Find common cell lines
  common_cells <- intersect(names(gene_expr), names(drug_sens))

  if (length(common_cells) < 10) {
    stop("Insufficient common cell lines for correlation analysis")
  }

  # Prepare data
  df <- data.frame(
    Gene = as.numeric(gene_expr[common_cells]),
    Drug = as.numeric(drug_sens[common_cells]),
    CellLine = common_cells,
    stringsAsFactors = FALSE
  )

  # Remove NA values
  df <- df[stats::complete.cases(df), ]

  if (nrow(df) < 10) {
    stop("Insufficient complete cases for correlation analysis")
  }

  # Calculate correlation
  cor_result <- analyze_correlation(df$Gene, df$Drug, method = cor_method)

  list(
    gene = gene,
    drug = drug,
    metric = metric,
    correlation = cor_result,
    data = df,
    n_samples = nrow(df)
  )
}

#' Batch Drug-Gene Correlation Analysis
#'
#' @description
#' Analyze correlation between multiple genes and a drug.
#'
#' @param genes Vector of gene symbols
#' @param drug Drug name
#' @param metric Sensitivity metric
#' @param cor_method Correlation method
#' @param n_workers Number of parallel workers
#' @param .progress Show progress bar
#' @return Data frame with correlation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze multiple genes with cisplatin
#' genes <- c("TP53", "BRCA1", "BRCA2", "ATM")
#' results <- analyze_drug_gene_batch(genes, "cisplatin")
#' head(results)
#' }
analyze_drug_gene_batch <- function(genes,
                                    drug,
                                    metric = "auc",
                                    cor_method = "spearman",
                                    n_workers = 4,
                                    .progress = TRUE) {
  # Query drug sensitivity once
  drug_sens <- query_drug_sensitivity(drug, metric = metric, source = "ccle")

  # Start mirai daemons
  mirai::daemons(n_workers)
  on.exit(mirai::daemons(0), add = TRUE)

  # Prepare argument list
  arg_list <- lapply(genes, function(gene) {
    list(
      gene = gene,
      drug_sens = drug_sens,
      metric = metric,
      cor_method = cor_method
    )
  })

  # Parallel analysis
  results <- mirai::mirai_map(
    arg_list,
    function(args) {
      tryCatch({
        # Query gene expression from CCLE
        xe <- UCSCXenaTools::XenaData
        host_url <- unique(xe$XenaHosts[xe$XenaHostNames == "publicHub"])
        dataset <- "ccle\xC3\x82\xC2\xB1CCLE_RNAseq_genes_rpkm_20180929"

        gene_result <- UCSCXenaTools::fetch_dense_values(
          host = host_url[1],
          dataset = dataset,
          identifiers = args$gene,
          check = FALSE,
          use_probeMap = TRUE
        )

        if (is.null(gene_result) || nrow(gene_result) == 0) {
          return(list(
            gene = args$gene,
            cor = NA,
            pvalue = NA,
            n = 0,
            error = "Gene not found"
          ))
        }

        gene_expr <- as.numeric(gene_result[1, ])
        names(gene_expr) <- colnames(gene_result)

        # Find common cell lines
        common_cells <- intersect(names(gene_expr), names(args$drug_sens))

        if (length(common_cells) < 5) {
          return(list(
            gene = args$gene,
            cor = NA,
            pvalue = NA,
            n = length(common_cells),
            error = "Insufficient common cell lines"
          ))
        }

        # Calculate correlation
        x <- as.numeric(gene_expr[common_cells])
        y <- as.numeric(args$drug_sens[common_cells])

        valid <- stats::complete.cases(x, y)
        if (sum(valid) < 5) {
          return(list(
            gene = args$gene,
            cor = NA,
            pvalue = NA,
            n = sum(valid),
            error = "Insufficient complete cases"
          ))
        }

        cor_result <- stats::cor.test(x[valid], y[valid], method = args$cor_method)

        list(
          gene = args$gene,
          cor = unname(cor_result$estimate),
          pvalue = cor_result$p.value,
          n = sum(valid),
          error = NULL
        )
      }, error = function(e) {
        list(
          gene = args$gene,
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
      gene = r$gene,
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

#' Get Available Drugs
#'
#' @description
#' Get list of available drugs for pharmacogenomic analysis.
#'
#' @param source Data source
#' @return Character vector of drug names
#' @export
#'
#' @examples
#' \dontrun{
#' drugs <- get_available_drugs()
#' head(drugs)
#' }
get_available_drugs <- function(source = c("ccle", "gdsc")) {
  source <- match.arg(source)

  # Common cancer drugs
  ccle_drugs <- c(
    "cisplatin",
    "carboplatin",
    "oxaliplatin",
    "paclitaxel",
    "docetaxel",
    "vinblastine",
    "vincristine",
    "doxorubicin",
    "epirubicin",
    "etoposide",
    "gemcitabine",
    "5-fluorouracil",
    "capecitabine",
    "methotrexate",
    "cyclophosphamide",
    "ifosfamide",
    "bleomycin",
    "mitomycin",
    "topotecan",
    "irinotecan",
    "erlotinib",
    "gefitinib",
    "afatinib",
    "osimertinib",
    "lapatinib",
    "trastuzumab",
    "cetuximab",
    "bevacizumab",
    "sunitinib",
    "sorafenib",
    "pazopanib",
    "axitinib",
    "regorafenib",
    "imatinib",
    "dasatinib",
    "nilotinib",
    "ponatinib",
    "bosutinib",
    "vemurafenib",
    "dabrafenib",
    "trametinib",
    "cobimetinib",
    "selumetinib",
    "crizotinib",
    "ceritinib",
    "alectinib",
    "brigatinib",
    "lorlatinib",
    "palbociclib",
    "ribociclib",
    "abemaciclib",
    "olaparib",
    "rucaparib",
    "niraparib",
    "talazoparib",
    "everolimus",
    "temsirolimus",
    "bortezomib",
    "carfilzomib",
    "ixazomib",
    "lenalidomide",
    "pomalidomide",
    "thalidomide",
    "venetoclax",
    "navitoclax",
    "azacitidine",
    "decitabine",
    "vorinostat",
    "romidepsin",
    "belinostat",
    "panobinostat",
    "entinostat",
    "tamoxifen",
    "raloxifene",
    "fulvestrant",
    "anastrozole",
    "letrozole",
    "exemestane",
    "bicalutamide",
    "enzalutamide",
    "abiraterone",
    "prednisone",
    "dexamethasone"
  )

  gdsc_drugs <- c(
    "Camptothecin",
    "Cisplatin",
    "Cyclophosphamide",
    "Docetaxel",
    "Doxorubicin",
    "Etoposide",
    "Gemcitabine",
    "Irinotecan",
    "Methotrexate",
    "Mitomycin",
    "Oxaliplatin",
    "Paclitaxel",
    "Topotecan",
    "Vinblastine",
    "Vincristine",
    "Vinorelbine",
    "5-FU",
    "Bleomycin",
    "Carboplatin",
    "Chlorambucil",
    "Dacarbazine",
    "Fludarabine",
    "Hydroxyurea",
    "Melphalan",
    "Mercaptopurine",
    "Thioguanine",
    "Trabectedin"
  )

  switch(source,
    "ccle" = ccle_drugs,
    "gdsc" = gdsc_drugs,
    ccle_drugs
  )
}

#' Analyze Drug Response by Mutation Status
#'
#' @description
#' Compare drug sensitivity between mutant and wild-type cell lines.
#'
#' @param drug Drug name
#' @param gene Gene to check mutation status
#' @param metric Sensitivity metric
#' @return List with comparison results
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare olaparib sensitivity in BRCA1 mutant vs WT cell lines
#' result <- analyze_drug_mutation("olaparib", "BRCA1")
#' print(result$test)
#' }
analyze_drug_mutation <- function(drug,
                                  gene,
                                  metric = "auc") {
  # Query drug sensitivity
  drug_sens <- query_drug_sensitivity(drug, metric = metric, source = "ccle")

  # Query mutation data
  mut_data <- query_mutation(gene, source = "ccle")

  # Find common cell lines
  common_cells <- intersect(names(drug_sens), names(mut_data))

  if (length(common_cells) < 10) {
    stop("Insufficient common cell lines")
  }

  # Prepare data
  df <- data.frame(
    CellLine = common_cells,
    Sensitivity = as.numeric(drug_sens[common_cells]),
    Mutation = as.character(mut_data[common_cells]),
    stringsAsFactors = FALSE
  )

  # Create binary mutation status
  df$Mutated <- ifelse(df$Mutation == "WT", "Wild-type", "Mutated")

  # Statistical test
  wt_values <- df$Sensitivity[df$Mutated == "Wild-type"]
  mut_values <- df$Sensitivity[df$Mutated == "Mutated"]

  if (length(wt_values) < 3 || length(mut_values) < 3) {
    stop("Insufficient samples in one or both groups")
  }

  test_result <- stats::wilcox.test(wt_values, mut_values)

  list(
    drug = drug,
    gene = gene,
    metric = metric,
    data = df,
    test = test_result,
    summary = list(
      wt_median = stats::median(wt_values, na.rm = TRUE),
      mut_median = stats::median(mut_values, na.rm = TRUE),
      wt_n = length(wt_values),
      mut_n = length(mut_values)
    )
  )
}

#' Create Drug Sensitivity Profile
#'
#' @description
#' Create a comprehensive drug sensitivity profile for a cell line or sample.
#'
#' @param sample_id Sample or cell line identifier
#' @param drugs Vector of drug names (NULL for all available)
#' @param source Data source
#' @return Data frame with sensitivity profile
#' @export
#'
#' @examples
#' \dontrun{
#' # Get drug sensitivity profile for MCF7 cell line
#' profile <- create_drug_profile("MCF7")
#' head(profile)
#' }
create_drug_profile <- function(sample_id,
                                drugs = NULL,
                                source = "ccle") {
  if (is.null(drugs)) {
    drugs <- get_available_drugs(source)
  }

  # Query sensitivity for each drug
  profile <- lapply(drugs, function(drug) {
    tryCatch({
      sens <- query_drug_sensitivity(drug, metric = "auc", source = source)
      if (sample_id %in% names(sens)) {
        data.frame(
          Drug = drug,
          AUC = as.numeric(sens[sample_id]),
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }, error = function(e) NULL)
  })

  # Combine results
  profile_df <- do.call(rbind, profile)

  if (is.null(profile_df) || nrow(profile_df) == 0) {
    stop("No drug sensitivity data found for sample: ", sample_id)
  }

  # Sort by AUC (lower = more sensitive)
  profile_df <- profile_df[order(profile_df$AUC), ]
  rownames(profile_df) <- NULL

  profile_df
}
