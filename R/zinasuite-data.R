#' ZinaSuite Data Query Functions
#'
#' @description
#' Core data query functions for ZinaSuite, based on UCSCXenaShiny implementation.
#' These functions provide reliable access to UCSC Xena data hubs.

#' Query Gene Expression from TCGA/GTEx
#'
#' @param identifier Gene symbol (e.g., "TP53", "BRCA1")
#' @param norm Normalization method: "tpm" (default), "fpkm", or "nc" (norm count)
#' @param host Xena hub host (default: "toilHub")
#' @return Named numeric vector of expression values
#' @export
#'
#' @examples
#' \donttest{
#' # Query TP53 expression
#' expr <- get_pancan_gene_value("TP53")
#' head(expr)
#'
#' # Query with FPKM normalization
#' expr <- get_pancan_gene_value("BRCA1", norm = "fpkm")
#' }
get_pancan_gene_value <- function(identifier, norm = c("tpm", "fpkm", "nc"), host = "toilHub") {
  norm <- match.arg(norm)

  # Determine dataset based on normalization
  dataset <- switch(norm,
    "tpm" = "TcgaTargetGtex_rsem_gene_tpm",
    "fpkm" = "TcgaTargetGtex_rsem_gene_fpkm",
    "nc" = "TcgaTargetGtex_RSEM_Hugo_norm_count"
  )

  # Query data
  result <- try_query_value(host, dataset, identifier)

  if (is.null(result)) {
    warning("No data returned for ", identifier)
    return(NULL)
  }

  # Convert to named vector
  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)

  # Add unit attribute
  unit <- switch(norm,
    "tpm" = "log2(tpm+0.001)",
    "fpkm" = "log2(fpkm+0.001)",
    "nc" = "log2(norm_count+1)"
  )
  attr(values, "unit") <- unit

  values
}

#' Query Protein Expression (RPPA)
#'
#' @param identifier Protein identifier
#' @param host Xena hub host (default: "toilHub")
#' @return Named numeric vector of protein expression values
#' @export
get_pancan_protein_value <- function(identifier, host = "toilHub") {
  dataset <- "TCGA-RPPA-pancan-clean.xena"
  result <- try_query_value(host, dataset, identifier)

  if (is.null(result)) {
    warning("No data returned for ", identifier)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "unit") <- "RPPA"

  values
}

#' Query Mutation Status
#'
#' @param identifier Gene symbol
#' @param host Xena hub host (default: "toilHub")
#' @return Named numeric vector (0 = wild-type, 1 = mutated)
#' @export
get_pancan_mutation_status <- function(identifier, host = "pancanAtlasHub") {
  dataset <- "mc3.v0.2.8.PUBLIC.nonsilentGene.xena"
  # Mutation data doesn't use probeMap
  result <- try_query_value(host, dataset, identifier, use_probeMap = FALSE)

  if (is.null(result)) {
    warning("No data returned for ", identifier)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)

  values
}

#' Query Copy Number Variation (CNV)
#'
#' @param identifier Gene symbol
#' @param use_thresholded Use thresholded data (default: TRUE)
#' @param host Xena hub host (default: "toilHub")
#' @return Named numeric vector of GISTIC2 scores
#' @export
get_pancan_cn_value <- function(identifier, use_thresholded = TRUE, host = "tcgaHub") {
  dataset <- if (use_thresholded) {
    "TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
  } else {
    "TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes"
  }

  result <- try_query_value(host, dataset, identifier)

  if (is.null(result)) {
    warning("No data returned for ", identifier)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "unit") <- "GISTIC2"

  values
}

#' Query Methylation Data
#'
#' @param identifier Gene symbol or probe ID
#' @param type Methylation array type: "450K" or "27K"
#' @param aggr Aggregation method for multiple probes
#' @param host Xena hub host (default: "toilHub")
#' @return Named numeric vector of beta values
#' @export
get_pancan_methylation_value <- function(identifier,
                                          type = c("450K", "27K"),
                                          aggr = c("mean", "median", "Q0", "Q25", "Q50", "Q75", "Q100"),
                                          host = "toilHub") {
  type <- match.arg(type)
  aggr <- match.arg(aggr)

  dataset <- if (type == "450K") {
    "TCGA.PANCAN.sampleMap_HumanMethylation450"
  } else {
    "TCGA.PANCAN.sampleMap_HumanMethylation27"
  }

  # For gene symbols, we need to aggregate multiple probes
  if (aggr != "mean" || !grepl("^cg", identifier)) {
    # Get probe mapping (simplified - would need actual probe map)
    result <- try_query_value(host, dataset, identifier)
  } else {
    result <- try_query_value(host, dataset, identifier)
  }

  if (is.null(result)) {
    warning("No data returned for ", identifier)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "unit") <- "beta"

  values
}

#' Query miRNA Expression
#'
#' @param identifier miRNA ID (e.g., "hsa-miR-21-5p")
#' @param host Xena hub host (default: "toilHub")
#' @return Named numeric vector of miRNA expression
#' @export
get_pancan_mirna_value <- function(identifier, host = "toilHub") {
  dataset <- "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
  result <- try_query_value(host, dataset, identifier)

  if (is.null(result)) {
    warning("No data returned for ", identifier)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "unit") <- "log2(norm_count+1)"

  values
}

#' Query Transcript Expression
#'
#' @param identifier Transcript ID or gene symbol
#' @param norm Normalization method: "tpm" or "fpkm"
#' @param host Xena hub host (default: "toilHub")
#' @return Named numeric vector of transcript expression
#' @export
get_pancan_transcript_value <- function(identifier, norm = c("tpm", "fpkm"), host = "toilHub") {
  norm <- match.arg(norm)

  dataset <- switch(norm,
    "tpm" = "TcgaTargetGtex_rsem_isoform_tpm",
    "fpkm" = "TcgaTargetGtex_rsem_isoform_fpkm"
  )

  result <- try_query_value(host, dataset, identifier)

  if (is.null(result)) {
    warning("No data returned for ", identifier)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "unit") <- paste0("log2(", norm, "+0.001)")

  values
}

#' Try Query Value with Retry Logic
#'
#' @param host Xena hub host
#' @param dataset Dataset name
#' @param identifier Gene/protein identifier
#' @param max_try Maximum retry attempts
#' @return Data matrix or NULL if failed
#' @keywords internal
try_query_value <- function(host, dataset, identifier, max_try = 5L, use_probeMap = TRUE) {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
    stop("Package 'UCSCXenaTools' is required.")
  }

  # Get host URL from UCSCXenaTools
  xe <- UCSCXenaTools::XenaData
  host_url <- unique(xe$XenaHosts[xe$XenaHostNames == host])

  if (length(host_url) == 0) {
    stop("Invalid host: ", host)
  }

  # Try to fetch data with retry
  # use_probeMap = TRUE allows querying by gene symbol for most datasets
  for (i in seq_len(max_try)) {
    result <- tryCatch({
      UCSCXenaTools::fetch_dense_values(
        host = host_url[1],
        dataset = dataset,
        identifiers = identifier,
        check = FALSE,
        use_probeMap = use_probeMap
      )
    }, error = function(e) {
      if (i == max_try) {
        warning("Failed after ", max_try, " attempts: ", conditionMessage(e))
        return(NULL)
      }
      Sys.sleep(0.5)
      NULL
    })

    if (!is.null(result) && nrow(result) > 0) {
      return(result)
    }
  }

  NULL
}

#' Unified Query Interface
#'
#' @param molecule Molecule identifier
#' @param data_type Data type: "mRNA", "protein", "mutation", "cnv", "methylation", "miRNA", "transcript"
#' @param database Database: "toil" (TCGA/GTEx), "ccle", "pcawg"
#' @param ... Additional parameters
#' @return Named numeric vector
#' @export
#'
#' @examples
#' \donttest{
#' # Query gene expression
#' expr <- query_pancan_value("TP53", data_type = "mRNA")
#'
#' # Query mutation
#' mut <- query_pancan_value("TP53", data_type = "mutation")
#'
#' # Query CNV
#' cnv <- query_pancan_value("MYC", data_type = "cnv")
#' }
query_pancan_value <- function(molecule,
                                data_type = c("mRNA", "protein", "mutation", "cnv", "methylation", "miRNA", "transcript"),
                                database = c("toil", "ccle", "pcawg"),
                                ...) {
  data_type <- match.arg(data_type)
  database <- match.arg(database)

  if (database == "toil") {
    result <- switch(data_type,
      "mRNA" = get_pancan_gene_value(molecule, ...),
      "protein" = get_pancan_protein_value(molecule, ...),
      "mutation" = get_pancan_mutation_status(molecule, ...),
      "cnv" = get_pancan_cn_value(molecule, ...),
      "methylation" = get_pancan_methylation_value(molecule, ...),
      "miRNA" = get_pancan_mirna_value(molecule, ...),
      "transcript" = get_pancan_transcript_value(molecule, ...),
      stop("Unknown data type: ", data_type)
    )
  } else if (database == "ccle") {
    # CCLE data functions would go here
    stop("CCLE database not yet implemented")
  } else if (database == "pcawg") {
    # PCAWG data functions would go here
    stop("PCAWG database not yet implemented")
  }

  result
}
