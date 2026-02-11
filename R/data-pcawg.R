#' PCAWG Data Source
#'
#' @description
#' Functions for querying data from PCAWG (Pan-Cancer Analysis of Whole Genomes).
#'
#' @name pcawg-data
NULL

#' Query PCAWG Gene Expression
#'
#' @param gene Gene symbol
#' @param host PCAWG hub host (default: "pcawgHub")
#' @return Named numeric vector of gene expression values
#' @export
#'
#' @examples
#' \dontrun{
#' tp53_expr <- query_pcawg_expression("TP53")
#' head(tp53_expr)
#' }
query_pcawg_expression <- function(gene, host = "pcawgHub") {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
    stop("Package 'UCSCXenaTools' is required")
  }

  # PCAWG gene expression dataset
  dataset <- "PCAWG.rnaseq.transcript.expr.FPKM.tsv"

  result <- try_query_value(host, dataset, gene)

  if (is.null(result)) {
    warning("No data returned for gene: ", gene)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "unit") <- "FPKM"
  attr(values, "source") <- "PCAWG"

  values
}

#' Query PCAWG Mutation Status
#'
#' @param gene Gene symbol
#' @param host PCAWG hub host (default: "pcawgHub")
#' @return Named numeric vector (0 = wild-type, 1 = mutated)
#' @export
#'
#' @examples
#' \dontrun{
#' tp53_mut <- query_pcawg_mutation("TP53")
#' table(tp53_mut)
#' }
query_pcawg_mutation <- function(gene, host = "pcawgHub") {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
    stop("Package 'UCSCXenaTools' is required")
  }

  # PCAWG mutation dataset
  dataset <- "PCAWG.consensus.2016.12.16.somaticSNV_subset_xena.tsv"

  result <- try_query_value(host, dataset, gene, use_probeMap = FALSE)

  if (is.null(result)) {
    warning("No data returned for gene: ", gene)
    return(NULL)
  }

  # Convert to binary mutation status
  values <- as.numeric(result[1, ])
  values <- ifelse(is.na(values), 0, 1)
  names(values) <- colnames(result)
  attr(values, "source") <- "PCAWG"

  values
}

#' Query PCAWG Copy Number Variation
#'
#' @param gene Gene symbol
#' @param host PCAWG hub host (default: "pcawgHub")
#' @return Named numeric vector of CNV values
#' @export
#'
#' @examples
#' \dontrun{
#' myc_cnv <- query_pcawg_cnv("MYC")
#' summary(myc_cnv)
#' }
query_pcawg_cnv <- function(gene, host = "pcawgHub") {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
    stop("Package 'UCSCXenaTools' is required")
  }

  # PCAWG CNV dataset
  dataset <- "PCAWG.consensus.2016.12.16.somaticCNV.tsv"

  result <- try_query_value(host, dataset, gene, use_probeMap = FALSE)

  if (is.null(result)) {
    warning("No data returned for gene: ", gene)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "source") <- "PCAWG"

  values
}

#' PCAWGData R6 Class
#'
#' @description
#' R6 class for accessing PCAWG (Pan-Cancer Analysis of Whole Genomes) data.
#'
#' @export
#' @examples
#' \dontrun{
#' # Create PCAWG data instance
#' pcawg <- PCAWGData$new()
#'
#' # Query gene expression
#' tp53_expr <- pcawg$get_gene_expression("TP53")
#'
#' # Query mutation status
#' tp53_mut <- pcawg$get_mutation_status("TP53")
#' }
PCAWGData <- R6::R6Class(
  "PCAWGData",
  inherit = DataSource,

  public = list(
    #' @description
    #' Initialize a new PCAWGData instance
    initialize = function() {
      super$initialize("PCAWG", "pcawgHub")
      private$init_dataset_map()
    },

    #' @description
    #' Get gene expression data
    #' @param gene Gene symbol
    #' @return Named numeric vector
    get_gene_expression = function(gene) {
      private$query_with_cache(gene, "expression", private$dataset_map$expression)
    },

    #' @description
    #' Get fusion data
    #' @param gene Gene symbol
    #' @return Named numeric vector (0/1)
    get_fusion = function(gene) {
      private$query_with_cache(gene, "fusion", private$dataset_map$fusion, use_probeMap = FALSE)
    },

    #' @description
    #' Get miRNA expression data
    #' @param mirna miRNA ID (e.g., "hsa-miR-21-5p")
    #' @param norm Normalization method: "TMM" (default) or "UQ"
    #' @return Named numeric vector
    get_mirna = function(mirna, norm = c("TMM", "UQ")) {
      norm <- match.arg(norm)
      dataset <- private$dataset_map[[paste0("mirna_", tolower(norm))]]
      private$query_with_cache(mirna, "miRNA", dataset)
    },

    #' @description
    #' Get promoter activity data
    #' @param promoter Promoter identifier (e.g., "prmtr.10000") or gene symbol
    #' @param type Promoter activity type: "raw" (default), "relative", or "outlier"
    #' @return Named numeric vector
    get_promoter = function(promoter, type = c("raw", "relative", "outlier")) {
      type <- match.arg(type)
      dataset <- private$dataset_map[[paste0("promoter_", type)]]
      private$query_with_cache(promoter, "promoter", dataset, use_probeMap = FALSE)
    },

    #' @description
    #' Get APOBEC mutagenesis data
    #' @param identifier APOBEC identifier
    #' @return Named numeric vector
    get_apobec = function(identifier = c("tCa_MutLoad_MinEstimate", "APOBECtCa_enrich",
                                         "A3A_or_A3B", "APOBEC_tCa_enrich_quartile")) {
      identifier <- match.arg(identifier)
      private$query_with_cache(identifier, "APOBEC", private$dataset_map$apobec, use_probeMap = FALSE)
    },

    #' @description
    #' Query data by type
    #' @param identifier Gene symbol or other identifier
    #' @param data_type Type of data: "expression", "fusion", "miRNA", "promoter", "APOBEC"
    #' @param ... Additional parameters (e.g., norm for miRNA)
    #' @return Named numeric vector
    query = function(identifier, data_type = "expression", ...) {
      switch(data_type,
        "expression" = self$get_gene_expression(identifier),
        "fusion" = self$get_fusion(identifier),
        "miRNA" = self$get_mirna(identifier, ...),
        "promoter" = self$get_promoter(identifier, ...),
        "APOBEC" = self$get_apobec(identifier),
        stop("Unknown data type: ", data_type)
      )
    }
  ),

  private = list(
    dataset_map = NULL,

    init_dataset_map = function() {
      # Use same dataset names as UCSCXenaShiny
      private$dataset_map <- list(
        expression = "tophat_star_fpkm_uq.v2_aliquot_gl.sp.log",
        fusion = "pcawg3_fusions_PKU_EBI.gene_centric.sp.xena",
        mirna_tmm = "x3t2m1.mature.TMM.mirna.matrix.log",
        mirna_uq = "x3t2m1.mature.UQ.mirna.matrix.log",
        promoter_raw = "rawPromoterActivity.sp",
        promoter_relative = "relativePromoterActivity.sp",
        promoter_outlier = "promoterCentricTable_0.2_1.0.sp",
        apobec = "MAF_Aug31_2016_sorted_A3A_A3B_comparePlus.sp"
      )
    },

    query_with_cache = function(identifier, data_type, dataset, use_probeMap = TRUE) {
      cache_key <- paste("pcawg", data_type, identifier, sep = "_")

      # Check cache
      cached <- private$cache_manager$get(cache_key)
      if (!is.null(cached)) {
        return(cached)
      }

      # Query from Xena with retry
      result <- NULL
      max_try <- 3
      for (i in seq_len(max_try)) {
        result <- tryCatch({
          private$fetch_from_xena(identifier, dataset, use_probeMap)
        }, error = function(e) {
          if (i == max_try) {
            warning("Failed to fetch ", identifier, " from PCAWG after ", max_try, " attempts: ", conditionMessage(e))
            return(NULL)
          }
          Sys.sleep(0.5)
          NULL
        })
        if (!is.null(result)) break
      }

      if (is.null(result)) {
        warning("No data returned for ", identifier, " from PCAWG ", data_type)
        return(NULL)
      }

      # Convert to named vector
      values <- as.numeric(result[1, ])
      names(values) <- colnames(result)
      attr(values, "source") <- "PCAWG"
      attr(values, "type") <- data_type

      # Cache result
      private$cache_manager$set(cache_key, values)

      values
    },

    fetch_from_xena = function(identifier, dataset, use_probeMap = TRUE) {
      if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
        stop("Package 'UCSCXenaTools' is required")
      }

      xe <- UCSCXenaTools::XenaData
      host_url <- unique(xe$XenaHosts[xe$XenaHostNames == private$host])

      if (length(host_url) == 0) {
        stop("Host not found: ", private$host)
      }

      result <- UCSCXenaTools::fetch_dense_values(
        host = host_url[1],
        dataset = dataset,
        identifiers = identifier,
        check = FALSE,
        use_probeMap = use_probeMap
      )

      if (is.null(result) || nrow(result) == 0) {
        return(NULL)
      }

      result
    }
  )
)
