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
    #' Get mutation status
    #' @param gene Gene symbol
    #' @return Named numeric vector (0/1)
    get_mutation_status = function(gene) {
      private$query_with_cache(gene, "mutation", private$dataset_map$mutation, use_probeMap = FALSE)
    },

    #' @description
    #' Get CNV data
    #' @param gene Gene symbol
    #' @return Named numeric vector
    get_cnv = function(gene) {
      private$query_with_cache(gene, "cnv", private$dataset_map$cnv, use_probeMap = FALSE)
    },

    #' @description
    #' Query data by type
    #' @param identifier Gene symbol
    #' @param data_type Type of data
    #' @return Named numeric vector
    query = function(identifier, data_type = "expression") {
      switch(data_type,
        "expression" = self$get_gene_expression(identifier),
        "mutation" = self$get_mutation_status(identifier),
        "cnv" = self$get_cnv(identifier),
        stop("Unknown data type: ", data_type)
      )
    }
  ),

  private = list(
    dataset_map = NULL,

    init_dataset_map = function() {
      private$dataset_map <- list(
        expression = "PCAWG.rnaseq.transcript.expr.FPKM.tsv",
        mutation = "PCAWG.consensus.2016.12.16.somaticSNV_subset_xena.tsv",
        cnv = "PCAWG.consensus.2016.12.16.somaticCNV.tsv"
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
