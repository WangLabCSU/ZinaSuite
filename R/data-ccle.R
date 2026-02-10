#' CCLE Data Source
#'
#' @description
#' Functions for querying data from CCLE (Cancer Cell Line Encyclopedia).
#'
#' @name ccle-data
NULL

#' Query CCLE Gene Expression
#'
#' @param gene Gene symbol
#' @param host CCLE hub host (default: "publicHub")
#' @return Named numeric vector of gene expression values
#' @export
#'
#' @examples
#' \dontrun{
#' tp53_expr <- query_ccle_expression("TP53")
#' head(tp53_expr)
#' }
query_ccle_expression <- function(gene, host = "publicHub") {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
    stop("Package 'UCSCXenaTools' is required")
  }

  # CCLE gene expression dataset
  dataset <- "ccle/CCLE_RNAseq_genes_rpkm_20180929"

  result <- try_query_value(host, dataset, gene)

  if (is.null(result)) {
    warning("No data returned for gene: ", gene)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "unit") <- "RPKM"
  attr(values, "source") <- "CCLE"

  values
}

#' Query CCLE Mutation Status
#'
#' @param gene Gene symbol
#' @param host CCLE hub host (default: "publicHub")
#' @return Named numeric vector (0 = wild-type, 1 = mutated)
#' @export
#'
#' @examples
#' \dontrun{
#' tp53_mut <- query_ccle_mutation("TP53")
#' table(tp53_mut)
#' }
query_ccle_mutation <- function(gene, host = "publicHub") {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
    stop("Package 'UCSCXenaTools' is required")
  }

  # CCLE mutation dataset
  dataset <- "ccle/CCLE_DepMap_18q3_maf_20180718"

  result <- try_query_value(host, dataset, gene, use_probeMap = FALSE)

  if (is.null(result)) {
    warning("No data returned for gene: ", gene)
    return(NULL)
  }

  # Convert to binary mutation status
  values <- as.numeric(result[1, ])
  values <- ifelse(is.na(values), 0, 1)
  names(values) <- colnames(result)
  attr(values, "source") <- "CCLE"

  values
}

#' Query CCLE Copy Number Variation
#'
#' @param gene Gene symbol
#' @param host CCLE hub host (default: "publicHub")
#' @return Named numeric vector of CNV values
#' @export
#'
#' @examples
#' \dontrun{
#' myc_cnv <- query_ccle_cnv("MYC")
#' summary(myc_cnv)
#' }
query_ccle_cnv <- function(gene, host = "publicHub") {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
    stop("Package 'UCSCXenaTools' is required")
  }

  # CCLE CNV dataset
  dataset <- "ccle/CCLE_gene_cn"

  result <- try_query_value(host, dataset, gene, use_probeMap = FALSE)

  if (is.null(result)) {
    warning("No data returned for gene: ", gene)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "source") <- "CCLE"

  values
}

#' Query CCLE Protein Expression (RPPA)
#'
#' @param protein Protein identifier
#' @param host CCLE hub host (default: "publicHub")
#' @return Named numeric vector of protein expression values
#' @export
#'
#' @examples
#' \dontrun{
#' tp53_prot <- query_ccle_protein("TP53")
#' head(tp53_prot)
#' }
query_ccle_protein <- function(protein, host = "publicHub") {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
    stop("Package 'UCSCXenaTools' is required")
  }

  # CCLE RPPA dataset
  dataset <- "ccle/CCLE_RPPA_20180123"

  result <- try_query_value(host, dataset, protein, use_probeMap = FALSE)

  if (is.null(result)) {
    warning("No data returned for protein: ", protein)
    return(NULL)
  }

  values <- as.numeric(result[1, ])
  names(values) <- colnames(result)
  attr(values, "unit") <- "RPPA"
  attr(values, "source") <- "CCLE"

  values
}

#' CCLEData R6 Class
#'
#' @description
#' R6 class for accessing CCLE (Cancer Cell Line Encyclopedia) data.
#'
#' @export
#' @examples
#' \dontrun{
#' # Create CCLE data instance
#' ccle <- CCLEData$new()
#'
#' # Query gene expression
#' tp53_expr <- ccle$get_gene_expression("TP53")
#'
#' # Query mutation status
#' tp53_mut <- ccle$get_mutation_status("TP53")
#' }
CCLEData <- R6::R6Class(
  "CCLEData",
  inherit = DataSource,

  public = list(
    #' @description
    #' Initialize a new CCLEData instance
    initialize = function() {
      super$initialize("CCLE", "publicHub")
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
    #' Get protein expression data
    #' @param protein Protein identifier
    #' @return Named numeric vector
    get_protein_expression = function(protein) {
      private$query_with_cache(protein, "protein", private$dataset_map$protein, use_probeMap = FALSE)
    },

    #' @description
    #' Query data by type
    #' @param identifier Gene or protein symbol
    #' @param data_type Type of data
    #' @return Named numeric vector
    query = function(identifier, data_type = "expression") {
      switch(data_type,
        "expression" = self$get_gene_expression(identifier),
        "mutation" = self$get_mutation_status(identifier),
        "cnv" = self$get_cnv(identifier),
        "protein" = self$get_protein_expression(identifier),
        stop("Unknown data type: ", data_type)
      )
    }
  ),

  private = list(
    dataset_map = NULL,

    init_dataset_map = function() {
      private$dataset_map <- list(
        expression = "ccle/CCLE_RNAseq_genes_rpkm_20180929",
        mutation = "ccle/CCLE_DepMap_18q3_maf_20180718",
        cnv = "ccle/CCLE_gene_cn",
        protein = "ccle/CCLE_RPPA_20180123"
      )
    },

    query_with_cache = function(identifier, data_type, dataset, use_probeMap = TRUE) {
      cache_key <- paste("ccle", data_type, identifier, sep = "_")

      # Check cache
      cached <- private$cache_manager$get(cache_key)
      if (!is.null(cached)) {
        return(cached)
      }

      # Query from Xena
      result <- private$fetch_from_xena(identifier, dataset, use_probeMap)

      if (is.null(result)) {
        warning("No data returned for ", identifier, " from CCLE ", data_type)
        return(NULL)
      }

      # Convert to named vector
      values <- as.numeric(result[1, ])
      names(values) <- colnames(result)
      attr(values, "source") <- "CCLE"
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
