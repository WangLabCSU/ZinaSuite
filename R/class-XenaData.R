#' XenaData R6 Class
#'
#' @description
#' Concrete implementation of DataSource for UCSC Xena data hubs.
#' Provides access to TCGA, PCAWG, CCLE, and other datasets.
#'
#' @examples
#' \dontrun{
#' # Create XenaData instance
#' xena <- XenaData$new(host = "toilHub")
#'
#' # Query gene expression
#' tp53_expr <- xena$get_gene_expression("TP53")
#'
#' # Query mutation status
#' tp53_mut <- xena$get_mutation_status("TP53")
#'
#' # Query CNV
#' myc_cnv <- xena$get_cnv("MYC")
#'
#' # Batch query with parallel processing
#' results <- xena$query_batch(
#'   identifiers = c("TP53", "BRCA1", "EGFR"),
#'   data_type = "mRNA",
#'   async = TRUE,
#'   n_workers = 4
#' )
#' }
#'
#' @export
XenaData <- R6::R6Class(
  "XenaData",
  inherit = DataSource,

  public = list(
    #' @description
    #' Initialize a new XenaData instance
    #' @param host Xena hub name (e.g., "toilHub", "tcgaHub", "pancanAtlasHub")
    #' @return XenaData instance
    initialize = function(host = "toilHub") {
      super$initialize("Xena", host)
      private$validate_host()
      private$init_dataset_map()
    },

    #' @description
    #' Query data for a single identifier
    #' @param identifier Gene symbol or other molecular identifier
    #' @param data_type Type of data: "mRNA", "protein", "mutation", "cnv", "methylation", "miRNA"
    #' @param ... Additional parameters (e.g., norm for mRNA)
    #' @return Named numeric vector or NULL if not found
    query = function(identifier, data_type = "mRNA", ...) {
      switch(data_type,
        "mRNA" = self$get_gene_expression(identifier, ...),
        "protein" = self$get_protein_expression(identifier),
        "mutation" = self$get_mutation_status(identifier),
        "cnv" = self$get_cnv(identifier, ...),
        "methylation" = self$get_methylation(identifier, ...),
        "miRNA" = self$get_mirna(identifier),
        "transcript" = self$get_transcript_expression(identifier, ...),
        stop("Unknown data type: ", data_type)
      )
    },

    #' @description
    #' Get gene expression data
    #' @param gene Gene symbol (e.g., "TP53")
    #' @param norm Normalization method: "tpm" (default), "fpkm", "nc"
    #' @return Named numeric vector of expression values
    get_gene_expression = function(gene, norm = c("tpm", "fpkm", "nc")) {
      norm <- match.arg(norm)
      dataset <- private$dataset_map[[paste0("gene_", norm)]]

      private$query_with_cache(gene, "mRNA", function(id, ...) {
        private$fetch_from_xena(id, dataset, use_probeMap = TRUE)
      }, norm = norm)
    },

    #' @description
    #' Get protein expression data (RPPA)
    #' @param protein Protein identifier
    #' @return Named numeric vector of protein expression
    get_protein_expression = function(protein) {
      dataset <- private$dataset_map[["protein"]]

      private$query_with_cache(protein, "protein", function(id, ...) {
        private$fetch_from_xena(id, dataset, use_probeMap = TRUE)
      })
    },

    #' @description
    #' Get mutation status
    #' @param gene Gene symbol
    #' @return Named numeric vector (0 = wild-type, 1 = mutated)
    get_mutation_status = function(gene) {
      # Mutation data is on pancanAtlasHub
      old_host <- private$host
      private$host <- "pancanAtlasHub"
      on.exit(private$host <- old_host)

      dataset <- "mc3.v0.2.8.PUBLIC.nonsilentGene.xena"

      private$query_with_cache(gene, "mutation", function(id, ...) {
        private$fetch_from_xena(id, dataset, use_probeMap = FALSE)
      })
    },

    #' @description
    #' Get copy number variation (CNV) data
    #' @param gene Gene symbol
    #' @param use_thresholded Use thresholded data (default: TRUE)
    #' @param ... Additional parameters
    #' @return Named numeric vector of GISTIC2 scores
    get_cnv = function(gene, use_thresholded = TRUE, ...) {
      # CNV data is on tcgaHub
      old_host <- private$host
      private$host <- "tcgaHub"
      on.exit(private$host <- old_host)

      # Ensure use_thresholded is a valid boolean
      if (is.null(use_thresholded) || length(use_thresholded) == 0 || is.na(use_thresholded)) {
        use_thresholded <- TRUE
      } else {
        use_thresholded <- as.logical(use_thresholded)[1]
        if (is.na(use_thresholded)) {
          use_thresholded <- TRUE
        }
      }

      dataset <- if (use_thresholded) {
        "TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
      } else {
        "TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes"
      }

      private$query_with_cache(gene, "cnv", function(id, ...) {
        private$fetch_from_xena(id, dataset, use_probeMap = TRUE)
      }, use_thresholded = use_thresholded)
    },

    #' @description
    #' Get methylation data
    #' @param gene Gene symbol
    #' @param type Methylation array type: "450K" (default) or "27K"
    #' @return Named numeric vector of beta values
    get_methylation = function(gene, type = c("450K", "27K")) {
      type <- match.arg(type)
      dataset <- private$dataset_map[[paste0("methylation_", type)]]

      private$query_with_cache(gene, "methylation", function(id, ...) {
        private$fetch_from_xena(id, dataset, use_probeMap = TRUE)
      }, type = type)
    },

    #' @description
    #' Get miRNA expression data
    #' @param mirna miRNA ID (e.g., "hsa-miR-21-5p")
    #' @return Named numeric vector of miRNA expression
    get_mirna = function(mirna) {
      dataset <- private$dataset_map[["mirna"]]

      private$query_with_cache(mirna, "miRNA", function(id, ...) {
        private$fetch_from_xena(id, dataset, use_probeMap = TRUE)
      })
    },

    #' @description
    #' Get transcript expression data
    #' @param transcript Transcript ID or gene symbol
    #' @param norm Normalization method: "tpm" (default) or "fpkm"
    #' @return Named numeric vector of transcript expression
    get_transcript_expression = function(transcript, norm = c("tpm", "fpkm")) {
      norm <- match.arg(norm)
      dataset <- private$dataset_map[[paste0("transcript_", norm)]]

      private$query_with_cache(transcript, "transcript", function(id, ...) {
        private$fetch_from_xena(id, dataset, use_probeMap = TRUE)
      }, norm = norm)
    },

    #' @description
    #' Batch query with parallel processing using mirai
    #' @param identifiers Vector of molecular identifiers
    #' @param data_type Type of data
    #' @param n_workers Number of parallel workers
    #' @param .progress Show progress bar
    #' @param ... Additional parameters
    #' @return Named list of query results
    query_batch_parallel = function(identifiers, data_type = "mRNA", n_workers = 4, .progress = TRUE, ...) {
      # Get dataset info for this data_type
      dataset_info <- private$get_dataset_info(data_type, ...)

      # Start mirai daemons
      mirai::daemons(n_workers)
      on.exit(mirai::daemons(0), add = TRUE)

      # Create a list of argument sets for each identifier
      arg_list <- lapply(identifiers, function(id) {
        list(id = id, host = dataset_info$host, dataset = dataset_info$dataset, use_probeMap = dataset_info$use_probeMap)
      })

      # Use mirai_map for parallel processing
      results <- mirai::mirai_map(
        arg_list,
        function(args) {
          # Directly use UCSCXenaTools in worker
          if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
            stop("UCSCXenaTools not available")
          }

          # Get host URL
          xe <- UCSCXenaTools::XenaData
          host_url <- unique(xe$XenaHosts[xe$XenaHostNames == args$host])

          if (length(host_url) == 0) {
            stop("Invalid host: ", args$host)
          }

          # Fetch data
          result <- UCSCXenaTools::fetch_dense_values(
            host = host_url[1],
            dataset = args$dataset,
            identifiers = args$id,
            check = FALSE,
            use_probeMap = args$use_probeMap
          )

          if (!is.null(result) && nrow(result) > 0) {
            values <- as.numeric(result[1, ])
            names(values) <- colnames(result)
            values
          } else {
            NULL
          }
        },
        .progress = .progress
      )

      # Collect and name results
      results_list <- results[]
      names(results_list) <- identifiers
      results_list
    }
  ),

  private = list(
    dataset_map = NULL,

    validate_host = function() {
      if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
        stop("Package 'UCSCXenaTools' is required")
      }

      xe <- UCSCXenaTools::XenaData
      valid_hosts <- unique(xe$XenaHostNames)

      if (!(private$host %in% valid_hosts)) {
        stop("Invalid host: ", private$host,
             ". Valid hosts: ", paste(valid_hosts, collapse = ", "))
      }
    },

    init_dataset_map = function() {
      private$dataset_map <- list(
        gene_tpm = "TcgaTargetGtex_rsem_gene_tpm",
        gene_fpkm = "TcgaTargetGtex_rsem_gene_fpkm",
        gene_nc = "TcgaTargetGtex_RSEM_Hugo_norm_count",
        protein = "TCGA-RPPA-pancan-clean.xena",
        methylation_450K = "TCGA.PANCAN.sampleMap_HumanMethylation450",
        methylation_27K = "TCGA.PANCAN.sampleMap_HumanMethylation27",
        mirna = "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena",
        transcript_tpm = "TcgaTargetGtex_rsem_isoform_tpm",
        transcript_fpkm = "TcgaTargetGtex_rsem_isoform_fpkm"
      )
    },

    get_dataset_info = function(data_type, ...) {
      args <- list(...)

      switch(data_type,
        "mRNA" = {
          norm <- args$norm %||% "tpm"
          list(
            host = private$host,
            dataset = private$dataset_map[[paste0("gene_", norm)]],
            use_probeMap = TRUE
          )
        },
        "protein" = list(
          host = private$host,
          dataset = private$dataset_map[["protein"]],
          use_probeMap = TRUE
        ),
        "mutation" = list(
          host = "pancanAtlasHub",
          dataset = "mc3.v0.2.8.PUBLIC.nonsilentGene.xena",
          use_probeMap = FALSE
        ),
        "cnv" = {
          # Default to TRUE if use_thresholded is not provided or invalid
          use_thresholded <- TRUE
          
          # Only process if use_thresholded is provided
          if (!is.null(args$use_thresholded)) {
            # Convert to logical
            logical_val <- as.logical(args$use_thresholded)
            # Check if conversion was successful
            if (!is.na(logical_val) && length(logical_val) > 0) {
              use_thresholded <- logical_val[1]
            }
          }
          
          list(
            host = "tcgaHub",
            dataset = if (use_thresholded) {
              "TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
            } else {
              "TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes"
            },
            use_probeMap = TRUE
          )
        },
        "methylation" = {
          type <- args$type %||% "450K"
          list(
            host = private$host,
            dataset = private$dataset_map[[paste0("methylation_", type)]],
            use_probeMap = TRUE
          )
        },
        "miRNA" = list(
          host = private$host,
          dataset = private$dataset_map[["mirna"]],
          use_probeMap = TRUE
        ),
        "transcript" = {
          norm <- args$norm %||% "tpm"
          list(
            host = private$host,
            dataset = private$dataset_map[[paste0("transcript_", norm)]],
            use_probeMap = TRUE
          )
        },
        stop("Unknown data type: ", data_type)
      )
    },

    fetch_from_xena = function(identifier, dataset, use_probeMap = TRUE, max_try = 5) {
      if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
        stop("Package 'UCSCXenaTools' is required")
      }

      # Get host URL
      xe <- UCSCXenaTools::XenaData
      host_url <- unique(xe$XenaHosts[xe$XenaHostNames == private$host])

      if (length(host_url) == 0) {
        stop("Invalid host: ", private$host)
      }

      # Try to fetch data with retry
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
          # Convert to named vector
          values <- as.numeric(result[1, ])
          names(values) <- colnames(result)

          # Standardize sample IDs for TCGA data
          if (grepl("TCGA", names(values)[1])) {
            names(values) <- standardize_sample_id(names(values), "tcga")
            # Deduplicate samples, preferring "A" suffix
            values <- deduplicate_samples(values, prefer_suffix = "A")
          }

          return(values)
        }
      }

      NULL
    },

    query_with_cache = function(identifier, data_type, query_fn, ...) {
      cache_key <- private$get_cache_key(identifier, data_type, ...)

      # Check cache
      cached <- private$cache_manager$get(cache_key)
      if (!is.null(cached)) {
        return(cached)
      }

      # Query and cache
      result <- query_fn(identifier, data_type, ...)
      if (!is.null(result)) {
        private$cache_manager$set(cache_key, result)
      }

      result
    },

    get_cache_key = function(identifier, data_type, ...) {
      args <- list(...)
      key_parts <- c(private$name, private$host, identifier, data_type, args)
      digest::digest(key_parts, algo = "xxhash32")
    }
  )
)
