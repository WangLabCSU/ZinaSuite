#' Query Multiple Molecules with Parallel Processing
#'
#' @description
#' Batch query multiple molecular identifiers using parallel processing via mirai.
#' This function provides significant performance improvements for batch operations.
#'
#' @param identifiers Vector of molecular identifiers
#' @param data_type Type of data (see \code{\link{query_molecule}})
#' @param source Data source: "tcga" (default), "pcawg", "ccle"
#' @param n_workers Number of parallel workers (default: 4)
#' @param .progress Show progress bar (default: TRUE)
#' @param ... Additional parameters passed to query functions
#' @return Named list of query results
#' @export
#'
#' @examples
#' \dontrun{
#' # Query multiple genes in parallel
#' genes <- c("TP53", "BRCA1", "EGFR", "KRAS", "MYC")
#' results <- query_molecules(
#'   identifiers = genes,
#'   data_type = "mRNA",
#'   n_workers = 4,
#'   .progress = TRUE
#' )
#'
#' # Access individual results
#' tp53_expr <- results$TP53
#' brca1_expr <- results$BRCA1
#'
#' # Query multiple data types
#' mut_results <- query_molecules(
#'   identifiers = genes,
#'   data_type = "mutation",
#'   n_workers = 4
#' )
#' }
query_molecules <- function(identifiers,
                            data_type = c("mRNA", "protein", "mutation", "cnv", "methylation", "miRNA", "transcript"),
                            source = c("tcga", "pcawg", "ccle"),
                            n_workers = 4,
                            .progress = TRUE,
                            ...) {
  data_type <- match.arg(data_type)
  source <- match.arg(source)

  # Map source to host
  host <- switch(source,
    "tcga" = "toilHub",
    "pcawg" = "pcawgHub",
    "ccle" = "publicHub",
    "toilHub"
  )

  # Use XenaData's parallel batch query
  xd <- XenaData$new(host = host)
  xd$query_batch_parallel(identifiers, data_type = data_type, n_workers = n_workers, .progress = .progress, ...)
}
