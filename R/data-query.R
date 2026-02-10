#' Query Molecular Data
#'
#' @description
#' Unified interface for querying molecular data from various sources (TCGA, PCAWG, CCLE).
#' This function provides a simple interface to fetch molecular data by gene symbol
#' or other identifiers.
#'
#' @param identifier Gene symbol, transcript ID, or other molecular identifier.
#'   For gene expression, use official gene symbols (e.g., "TP53", "BRCA1").
#' @param data_type Type of molecular data to query. Options:
#'   - "mRNA": Gene expression (TPM/FPKM)
#'   - "protein": Protein expression (RPPA)
#'   - "mutation": Mutation status
#'   - "cnv": Copy number variation
#'   - "methylation": DNA methylation
#'   - "miRNA": miRNA expression
#'   - "transcript": Transcript-level expression
#' @param source Data source. One of "tcga", "pcawg", or "ccle".
#' @param dataset Specific dataset name (optional). If NULL, uses default dataset.
#' @return A named vector of values with sample IDs as names.
#' @export
#'
#' @examples
#' \donttest{
#' # Query TP53 gene expression from TCGA
#' tp53_expr <- query_molecule("TP53", data_type = "mRNA", source = "tcga")
#' head(tp53_expr)
#'
#' # Query BRCA1 mutation status
#' brca1_mut <- query_molecule("BRCA1", data_type = "mutation", source = "tcga")
#' table(brca1_mut)
#'
#' # Query EGFR CNV from PCAWG
#' egfr_cnv <- query_molecule("EGFR", data_type = "cnv", source = "pcawg")
#'
#' # Query miR-21 expression
#' mir21 <- query_molecule("hsa-miR-21-5p", data_type = "miRNA", source = "tcga")
#' }
query_molecule <- function(identifier,
                           data_type = c("mRNA", "protein", "mutation", "cnv", "methylation", "miRNA", "transcript"),
                           source = c("tcga", "pcawg", "ccle"),
                           dataset = NULL) {
  data_type <- match.arg(data_type)
  source <- match.arg(source)

  # Validate identifier
  if (!is_valid_gene(identifier)) {
    stop("Invalid identifier: ", identifier)
  }

  # Get data source instance
  ds <- get_data_source(source)

  # Query data
  ds$query(identifier, data_type, dataset)
}

#' Query Gene Expression
#'
#' @description
#' Convenience function to query gene expression data.
#'
#' @param gene Gene symbol (e.g., "TP53", "BRCA1")
#' @param source Data source. One of "tcga", "pcawg", or "ccle".
#' @param dataset Specific dataset name (optional). If NULL, uses default dataset.
#' @return Named vector of expression values
#' @export
#'
#' @examples
#' \donttest{
#' # Query TP53 expression from TCGA
#' tp53 <- query_gene_expression("TP53", source = "tcga")
#' summary(tp53)
#'
#' # Query multiple genes
#' genes <- c("TP53", "BRCA1", "EGFR")
#' expr_list <- lapply(genes, query_gene_expression, source = "tcga")
#' names(expr_list) <- genes
#' }
query_gene_expression <- function(gene,
                                  source = c("tcga", "pcawg", "ccle"),
                                  dataset = NULL) {
  source <- match.arg(source)
  query_molecule(gene, data_type = "mRNA", source = source, dataset = dataset)
}

#' Query Mutation Status
#'
#' @description
#' Convenience function to query mutation status.
#'
#' @param gene Gene symbol
#' @param source Data source. One of "tcga", "pcawg", or "ccle".
#' @param dataset Specific dataset name (optional). If NULL, uses default dataset.
#' @return Named vector with mutation status
#' @export
#'
#' @examples
#' \donttest{
#' # Query TP53 mutation status
#' tp53_mut <- query_mutation("TP53", source = "tcga")
#' table(tp53_mut)
#'
#' # Query multiple genes
#' genes <- c("TP53", "BRCA1", "PIK3CA")
#' mut_list <- lapply(genes, query_mutation, source = "tcga")
#' names(mut_list) <- genes
#' }
query_mutation <- function(gene,
                           source = c("tcga", "pcawg", "ccle"),
                           dataset = NULL) {
  source <- match.arg(source)
  query_molecule(gene, data_type = "mutation", source = source, dataset = dataset)
}

#' Query Copy Number Variation
#'
#' @description
#' Convenience function to query CNV data.
#'
#' @param gene Gene symbol
#' @param source Data source. One of "tcga", "pcawg", or "ccle".
#' @param dataset Specific dataset name (optional). If NULL, uses default dataset.
#' @return Named vector of CNV values (GISTIC2 scores)
#' @export
#'
#' @examples
#' \donttest{
#' # Query TP53 CNV
#' tp53_cnv <- query_cnv("TP53", source = "tcga")
#' table(tp53_cnv)
#' }
query_cnv <- function(gene,
                      source = c("tcga", "pcawg", "ccle"),
                      dataset = NULL) {
  source <- match.arg(source)
  query_molecule(gene, data_type = "cnv", source = source, dataset = dataset)
}

#' Query Protein Expression
#'
#' @description
#' Convenience function to query protein expression (RPPA) data.
#'
#' @param protein Protein name or gene symbol
#' @param source Data source. One of "tcga", "pcawg", or "ccle".
#' @param dataset Specific dataset name (optional). If NULL, uses default dataset.
#' @return Named vector of protein expression values
#' @export
#'
#' @examples
#' \donttest{
#' # Query AKT protein expression
#' akt_prot <- query_protein("AKT", source = "tcga")
#' summary(akt_prot)
#' }
query_protein <- function(protein,
                          source = c("tcga", "pcawg", "ccle"),
                          dataset = NULL) {
  source <- match.arg(source)
  query_molecule(protein, data_type = "protein", source = source, dataset = dataset)
}

#' Query DNA Methylation
#'
#' @description
#' Convenience function to query DNA methylation data.
#'
#' @param gene Gene symbol or probe ID
#' @param type Methylation array type: "450K" or "27K"
#' @param source Data source. One of "tcga", "pcawg", or "ccle".
#' @param dataset Specific dataset name (optional). If NULL, uses default dataset.
#' @return Named vector of methylation beta values
#' @export
#'
#' @examples
#' \donttest{
#' # Query TP53 methylation
#' tp53_meth <- query_methylation("TP53", source = "tcga")
#' summary(tp53_meth)
#' }
query_methylation <- function(gene,
                              type = c("450K", "27K"),
                              source = c("tcga", "pcawg", "ccle"),
                              dataset = NULL) {
  source <- match.arg(source)
  type <- match.arg(type)

  ds <- get_data_source(source)
  ds$get_methylation(gene, dataset = dataset)
}

#' Query miRNA Expression
#'
#' @description
#' Convenience function to query miRNA expression data.
#'
#' @param mirna miRNA ID (e.g., "hsa-miR-21-5p")
#' @param source Data source. One of "tcga", "pcawg", or "ccle".
#' @param dataset Specific dataset name (optional). If NULL, uses default dataset.
#' @return Named vector of miRNA expression values
#' @export
#'
#' @examples
#' \donttest{
#' # Query miR-21 expression
#' mir21 <- query_mirna("hsa-miR-21-5p", source = "tcga")
#' summary(mir21)
#' }
query_mirna <- function(mirna,
                        source = c("tcga", "pcawg", "ccle"),
                        dataset = NULL) {
  source <- match.arg(source)
  query_molecule(mirna, data_type = "miRNA", source = source, dataset = dataset)
}
