#' Sample ID Utilities for ZinaSuite
#'
#' @description
#' Functions for handling sample IDs across different data sources (TCGA, CCLE, PCAWG).
#' These functions ensure consistent sample ID processing and matching.
#'
#' @name sample_id_utils
NULL

#' Standardize Sample ID
#'
#' @description
#' Standardize sample ID based on data source.
#' For TCGA: extract first 15 characters (TCGA-XX-XXXX-XX)
#' For CCLE: keep as-is (cell line names)
#' For PCAWG: keep as-is (ICGC sample IDs)
#'
#' @param id Sample ID or vector of sample IDs
#' @param source Data source: "tcga", "ccle", "pcawg", or "auto" (default)
#' @return Standardized sample ID(s)
#' @export
#'
#' @examples
#' \dontrun{
#' # TCGA sample ID
#' standardize_sample_id("TCGA-19-1787-01A-01R-1564-08")
#' # Returns: "TCGA-19-1787-01"
#'
#' # CCLE sample ID (unchanged)
#' standardize_sample_id("MCF7_BREAST", source = "ccle")
#' # Returns: "MCF7_BREAST"
#' }
standardize_sample_id <- function(id, source = "auto") {
  if (source == "auto") {
    # Auto-detect source based on ID pattern
    if (grepl("^TCGA", id[1])) {
      source <- "tcga"
    } else if (grepl("^GTEX", id[1])) {
      source <- "gtex"
    } else if (grepl("^DO", id[1]) || grepl("^SP", id[1])) {
      source <- "pcawg"
    } else {
      source <- "ccle"
    }
  }

  switch(source,
    "tcga" = substr(id, 1, 15),
    "gtex" = id,  # GTEx IDs don't need standardization
    "ccle" = id,  # CCLE uses cell line names
    "pcawg" = id, # PCAWG uses ICGC IDs
    id # default: return as-is
  )
}

#' Extract Barcode from Sample ID
#'
#' @description
#' Extract the first 12 characters (TCGA-XX-XXXX) from TCGA sample ID.
#' This is used for matching with clinical data.
#'
#' @param id Sample ID or vector of sample IDs
#' @return Barcode(s) (first 12 characters)
#' @export
#'
#' @examples
#' \dontrun{
#' extract_barcode("TCGA-19-1787-01A-01R-1564-08")
#' # Returns: "TCGA-19-1787"
#' }
extract_barcode <- function(id) {
  substr(id, 1, 12)
}

#' Determine TCGA Sample Type
#'
#' @description
#' Determine if a TCGA sample is tumor or normal based on sample type code.
#' Sample type code is at positions 14-15 of the standardized sample ID.
#' 01-09: Tumor samples
#' 10-19: Normal samples
#' 20-29: Control samples
#'
#' @param id Sample ID or vector of sample IDs
#' @return Character vector: "tumor", "normal", or NA
#' @export
#'
#' @examples
#' \dontrun{
#' determine_tcga_type("TCGA-19-1787-01")  # Returns: "tumor"
#' determine_tcga_type("TCGA-B2-5641-11")  # Returns: "normal"
#' }
determine_tcga_type <- function(id) {
  # Extract sample type code (positions 14-15)
  type_code <- substr(id, 14, 15)
  code_num <- as.numeric(type_code)

  ifelse(is.na(code_num), NA,
    ifelse(code_num >= 1 & code_num <= 9, "tumor",
      ifelse(code_num >= 10 & code_num <= 19, "normal", NA)
    )
  )
}

#' Determine TCGA Cancer Type
#'
#' @description
#' Extract cancer type from TCGA sample ID.
#' Uses the tissue source site (TSS) code to map to cancer type.
#'
#' @param id Sample ID or vector of sample IDs
#' @return Cancer type abbreviation (e.g., "BRCA", "LUAD")
#' @export
#'
#' @examples
#' \dontrun{
#' determine_tcga_cancer("TCGA-19-1787-01")  # Returns: "GBM"
#' }
determine_tcga_cancer <- function(id) {
  # Extract TSS code (part 2 of TCGA ID)
  parts <- strsplit(id, "-")
  tss_codes <- sapply(parts, function(x) {
    if (length(x) >= 2) x[2] else NA
  })

  # TSS to cancer type mapping (simplified)
  # Based on TCGA Tissue Source Site Codes
  tss_map <- c(
    # BRCA
    "A1" = "BRCA", "A2" = "BRCA", "A3" = "BRCA", "A4" = "BRCA",
    "A5" = "BRCA", "A6" = "BRCA", "A7" = "BRCA", "A8" = "BRCA",
    # GBM
    "02" = "GBM", "06" = "GBM", "19" = "GBM", "26" = "GBM",
    # LUAD
    "18" = "LUAD", "38" = "LUAD", "44" = "LUAD", "46" = "LUAD",
    "49" = "LUAD", "4B" = "LUAD", "4C" = "LUAD", "4J" = "LUAD",
    # LUSC
    "18" = "LUSC", "22" = "LUSC", "33" = "LUSC", "37" = "LUSC",
    "39" = "LUSC", "43" = "LUSC", "51" = "LUSC", "52" = "LUSC",
    # COAD/READ
    "01" = "COAD", "05" = "COAD", "09" = "COAD", "0A" = "COAD",
    "0D" = "COAD", "0E" = "COAD", "0F" = "COAD", "14" = "COAD",
    "15" = "COAD", "17" = "COAD", "1A" = "COAD", "1B" = "COAD",
    "20" = "COAD", "21" = "COAD", "23" = "COAD", "24" = "COAD",
    "25" = "COAD", "27" = "COAD", "28" = "COAD", "29" = "COAD",
    "2A" = "COAD", "2C" = "COAD", "2D" = "COAD", "2E" = "COAD",
    "2F" = "COAD", "2G" = "COAD", "2H" = "COAD", "2J" = "COAD",
    "2K" = "COAD", "2L" = "COAD", "2M" = "COAD", "2N" = "COAD",
    "2P" = "COAD", "2Q" = "COAD", "2R" = "COAD", "2T" = "COAD",
    "2U" = "COAD", "2V" = "COAD", "2W" = "COAD", "2X" = "COAD",
    "2Y" = "COAD", "2Z" = "COAD", "30" = "COAD", "31" = "COAD",
    "32" = "COAD", "34" = "COAD", "35" = "COAD", "36" = "COAD",
    "37" = "COAD", "3C" = "COAD", "3E" = "COAD", "3G" = "COAD",
    "3H" = "COAD", "3K" = "COAD", "3L" = "COAD", "3M" = "COAD",
    "3N" = "COAD", "3P" = "COAD", "3Q" = "COAD", "3R" = "COAD",
    "3S" = "COAD", "3T" = "COAD", "3U" = "COAD", "3X" = "COAD",
    "3Z" = "COAD", "41" = "COAD", "42" = "COAD", "43" = "COAD",
    "44" = "COAD", "46" = "COAD", "47" = "COAD", "48" = "COAD",
    "49" = "COAD", "4A" = "COAD", "4B" = "COAD", "4C" = "COAD",
    "4E" = "COAD", "4G" = "COAD", "4H" = "COAD", "4J" = "COAD",
    "4K" = "COAD", "4L" = "COAD", "4N" = "COAD", "4P" = "COAD",
    "4Q" = "COAD", "4R" = "COAD", "4T" = "COAD", "4V" = "COAD",
    "4W" = "COAD", "4X" = "COAD", "4Z" = "COAD", "50" = "COAD",
    "51" = "COAD", "52" = "COAD", "53" = "COAD", "55" = "COAD",
    "56" = "COAD", "57" = "COAD", "58" = "COAD", "59" = "COAD",
    "5B" = "COAD", "5C" = "COAD", "5D" = "COAD", "5E" = "COAD",
    "5F" = "COAD", "5G" = "COAD", "5H" = "COAD", "5J" = "COAD",
    "5K" = "COAD", "5L" = "COAD", "5M" = "COAD", "5N" = "COAD",
    "5P" = "COAD", "5Q" = "COAD", "5R" = "COAD", "5S" = "COAD",
    "5T" = "COAD", "5U" = "COAD", "5V" = "COAD", "5W" = "COAD",
    "5X" = "COAD", "5Y" = "COAD", "5Z" = "COAD"
  )

  # Return mapped cancer type or "Unknown"
  unname(tss_map[tss_codes]) %||% "Unknown"
}

#' Match Samples Between Different Data Sources
#'
#' @description
#' Match samples between two datasets with potentially different ID formats.
#' Returns the intersection of standardized sample IDs.
#'
#' @param ids1 First set of sample IDs
#' @param ids2 Second set of sample IDs
#' @param source1 Source of first set ("tcga", "ccle", "pcawg", "auto")
#' @param source2 Source of second set ("tcga", "ccle", "pcawg", "auto")
#' @return List with matched IDs and indices
#' @export
#'
#' @examples
#' \dontrun{
#' ids1 <- c("TCGA-19-1787-01A", "TCGA-B2-5641-11A")
#' ids2 <- c("TCGA-19-1787-01", "TCGA-B2-5641-11")
#' match <- match_samples(ids1, ids2)
#' }
match_samples <- function(ids1, ids2, source1 = "auto", source2 = "auto") {
  # Standardize both sets
  std1 <- standardize_sample_id(ids1, source1)
  std2 <- standardize_sample_id(ids2, source2)

  # Find intersection
  common <- intersect(std1, std2)

  # Get indices
  idx1 <- match(common, std1)
  idx2 <- match(common, std2)

  list(
    common_ids = common,
    idx1 = idx1,
    idx2 = idx2,
    n_matched = length(common)
  )
}

#' Deduplicate Samples
#'
#' @description
#' Remove duplicate samples, preferring samples with specific suffix.
#' Default prefers samples with "A" suffix (e.g., TCGA-XX-XXXX-01A).
#'
#' @param data Named vector or data frame with sample IDs as names/rownames
#' @param prefer_suffix Suffix to prefer when deduplicating (default: "A")
#' @return Deduplicated data
#' @export
#'
#' @examples
#' \dontrun{
#' data <- c(5.2, 6.1, 7.3)
#' names(data) <- c("TCGA-19-1787-01", "TCGA-19-1787-01A", "TCGA-B2-5641-11")
#' dedup <- deduplicate_samples(data)
#' }
deduplicate_samples <- function(data, prefer_suffix = "A") {
  if (is.data.frame(data)) {
    ids <- rownames(data)
  } else {
    ids <- names(data)
  }

  # Standardize IDs for deduplication
  std_ids <- standardize_sample_id(ids, "tcga")

  # Find duplicates
  dup_mask <- duplicated(std_ids, fromLast = FALSE)
  dup_ids <- std_ids[dup_mask]

  if (length(dup_ids) == 0) {
    return(data) # No duplicates
  }

  # For each duplicate group, prefer the one with specified suffix
  keep_idx <- seq_along(ids)

  for (dup_id in unique(dup_ids)) {
    group_idx <- which(std_ids == dup_id)
    group_ids <- ids[group_idx]

    # Check if any has preferred suffix
    has_suffix <- grepl(paste0(prefer_suffix, "$"), group_ids)

    if (any(has_suffix)) {
      # Keep the first one with preferred suffix
      keep <- group_idx[has_suffix][1]
    } else {
      # Keep the first one
      keep <- group_idx[1]
    }

    # Mark others for removal
    keep_idx <- setdiff(keep_idx, setdiff(group_idx, keep))
  }

  # Return deduplicated data
  if (is.data.frame(data)) {
    data[keep_idx, , drop = FALSE]
  } else {
    data[keep_idx]
  }
}

#' Filter TCGA Tumor Samples
#'
#' @description
#' Filter TCGA data to keep only tumor samples (sample type 01-09).
#'
#' @param data Named vector or data frame with TCGA sample IDs
#' @return Filtered data with only tumor samples
#' @export
#'
#' @examples
#' \dontrun{
#' data <- query_gene_expression("TP53")
#' tumor_data <- filter_tcga_tumor(data)
#' }
filter_tcga_tumor <- function(data) {
  if (is.data.frame(data)) {
    ids <- rownames(data)
  } else {
    ids <- names(data)
  }

  is_tumor <- determine_tcga_type(ids) == "tumor"

  if (is.data.frame(data)) {
    data[is_tumor, , drop = FALSE]
  } else {
    data[is_tumor]
  }
}

#' Filter TCGA Normal Samples
#'
#' @description
#' Filter TCGA data to keep only normal samples (sample type 10-19).
#'
#' @param data Named vector or data frame with TCGA sample IDs
#' @return Filtered data with only normal samples
#' @export
#'
#' @examples
#' \dontrun{
#' data <- query_gene_expression("TP53")
#' normal_data <- filter_tcga_normal(data)
#' }
filter_tcga_normal <- function(data) {
  if (is.data.frame(data)) {
    ids <- rownames(data)
  } else {
    ids <- names(data)
  }

  is_normal <- determine_tcga_type(ids) == "normal"

  if (is.data.frame(data)) {
    data[is_normal, , drop = FALSE]
  } else {
    data[is_normal]
  }
}

#' Prepare Data for Analysis
#'
#' @description
#' Prepare and merge multiple datasets for analysis.
#' Handles sample ID standardization and matching automatically.
#'
#' @param ... Named vectors or data frames to merge
#' @param source Source type ("tcga", "ccle", "pcawg", "auto")
#' @return Data frame with merged data
#' @export
#'
#' @examples
#' \dontrun{
#' g1 <- query_gene_expression("TP53")
#' g2 <- query_gene_expression("BRCA1")
#' merged <- prepare_analysis_data(gene1 = g1, gene2 = g2)
#' }
prepare_analysis_data <- function(..., source = "auto") {
  data_list <- list(...)

  if (length(data_list) == 0) {
    return(NULL)
  }

  # Standardize all sample IDs
  std_ids_list <- lapply(data_list, function(x) {
    if (is.data.frame(x)) {
      standardize_sample_id(rownames(x), source)
    } else {
      standardize_sample_id(names(x), source)
    }
  })

  # Find common samples
  common_ids <- std_ids_list[[1]]
  if (length(std_ids_list) > 1) {
    for (i in 2:length(std_ids_list)) {
      common_ids <- intersect(common_ids, std_ids_list[[i]])
    }
  }

  if (length(common_ids) == 0) {
    warning("No common samples found between datasets")
    return(data.frame())
  }

  # Create merged data frame
  result <- data.frame(Sample = common_ids)

  for (i in seq_along(data_list)) {
    data <- data_list[[i]]
    name <- names(data_list)[i]

    if (is.data.frame(data)) {
      values <- data[match(common_ids, std_ids_list[[i]]), , drop = FALSE]
      result <- cbind(result, values)
    } else {
      values <- data[match(common_ids, std_ids_list[[i]])]
      result[[name]] <- values
    }
  }

  result
}
