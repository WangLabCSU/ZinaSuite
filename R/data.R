#' TCGA Clinical Data
#'
#' Clinical phenotype data for TCGA samples.
#'
#' @docType data
#' @name tcga_clinical
#' @format A data frame with TCGA clinical information
#' @source UCSC Xena Toil Hub
#' @examples
#' data("tcga_clinical")
#' head(tcga_clinical)
NULL

#' TCGA Clinical Data (Fine-grained)
#'
#' Cleaned TCGA clinical data optimized for sample grouping.
#'
#' @docType data
#' @name tcga_clinical_fine
#' @format A data frame with cleaned clinical variables
#' @source UCSC Xena Toil Hub
#' @examples
#' data("tcga_clinical_fine")
#' head(tcga_clinical_fine)
NULL

#' TCGA Survival Data
#'
#' Overall survival data for TCGA samples.
#'
#' @docType data
#' @name tcga_surv
#' @format A data frame with OS and PFS information
#' @source UCSC Xena Toil Hub
#' @examples
#' data("tcga_surv")
#' head(tcga_surv)
NULL

#' TCGA GTEx Sample Information
#'
#' Sample information for TCGA and GTEx samples.
#'
#' @docType data
#' @name tcga_gtex
#' @format A data frame with sample metadata
#' @source UCSC Xena Toil Hub
#' @examples
#' data("tcga_gtex")
#' head(tcga_gtex)
NULL

#' TCGA Purity Data
#'
#' Tumor purity estimates for TCGA samples.
#'
#' @docType data
#' @name tcga_purity
#' @format A data frame with purity scores
#' @source Aran et al., Nat Commun 2015
#' @examples
#' data("tcga_purity")
#' head(tcga_purity)
NULL

#' TCGA Subtype Data
#'
#' Molecular subtype classifications for TCGA samples.
#'
#' @docType data
#' @name tcga_subtypes
#' @format A data frame with subtype assignments
#' @source UCSC Xena
#' @examples
#' data("tcga_subtypes")
#' head(tcga_subtypes)
NULL

#' TCGA Organ Data
#'
#' Anatomical location mapping for TCGA cancer types.
#'
#' @docType data
#' @name TCGA.organ
#' @format A data frame with organ mappings
#' @source UCSC Xena
#' @examples
#' data("TCGA.organ")
#' head(TCGA.organ)
NULL

#' PCAWG Sample Information
#'
#' Phenotype data for PCAWG samples.
#'
#' @docType data
#' @name pcawg_info
#' @format A data frame with PCAWG sample metadata
#' @source UCSC Xena PCAWG Hub
#' @examples
#' data("pcawg_info")
#' head(pcawg_info)
NULL

#' PCAWG Sample Information (Fine-grained)
#'
#' Cleaned PCAWG phenotype data for grouping.
#'
#' @docType data
#' @name pcawg_info_fine
#' @format A data frame with cleaned PCAWG variables
#' @source UCSC Xena PCAWG Hub
#' @examples
#' data("pcawg_info_fine")
#' head(pcawg_info_fine)
NULL

#' PCAWG Purity Data
#'
#' Tumor purity estimates for PCAWG samples.
#'
#' @docType data
#' @name pcawg_purity
#' @format A data frame with purity scores
#' @source UCSC Xena PCAWG Hub
#' @examples
#' data("pcawg_purity")
#' head(pcawg_purity)
NULL

#' CCLE Sample Information
#'
#' Phenotype data for CCLE cell lines.
#'
#' @docType data
#' @name ccle_info
#' @format A data frame with CCLE cell line metadata
#' @source UCSC Xena Public Hub
#' @examples
#' data("ccle_info")
#' head(ccle_info)
NULL

#' CCLE Sample Information (Fine-grained)
#'
#' Cleaned CCLE phenotype data for grouping.
#'
#' @docType data
#' @name ccle_info_fine
#' @format A data frame with cleaned CCLE variables
#' @source UCSC Xena Public Hub
#' @examples
#' data("ccle_info_fine")
#' head(ccle_info_fine)
NULL

#' CCLE ABSOLUTE Data
#'
#' ABSOLUTE purity and ploidy estimates for CCLE cell lines.
#'
#' @docType data
#' @name ccle_absolute
#' @format A data frame with ABSOLUTE results
#' @source UCSC Xena Public Hub
#' @examples
#' data("ccle_absolute")
#' head(ccle_absolute)
NULL

#' TCGA TIL Data
#'
#' Tumor-infiltrating lymphocyte (TIL) fractions for TCGA samples.
#'
#' @docType data
#' @name tcga_TIL
#' @format A data frame with TIL fractions
#' @source Thorsson et al., Immunity 2018
#' @examples
#' data("tcga_TIL")
#' head(tcga_TIL)
NULL

#' TCGA TMB Data
#'
#' Tumor mutation burden (TMB) for TCGA samples.
#'
#' @docType data
#' @name tcga_tmb
#' @format A data frame with TMB scores
#' @source GDC Pan-Cancer Atlas
#' @examples
#' data("tcga_tmb")
#' head(tcga_tmb)
NULL

#' TCGA MSI Data
#'
#' Microsatellite instability (MSI) scores for TCGA samples.
#'
#' @docType data
#' @name tcga_MSI
#' @format A data frame with MSI scores
#' @source GDC Pan-Cancer Atlas
#' @examples
#' data("tcga_MSI")
#' head(tcga_MSI)
NULL

#' TCGA Stemness Data
#'
#' Tumor stemness scores for TCGA samples.
#'
#' @docType data
#' @name tcga_stemness
#' @format A data frame with stemness scores
#' @source Malta et al., Cell 2018
#' @examples
#' data("tcga_stemness")
#' head(tcga_stemness)
NULL

#' TCGA Genome Instability Data
#'
#' Genome instability metrics for TCGA samples.
#'
#' @docType data
#' @name tcga_genome_instability
#' @format A data frame with instability scores
#' @source GDC Pan-Cancer Atlas
#' @examples
#' data("tcga_genome_instability")
#' head(tcga_genome_instability)
NULL
