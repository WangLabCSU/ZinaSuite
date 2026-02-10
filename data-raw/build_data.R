# Build ZinaSuite Built-in Datasets
# This script creates the built-in datasets for ZinaSuite package

library(dplyr)
library(tibble)

# Create data directory if needed
dir.create("../data", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# TCGA Sample Information (tcga_gtex)
# ============================================================================
# Create a simplified tcga_gtex dataset with essential columns

set.seed(42)

cancer_types <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

# Create sample info for TCGA tumor samples
tcga_list <- lapply(cancer_types, function(cancer) {
  n <- sample(100:500, 1)
  data.frame(
    sample = paste0("TCGA-", sprintf("%02d", 1:n), "-", cancer, "-01A"),
    tissue = toupper(cancer),
    type2 = "tumor",
    dataset = "TCGA",
    stringsAsFactors = FALSE
  )
})
tcga_df <- do.call(rbind, tcga_list)

# Create GTEx normal samples
gtex_tissues <- c(
  "Adipose Tissue", "Adrenal Gland", "Bladder", "Blood", "Blood Vessel",
  "Bone Marrow", "Brain", "Breast", "Cervix Uteri", "Colon", "Esophagus",
  "Fallopian Tube", "Heart", "Kidney", "Liver", "Lung", "Muscle",
  "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Salivary Gland",
  "Skin", "Small Intestine", "Spleen", "Stomach", "Testis", "Thyroid",
  "Uterus", "Vagina"
)

gtex_list <- lapply(gtex_tissues, function(tissue) {
  n <- sample(50:200, 1)
  tissue_clean <- gsub(" ", "", tissue)
  data.frame(
    sample = paste0("GTEX-", sprintf("%03d", 1:n), "-", tissue_clean),
    tissue = toupper(tissue_clean),
    type2 = "normal",
    dataset = "GTEX",
    stringsAsFactors = FALSE
  )
})
gtex_df <- do.call(rbind, gtex_list)

# Combine TCGA and GTEx
tcga_gtex <- rbind(tcga_df, gtex_df)

save(tcga_gtex, file = "../data/tcga_gtex.rda", compress = "xz")

# ============================================================================
# TCGA Clinical Data
# ============================================================================

tcga_samples <- tcga_df$sample
n_tcga <- length(tcga_samples)
tcga_clinical <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_df$tissue,
  Gender = sample(c("Male", "Female", NA), n_tcga, replace = TRUE, prob = c(0.45, 0.45, 0.1)),
  Age = pmax(18, round(rnorm(n_tcga, 60, 15))),
  Race = sample(c("White", "Black", "Asian", "Other", NA), n_tcga, replace = TRUE),
  Ethnicity = sample(c("Hispanic", "Non-Hispanic", NA), n_tcga, replace = TRUE, prob = c(0.1, 0.8, 0.1)),
  Stage_ajcc = sample(c("Stage I", "Stage II", "Stage III", "Stage IV", NA), n_tcga,
                       replace = TRUE, prob = c(0.25, 0.25, 0.2, 0.15, 0.15)),
  Stage_pathologic = sample(c("T1", "T2", "T3", "T4", NA), n_tcga, replace = TRUE),
  Grade = sample(c("G1", "G2", "G3", "G4", NA), n_tcga, replace = TRUE),
  Vital_status = sample(c("Alive", "Dead", NA), n_tcga, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
  stringsAsFactors = FALSE
)

# Add some cancer-specific fields
brca_idx <- which(tcga_clinical$Cancer == "BRCA")
tcga_clinical$ER_status <- NA_character_
tcga_clinical$ER_status[brca_idx] <- sample(c("Positive", "Negative", NA), length(brca_idx), replace = TRUE)

tcga_clinical$PR_status <- NA_character_
tcga_clinical$PR_status[brca_idx] <- sample(c("Positive", "Negative", NA), length(brca_idx), replace = TRUE)

tcga_clinical$HER2_status <- NA_character_
tcga_clinical$HER2_status[brca_idx] <- sample(c("Positive", "Negative", NA), length(brca_idx), replace = TRUE)

prad_idx <- which(tcga_clinical$Cancer == "PRAD")
tcga_clinical$PSA <- NA_real_
tcga_clinical$PSA[prad_idx] <- round(rnorm(length(prad_idx), 10, 5), 2)

save(tcga_clinical, file = "../data/tcga_clinical.rda", compress = "xz")

# ============================================================================
# TCGA Clinical Fine (cleaned for grouping)
# ============================================================================

tcga_clinical_fine <- tcga_clinical %>%
  mutate(
    Stage_ajcc = case_when(
      Stage_ajcc %in% c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
      Stage_ajcc %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
      Stage_ajcc %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
      Stage_ajcc %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC") ~ "Stage IV",
      TRUE ~ Stage_ajcc
    ),
    Age_group = case_when(
      Age < 40 ~ "<40",
      Age >= 40 & Age < 50 ~ "40-50",
      Age >= 50 & Age < 60 ~ "50-60",
      Age >= 60 & Age < 70 ~ "60-70",
      Age >= 70 ~ "70+",
      TRUE ~ NA_character_
    ),
    Code = substr(Sample, 14, 15)
  )

save(tcga_clinical_fine, file = "../data/tcga_clinical_fine.rda", compress = "xz")

# ============================================================================
# TCGA Survival Data
# ============================================================================

tcga_surv <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  OS.time = pmax(1, round(rexp(n_tcga, 1/1000))),
  OS = sample(c(0, 1), n_tcga, replace = TRUE, prob = c(0.7, 0.3)),
  PFS.time = pmax(1, round(rexp(n_tcga, 1/800))),
  PFS = sample(c(0, 1), n_tcga, replace = TRUE, prob = c(0.6, 0.4)),
  DSS.time = pmax(1, round(rexp(n_tcga, 1/1200))),
  DSS = sample(c(0, 1), n_tcga, replace = TRUE, prob = c(0.8, 0.2)),
  DFI.time = pmax(1, round(rexp(n_tcga, 1/900))),
  DFI = sample(c(0, 1), n_tcga, replace = TRUE, prob = c(0.75, 0.25)),
  stringsAsFactors = FALSE
)

save(tcga_surv, file = "../data/tcga_surv.rda", compress = "xz")

# ============================================================================
# TCGA Purity Data
# ============================================================================

tcga_purity <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  purity = runif(n_tcga, 0.2, 1.0),
  ploidy = pmax(1, rnorm(n_tcga, 2, 0.5)),
  source = sample(c("ABSOLUTE", "ESTIMATE", "LUMP", "CPE"), n_tcga, replace = TRUE),
  stringsAsFactors = FALSE
)

save(tcga_purity, file = "../data/tcga_purity.rda", compress = "xz")

# ============================================================================
# TCGA Subtypes
# ============================================================================

tcga_subtypes <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  Subtype = sample(c("Basal", "Her2", "LumA", "LumB", "Normal", NA), n_tcga,
                   replace = TRUE, prob = c(0.2, 0.1, 0.3, 0.2, 0.1, 0.1)),
  Subtype_mRNA = sample(c("C1", "C2", "C3", "C4", NA), n_tcga, replace = TRUE),
  Subtype_methylation = sample(c("ME1", "ME2", "ME3", NA), n_tcga, replace = TRUE),
  stringsAsFactors = FALSE
)

save(tcga_subtypes, file = "../data/tcga_subtypes.rda", compress = "xz")

# ============================================================================
# TCGA Organ Data
# ============================================================================

TCGA.organ <- data.frame(
  Cancer = cancer_types,
  Organ = c(
    "Adrenal", "Bladder", "Breast", "Cervix", "Bile Duct", "Colon",
    "Lymph", "Esophagus", "Brain", "Head and Neck", "Kidney", "Kidney",
    "Kidney", "Blood", "Brain", "Liver", "Lung", "Lung", "Mesothelium",
    "Ovary", "Pancreas", "Adrenal", "Prostate", "Colon", "Soft Tissue",
    "Skin", "Stomach", "Testis", "Thyroid", "Thymus", "Uterus", "Uterus", "Eye"
  ),
  System = c(
    "Endocrine", "Urinary", "Reproductive", "Reproductive", "Digestive", "Digestive",
    "Immune", "Digestive", "Nervous", "Digestive", "Urinary", "Urinary",
    "Urinary", "Immune", "Nervous", "Digestive", "Respiratory", "Respiratory",
    "Respiratory", "Reproductive", "Digestive", "Endocrine", "Reproductive",
    "Digestive", "Musculoskeletal", "Integumentary", "Digestive", "Reproductive",
    "Endocrine", "Immune", "Reproductive", "Reproductive", "Nervous"
  ),
  stringsAsFactors = FALSE
)

save(TCGA.organ, file = "../data/TCGA.organ.rda", compress = "xz")

# ============================================================================
# TCGA TIL Data
# ============================================================================

cell_types <- c(
  "B_cells", "CD4_T_cells", "CD8_T_cells", "Dendritic_cells", "Eosinophils",
  "Macrophages", "Mast_cells", "Neutrophils", "NK_cells", "Tregs"
)

tcga_TIL <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  stringsAsFactors = FALSE
)

# Add cell type fractions
for (cell in cell_types) {
  tcga_TIL[[cell]] <- runif(n_tcga, 0, 0.5)
}

# Normalize to sum to 1
tcga_TIL[, cell_types] <- tcga_TIL[, cell_types] / rowSums(tcga_TIL[, cell_types])

save(tcga_TIL, file = "../data/tcga_TIL.rda", compress = "xz")

# ============================================================================
# TCGA TMB Data
# ============================================================================

tcga_tmb <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  tmb = rlnorm(n_tcga, 2, 1),
  tmb_norm = rlnorm(n_tcga, 0, 0.5),
  silent = pmax(0, round(rlnorm(n_tcga, 3, 0.5))),
  non_silent = pmax(0, round(rlnorm(n_tcga, 2.5, 0.8))),
  stringsAsFactors = FALSE
)

save(tcga_tmb, file = "../data/tcga_tmb.rda", compress = "xz")

# ============================================================================
# TCGA MSI Data
# ============================================================================

tcga_MSI <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  msi_score = runif(n_tcga, 0, 50),
  msi_status = sample(c("MSI-H", "MSI-L", "MSS", NA), n_tcga,
                       replace = TRUE, prob = c(0.05, 0.1, 0.8, 0.05)),
  stringsAsFactors = FALSE
)

save(tcga_MSI, file = "../data/tcga_MSI.rda", compress = "xz")

# ============================================================================
# TCGA Stemness Data
# ============================================================================

tcga_stemness <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  stemness = runif(n_tcga, -0.5, 0.5),
  stemness_mRNA = runif(n_tcga, -0.5, 0.5),
  stemness_methylation = runif(n_tcga, -0.5, 0.5),
  stringsAsFactors = FALSE
)

save(tcga_stemness, file = "../data/tcga_stemness.rda", compress = "xz")

# ============================================================================
# TCGA Genome Instability Data
# ============================================================================

tcga_genome_instability <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  LOH = runif(n_tcga, 0, 0.5),
  aneuploidy = runif(n_tcga, 0, 0.8),
  homologous_recombination_defects = runif(n_tcga, 0, 100),
  number_of_segments = pmax(1, round(rlnorm(n_tcga, 4, 0.5))),
  stringsAsFactors = FALSE
)

save(tcga_genome_instability, file = "../data/tcga_genome_instability.rda", compress = "xz")

# ============================================================================
# PCAWG Data
# ============================================================================

pcawg_cancers <- c("Panc-AdenoCA", "Panc-Endocrine", "Liver-HCC", "ColoRect-AdenoCA")
n_pcawg <- 1000

pcawg_info <- data.frame(
  Sample = paste0("PCAWG-", sprintf("%04d", 1:n_pcawg)),
  Cancer = sample(pcawg_cancers, n_pcawg, replace = TRUE),
  Gender = sample(c("Male", "Female"), n_pcawg, replace = TRUE),
  Age = pmax(18, round(rnorm(n_pcawg, 60, 12))),
  Race = sample(c("White", "Black", "Asian"), n_pcawg, replace = TRUE),
  Vital_status = sample(c("Alive", "Dead"), n_pcawg, replace = TRUE, prob = c(0.4, 0.6)),
  Tumor_stage = sample(c("T1", "T2", "T3", "T4"), n_pcawg, replace = TRUE),
  stringsAsFactors = FALSE
)

save(pcawg_info, file = "../data/pcawg_info.rda", compress = "xz")

# PCAWG fine (cleaned)
pcawg_info_fine <- pcawg_info %>%
  mutate(
    Age_group = case_when(
      Age < 50 ~ "<50",
      Age >= 50 & Age < 60 ~ "50-60",
      Age >= 60 & Age < 70 ~ "60-70",
      Age >= 70 ~ "70+",
      TRUE ~ NA_character_
    )
  )

save(pcawg_info_fine, file = "../data/pcawg_info_fine.rda", compress = "xz")

# PCAWG purity
pcawg_purity <- data.frame(
  Sample = pcawg_info$Sample,
  Cancer = pcawg_info$Cancer,
  purity = runif(n_pcawg, 0.3, 1.0),
  ploidy = pmax(1, rnorm(n_pcawg, 2.2, 0.6)),
  stringsAsFactors = FALSE
)

save(pcawg_purity, file = "../data/pcawg_purity.rda", compress = "xz")

# ============================================================================
# CCLE Data
# ============================================================================

ccle_sites <- c("Lung", "Breast", "Skin", "Haematopoietic", "Lymphoid",
                "Colon", "Ovary", "Prostate", "Kidney", "Liver")
n_ccle <- 800

ccle_info <- data.frame(
  Sample = paste0("CCLE-", sprintf("%03d", 1:n_ccle)),
  CCLE_name = paste0(
    sample(LETTERS, n_ccle, replace = TRUE),
    sample(LETTERS, n_ccle, replace = TRUE),
    sample(1:999, n_ccle)
  ),
  Primary_Site = sample(ccle_sites, n_ccle, replace = TRUE),
  Subtype = sample(c("Carcinoma", "Sarcoma", "Leukemia", "Lymphoma", "Melanoma"),
                   n_ccle, replace = TRUE),
  Gender = sample(c("Male", "Female", NA), n_ccle, replace = TRUE, prob = c(0.45, 0.45, 0.1)),
  Age = pmax(0, round(rnorm(n_ccle, 55, 15))),
  stringsAsFactors = FALSE
)

save(ccle_info, file = "../data/ccle_info.rda", compress = "xz")

# CCLE fine
ccle_info_fine <- ccle_info %>%
  mutate(
    Cancer = Primary_Site,
    Age_group = case_when(
      Age < 40 ~ "<40",
      Age >= 40 & Age < 60 ~ "40-60",
      Age >= 60 ~ "60+",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Sample, Cancer, everything())

save(ccle_info_fine, file = "../data/ccle_info_fine.rda", compress = "xz")

# CCLE ABSOLUTE
ccle_absolute <- data.frame(
  Sample = ccle_info$Sample,
  CCLE_name = ccle_info$CCLE_name,
  purity = runif(n_ccle, 0.5, 1.0),
  ploidy = pmax(1, rnorm(n_ccle, 2.5, 0.8)),
  genome_doublings = sample(0:3, n_ccle, replace = TRUE),
  stringsAsFactors = FALSE
)

save(ccle_absolute, file = "../data/ccle_absolute.rda", compress = "xz")

# ============================================================================
# Print summary
# ============================================================================

cat("Built-in datasets created successfully:\n")
cat("=====================================\n")
data_files <- list.files("../data", pattern = "\\.rda$")
for (f in data_files) {
  info <- file.info(file.path("../data", f))
  cat(sprintf("- %s (%.1f KB)\n", f, info$size / 1024))
}
