# Build Correct TCGA Clinical Data for ZinaSuite
# This script creates tcga_clinical and tcga_clinical_fine from UCSCXenaShiny real data

library(dplyr)

# Load UCSCXenaShiny data
load("../../data/tcga_clinical.rda")
load("../../data/tcga_surv.rda")

cat("UCSCXenaShiny tcga_clinical:", nrow(tcga_clinical), "samples\n")
cat("UCSCXenaShiny tcga_surv:", nrow(tcga_surv), "samples\n")

# ============================================================================
# Create tcga_clinical for ZinaSuite (simplified format)
# ============================================================================

# Map UCSCXenaShiny columns to ZinaSuite format
tcga_clinical_zina <- data.frame(
  Sample = tcga_clinical$sample,
  Cancer = tcga_clinical$type,
  Gender = tcga_clinical$gender,
  Age = tcga_clinical$age_at_initial_pathologic_diagnosis,
  Race = tcga_clinical$race,
  Ethnicity = NA_character_,  # Not available in UCSCXenaShiny data
  Stage_ajcc = tcga_clinical$ajcc_pathologic_tumor_stage,
  Stage_pathologic = NA_character_,  # Not directly available
  Grade = tcga_clinical$histological_grade,
  Vital_status = tcga_clinical$vital_status,
  stringsAsFactors = FALSE
)

# Clean up Stage_ajcc
tcga_clinical_zina$Stage_ajcc <- gsub("^Stage ", "Stage ", tcga_clinical_zina$Stage_ajcc)
tcga_clinical_zina$Stage_ajcc <- gsub("^[ABC]$", "", tcga_clinical_zina$Stage_ajcc)

# Remove rows with NA Sample
 tcga_clinical_zina <- tcga_clinical_zina[!is.na(tcga_clinical_zina$Sample), ]

# Keep only TCGA samples (15 character IDs)
tcga_clinical_zina <- tcga_clinical_zina[nchar(tcga_clinical_zina$Sample) == 15, ]

cat("ZinaSuite tcga_clinical:", nrow(tcga_clinical_zina), "samples\n")

# Save
tcga_clinical <- tcga_clinical_zina
save(tcga_clinical, file = "../data/tcga_clinical.rda", compress = "xz")

# ============================================================================
# Create tcga_clinical_fine (cleaned for grouping)
# ============================================================================

tcga_clinical_fine <- tcga_clinical %>%
  mutate(
    Stage_ajcc = case_when(
      grepl("Stage I[^V]", Stage_ajcc) ~ "Stage I",
      grepl("Stage II", Stage_ajcc) ~ "Stage II",
      grepl("Stage III", Stage_ajcc) ~ "Stage III",
      grepl("Stage IV", Stage_ajcc) ~ "Stage IV",
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

cat("ZinaSuite tcga_clinical_fine:", nrow(tcga_clinical_fine), "samples\n")

# Save
save(tcga_clinical_fine, file = "../data/tcga_clinical_fine.rda", compress = "xz")

# ============================================================================
# Create tcga_surv (survival data)
# ============================================================================

# Merge survival data with cancer type
tcga_surv_zina <- tcga_surv %>%
  inner_join(
    tcga_clinical %>% select(Sample, Cancer),
    by = c("sample" = "Sample")
  ) %>%
  select(Sample = sample, Cancer, everything())

cat("ZinaSuite tcga_surv:", nrow(tcga_surv_zina), "samples\n")

# Save
tcga_surv <- tcga_surv_zina
save(tcga_surv, file = "../data/tcga_surv.rda", compress = "xz")

cat("\nDone! Clinical data created successfully.\n")
cat("Sample ID format check:\n")
cat("  tcga_clinical:", head(tcga_clinical$Sample, 1), "\n")
cat("  tcga_clinical_fine:", head(tcga_clinical_fine$Sample, 1), "\n")
cat("  tcga_surv:", head(tcga_surv$Sample, 1), "\n")
