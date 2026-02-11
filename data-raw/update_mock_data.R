# Update Mock Data with Correct Sample IDs
# This script updates all mock datasets to use the same sample ID format as real data

library(dplyr)

# Load the real tcga_clinical to get correct sample IDs
load("../data/tcga_clinical.rda")
tcga_samples <- tcga_clinical$Sample
n_tcga <- length(tcga_samples)

cat("Using", n_tcga, "TCGA samples from real clinical data\n")
cat("First sample:", head(tcga_samples, 1), "\n")

# ============================================================================
# Update TCGA TIL Data
# ============================================================================
cat("\n=== Updating tcga_TIL ===\n")

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
set.seed(42)
for (cell in cell_types) {
  tcga_TIL[[cell]] <- runif(n_tcga, 0, 0.5)
}

# Normalize to sum to 1
tcga_TIL[, cell_types] <- tcga_TIL[, cell_types] / rowSums(tcga_TIL[, cell_types])

save(tcga_TIL, file = "../data/tcga_TIL.rda", compress = "xz")
cat("Saved tcga_TIL with", nrow(tcga_TIL), "samples\n")

# ============================================================================
# Update TCGA TMB Data
# ============================================================================
cat("\n=== Updating tcga_tmb ===\n")

set.seed(42)
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
cat("Saved tcga_tmb with", nrow(tcga_tmb), "samples\n")

# ============================================================================
# Update TCGA MSI Data
# ============================================================================
cat("\n=== Updating tcga_MSI ===\n")

set.seed(42)
tcga_MSI <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  msi_score = runif(n_tcga, 0, 1),
  msi_status = sample(c("MSI-H", "MSI-L", "MSS"), n_tcga, replace = TRUE, prob = c(0.1, 0.2, 0.7)),
  stringsAsFactors = FALSE
)

save(tcga_MSI, file = "../data/tcga_MSI.rda", compress = "xz")
cat("Saved tcga_MSI with", nrow(tcga_MSI), "samples\n")

# ============================================================================
# Update TCGA Stemness Data
# ============================================================================
cat("\n=== Updating tcga_stemness ===\n")

set.seed(42)
tcga_stemness <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  rnaseq = rnorm(n_tcga, 0, 1),
  methylation = rnorm(n_tcga, 0, 1),
  stringsAsFactors = FALSE
)

save(tcga_stemness, file = "../data/tcga_stemness.rda", compress = "xz")
cat("Saved tcga_stemness with", nrow(tcga_stemness), "samples\n")

# ============================================================================
# Update TCGA Genome Instability Data
# ============================================================================
cat("\n=== Updating tcga_genome_instability ===\n")

set.seed(42)
tcga_genome_instability <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  LOH = runif(n_tcga, 0, 0.5),
  HRD = runif(n_tcga, 0, 0.5),
  LST = runif(n_tcga, 0, 0.5),
  ntel = runif(n_tcga, 0, 0.5),
  stringsAsFactors = FALSE
)

save(tcga_genome_instability, file = "../data/tcga_genome_instability.rda", compress = "xz")
cat("Saved tcga_genome_instability with", nrow(tcga_genome_instability), "samples\n")

# ============================================================================
# Update TCGA Purity Data
# ============================================================================
cat("\n=== Updating tcga_purity ===\n")

set.seed(42)
tcga_purity <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  purity = runif(n_tcga, 0.2, 1.0),
  ploidy = pmax(1, rnorm(n_tcga, 2, 0.5)),
  source = sample(c("ABSOLUTE", "ESTIMATE", "LUMP", "CPE"), n_tcga, replace = TRUE),
  stringsAsFactors = FALSE
)

save(tcga_purity, file = "../data/tcga_purity.rda", compress = "xz")
cat("Saved tcga_purity with", nrow(tcga_purity), "samples\n")

# ============================================================================
# Update TCGA Subtypes Data
# ============================================================================
cat("\n=== Updating tcga_subtypes ===\n")

set.seed(42)
tcga_subtypes <- data.frame(
  Sample = tcga_samples,
  Cancer = tcga_clinical$Cancer,
  subtype = sample(c("Basal", "Classical", "Mesenchymal", "Proneural", NA), n_tcga, replace = TRUE),
  stringsAsFactors = FALSE
)

save(tcga_subtypes, file = "../data/tcga_subtypes.rda", compress = "xz")
cat("Saved tcga_subtypes with", nrow(tcga_subtypes), "samples\n")

cat("\n=== All mock data updated successfully! ===\n")
