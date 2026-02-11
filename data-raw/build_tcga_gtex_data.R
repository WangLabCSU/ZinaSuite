# Build Correct TCGA-GTEx Sample Info for ZinaSuite
# This script creates tcga_gtex from UCSCXenaShiny real data

library(dplyr)

# Load UCSCXenaShiny data
load("../../data/tcga_gtex.rda")

cat("UCSCXenaShiny tcga_gtex:", nrow(tcga_gtex), "samples\n")

# Convert to ZinaSuite format
tcga_gtex_zina <- data.frame(
  Sample = as.character(tcga_gtex$sample),
  Tissue = tcga_gtex$tissue,
  Type = tcga_gtex$type,
  Dataset = tcga_gtex$type2,
  stringsAsFactors = FALSE
)

# Remove rows with NA Sample
tcga_gtex_zina <- tcga_gtex_zina[!is.na(tcga_gtex_zina$Sample), ]

cat("ZinaSuite tcga_gtex:", nrow(tcga_gtex_zina), "samples\n")
cat("First sample ID:", head(tcga_gtex_zina$Sample, 1), "\n")
cat("ID length:", nchar(head(tcga_gtex_zina$Sample, 1)), "\n")

# Check TCGA samples
tcga_samples <- tcga_gtex_zina$Sample[grep("^TCGA", tcga_gtex_zina$Sample)]
cat("TCGA samples:", length(tcga_samples), "\n")

# Check GTEx samples
gtex_samples <- tcga_gtex_zina$Sample[grep("^GTEX", tcga_gtex_zina$Sample)]
cat("GTEx samples:", length(gtex_samples), "\n")

# Save
tcga_gtex <- tcga_gtex_zina
save(tcga_gtex, file = "../data/tcga_gtex.rda", compress = "xz")

cat("\nDone! tcga_gtex created successfully.\n")
