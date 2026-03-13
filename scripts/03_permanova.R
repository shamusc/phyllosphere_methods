# PERMANOVA Analysis for Community Composition
# Based on Boutin et al. (2024) doi:10.1139/cjm-2023-0215

library(phyloseq)
library(vegan)
library(DESeq2)
library(tidyverse)

# --------------------------------------------------
# Load data and build phyloseq objects
# --------------------------------------------------
seqtab_16S <- readRDS("data/processed/16S/seqtab_16S.rds")
taxa_16S <- readRDS("data/processed/16S/taxa_16S.rds")
seqtab_ITS <- readRDS("data/processed/ITS/seqtab_ITS.rds")
taxa_ITS <- readRDS("data/processed/ITS/taxa_ITS.rds")
metadata <- read.csv("data/metadata/sample_metadata.csv")

OUT_PATH <- "results/tables"
dir.create(OUT_PATH, recursive = TRUE, showWarnings = FALSE)

# Build metadata keyed by SRA_Accession (matches seqtab rownames)
build_sample_data <- function(meta, amplicon_type) {
  meta_amp <- meta[meta$Amplicon == amplicon_type, ]
  rownames(meta_amp) <- meta_amp$SRA_Accession
  return(meta_amp)
}

meta_16S <- build_sample_data(metadata, "16S")
meta_ITS <- build_sample_data(metadata, "ITS")

# Subset seqtabs to only samples present in metadata (excludes controls)
shared_16S <- intersect(rownames(seqtab_16S), rownames(meta_16S))
shared_ITS <- intersect(rownames(seqtab_ITS), rownames(meta_ITS))

if (length(shared_16S) != nrow(meta_16S)) {
  stop(paste("16S sample mismatch: seqtab has", length(shared_16S),
             "of", nrow(meta_16S), "metadata samples"))
}
if (length(shared_ITS) != nrow(meta_ITS)) {
  stop(paste("ITS sample mismatch: seqtab has", length(shared_ITS),
             "of", nrow(meta_ITS), "metadata samples"))
}

message("16S: ", length(shared_16S), " samples (", 
        nrow(seqtab_16S) - length(shared_16S), " controls excluded)")
message("ITS: ", length(shared_ITS), " samples (",
        nrow(seqtab_ITS) - length(shared_ITS), " controls excluded)")

ps_16S <- phyloseq(
  otu_table(seqtab_16S[shared_16S, ], taxa_are_rows = FALSE),
  sample_data(meta_16S[shared_16S, ]),
  tax_table(taxa_16S)
)

ps_ITS <- phyloseq(
  otu_table(seqtab_ITS[shared_ITS, ], taxa_are_rows = FALSE),
  sample_data(meta_ITS[shared_ITS, ]),
  tax_table(taxa_ITS)
)

# Save phyloseq objects for downstream scripts
saveRDS(ps_16S, file.path(OUT_PATH, "ps_16S.rds"))
saveRDS(ps_ITS, file.path(OUT_PATH, "ps_ITS.rds"))
message("Saved phyloseq objects: ps_16S.rds, ps_ITS.rds")

# --------------------------------------------------
# Variance stabilizing transformation (DESeq2)
# --------------------------------------------------
vst_transform <- function(ps) {
  dds <- phyloseq_to_deseq2(ps, ~ 1)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- estimateDispersions(dds, fitType = "local")
  vst_counts <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
  vst_counts[vst_counts < 0] <- 0
  return(vst_counts)
}

message("Applying variance stabilizing transformation...")
vst_16S <- vst_transform(ps_16S)
vst_ITS <- vst_transform(ps_ITS)

# Calculate Bray-Curtis distances
bc_16S <- vegdist(t(vst_16S), method = "bray")
bc_ITS <- vegdist(t(vst_ITS), method = "bray")

# --------------------------------------------------
# PERMANOVA: May samples (flowers + leaves)
# --------------------------------------------------
run_permanova_may <- function(dist_mat, ps, label) {
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  may_idx <- which(meta$Time == "May")
  dist_may <- as.dist(as.matrix(dist_mat)[may_idx, may_idx])
  meta_may <- data.frame(meta[may_idx, ], stringsAsFactors = FALSE)
  
  message("\n=== PERMANOVA: ", label, " (May - Flowers & Leaves) ===")
  
  # Full model
  perm_full <- adonis2(
    dist_may ~ Cultivar * Site * Compartment,
    data = meta_may,
    permutations = 999
  )
  print(perm_full)
  
  # Flowers only
  flower_idx <- which(meta_may$Compartment == "Flower")
  if (length(flower_idx) > 0) {
    dist_flower <- as.dist(as.matrix(dist_may)[flower_idx, flower_idx])
    meta_flower <- data.frame(meta_may[flower_idx, ], stringsAsFactors = FALSE)
    
    perm_flower <- adonis2(
      dist_flower ~ Cultivar * Site,
      data = meta_flower,
      permutations = 999
    )
    message("\n--- Flowers only ---")
    print(perm_flower)
  }
  
  # Leaves only
  leaf_idx <- which(meta_may$Compartment == "Leaf")
  if (length(leaf_idx) > 0) {
    dist_leaf <- as.dist(as.matrix(dist_may)[leaf_idx, leaf_idx])
    meta_leaf <- data.frame(meta_may[leaf_idx, ], stringsAsFactors = FALSE)
    
    perm_leaf <- adonis2(
      dist_leaf ~ Cultivar * Site,
      data = meta_leaf,
      permutations = 999
    )
    message("\n--- Leaves only (May) ---")
    print(perm_leaf)
  }
  
  return(list(full = perm_full))
}

# --------------------------------------------------
# PERMANOVA: Leaf samples across season
# --------------------------------------------------
run_permanova_leaves <- function(dist_mat, ps, label) {
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  leaf_idx <- which(meta$Compartment == "Leaf")
  dist_leaf <- as.dist(as.matrix(dist_mat)[leaf_idx, leaf_idx])
  meta_leaf <- data.frame(meta[leaf_idx, ], stringsAsFactors = FALSE)
  
  message("\n=== PERMANOVA: ", label, " (Leaves across season) ===")
  
  perm_leaf <- adonis2(
    dist_leaf ~ Time * Cultivar * Site,
    data = meta_leaf,
    permutations = 999,
    strata = meta_leaf$TreeID
  )
  print(perm_leaf)
  
  return(perm_leaf)
}

# Run analyses
results_16S_may <- run_permanova_may(bc_16S, ps_16S, "Bacteria")
results_ITS_may <- run_permanova_may(bc_ITS, ps_ITS, "Fungi")
results_16S_leaves <- run_permanova_leaves(bc_16S, ps_16S, "Bacteria")
results_ITS_leaves <- run_permanova_leaves(bc_ITS, ps_ITS, "Fungi")

# Save results
saveRDS(
  list(
    bacteria_may = results_16S_may,
    fungi_may = results_ITS_may,
    bacteria_leaves = results_16S_leaves,
    fungi_leaves = results_ITS_leaves
  ),
  file.path(OUT_PATH, "permanova_results.rds")
)

message("\nResults saved to ", OUT_PATH)
