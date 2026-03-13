# Alpha Diversity Analysis
# Based on Boutin et al. (2024) doi:10.1139/cjm-2023-0215

library(phyloseq)
library(vegan)
library(lme4)
library(emmeans)
library(rstatix)
library(tidyverse)

# --------------------------------------------------
# Load phyloseq objects (built by 03_permanova.R)
# --------------------------------------------------
ps_16S <- readRDS("results/tables/ps_16S.rds")
ps_ITS <- readRDS("results/tables/ps_ITS.rds")

OUT_PATH <- "results/tables"
dir.create(OUT_PATH, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------
# Rarefied Shannon diversity (100 iterations)
# --------------------------------------------------
calculate_rarefied_shannon <- function(ps, depth, n_iter = 100) {
  message("Calculating rarefied Shannon diversity (", n_iter, " iterations, depth=", depth, ")...")
  
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  
  # Filter samples with sufficient depth
  sufficient <- rowSums(otu) >= depth
  n_dropped <- sum(!sufficient)
  if (n_dropped > 0) {
    message("  Dropping ", n_dropped, " samples with <", depth, " reads")
  }
  otu <- otu[sufficient, ]
  message("  Rarefying ", nrow(otu), " samples...")
  
  shannon_list <- vector("list", n_iter)
  
  for (i in seq_len(n_iter)) {
    rarefied <- rrarefy(otu, sample = depth)
    shannon_list[[i]] <- diversity(rarefied, index = "shannon")
  }
  
  # Calculate mean Shannon across iterations
  shannon_mat <- do.call(rbind, shannon_list)
  shannon_mean <- colMeans(shannon_mat)
  shannon_sd <- apply(shannon_mat, 2, sd)
  
  return(data.frame(
    SampleID = names(shannon_mean),
    Shannon = shannon_mean,
    Shannon_SD = shannon_sd,
    row.names = NULL
  ))
}

# Rarefaction depths from paper
# Leaves: 16S = 3900, ITS = 4900
# Flowers: 16S = 6900, ITS = 3600

ps_16S_leaves <- subset_samples(ps_16S, Compartment == "Leaf")
ps_ITS_leaves <- subset_samples(ps_ITS, Compartment == "Leaf")
ps_16S_flowers <- subset_samples(ps_16S, Compartment == "Flower")
ps_ITS_flowers <- subset_samples(ps_ITS, Compartment == "Flower")

shannon_16S_leaves <- calculate_rarefied_shannon(ps_16S_leaves, depth = 3900)
shannon_ITS_leaves <- calculate_rarefied_shannon(ps_ITS_leaves, depth = 4900)
shannon_16S_flowers <- calculate_rarefied_shannon(ps_16S_flowers, depth = 6900)
shannon_ITS_flowers <- calculate_rarefied_shannon(ps_ITS_flowers, depth = 3600)

# Merge diversity with metadata
# Extract all metadata once as a plain data.frame, then merge by SRA accession
all_meta <- data.frame(sample_data(ps_16S), stringsAsFactors = FALSE)
all_meta$SampleID <- rownames(all_meta)
all_meta_ITS <- data.frame(sample_data(ps_ITS), stringsAsFactors = FALSE)
all_meta_ITS$SampleID <- rownames(all_meta_ITS)

div_16S_leaves <- merge(shannon_16S_leaves, all_meta, by = "SampleID")
div_ITS_leaves <- merge(shannon_ITS_leaves, all_meta_ITS, by = "SampleID")
div_16S_flowers <- merge(shannon_16S_flowers, all_meta, by = "SampleID")
div_ITS_flowers <- merge(shannon_ITS_flowers, all_meta_ITS, by = "SampleID")

# --------------------------------------------------
# Linear mixed models for leaves
# --------------------------------------------------
message("\n=== Linear Mixed Models: Bacterial Diversity ===")

# Bacteria: Shannon ~ Site * Time + (1|TreeID)
lmm_16S <- lmer(Shannon ~ Site * Time + (1 | TreeID), data = div_16S_leaves)
summary(lmm_16S)

# Pairwise comparisons
emm_16S <- emmeans(lmm_16S, ~ Site * Time)
pairs_16S <- pairs(emm_16S, adjust = "tukey")
print(pairs_16S)

message("\n=== Linear Mixed Models: Fungal Diversity ===")

# Fungi: Shannon ~ Site * Time + Site * Cultivar + (1|TreeID)
lmm_ITS <- lmer(Shannon ~ Site * Time + Site * Cultivar + (1 | TreeID), data = div_ITS_leaves)
summary(lmm_ITS)

emm_ITS <- emmeans(lmm_ITS, ~ Site * Time)
pairs_ITS <- pairs(emm_ITS, adjust = "tukey")
print(pairs_ITS)

# --------------------------------------------------
# Kruskal-Wallis + Dunn's test (non-parametric)
# --------------------------------------------------
message("\n=== Non-parametric tests: Cultivar differences ===")

run_kruskal_dunn <- function(div_df, var_name) {
  results_list <- list()
  
  for (site in unique(div_df$Site)) {
    for (time in unique(div_df$Time)) {
      subset_df <- div_df %>%
        filter(Site == site, Time == time)
      
      if (nrow(subset_df) >= 6 && length(unique(subset_df$Cultivar)) >= 2) {
        kw <- kruskal_test(subset_df, Shannon ~ Cultivar)
        
        if (kw$p < 0.05) {
          dunn <- dunn_test(subset_df, Shannon ~ Cultivar, p.adjust.method = "holm")
          dunn$Site <- site
          dunn$Time <- time
          dunn$Variable <- var_name
          results_list[[paste(site, time, sep = "_")]] <- dunn
        }
      }
    }
  }
  
  if (length(results_list) > 0) {
    return(bind_rows(results_list))
  } else {
    message("  No significant Kruskal-Wallis results for ", var_name)
    return(NULL)
  }
}

dunn_16S <- run_kruskal_dunn(div_16S_leaves, "Bacteria")
dunn_ITS <- run_kruskal_dunn(div_ITS_leaves, "Fungi")

# Save results
write.csv(div_16S_leaves, file.path(OUT_PATH, "alpha_diversity_16S_leaves.csv"), row.names = FALSE)
write.csv(div_ITS_leaves, file.path(OUT_PATH, "alpha_diversity_ITS_leaves.csv"), row.names = FALSE)
write.csv(div_16S_flowers, file.path(OUT_PATH, "alpha_diversity_16S_flowers.csv"), row.names = FALSE)
write.csv(div_ITS_flowers, file.path(OUT_PATH, "alpha_diversity_ITS_flowers.csv"), row.names = FALSE)

saveRDS(
  list(
    lmm_bacteria = lmm_16S,
    lmm_fungi = lmm_ITS,
    emmeans_bacteria = emm_16S,
    emmeans_fungi = emm_ITS,
    dunn_bacteria = dunn_16S,
    dunn_fungi = dunn_ITS
  ),
  file.path(OUT_PATH, "alpha_diversity_models.rds")
)

message("\nResults saved to ", OUT_PATH)
