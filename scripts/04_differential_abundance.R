# Differential Abundance Analysis with ANCOM-BC2
# Based on Boutin et al. (2024) doi:10.1139/cjm-2023-0215

library(phyloseq)
library(ANCOMBC)
library(tidyverse)

# --------------------------------------------------
# Load phyloseq objects (built by 03_permanova.R)
# --------------------------------------------------
ps_16S <- readRDS("results/tables/ps_16S.rds")
ps_ITS <- readRDS("results/tables/ps_ITS.rds")

OUT_PATH <- "results/tables"
FIG_PATH <- "results/figures"
dir.create(OUT_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_PATH, recursive = TRUE, showWarnings = FALSE)

# Aggregate to genus level
ps_16S_genus <- tax_glom(ps_16S, taxrank = "Genus", NArm = FALSE)
ps_ITS_genus <- tax_glom(ps_ITS, taxrank = "Genus", NArm = FALSE)

# Filter to leaf samples only
ps_16S_leaves <- subset_samples(ps_16S_genus, Compartment == "Leaf")
ps_ITS_leaves <- subset_samples(ps_ITS_genus, Compartment == "Leaf")

# --------------------------------------------------
# ANCOM-BC2: Time comparisons
# --------------------------------------------------
run_ancombc_time <- function(ps, label) {
  message("\n=== ANCOM-BC2 Time Comparison: ", label, " ===")
  
  result <- ancombc2(
    data = ps,
    fix_formula = "Time",
    p_adj_method = "holm",
    prv_cut = 0.1,
    lib_cut = 1000,
    group = "Time",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    global = TRUE,
    pairwise = TRUE,
    dunnet = FALSE,
    trend = FALSE
  )
  
  return(result)
}

# --------------------------------------------------
# ANCOM-BC2: Site comparisons
# --------------------------------------------------
run_ancombc_site <- function(ps, label) {
  message("\n=== ANCOM-BC2 Site Comparison: ", label, " ===")
  
  result <- ancombc2(
    data = ps,
    fix_formula = "Site",
    p_adj_method = "holm",
    prv_cut = 0.1,
    lib_cut = 1000,
    group = "Site",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    global = TRUE,
    pairwise = TRUE,
    dunnet = FALSE,
    trend = FALSE
  )
  
  return(result)
}

# Run analyses
ancom_16S_time <- run_ancombc_time(ps_16S_leaves, "Bacteria")
ancom_ITS_time <- run_ancombc_time(ps_ITS_leaves, "Fungi")
ancom_16S_site <- run_ancombc_site(ps_16S_leaves, "Bacteria")
ancom_ITS_site <- run_ancombc_site(ps_ITS_leaves, "Fungi")

# Extract significant results
extract_significant <- function(result, comparison_type) {
  res_df <- result$res
  
  # Get log fold changes and p-values
  lfc_cols <- grep("^lfc_", colnames(res_df), value = TRUE)
  pval_cols <- grep("^q_", colnames(res_df), value = TRUE)
  
  sig_results <- res_df %>%
    select(taxon, all_of(lfc_cols), all_of(pval_cols)) %>%
    pivot_longer(
      cols = -taxon,
      names_to = c("metric", "comparison"),
      names_pattern = "(.+)_(.+)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    filter(q < 0.05) %>%
    arrange(q)
  
  return(sig_results)
}

sig_16S_time <- extract_significant(ancom_16S_time, "time")
sig_ITS_time <- extract_significant(ancom_ITS_time, "time")
sig_16S_site <- extract_significant(ancom_16S_site, "site")
sig_ITS_site <- extract_significant(ancom_ITS_site, "site")

# Save results
write.csv(sig_16S_time, file.path(OUT_PATH, "diff_abundance_16S_time.csv"), row.names = FALSE)
write.csv(sig_ITS_time, file.path(OUT_PATH, "diff_abundance_ITS_time.csv"), row.names = FALSE)
write.csv(sig_16S_site, file.path(OUT_PATH, "diff_abundance_16S_site.csv"), row.names = FALSE)
write.csv(sig_ITS_site, file.path(OUT_PATH, "diff_abundance_ITS_site.csv"), row.names = FALSE)

saveRDS(
  list(
    bacteria_time = ancom_16S_time,
    fungi_time = ancom_ITS_time,
    bacteria_site = ancom_16S_site,
    fungi_site = ancom_ITS_site
  ),
  file.path(OUT_PATH, "ancombc_results.rds")
)

message("\nResults saved to ", OUT_PATH)
