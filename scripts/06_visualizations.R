# Visualizations
# Based on Boutin et al. (2024) doi:10.1139/cjm-2023-0215

library(phyloseq)
library(vegan)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# --------------------------------------------------
# Load phyloseq objects (built by 03_permanova.R)
# --------------------------------------------------
ps_16S <- readRDS("results/tables/ps_16S.rds")
ps_ITS <- readRDS("results/tables/ps_ITS.rds")

FIG_PATH <- "results/figures"
dir.create(FIG_PATH, recursive = TRUE, showWarnings = FALSE)

# Theme for plots
theme_phyllosphere <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Color palettes
site_colors <- c("A" = "#E69F00", "B" = "#56B4E9", "C" = "#009E73")
time_colors <- c("May" = "#66c2a5", "July" = "#fc8d62", "August" = "#8da0cb")
compartment_colors <- c("Flower" = "#e78ac3", "Leaf" = "#a6d854")

# --------------------------------------------------
# PCoA ordinations
# --------------------------------------------------
plot_pcoa <- function(ps, color_var, shape_var = NULL, title) {
  # VST transformation
  dds <- phyloseq_to_deseq2(ps, ~ 1)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- estimateDispersions(dds, fitType = "local")
  vst_counts <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
  vst_counts[vst_counts < 0] <- 0
  
  # PCoA
  bc_dist <- vegdist(t(vst_counts), method = "bray")
  pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
  
  # Calculate variance explained
  eig <- pcoa_res$eig
  var_explained <- round(eig / sum(eig[eig > 0]) * 100, 1)
  
  # Build plot data
  plot_df <- data.frame(
    PCo1 = pcoa_res$points[, 1],
    PCo2 = pcoa_res$points[, 2],
    as.data.frame(sample_data(ps))
  )
  
  p <- ggplot(plot_df, aes(x = .data$PCo1, y = .data$PCo2, color = .data[[color_var]]))
  
  if (!is.null(shape_var)) {
    p <- p + geom_point(aes(shape = .data[[shape_var]]), size = 3, alpha = 0.8)
  } else {
    p <- p + geom_point(size = 3, alpha = 0.8)
  }
  
  p <- p +
    stat_ellipse(level = 0.95, linetype = "dashed") +
    labs(
      x = paste0("PCo1 (", var_explained[1], "%)"),
      y = paste0("PCo2 (", var_explained[2], "%)"),
      title = title
    ) +
    theme_phyllosphere
  
  return(p)
}

# May samples: Flowers vs Leaves
ps_16S_may <- subset_samples(ps_16S, Time == "May")
ps_ITS_may <- subset_samples(ps_ITS, Time == "May")

p_16S_may <- plot_pcoa(ps_16S_may, "Compartment", "Site", "Bacteria (May)")
p_ITS_may <- plot_pcoa(ps_ITS_may, "Compartment", "Site", "Fungi (May)")

# Leaves by time
ps_16S_leaves <- subset_samples(ps_16S, Compartment == "Leaf")
ps_ITS_leaves <- subset_samples(ps_ITS, Compartment == "Leaf")

p_16S_time <- plot_pcoa(ps_16S_leaves, "Time", "Site", "Bacteria (Leaves)")
p_ITS_time <- plot_pcoa(ps_ITS_leaves, "Time", "Site", "Fungi (Leaves)")

# Save PCoA plots
ggsave(file.path(FIG_PATH, "pcoa_may_bacteria.pdf"), p_16S_may, width = 8, height = 6)
ggsave(file.path(FIG_PATH, "pcoa_may_fungi.pdf"), p_ITS_may, width = 8, height = 6)
ggsave(file.path(FIG_PATH, "pcoa_leaves_bacteria.pdf"), p_16S_time, width = 8, height = 6)
ggsave(file.path(FIG_PATH, "pcoa_leaves_fungi.pdf"), p_ITS_time, width = 8, height = 6)

# --------------------------------------------------
# Relative abundance barplots (top 20 genera)
# --------------------------------------------------
plot_relative_abundance <- function(ps, title) {
  # Aggregate to genus
  ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
  
  # Transform to relative abundance
  ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x) * 100)
  
  # Melt to long format
  melt_df <- psmelt(ps_rel)
  
  # Get top 20 genera
  top20 <- melt_df %>%
    group_by(Genus) %>%
    summarize(mean_abund = mean(Abundance)) %>%
    arrange(desc(mean_abund)) %>%
    head(20) %>%
    pull(Genus)
  
  # Assign "Other" to non-top20
  melt_df <- melt_df %>%
    mutate(Genus_plot = ifelse(Genus %in% top20, Genus, "Other"))
  
  # Summarize by Site, Time, Cultivar
  summary_df <- melt_df %>%
    group_by(Site, Time, Cultivar, Genus_plot) %>%
    summarize(Abundance = mean(Abundance), .groups = "drop")
  
  # Plot
  p <- ggplot(summary_df, aes(x = Cultivar, y = Abundance, fill = Genus_plot)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(Site ~ Time) +
    labs(
      x = "Cultivar",
      y = "Relative Abundance (%)",
      fill = "Genus",
      title = title
    ) +
    theme_phyllosphere +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

p_abund_16S <- plot_relative_abundance(ps_16S_leaves, "Bacterial Community Composition")
p_abund_ITS <- plot_relative_abundance(ps_ITS_leaves, "Fungal Community Composition")

ggsave(file.path(FIG_PATH, "abundance_bacteria.pdf"), p_abund_16S, width = 12, height = 10)
ggsave(file.path(FIG_PATH, "abundance_fungi.pdf"), p_abund_ITS, width = 12, height = 10)

# --------------------------------------------------
# Alpha diversity plots (requires 05_alpha_diversity.R to have run first)
# --------------------------------------------------
alpha_16S_file <- "results/tables/alpha_diversity_16S_leaves.csv"
alpha_ITS_file <- "results/tables/alpha_diversity_ITS_leaves.csv"

if (file.exists(alpha_16S_file) && file.exists(alpha_ITS_file)) {
  div_16S <- read.csv(alpha_16S_file)
  div_ITS <- read.csv(alpha_ITS_file)
  
  plot_alpha <- function(div_df, title) {
    ggplot(div_df, aes(x = Time, y = Shannon, color = Cultivar)) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
      facet_wrap(~ Site) +
      labs(
        x = "Time",
        y = "Shannon Diversity",
        title = title
      ) +
      theme_phyllosphere
  }
  
  p_alpha_16S <- plot_alpha(div_16S, "Bacterial Alpha Diversity")
  p_alpha_ITS <- plot_alpha(div_ITS, "Fungal Alpha Diversity")
  
  ggsave(file.path(FIG_PATH, "alpha_diversity_bacteria.pdf"), p_alpha_16S, width = 10, height = 6)
  ggsave(file.path(FIG_PATH, "alpha_diversity_fungi.pdf"), p_alpha_ITS, width = 10, height = 6)
} else {
  message("Alpha diversity CSVs not found — run 05_alpha_diversity.R first for those plots")
}

message("All figures saved to ", FIG_PATH)
