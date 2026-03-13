# Install required packages for phyllosphere microbiome analysis
# Based on Boutin et al. (2024) doi:10.1139/cjm-2023-0215

# CRAN packages
cran_packages <- c(
  "tidyverse",
  "vegan",
  "ggplot2",
  "rstatix",
  "lme4",
  "emmeans",
  "ape",
  "gridExtra"
)

# Bioconductor packages
bioc_packages <- c(
  "dada2",
  "phyloseq",
  "DESeq2",
  "ANCOMBC"
)

# Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

message("All packages installed successfully.")
