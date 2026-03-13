# Phyllosphere Microbiome Methods Replication

Replication of the computational methods from Tembrock et al. (2024), with an investigation into how reference database updates affect genus-level taxonomy assignments.

**Original paper:** Tembrock et al. (2024). "Investigating the spatiotemporal dynamics of apple tree phyllosphere bacterial and fungal communities across cultivars in orchards." *Canadian Journal of Microbiology*. [DOI: 10.1139/cjm-2023-0215](https://doi.org/10.1139/cjm-2023-0215)

## What this repo contains

1. **Full DADA2 pipelines** for 16S (bacterial) and ITS (fungal) amplicon sequences
2. **Statistical analyses**: PERMANOVA, ANCOM-BC2 differential abundance, alpha diversity (Shannon index with linear mixed models)
3. **Taxonomy comparison**: ASV-level sequence matching between our results and the authors' original phyloseq objects, with rank-by-rank taxonomy comparison across database versions

## Key finding

The pipeline replicates cleanly. Per-sample sequence counts correlate at r > 0.999, ASV counts are within 1.5%, and all statistical conclusions hold.

However, running taxonomy assignment against current reference databases (SILVA v138.1, UNITE 2025-02-19) instead of the versions used in the paper (SILVA v132, UNITE v7.2) causes 43% of differentially abundant genera to carry different names. The underlying biology is identical; only the labels change. Notable reclassifications include *Cryptococcus* (split into 4 genera), *Lactobacillus* (split into 11 genera), and *Mycosphaerella* (partially absorbed into *Cladosporium*).

See `results/figures/taxonomy_sankey.png` for a visual summary.

## Original data sources

| Resource | DOI/Accession |
|----------|---------------|
| Raw sequences (NCBI SRA) | [PRJNA1081804](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1081804) |
| 16S R scripts | [10.6084/m9.figshare.24570064](https://doi.org/10.6084/m9.figshare.24570064) |
| ITS R scripts | [10.6084/m9.figshare.24570055](https://doi.org/10.6084/m9.figshare.24570055) |
| 16S data (phyloseq object) | [10.6084/m9.figshare.24651720](https://doi.org/10.6084/m9.figshare.24651720) |
| ITS data (phyloseq object) | [10.6084/m9.figshare.24651945](https://doi.org/10.6084/m9.figshare.24651945) |
| Sample metadata | [10.6084/m9.figshare.25300561](https://doi.org/10.6084/m9.figshare.25300561) |

## Reproducing the analysis

### Prerequisites

- R 4.5+ with [renv](https://rstudio.github.io/renv/)
- ~10 GB disk space for raw data and reference databases
- cutadapt (for ITS primer removal)
- SRA Toolkit (for downloading raw sequences)

### Setup

```bash
# Clone the repo
git clone https://github.com/YOUR_USERNAME/phyllosphere_methods.git
cd phyllosphere_methods

# Restore R dependencies
Rscript -e "renv::restore()"

# Download raw sequences from SRA
bash scripts/00_download_data.sh

# Download reference databases:
# - SILVA v138.1: https://zenodo.org/record/4587955
# - UNITE general release: https://unite.ut.ee/repository.php
# Place in data/reference/
```

### Running the pipeline

Scripts are numbered and should be run in order:

```bash
Rscript scripts/01_dada2_16S.R        # 16S ASV inference + taxonomy
Rscript scripts/02_dada2_ITS.R        # ITS ASV inference + taxonomy
Rscript scripts/03_permanova.R        # Community composition (PERMANOVA)
Rscript scripts/04_differential_abundance.R  # ANCOM-BC2
Rscript scripts/05_alpha_diversity.R   # Shannon diversity + LMMs
Rscript scripts/06_visualizations.R    # PCoA, abundance, diversity plots
```

## Directory structure

```
phyllosphere_methods/
├── scripts/                # Numbered analysis scripts
├── data/
│   ├── metadata/           # Sample metadata and SRA mappings
│   ├── raw/                # Raw FASTQ files (not tracked, download from SRA)
│   ├── processed/          # DADA2 intermediate outputs (generated)
│   └── reference/          # SILVA and UNITE databases (not tracked)
├── results/
│   ├── figures/            # PCoA, abundance, diversity, and Sankey plots
│   └── tables/             # Phyloseq objects, ANCOM-BC2 results, diversity tables
├── renv.lock               # R package versions for reproducibility
└── README.md
```

## Pipeline details

### Sequence processing (DADA2)

- **16S**: Primers 799F/1115R, trimming forward 210 bp / reverse 160 bp, taxonomy via SILVA v138.1
- **ITS**: Primers ITS1F/ITS2, cutadapt primer removal, taxonomy via UNITE general release (2025-02-19), `tryRC = TRUE`
- Both: chimera removal, ASV filtering (>=10 sequences across all samples)

### Statistical analyses

- **PERMANOVA**: DESeq2 variance-stabilizing transformation, Bray-Curtis dissimilarity, `adonis2` with 999 permutations
- **Differential abundance**: ANCOM-BC2 at genus level, testing Time and Site effects on leaf samples
- **Alpha diversity**: Shannon index via 100x rarefaction, linear mixed models (`lmer`) with TreeID as random effect, `emmeans` for pairwise contrasts

## R dependencies

All dependencies are pinned in `renv.lock`. Key packages:

- dada2, phyloseq, Biostrings
- DESeq2, ANCOMBC, vegan
- lme4, emmeans, rstatix
- tidyverse, ggplot2

## License

This repo contains analysis code only. Raw sequence data is available from NCBI SRA under BioProject [PRJNA1081804](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1081804). The original study data and scripts are available from Figshare (see links above).
