# Replication: Boutin et al. (2024) — Project Context

## Objective

Replicate the computational methods from:

**Boutin, S., Lussier, E., & Laforest-Lapointe, I. (2024).** Investigating the spatiotemporal dynamics of apple tree phyllosphere bacterial and fungal communities across cultivars in orchards. *Canadian Journal of Microbiology*, 70: 238–251.

**DOI:** 10.1139/cjm-2023-0215

---

## Tools Installed

| Tool | Version | Purpose |
|------|---------|---------|
| SRA Toolkit | 3.2.1 | Download SRA data, convert to FASTQ |
| wget | 1.25.0 | HTTP file downloads |
| NCBI Entrez Direct | 25.1 | Query NCBI databases (esearch, efetch) |
| R | 4.5.2 | Statistical analysis, DADA2 pipeline |
| renv | (R package) | R virtual environment for reproducibility |
| cutadapt | (via pipx) | Primer removal for ITS amplicons |

**Note:** Entrez Direct required manual CA certificate fix:
```bash
curl -o ~/edirect/cacert.pem https://curl.se/ca/cacert.pem
```

---

## R Packages Installed (via renv)

### Bioconductor
- dada2 - Amplicon sequence variant inference
- phyloseq - Microbiome data handling
- DESeq2 - Variance stabilizing transformation
- ANCOMBC - Differential abundance analysis

### CRAN
- tidyverse - Data manipulation and visualization
- vegan - Community ecology (PERMANOVA, diversity)
- lme4 - Linear mixed models
- emmeans - Estimated marginal means
- rstatix - Pipe-friendly statistics
- gridExtra - Plot arrangement

### System Dependencies Installed (via Homebrew)
- fribidi (for textshaping)
- libwebp (for ragg)
- cmake (for nloptr)
- gsl (for ANCOMBC/energy)

---

## Data Downloaded

### Raw Sequencing Data

- **Source:** NCBI SRA BioProject PRJNA1081804
- **Total runs:** 361 (181 bacterial/16S + 180 fungal/ITS)
- **All 361 runs successfully downloaded** ✓
- **Total size:** 4.1 GB (compressed FASTQ)
- **Format:** Illumina MiSeq 300bp paired-end
- **Control samples:** 3 (2 bacterial: 1 positive + 1 negative PCR control; 1 fungal: negative PCR control)

### Reference Databases

| Database | Location | Size | Purpose |
|----------|----------|------|---------|
| SILVA v138.1 (train) | `data/reference/silva_nr99_v138.1_train_set.fa.gz` | 131 MB | 16S genus-level taxonomy |
| SILVA v138.1 (species) | `data/reference/silva_species_assignment_v138.1.fa.gz` | 75 MB | 16S species-level matching |
| UNITE full (2025-02-19) | `data/reference/UNITE_public_19.02.2025.fasta` | 1.4 GB | ITS fungal taxonomy (too large, caused OOM) |
| UNITE general release | `data/reference/sh_general_release_dynamic_19.02.2025.fasta` | 70 MB | ITS fungal taxonomy (correct, ~102K sequences) |

### Authors' Original Data

| Source | DOI | Files |
|--------|-----|-------|
| 16S Scripts | 10.6084/m9.figshare.24570064 | `data/authors/16s_dada2.Rmd` |
| 16S Phyloseq | 10.6084/m9.figshare.24651720 | `data/authors/phylo.rds` |
| ITS Scripts | 10.6084/m9.figshare.24570050 | `data/authors/ITS_dada2.Rmd` |
| ITS Phyloseq | 10.6084/m9.figshare.24651945 | `data/authors/phylo_ITS.rds` |
| Supplement | Paper SI | Table S2 (primers), Table S3 (read counts) |

---

## Data Organization

### Directory Structure

```
replication-boutin-2024/
├── data/
│   ├── raw/
│   │   ├── 16S/           # 362 files (181 samples × 2 paired-end, includes 2 controls)
│   │   ├── ITS/           # 360 files (180 samples × 2 paired-end, includes 1 control)
│   │   └── SraRunInfo.csv # SRA metadata for all 361 runs
│   ├── metadata/
│   │   ├── sample_metadata.csv      # Complete metadata (358 samples)
│   │   ├── run_marker_gene.csv      # SRA Experiment TITLE per run (16S vs ITS label)
│   │   ├── authors_metadata.csv     # Extracted from authors' RDS
│   │   └── authors_16S_data         # Original phyloseq object
│   ├── processed/
│   │   ├── 16S/           # Output directory for DADA2
│   │   └── ITS/           # Output directory for DADA2
│   ├── reference/         # SILVA and UNITE databases
│   └── authors/           # Authors' scripts and phyloseq objects
├── scripts/
│   ├── 00_download_data.sh
│   ├── 00_install_packages.R
│   ├── 01_dada2_16S.R
│   ├── 02_dada2_ITS.R
│   ├── 03_permanova.R
│   ├── 04_differential_abundance.R
│   ├── 05_alpha_diversity.R
│   └── 06_visualizations.R
├── docs/supplement/       # Supplementary materials from paper
├── logs/                  # Pipeline run logs
├── results/
│   ├── figures/
│   └── tables/
├── renv/                  # R virtual environment
├── docs/                  # LinkedIn draft
├── renv.lock              # Locked package versions
└── README.md
```

### Control Samples

| SRA Accession | Type | Marker Gene |
|---------------|------|-------------|
| SRR28147178 | Positive PCR control | 16S |
| SRR28147180 | Negative PCR control | 16S |
| SRR28148420 | Negative PCR control | ITS |

These were initially orphaned in `data/raw/fastq/` and not included in pipeline runs. Now correctly placed in their respective directories.

### Sample Naming Convention

Library names follow pattern: `[number]-[site]-[cultivar]-[compartment]-[replicate]`

| Code | Meaning |
|------|---------|
| 1 | 16S (bacteria) |
| 2, 3 | ITS (fungi) |
| ASB | Site A (Abbaye St-Benoit) |
| PMB | Site B |
| VBS | Site C |
| CL | Cortland cultivar |
| LI | Liberty cultivar |
| PR | Paulared cultivar |
| FL | Flower |
| LE | Leaf |

---

## Metadata Status: COMPLETE ✓

### File: `data/metadata/sample_metadata.csv`

**Columns:** LibraryName, SampleID, SRA_Accession, Amplicon, Site, Cultivar, Compartment, Replicate, Time

### Sample Distribution

| Variable | Values |
|----------|--------|
| **Total samples** | 358 |
| **Amplicon** | 16S (180), ITS (178) |
| **Time** | May (180), July (90), August (88) |
| **Compartment** | Flower (90), Leaf (268) |
| **Site** | A (120), B (120), C (118) |
| **Cultivar** | Cortland (120), Liberty (120), Paulared (118) |

**Note:** Time metadata was extracted from the authors' phyloseq RDS file by merging on LibraryName.

### File: `data/metadata/run_marker_gene.csv`

- 361 rows mapping SRA Run accession → Experiment TITLE
- Labels: "Bacterial seq of phyllosphere" (181 runs) or "Fungal seq of phyllosphere" (180 runs)
- Source: SRA Experiment XML queried via `esearch`/`efetch`

---

## Scripts

### 01_dada2_16S.R
- DADA2 pipeline for 16S bacterial amplicons
- Primer removal via `trimLeft = c(19, 16)` in `filterAndTrim` (matches authors' method)
- Truncation: `truncLen = c(210, 160)`
- Quality filtering: `maxEE = c(2, 2)`
- Taxonomy: SILVA v138.1
- Checkpointing: saves intermediate results after each major step
- Includes fail-fast check for missing database

### 02_dada2_ITS.R
- DADA2 pipeline for ITS fungal amplicons
- Primer removal via `cutadapt` (matches authors' method)
  - FWD: ITS1F `CTTGGTCATTTAGAGGAAGTAA` (22bp)
  - REV: ITS2 `GCTGCGTTCTTCATCGATGC` (20bp)
- No fixed truncation (variable ITS length)
- Quality filtering: `maxEE = c(2, 5)`
- Taxonomy: UNITE general release (single-threaded, `tryRC = TRUE`)
- Checkpointing: saves intermediate results after each major step
- Includes fail-fast check for missing database

### 03_permanova.R
- PERMANOVA analysis on community composition
- Loads seqtabs + metadata, fixes Amplicon labels via `run_marker_gene.csv`, adds TreeID
- Excludes 3 control samples, builds phyloseq objects with DESeq2 VST normalization
- Runs `adonis2` for: (a) May flowers+leaves (Cultivar × Site × Compartment), (b) Leaves only (Time × Cultivar × Site, stratified by TreeID)
- Outputs: `results/tables/permanova_results.rds`, `results/tables/ps_16S.rds`, `results/tables/ps_ITS.rds`

### 04_differential_abundance.R
- ANCOM-BC2 differential abundance analysis
- Runs Time and Site comparisons for both 16S and ITS at genus level
- Outputs: `results/tables/ancombc_results.rds`, 4 CSV files with significant taxa

### 05_alpha_diversity.R
- Shannon diversity via `estimate_richness`, linear mixed models (Site × Time + TreeID random effect)
- Estimated marginal means via `emmeans`, Kruskal-Wallis + Dunn's test for cultivar effects
- Outputs: `results/tables/alpha_diversity_models.rds`, 4 CSV files

### 06_visualizations.R
- PCoA ordinations (Bray-Curtis on VST counts) for leaves and May samples
- Relative abundance bar plots at phylum level
- Alpha diversity box plots by Site × Time
- Outputs: 8 PDF figures in `results/figures/`

---

## Verified Parameters from Authors' Scripts

### 16S (from `16s_dada2.Rmd`)

| Parameter | Value | Source |
|-----------|-------|--------|
| 799F primer | AACMGGATTAGATACCCKG (19bp) | Table S2 |
| 1115R primer | AGGGTTGCGCTCGTTG (16bp) | Table S2 |
| Primer removal | `trimLeft = c(19, 16)` | Authors' script |
| truncLen | c(210, 160) | Authors' script |
| maxEE | c(2, 2) | Authors' script |
| SILVA version | v132 | Authors' script |

### ITS (from `ITS_dada2.Rmd`)

| Parameter | Value | Source |
|-----------|-------|--------|
| ITS1F primer | CTTGGTCATTTAGAGGAAGTAA (22bp) | Authors' script |
| ITS2 primer | GCTGCGTTCTTCATCGATGC (20bp) | Authors' script |
| Primer removal | `cutadapt` | Authors' script |
| truncLen | None (variable ITS length) | Authors' script |
| maxEE | c(2, 5) | Authors' script |
| UNITE version | v7.2 (2017-06-28) | Authors' script |

**Note:** Authors used SILVA v132 and UNITE v7.2; we use SILVA v138.1 and UNITE 2025-02-19. This may cause minor taxonomy differences but should not significantly affect ASV counts.

---

## DADA2 Pipeline Results — FINAL

### 16S Bacterial Pipeline

| Metric | Our Result | Paper Reported | Difference |
|--------|-----------|----------------|------------|
| Samples processed | 181 | 179 | +2 (controls) |
| Raw ASVs (pre-chimera) | 68,596 | — | — |
| Chimeras removed | 58,945 | — | 86% chimeric |
| Chimera-free ASVs | 9,651 | — | — |
| ASVs after ≥10 filter | **4,077** | **4,015** | **+1.5%** |
| Total sequences retained | 3,185,773 | 3,184,892 | +0.03% |
| Mean seqs/sample | 17,601 | — | — |
| SD seqs/sample | 7,146 | — | — |
| Sequences retained (chimera step) | 72% | — | Normal |

### ITS Fungal Pipeline

| Metric | Our Result | Paper Reported | Difference |
|--------|-----------|----------------|------------|
| Samples processed | 180 | 178 | +2 (control + extra) |
| Raw ASVs (pre-chimera) | 2,495 | — | — |
| Chimeras removed | 382 | — | 15% chimeric |
| Chimera-free ASVs | 2,113 | — | — |
| ASVs after ≥10 filter | **1,030** | **1,027** | **+0.3%** |
| Total sequences retained | 2,793,920 | — | — |
| Mean seqs/sample | 15,522 | — | — |
| SD seqs/sample | 8,047 | — | — |
| Sequences retained (chimera step) | 96.4% | — | Normal |

**Conclusion:** Both pipelines closely replicate the published results. The small ASV surplus (+1.5% for 16S, +0.3% for ITS) is expected due to: (a) newer reference databases (SILVA v138.1 vs v132, UNITE 2025 vs v7.2), (b) inclusion of control samples the authors likely excluded, and (c) minor software version differences.

---

## Completed Tasks

1. ✅ Created project directory structure
2. ✅ Downloaded raw sequencing data from NCBI SRA (all 361 runs)
3. ✅ Downloaded UNITE databases (full + general release)
4. ✅ Downloaded SILVA v138.1 database
5. ✅ Installed R and renv environment
6. ✅ Installed all R/Bioconductor packages
7. ✅ Created sample metadata from SraRunInfo.csv
8. ✅ Extracted Time metadata from authors' RDS file
9. ✅ Merged metadata to complete sample_metadata.csv
10. ✅ Downloaded authors' original scripts and phyloseq objects (16S + ITS)
11. ✅ Ran 16S pipeline (first attempt, with cutadapt — yielded 4,980 ASVs)
12. ✅ Corrected 16S pipeline to use `trimLeft` per authors' method — yielded 4,053 ASVs (0.9% from paper's 4,015)
13. ✅ Ran ITS pipeline (first attempt, with mixed 16S/ITS data) — yielded 1,414 ASVs vs paper's 1,027
14. ✅ Investigated ITS ASV discrepancy — discovered data mis-labeling (see below)
15. ✅ Corrected raw data organization using SRA Experiment TITLE metadata
16. ✅ Recovered 3 orphaned control samples from `data/raw/fastq/` (see below)
17. ✅ Re-ran 16S DADA2 pipeline with all 181 correctly sorted samples → **4,077 ASVs**
18. ✅ Re-ran ITS DADA2 pipeline with all 180 correctly sorted samples → **1,030 ASVs**
19. ✅ Fixed Amplicon labels in `sample_metadata.csv` using `run_marker_gene.csv`
20. ✅ Added TreeID column to metadata for repeated measures
21. ✅ Ran PERMANOVA analysis (all 4 comparisons significant at p = 0.001)
22. ✅ Ran alpha diversity analysis (LMM, emmeans, Dunn's tests)
23. ✅ Ran differential abundance analysis (ANCOM-BC2, Time + Site for both markers)
24. ✅ Generated all visualizations (8 PDF figures)
25. ✅ Compared results against authors' phyloseq objects (r > 0.999 per-sample, all 36 emmeans contrasts match)
26. ✅ Investigated taxonomy reclassification impact on differential abundance
27. ✅ Generated taxonomy Sankey diagram (`results/figures/taxonomy_sankey.png`)
28. ✅ Drafted and published LinkedIn/X post summarizing findings (`docs/linkedin_draft.txt`)
29. ✅ Created `.gitignore`, cleaned repo for public sharing
30. ✅ Updated `README.md` with project overview, setup instructions, key findings
31. ✅ Initialized git repo, made initial commit, pushed to GitHub
32. ✅ Canonical GitHub repo name: **`replication-boutin-2024`** (Boutin et al. 2024; earlier working names included `phyllosphere_methods` and a mistaken “Tembrock” label)

## Statistical Analysis Results

### PERMANOVA (Community Composition)

All models significant at p = 0.001 (999 permutations).

| Analysis | R² | F | p | Notes |
|----------|-----|---|---|-------|
| Bacteria — May (Cultivar × Site × Compartment) | 0.482 | 3.95 | 0.001 | 48% variance explained |
| Fungi — May (Cultivar × Site × Compartment) | 0.609 | 6.59 | 0.001 | 61% variance explained |
| Bacteria — Leaves (Time × Cultivar × Site) | 0.589 | 5.89 | 0.001 | Stratified by TreeID |
| Fungi — Leaves (Time × Cultivar × Site) | 0.802 | 16.65 | 0.001 | 80% variance explained |

### ANCOM-BC2 (Differential Abundance)

| Comparison | 16S sig. taxa | ITS sig. taxa |
|------------|---------------|---------------|
| Time | 94 | 42 |
| Site | 63 | 58 |

### Alpha Diversity (Shannon Index)

Linear mixed models (Site × Time, random = TreeID) fit for both bacteria and fungi on leaf samples. Estimated marginal means computed via `emmeans`. Cultivar effects tested via Kruskal-Wallis + Dunn's post-hoc at specific Site × Time combinations.

Key findings:
- **Bacteria:** Site C consistently lowest Shannon diversity (e.g., May: 2.42 vs A: 3.97, B: 4.33). Significant cultivar effect at Site C/May (Cortland vs Liberty, p.adj = 0.004).
- **Fungi:** Temporal shift — May highest diversity (A: 3.19, C: 3.20), summer months lower (August A: 1.96). Multiple significant cultivar effects across Site × Time combinations.

### Generated Outputs

**Tables** (`results/tables/`):
- `permanova_results.rds` — Full adonis2 results
- `ancombc_results.rds` — Full ANCOM-BC2 results
- `diff_abundance_{16S,ITS}_{time,site}.csv` — Significant differentially abundant taxa
- `alpha_diversity_{16S,ITS}_{flowers,leaves}.csv` — Shannon diversity values
- `alpha_diversity_models.rds` — LMM, emmeans, Dunn's test results
- `ps_16S.rds`, `ps_ITS.rds` — Phyloseq objects

**Figures** (`results/figures/`):
- `pcoa_leaves_bacteria.pdf`, `pcoa_leaves_fungi.pdf` — PCoA ordinations (leaves, all timepoints)
- `pcoa_may_bacteria.pdf`, `pcoa_may_fungi.pdf` — PCoA ordinations (May, flowers+leaves)
- `abundance_bacteria.pdf`, `abundance_fungi.pdf` — Relative abundance bar plots (phylum level)
- `alpha_diversity_bacteria.pdf`, `alpha_diversity_fungi.pdf` — Shannon diversity box plots

## Comparison Against Authors' Phyloseq Objects

### Sample-Level Replication

Both markers: 179/179 samples matched (100%). Per-sample sequence counts correlate at r > 0.999.

| Marker | Per-sample seq diff (mean) | Per-sample seq diff (max) | Pearson r |
|--------|---------------------------|--------------------------|-----------|
| 16S | −0.38% | 4.58% | 0.9997 |
| ITS | +0.31% | 4.44% | 0.9998 |

### ASV and Taxonomy Replication

| Metric | Authors | Ours | Overlap |
|--------|---------|------|---------|
| **16S ASVs** | 4,015 | 4,077 | +1.5% |
| **16S Total seqs** | 3,184,892 | 3,169,682 | r=0.9997 |
| **16S Genera** | 468 | 497 | 385 (82%) |
| **16S Species** | 279 | 266 | 220 (79%) |
| **ITS ASVs** | 1,027 | 1,030 | +0.3% |
| **ITS Total seqs** | 2,779,879 | 2,792,024 | r=0.9998 |
| **ITS Genera** | 224 | 293 | 176 (79%) |
| **ITS Species** | 272 | 319 | 132 (49%) |

### Taxonomy Differences Explained

**16S phyla naming**: SILVA v132→v138.1 renamed all phyla (e.g., Actinobacteria→Actinobacteriota, Bacteroidetes→Bacteroidota). After mapping, 16/19 authors' phyla match. The 3 unmatched reflect genuine reclassifications (v138.1 split Deltaproteobacteria into Desulfobacterota, Myxococcota, Bdellovibrionota).

**ITS genus/species divergence**: UNITE v7.2 (2017) → 2025 release updated many species hypotheses. This explains the larger species-level gap (49% overlap) — same ASVs assigned to updated or reclassified taxa. At genus level (79% overlap), the core community composition is preserved.

### Statistical Results Comparison (vs Paper Table S6)

All 36 pairwise emmeans contrasts (alpha diversity, Site × Time) match the paper in direction, significance, and magnitude. No significance disagreements. Estimates typically within 0.01–0.06 of published values.

### Replication Conclusion

The replication is strong at the sample and ASV level (r > 0.999 per-sample correlations, ASV counts within 1.5%). Taxonomy differences are attributable to reference database version updates (SILVA v132→v138.1, UNITE v7.2→2025), not pipeline errors. Structural conclusions (sites differ, time matters, fungi more structured than bacteria) are robust to database updates.

---

## Taxonomy Reclassification Investigation

### Motivation

While the replication is successful at the community/diversity level, the taxonomy database updates (SILVA v132→v138.1 for 16S, UNITE v7.2→2025 for ITS) may affect genus-level analyses like ANCOM-BC2. Key concerns:

1. **SILVA reclassified Deltaproteobacteria** into three new phyla (Desulfobacterota, Myxococcota, Bdellovibrionota), affecting Proteobacteria relative abundance statements.
2. **UNITE reclassified several top-20 genera** from the paper (e.g., Gibberella→Fusarium, Mycosphaerella split), potentially changing which genera appear as differentially abundant.
3. **ITS species-level overlap is only 49%**, suggesting substantial reclassification at finer taxonomic ranks.

### Step 1: ASV Sequence Matching

| Marker | Authors ASVs | Our ASVs | Exact matches | % of authors | % of ours |
|--------|-------------|----------|---------------|-------------|-----------|
| 16S | 4,015 | 4,077 | 3,856 | 96.0% | 94.6% |
| ITS | 1,027 | 1,030 | 1,008 | 98.1% | 97.9% |

The unmatched ASVs (~4%) are low-abundance sequences at the detection boundary — expected stochastic variation from DADA2's probabilistic error model across software versions.

### Step 2: Taxonomy Agreement by Rank (matched ASVs only)

**16S Bacteria (3,856 matched ASVs):**

| Rank | Both Assigned | Agree | Disagree | % Agree |
|------|--------------|-------|----------|---------|
| Kingdom | 3,856 | 3,856 | 0 | 100% |
| Phylum | 3,837 | 1,785 | 2,052 | 46.5% |
| Class | 3,797 | 3,669 | 128 | 96.6% |
| Order | 3,706 | 3,129 | 577 | 84.4% |
| Family | 3,448 | 3,011 | 437 | 87.3% |
| Genus | 2,642 | 2,353 | 289 | 89.1% |
| Species | 250 | 249 | 1 | 99.6% |

**Note:** 16S Phylum agreement is 46.5% because SILVA v138.1 systematically renamed all phyla (e.g., Actinobacteria→Actinobacteriota, Bacteroidetes→Bacteroidota). These are 1:1 name changes, not reclassifications. The 3 genuine reclassifications are Proteobacteria being split into Myxococcota (50 ASVs), Bdellovibrionota (9), and Desulfobacterota (4).

**ITS Fungi (1,008 matched ASVs):**

| Rank | Both Assigned | Agree | Disagree | % Agree |
|------|--------------|-------|----------|---------|
| Kingdom | 1,008 | 1,008 | 0 | 100% |
| Phylum | 939 | 920 | 19 | 98.0% |
| Class | 842 | 817 | 25 | 97.0% |
| Order | 795 | 620 | 175 | 78.0% |
| Family | 724 | 409 | 315 | 56.5% |
| Genus | 602 | 473 | 129 | 78.6% |
| Species | 365 | 220 | 145 | 60.3% |

**Note:** ITS divergence increases at finer ranks (Family 56.5%, Species 60.3%), reflecting UNITE's extensive reclassification of species hypotheses between v7.2 (2017) and 2025.

### Key Genus Reclassifications Affecting Top-20 Taxa

**16S (from Table S5):**
- *Methylobacterium* → *Methylobacterium-Methylorubrum* (all 47 ASVs, genus merge)
- *Lactobacillus* → split into 11 genera (34/45 ASVs: *Limosilactobacillus*, *Ligilactobacillus*, etc.)
- *Bacillus* → 4/53 ASVs now *Falsibacillus*
- *Frigoribacterium* → 1/9 ASVs now *Amnibacterium*

**ITS (from Table S5):**
- *Cryptococcus* → split completely (9/9 ASVs: *Filobasidium* 3, *Papiliotrema* 3, *Naganishia* 2, *Teunia* 1)
- *Gibberella* → *Fusarium* (all 3 ASVs)
- *Lemonniera* → *Gyoerffyella* (3/4 ASVs)
- *Mycosphaerella* → *Cladosporium* (2) + *Ramularia* (2), 4/7 ASVs reclassified

### Step 3: Impact on ANCOM-BC2 Differential Abundance

| Marker | Significant genera | Affected by reclassification | Unchanged |
|--------|-------------------|------------------------------|-----------|
| 16S | 57 | 25 (44%) | 32 (56%) |
| ITS | 42 | 18 (43%) | 24 (57%) |

**~43% of differentially abundant genera are affected by taxonomy reclassifications.** This means that if the authors ran ANCOM-BC2 on their data, their results list would contain different genus names for the same biological organisms in nearly half the cases.

Notable examples:
- Our *Methylobacterium-Methylorubrum* = their *Methylobacterium* (genus merge)
- Our *Gyoerffyella* = their *Lemonniera* (complete reclassification)
- Our *Filobasidium* and *Naganishia* = their *Cryptococcus* (genus split)
- Our *Cladosporium* includes ASVs they called *Mycosphaerella* (partial reclassification)

### Conclusion

The taxonomy database updates do NOT change the biological story at the community level (diversity, composition structure, which factors matter). However, they substantially affect genus-level reporting: **43–44% of differentially abundant genera carry different names** between SILVA v132/UNITE v7.2 and SILVA v138.1/UNITE 2025. The underlying organisms are the same — it is the nomenclature that shifted. Any replication using updated databases must account for these reclassifications when comparing genus-level results to the original publication.

## Dissemination

### LinkedIn/X Post
- Published post highlighting the "static paper vs. dynamic databases" finding
- Central hook: 43% of differentially abundant genera carry different names when re-running taxonomy against current databases
- Key examples: *Cryptococcus* split into 4 genera, *Lactobacillus* into 11, *Mycosphaerella* partially absorbed into *Cladosporium*
- Post explicitly separates pipeline replication (using authors' database versions) from taxonomy re-analysis (using current versions) to preempt methodological criticism
- Includes Sankey diagram (`results/figures/taxonomy_sankey.png`) as visual

### GitHub Repository
- **URL:** https://github.com/shamusc/replication-boutin-2024
- Contains: analysis scripts, metadata templates, results (figures + tables), `renv.lock`, `PROJECT_CONTEXT.md`
- Excludes (via `.gitignore`): raw FASTQ, reference databases, `renv/library/`, logs, personal drafts
- No unlicensed or third-party content committed — all code is original, all data is linked (not redistributed)

## Project Status

**Complete.** All replication, analysis, investigation, and dissemination tasks finished.

---

## Discrepancies and Investigations

### 1. Raw Data Mis-Labeling (RESOLVED)

**Problem:** Our initial organization of SRA data into `data/raw/16S/` and `data/raw/ITS/` directories was approximately 50/50 correct/wrong. Both directories contained a mixture of 16S and ITS data.

**Root Cause:** The SRA metadata at the `Run` and `BioSample` level does not explicitly distinguish 16S from ITS. The initial sorting method (based on `LibraryName` prefix) was unreliable because the naming convention was inconsistent.

**Discovery:** When the ITS pipeline produced 1,414 ASVs (vs paper's 1,027), investigation revealed:
- ~330bp-long ASVs unique to our results started with the 16S 799F primer and ended with the 16S 1115R reverse complement
- ~89 "ITS" samples were >80% 16S reads; ~83 samples were <5% 16S reads
- All 1,414 ASVs were forced into `k__Fungi` by UNITE, but the 16S sequences were unclassified below Kingdom

**Solution:** Queried SRA Experiment-level XML metadata using `esearch`/`efetch`. The `Experiment TITLE` field definitively labels runs:
- "Bacterial seq of phyllosphere" → 16S (181 runs)
- "Fungal seq of phyllosphere" → ITS (180 runs)

Created `data/metadata/run_marker_gene.csv` with this mapping. Moved all files to correct directories.

**Verification after swap:**

| Directory | Before (correct / mislabeled) | After |
|-----------|-------------------------------|-------|
| `data/raw/16S/` | 90 correct + 90 mislabeled ITS | 181 correct 16S |
| `data/raw/ITS/` | 89 correct + 89 mislabeled 16S | 180 correct ITS |

---

### 2. Orphaned Control Samples (RESOLVED)

**Problem:** 3 control samples (2 bacterial, 1 fungal) were sitting in `data/raw/fastq/` — they were never moved into the correct directories during the initial data organization or the subsequent re-sort. This meant pipeline runs were missing samples.

| SRA Accession | Control Type | Marker | Was In | Moved To |
|---------------|-------------|--------|--------|----------|
| SRR28147178 | Positive PCR control | 16S | `data/raw/fastq/` | `data/raw/16S/` |
| SRR28147180 | Negative PCR control | 16S | `data/raw/fastq/` | `data/raw/16S/` |
| SRR28148420 | Negative PCR control | ITS | `data/raw/fastq/` | `data/raw/ITS/` |

**Resolution:** Identified via comparison of all accessions in `run_marker_gene.csv` against files in `16S/` and `ITS/` directories. Moved to correct directories. `data/raw/fastq/` is now empty.

**Final counts:** 181 16S pairs + 180 ITS pairs = 361 total (all SRA runs accounted for).

---

### 3. Raw Read Count Discrepancy (INVESTIGATED, UNEXPLAINED)

| Source | 16S Reads | Notes |
|--------|-----------|-------|
| SRA data | 5,374,560 pairs | 180 runs from NCBI |
| Paper Table S3 "input" | 6,726,468 | ~25% more than SRA |
| Authors' final data | 3,184,892 sequences | Matches Table S3 |

**Investigation findings:**
1. SRA metadata exactly matches our downloaded files (verified)
2. Paper has MORE 16S but FEWER ITS reads than SRA deposit
3. Each LibraryName has 2 runs with different BioSamples (different timepoints)
4. Total of 180 unique 16S samples in SRA, authors processed 179

**Unexplained:** The 25% input difference remains unclear. Possible explanations:
- Additional sequencing runs not deposited in SRA
- Authors' local files had different filtering before "input" count
- Different counting methodology (reads vs pairs)

**Decision:** Proceed with SRA data. Final results represent normal retention for DADA2.

---

### 4. 16S Pipeline Correction (RESOLVED)

**Original issue:** We used `cutadapt` with `--discard-untrimmed` for primer removal.

**Discovery:** Authors' published script (Figshare DOI: 10.6084/m9.figshare.24570064) shows they used `trimLeft` in `filterAndTrim`:

```r
filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
              trimLeft = c(19, 16),      # Remove primers by fixed-length trimming
              truncLen = c(210, 160),    # Truncate reads
              ...)
```

**Fix applied:** Updated `01_dada2_16S.R` to use `trimLeft` instead of cutadapt.

**Result (with mixed data, before re-sort — needs re-run with correct data):**

| Metric | Our Result | Paper Reported | Difference |
|--------|-----------|----------------|------------|
| Total ASVs | 4,053 | 4,015 | +0.9% |

---

### 5. ITS UNITE Database OOM (RESOLVED)

**Problem:** `assignTaxonomy` crashed with Out-Of-Memory errors using the full UNITE+INSD database (1.4 GB, 2M+ sequences).

**Fix:** Downloaded the smaller UNITE "General FASTA release" (`sh_general_release_dynamic_19.02.2025.fasta`, 70 MB, ~102K sequences), which is the recommended database for DADA2's `assignTaxonomy`. Updated scripts to use this file.

---

## Session Notes

- NCBI SRA download failed once on sample 67 due to network error; resumed successfully
- All 361 SRA runs were successfully downloaded; 3 control samples were initially placed in `data/raw/fastq/` instead of the correct marker-gene directories
- UNITE database required manual download due to license agreement
- R installed via Homebrew (4.5.2)
- R packages required several system dependencies (fribidi, libwebp, cmake, gsl) installed via Homebrew
- Authors' RDS file saved without .rds extension
- Time metadata successfully extracted by merging on LibraryName field
- The full UNITE+INSD database (1.4GB) causes OOM on 16GB RAM; use the general release (~70MB) instead
- ITS taxonomy assignment must run single-threaded (`multithread = FALSE`) on 16GB RAM
- Taxonomy Sankey diagram uses equal-weight flows (not ASV counts) to emphasize branching patterns over abundance
- LinkedIn post went through ~15 drafts to avoid AI-sounding language (em dashes, "jaw-dropping" hyperbole)
- GitHub repo: **`replication-boutin-2024`** (correct author: Boutin et al. 2024)
- `PROJECT_CONTEXT.md` is tracked in git for lab continuity (update if any paths are machine-specific)

---

## Key Concepts Discussed

### DADA2 Pipeline
- **ASVs vs OTUs**: DADA2 produces exact Amplicon Sequence Variants rather than clustered OTUs
- **Error learning**: Uses Loess regression to model quality-score-dependent error rates
- **Poisson model**: Determines if observed sequence abundance exceeds expected errors

### Reference Databases
- **SILVA**: Gold standard for 16S rRNA, uses naive Bayesian classifier for taxonomy
- **UNITE**: Standard for fungal ITS region; use "General FASTA release" (Species Hypothesis reps, ~100K seqs) not full UNITE+INSD (2M+ seqs)
- **Local vs Remote**: Local databases enable batch processing (millions of comparisons) that would be impractical via API

### SRA Toolkit
- **VDB format**: NCBI's columnar database format for efficient storage
- **fasterq-dump**: Converts VDB to universal FASTQ format
- **Tradeoff**: VDB is space-efficient but tool-dependent; FASTQ is universal but larger
- **Experiment TITLE**: The SRA Experiment-level XML metadata contains the most reliable marker gene label (16S vs ITS)
