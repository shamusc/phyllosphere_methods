# DADA2 Pipeline for 16S rRNA (Bacteria)
# Based on Boutin et al. (2024) doi:10.1139/cjm-2023-0215
# Primers: 799F-1115R (v5-v6 region, chloroplast-excluding)
#
# MATCHES AUTHORS' SCRIPT: Uses trimLeft for primer removal (not cutadapt)
# Source: https://doi.org/10.6084/m9.figshare.24570064 (16s_2022_dada2_21avril.Rmd)
# CHECKPOINT SAVES: Intermediate results saved after each expensive step

library(dada2)
library(tidyverse)

# Configuration
RAW_PATH <- "data/raw/16S"
OUT_PATH <- "data/processed/16S"
CHECKPOINT_PATH <- file.path(OUT_PATH, "checkpoints")
SILVA_DB <- "data/reference/silva_nr99_v138.1_train_set.fa.gz"
SILVA_SPECIES <- "data/reference/silva_species_assignment_v138.1.fa.gz"

# Primer lengths (from paper Table S2)
# 799F: AACMGGATTAGATACCCKG (19 bp) - Chelius & Triplett 2001
# 1115R: AGGGTTGCGCTCGTTG (16 bp) - Redford & Fierer 2009
PRIMER_F_LEN <- 19  # trimLeft for forward reads
PRIMER_R_LEN <- 16  # trimLeft for reverse reads

# Trimming parameters (from authors' script)
TRUNC_F <- 210  # Forward read truncation
TRUNC_R <- 160  # Reverse read truncation
MIN_SEQ_THRESHOLD <- 10  # Minimum sequences per ASV across dataset

# Helper function for tracking
getN <- function(x) sum(getUniques(x))

# Fail fast if SILVA database is missing
if (!file.exists(SILVA_DB)) {
  stop(
    "SILVA database not found at: ", SILVA_DB, "\n",
    "Run the download script or manually download from:\n",
    "https://zenodo.org/record/4587955"
  )
}

dir.create(OUT_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_PATH, "filtered"), recursive = TRUE, showWarnings = FALSE)

# List input files (SRA naming: _1.fastq.gz, _2.fastq.gz)
fnFs <- sort(list.files(RAW_PATH, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(RAW_PATH, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))

if (length(fnFs) == 0) {
  stop("No FASTQ files found in ", RAW_PATH)
}

# Extract sample names (SRR accessions)
sample_names <- sub("_1\\.fastq\\.gz$", "", basename(fnFs))
message("Found ", length(fnFs), " samples")

# Quality filtered output paths
filtFs <- file.path(OUT_PATH, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(OUT_PATH, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# ============================================================================
# STEP 1: Filter and trim (including primer removal via trimLeft)
# ============================================================================
# NOTE: trimLeft removes primers BEFORE truncLen is applied
# This matches the authors' approach exactly
checkpoint_filter <- file.path(CHECKPOINT_PATH, "01_filter_output.rds")

if (file.exists(checkpoint_filter)) {
  message("CHECKPOINT: Loading filter results from ", checkpoint_filter)
  out <- readRDS(checkpoint_filter)
} else {
  message("STEP 1: Filtering, trimming, and removing primers...")
  message("  trimLeft = c(", PRIMER_F_LEN, ", ", PRIMER_R_LEN, ") - removes primer sequences")
  message("  truncLen = c(", TRUNC_F, ", ", TRUNC_R, ") - truncates reads")
  
  out <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    trimLeft = c(PRIMER_F_LEN, PRIMER_R_LEN),  # Remove primers (authors' method)
    truncLen = c(TRUNC_F, TRUNC_R),
    maxN = 0,
    maxEE = c(2, 2),
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
  )
  saveRDS(out, checkpoint_filter)
  message("CHECKPOINT SAVED: ", checkpoint_filter)
}

# ============================================================================
# STEP 2: Learn error rates
# ============================================================================
checkpoint_errF <- file.path(CHECKPOINT_PATH, "02_errF.rds")
checkpoint_errR <- file.path(CHECKPOINT_PATH, "02_errR.rds")

if (file.exists(checkpoint_errF) && file.exists(checkpoint_errR)) {
  message("CHECKPOINT: Loading error models...")
  errF <- readRDS(checkpoint_errF)
  errR <- readRDS(checkpoint_errR)
} else {
  message("STEP 2: Learning error rates...")
  
  # Only use filtered files that exist
  filtFs_exist <- filtFs[file.exists(filtFs)]
  filtRs_exist <- filtRs[file.exists(filtRs)]
  
  errF <- learnErrors(filtFs_exist, multithread = TRUE)
  saveRDS(errF, checkpoint_errF)
  message("CHECKPOINT SAVED: ", checkpoint_errF)
  
  errR <- learnErrors(filtRs_exist, multithread = TRUE)
  saveRDS(errR, checkpoint_errR)
  message("CHECKPOINT SAVED: ", checkpoint_errR)
}

# Save error plots
pdf(file.path(OUT_PATH, "error_rates_forward.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf(file.path(OUT_PATH, "error_rates_reverse.pdf"))
plotErrors(errR, nominalQ = TRUE)
dev.off()

# ============================================================================
# STEP 3: Dereplicate sequences
# ============================================================================
checkpoint_derep <- file.path(CHECKPOINT_PATH, "03_derep.rds")

# Only use filtered files that exist
filtFs_exist <- filtFs[file.exists(filtFs)]
filtRs_exist <- filtRs[file.exists(filtRs)]
sample_names_exist <- names(filtFs_exist)

if (file.exists(checkpoint_derep)) {
  message("CHECKPOINT: Loading dereplicated sequences...")
  derep_data <- readRDS(checkpoint_derep)
  derepFs <- derep_data$derepFs
  derepRs <- derep_data$derepRs
} else {
  message("STEP 3: Dereplicating sequences...")
  derepFs <- derepFastq(filtFs_exist)
  derepRs <- derepFastq(filtRs_exist)
  names(derepFs) <- sample_names_exist
  names(derepRs) <- sample_names_exist
  saveRDS(list(derepFs = derepFs, derepRs = derepRs), checkpoint_derep)
  message("CHECKPOINT SAVED: ", checkpoint_derep)
}

# ============================================================================
# STEP 4: Sample inference (DADA2 core algorithm)
# ============================================================================
checkpoint_dada <- file.path(CHECKPOINT_PATH, "04_dada.rds")

if (file.exists(checkpoint_dada)) {
  message("CHECKPOINT: Loading DADA2 inference results...")
  dada_data <- readRDS(checkpoint_dada)
  dadaFs <- dada_data$dadaFs
  dadaRs <- dada_data$dadaRs
} else {
  message("STEP 4: Running DADA2 sample inference...")
  dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
  dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
  saveRDS(list(dadaFs = dadaFs, dadaRs = dadaRs), checkpoint_dada)
  message("CHECKPOINT SAVED: ", checkpoint_dada)
}

# ============================================================================
# STEP 5: Merge paired reads
# ============================================================================
checkpoint_merged <- file.path(CHECKPOINT_PATH, "05_merged.rds")

if (file.exists(checkpoint_merged)) {
  message("CHECKPOINT: Loading merged reads...")
  mergers <- readRDS(checkpoint_merged)
} else {
  message("STEP 5: Merging paired reads...")
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
  saveRDS(mergers, checkpoint_merged)
  message("CHECKPOINT SAVED: ", checkpoint_merged)
}

# ============================================================================
# STEP 6: Construct sequence table
# ============================================================================
checkpoint_seqtab <- file.path(CHECKPOINT_PATH, "06_seqtab.rds")

if (file.exists(checkpoint_seqtab)) {
  message("CHECKPOINT: Loading sequence table...")
  seqtab <- readRDS(checkpoint_seqtab)
} else {
  message("STEP 6: Constructing sequence table...")
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, checkpoint_seqtab)
  message("CHECKPOINT SAVED: ", checkpoint_seqtab)
}
message("Sequence table dimensions: ", dim(seqtab)[1], " samples, ", dim(seqtab)[2], " ASVs")

# ============================================================================
# STEP 7: Remove chimeras
# ============================================================================
checkpoint_nochim <- file.path(CHECKPOINT_PATH, "07_seqtab_nochim.rds")

if (file.exists(checkpoint_nochim)) {
  message("CHECKPOINT: Loading chimera-free table...")
  seqtab_nochim <- readRDS(checkpoint_nochim)
} else {
  message("STEP 7: Removing chimeras...")
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
  saveRDS(seqtab_nochim, checkpoint_nochim)
  message("CHECKPOINT SAVED: ", checkpoint_nochim)
}
message("Chimera-free ASVs: ", ncol(seqtab_nochim))
message("Sequences retained: ", round(sum(seqtab_nochim) / sum(seqtab) * 100, 1), "%")

# ============================================================================
# STEP 8: Filter low-abundance ASVs
# ============================================================================
message("STEP 8: Filtering ASVs with <", MIN_SEQ_THRESHOLD, " total sequences...")
asv_sums <- colSums(seqtab_nochim)
seqtab_filtered <- seqtab_nochim[, asv_sums >= MIN_SEQ_THRESHOLD]
message("ASVs after filtering: ", ncol(seqtab_filtered))

# Save filtered sequence table immediately
saveRDS(seqtab_filtered, file.path(OUT_PATH, "seqtab_16S.rds"))
message("SAVED: seqtab_16S.rds")

# ============================================================================
# STEP 9: Track reads through pipeline
# ============================================================================
message("STEP 9: Creating read tracking table...")

# Match sample names between filter output and downstream objects
track_samples <- intersect(rownames(out), names(dadaFs))

track <- cbind(
  out[track_samples, ],
  sapply(dadaFs[track_samples], getN),
  sapply(dadaRs[track_samples], getN),
  sapply(mergers[track_samples], getN),
  rowSums(seqtab_nochim[track_samples, , drop = FALSE]),
  rowSums(seqtab_filtered[track_samples, , drop = FALSE])
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "final")
write.csv(track, file.path(OUT_PATH, "read_tracking.csv"))
message("SAVED: read_tracking.csv")

# ============================================================================
# STEP 10: Assign taxonomy with SILVA
# ============================================================================
checkpoint_taxa <- file.path(CHECKPOINT_PATH, "10_taxa.rds")

if (file.exists(checkpoint_taxa)) {
  message("CHECKPOINT: Loading taxonomy...")
  taxa <- readRDS(checkpoint_taxa)
} else {
  message("STEP 10: Assigning taxonomy with SILVA (this may take a while)...")
  taxa <- assignTaxonomy(seqtab_filtered, SILVA_DB, multithread = TRUE)
  saveRDS(taxa, checkpoint_taxa)
  message("CHECKPOINT SAVED: ", checkpoint_taxa)
  
  message("Adding species-level assignments...")
  taxa <- addSpecies(taxa, SILVA_SPECIES)
}

# Save final taxonomy
saveRDS(taxa, file.path(OUT_PATH, "taxa_16S.rds"))
message("SAVED: taxa_16S.rds")

# ============================================================================
# Summary
# ============================================================================
message("\n=== 16S Processing Summary ===")
message("Total sequences: ", sum(seqtab_filtered))
message("Mean sequences per sample: ", round(mean(rowSums(seqtab_filtered))))
message("SD sequences per sample: ", round(sd(rowSums(seqtab_filtered))))
message("Total ASVs: ", ncol(seqtab_filtered))
message("Paper reported: 4,015 ASVs")

message("\n=== Checkpoint files saved in ", CHECKPOINT_PATH, " ===")
message("To re-run from scratch, delete the checkpoints directory.")
message("\nDone. Final outputs saved to ", OUT_PATH)
