# DADA2 Pipeline for ITS (Fungi)
# Based on Boutin et al. (2024) doi:10.1139/cjm-2023-0215
# Primers: ITS-1F / ITS2
#
# MATCHES AUTHORS' SCRIPT: Uses cutadapt for primer removal (ITS has variable length)
# Source: https://doi.org/10.6084/m9.figshare.24570052 (ITS_2022_dada2_18avril.Rmd)
# CHECKPOINT SAVES: Intermediate results saved after each expensive step

library(dada2)
library(tidyverse)

# Configuration
RAW_PATH <- "data/raw/ITS"
OUT_PATH <- "data/processed/ITS"
CHECKPOINT_PATH <- file.path(OUT_PATH, "checkpoints")
CUTADAPT_PATH <- file.path(OUT_PATH, "cutadapt")
UNITE_DB <- "data/reference/sh_general_release_dynamic_19.02.2025.fasta"

MIN_SEQ_THRESHOLD <- 10  # Minimum sequences per ASV across dataset

# Primer sequences (from paper Table S2)
# ITS-1F: Gardes & Bruns 1993
# ITS2: White et al. 1990
PRIMER_F <- "CTTGGTCATTTAGAGGAAGTAA"  # ITS-1F (22 bp)
PRIMER_R <- "GCTGCGTTCTTCATCGATGC"    # ITS2 (20 bp)

# Helper function for tracking
getN <- function(x) sum(getUniques(x))

# Fail fast if UNITE database is missing
if (!file.exists(UNITE_DB)) {
  stop(
    "UNITE database not found at: ", UNITE_DB, "\n\n",
    "Download from: https://doi.plutof.ut.ee/doi/10.15156/BIO/2959336\n",
    "  OR: https://unite.ut.ee/repository.php (select UNITE+INSD, DADA2 format)\n\n",
    "After download:\n",
    "  1. Extract the archive\n",
    "  2. Rename the .fasta file to match: ", basename(UNITE_DB), "\n",
    "  3. Place in data/reference/"
  )
}

# Check if cutadapt is available
Sys.setenv(PATH = paste0(Sys.getenv("HOME"), "/.local/bin:", Sys.getenv("PATH")))
cutadapt_path <- Sys.which("cutadapt")
if (cutadapt_path == "") {
  stop("cutadapt not found. Install with: pipx install cutadapt")
}
message("Using cutadapt: ", cutadapt_path)

dir.create(OUT_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(CUTADAPT_PATH, recursive = TRUE, showWarnings = FALSE)
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

# ============================================================================
# STEP 0: Remove primers with cutadapt (authors' method for ITS)
# ============================================================================
# ITS requires cutadapt (not trimLeft) because ITS region has variable length.
# Authors used: -g FWD -a REV_RC -G REV -A FWD_RC -n 2 -m 20
# NO --discard-untrimmed (keeps all reads)
checkpoint_cutadapt <- file.path(CHECKPOINT_PATH, "00_cutadapt.rds")

# Cutadapt output paths
cutFs <- file.path(CUTADAPT_PATH, basename(fnFs))
cutRs <- file.path(CUTADAPT_PATH, basename(fnRs))

if (file.exists(checkpoint_cutadapt)) {
  message("CHECKPOINT: Cutadapt already completed")
} else {
  message("STEP 0: Removing primers with cutadapt...")
  
  # Reverse complements (using dada2 internal function)
  FWD_RC <- dada2:::rc(PRIMER_F)
  REV_RC <- dada2:::rc(PRIMER_R)
  
  message("  FWD: ", PRIMER_F)
  message("  REV: ", PRIMER_R)
  message("  FWD_RC: ", FWD_RC)
  message("  REV_RC: ", REV_RC)
  
  # R1 flags: trim FWD from 5' and REV_RC from 3'
  R1_flags <- paste("-g", PRIMER_F, "-a", REV_RC)
  # R2 flags: trim REV from 5' and FWD_RC from 3'
  R2_flags <- paste("-G", PRIMER_R, "-A", FWD_RC)
  
  for (i in seq_along(fnFs)) {
    if (file.exists(cutFs[i]) && file.size(cutFs[i]) > 0) {
      message("  ", sample_names[i], " - already processed, skipping")
      next
    }
    
    message("  Processing ", sample_names[i], " (", i, "/", length(fnFs), ")...")
    system2(cutadapt_path, args = c(
      R1_flags, R2_flags,
      "-m", 20,   # minimum length 20bp
      "-n", 2,    # remove up to 2 adapters per read (needed for FWD+REV)
      "-o", cutFs[i], "-p", cutRs[i],
      fnFs[i], fnRs[i]
    ), stdout = "", stderr = "")
  }
  
  saveRDS(TRUE, checkpoint_cutadapt)
  message("CHECKPOINT SAVED: ", checkpoint_cutadapt)
}

# Verify cutadapt output
cutFs <- sort(list.files(CUTADAPT_PATH, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
cutRs <- sort(list.files(CUTADAPT_PATH, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
sample_names <- sub("_1\\.fastq\\.gz$", "", basename(cutFs))
message("Cutadapt output files: ", length(cutFs))

if (length(cutFs) == 0) {
  stop("No cutadapt output files found. Check cutadapt installation.")
}

# Quality filtered output paths (from cutadapt output)
filtFs <- file.path(OUT_PATH, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(OUT_PATH, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# ============================================================================
# STEP 1: Filter and trim (no truncLen for ITS — variable length region)
# ============================================================================
# Authors used: maxN=0, maxEE=c(2,5), truncQ=2, minLen=50, rm.phix=TRUE
checkpoint_filter <- file.path(CHECKPOINT_PATH, "01_filter_output.rds")

if (file.exists(checkpoint_filter)) {
  message("CHECKPOINT: Loading filter results from ", checkpoint_filter)
  out <- readRDS(checkpoint_filter)
} else {
  message("STEP 1: Filtering reads (from cutadapt output)...")
  out <- filterAndTrim(
    cutFs, filtFs, cutRs, filtRs,
    maxN = 0,
    maxEE = c(2, 5),   # Authors' values: lenient on reverse reads
    truncQ = 2,
    minLen = 50,
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

# Only use filtered files that exist
filtFs_exist <- filtFs[file.exists(filtFs)]
filtRs_exist <- filtRs[file.exists(filtRs)]

if (file.exists(checkpoint_errF) && file.exists(checkpoint_errR)) {
  message("CHECKPOINT: Loading error models...")
  errF <- readRDS(checkpoint_errF)
  errR <- readRDS(checkpoint_errR)
} else {
  message("STEP 2: Learning error rates...")
  
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
# STEP 3+4: Dereplicate + Sample inference (combined to save memory)
# ============================================================================
# NOTE: dada() handles dereplication internally when given file paths.
# This avoids loading all samples into memory at once (OOM on 16GB).
filtFs_exist <- filtFs[file.exists(filtFs)]
filtRs_exist <- filtRs[file.exists(filtRs)]
sample_names_exist <- names(filtFs_exist)

checkpoint_dada <- file.path(CHECKPOINT_PATH, "04_dada.rds")

if (file.exists(checkpoint_dada)) {
  message("CHECKPOINT: Loading DADA2 inference results...")
  dada_data <- readRDS(checkpoint_dada)
  dadaFs <- dada_data$dadaFs
  dadaRs <- dada_data$dadaRs
} else {
  message("STEP 3+4: Dereplicating and running DADA2 sample inference...")
  dadaFs <- dada(filtFs_exist, err = errF, multithread = TRUE)
  dadaRs <- dada(filtRs_exist, err = errR, multithread = TRUE)
  names(dadaFs) <- sample_names_exist
  names(dadaRs) <- sample_names_exist
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
  mergers <- mergePairs(dadaFs, filtFs_exist, dadaRs, filtRs_exist, verbose = TRUE)
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
saveRDS(seqtab_filtered, file.path(OUT_PATH, "seqtab_ITS.rds"))
message("SAVED: seqtab_ITS.rds")

# ============================================================================
# STEP 9: Track reads through pipeline
# ============================================================================
message("STEP 9: Creating read tracking table...")

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
# STEP 10: Assign taxonomy with UNITE
# ============================================================================
checkpoint_taxa <- file.path(CHECKPOINT_PATH, "10_taxa.rds")

if (file.exists(checkpoint_taxa)) {
  message("CHECKPOINT: Loading taxonomy...")
  taxa <- readRDS(checkpoint_taxa)
} else {
  message("STEP 10: Assigning taxonomy with UNITE...")
  taxa <- assignTaxonomy(seqtab_filtered, UNITE_DB, multithread = TRUE, tryRC = TRUE)
  saveRDS(taxa, checkpoint_taxa)
  message("CHECKPOINT SAVED: ", checkpoint_taxa)
}

# Save final taxonomy
saveRDS(taxa, file.path(OUT_PATH, "taxa_ITS.rds"))
message("SAVED: taxa_ITS.rds")

# ============================================================================
# Summary
# ============================================================================
message("\n=== ITS Processing Summary ===")
message("Total sequences: ", sum(seqtab_filtered))
message("Mean sequences per sample: ", round(mean(rowSums(seqtab_filtered))))
message("SD sequences per sample: ", round(sd(rowSums(seqtab_filtered))))
message("Total ASVs: ", ncol(seqtab_filtered))
message("Paper reported: 1,027 ASVs")

message("\n=== Checkpoint files saved in ", CHECKPOINT_PATH, " ===")
message("To re-run from scratch, delete the checkpoints directory.")
message("\nDone. Final outputs saved to ", OUT_PATH)
