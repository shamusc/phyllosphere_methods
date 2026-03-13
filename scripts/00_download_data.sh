#!/bin/bash
# Download raw data from NCBI SRA and reference databases
# Based on Boutin et al. (2024) doi:10.1139/cjm-2023-0215

set -e

# Create directories
mkdir -p data/raw/16S data/raw/ITS data/reference data/metadata

echo "=== Downloading raw sequences from NCBI SRA ==="
echo "BioProject: PRJNA1081804"

# Install SRA Toolkit if not available
if ! command -v fasterq-dump &> /dev/null; then
    echo "SRA Toolkit not found. Please install it first:"
    echo "  conda install -c bioconda sra-tools"
    echo "  OR"
    echo "  brew install sratoolkit"
    exit 1
fi

# Download SRA accession list
# You may need to get the full list from NCBI
echo "Fetching SRA run info..."
esearch -db sra -query "PRJNA1081804" | efetch -format runinfo > data/raw/SraRunInfo.csv

# Download FASTQ files (this will take time)
echo "Downloading FASTQ files..."
cut -d',' -f1 data/raw/SraRunInfo.csv | tail -n +2 | while read acc; do
    echo "Downloading $acc..."
    fasterq-dump --split-files --outdir data/raw "$acc"
    gzip data/raw/${acc}*.fastq
done

# Separate 16S and ITS files based on library strategy (manual step may be needed)
echo "NOTE: You may need to manually separate 16S and ITS files based on metadata"

echo ""
echo "=== Downloading reference databases ==="

# SILVA database for 16S
echo "Downloading SILVA v138.1..."
wget -P data/reference/ https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
wget -P data/reference/ https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz

# UNITE database for ITS (manual download required)
echo ""
echo "=== UNITE Database (manual download required) ==="
echo "UNITE requires accepting license terms, so cannot be auto-downloaded."
echo ""
echo "Option 1 - Current release via DOI:"
echo "  https://doi.plutof.ut.ee/doi/10.15156/BIO/2959336"
echo "  Look for: sh_general_release_dynamic_*.fasta (DADA2 format)"
echo ""
echo "Option 2 - Via UNITE website:"
echo "  https://unite.ut.ee/repository.php"
echo "  Download: UNITE+INSD release, FASTA format for DADA2"
echo ""
echo "After download, extract and place .fasta file in: data/reference/"
echo ""
echo "NOTE: You can skip ITS/fungi and run 16S/bacteria analysis first."

echo ""
echo "=== Downloading original scripts from Figshare ==="

# 16S scripts
wget -O data/metadata/16S_scripts.zip "https://figshare.com/ndownloader/files/$(curl -s https://api.figshare.com/v2/articles/24570064 | jq -r '.files[0].id')" || echo "Manual download needed: https://doi.org/10.6084/m9.figshare.24570064"

# ITS scripts  
wget -O data/metadata/ITS_scripts.zip "https://figshare.com/ndownloader/files/$(curl -s https://api.figshare.com/v2/articles/24570055 | jq -r '.files[0].id')" || echo "Manual download needed: https://doi.org/10.6084/m9.figshare.24570055"

echo ""
echo "=== Download complete ==="
echo "Next steps:"
echo "1. Verify FASTQ files in data/raw/"
echo "2. Download UNITE database manually if needed"
echo "3. Create sample_metadata.csv based on SraRunInfo.csv"
echo "4. Run R scripts in order: 01_dada2_16S.R, 02_dada2_ITS.R, etc."
