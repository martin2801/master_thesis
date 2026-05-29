#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Functional/04_clean_BLAST'
output_base="${base_dir}/output"
annotations="${base_dir}/annotations"

# Create output directories
mkdir -p "${base_dir}"
mkdir -p "${output_base}"
mkdir -p "${annotations}"
cd "${base_dir}"

# Copy annotations from previous step
cp /home/senekowitsch/Thesis/Functional/02_eggnog/output/*/*/*.annotations "${annotations}/"

cd "${base_dir}"

# Run R script to clean BLAST results
Rscript /home/senekowitsch/Thesis/Functional/04_prepare_BLAST_data.r