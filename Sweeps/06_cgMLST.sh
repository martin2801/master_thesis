#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Sweeps/06_cgMLST'
output_base="${base_dir}/output"
output_cgMLST="${output_base}/cgMLST"

GENOMES='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/genomes'
SCHEMA="${base_dir}/salmonella_cgmlst_schema"
THREADS=30

# Create output directory
cd "${base_dir}"
mkdir -p "${output_base}"
mkdir -p "${output_cgMLST}"
mkdir -p "${SCHEMA}"

# Download and prepare the chewBBACA schema
conda activate chewBBACA
echo "Downloading and preparing the chewBBACA schema..."
chewBBACA.py DownloadSchema -sp 14 -sc 1 -o "${SCHEMA}"
# Variable to the schema file


# Run chewBBACA allele calling
echo "Running chewBBACA schema creation..."
echo "on genomes in ${GENOMES}"
chewBBACA.py AlleleCall \
  -i "${GENOMES}" \
  -g "${SCHEMA}/Salmonella_enterica_wgMLST" \
  -o "${output_cgMLST}" \
  --cpu "${THREADS}"

# Check what the output folder is named and set it as a variable
output_allelecall=$(ls -d ${output_cgMLST}/results*/ | head -n 1)
echo "Allele calling results are in: ${output_allelecall}"

# Run Evaluator to check the quality of the allele calling results
chewBBACA.py AlleleCallEvaluator \
  -i "${output_allelecall}" \
  -g "${SCHEMA}/Salmonella_enterica_wgMLST" \
  -o "${output_cgMLST}/EvaluatorResults" \
  --cpu "${THREADS}"

# Extract the cgMLST matrix from the allele calling results
chewBBACA.py ExtractCgMLST \
  -i "${output_allelecall}/results_alleles.tsv" \
  -o "${output_cgMLST}/cgmlst_matrix/" \
  --t 0.95    # this threshold is what converts wgMLST → effective cgMLST

conda activate phylo_pipeline
python cgmlst_sweep_comparison.py \
    --profiles ${output_cgMLST}/EvaluatorResults/cgMLST_profiles.tsv \
    --sweeps /home/senekowitsch/Thesis/Sweeps/05_check_distance/output/sweeps_bottomup_clonal_5x.txt \
    --sweep_tree /home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/validate_sweeps/full_tree.treefile \
    --outdir ${output_cgMLST}/cgmlst_sweep_results/

