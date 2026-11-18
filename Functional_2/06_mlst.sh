#!/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh
THREADS=30

# Variables
base_dir='/home/senekowitsch/Thesis/Functional_2/06_mlst'
output_base="${base_dir}/output"
output_mlst="${output_base}/mlst"

# Inputs
assemblies_infantis="/home/senekowitsch/Thesis/Functional_2/00_data"  # *.fna
# Group/lineage lists
LIST_SWEEP7="/home/senekowitsch/Thesis/Functional_2/sweep7_genomes.txt"
# for tree visualization
infantis_tree="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/tree/full_tree_renamed.treefile"

# ------------------------------------------
# Create directories
# ------------------------------------------
mkdir -p "${base_dir}"
mkdir -p "${output_base}"
mkdir -p "${output_mlst}"

cd "${base_dir}"

# ------------------------------------------
# Check database for Salmonella - Achtman scheme
# ------------------------------------------
conda activate pESI
mlst --list | tr ' ' '\n' | grep -i enter

# ------------------------------------------
# Run MLST on all infantis assemblies
# ------------------------------------------
cd "${output_mlst}"
mlst --scheme senterica_achtman_2 --threads ${THREADS} \
     ${assemblies_infantis}/*.fna > mlst_results.tsv

# ------------------------------------------
# Plot tree with MLST results
# ------------------------------------------
conda activate phylo_pipeline

python3 "${base_dir}/plot_mlst_tree.py" \
    --tree   "${infantis_tree}" \
    --mlst   mlst_results.tsv \
    --sweep7 "${LIST_SWEEP7}" \
    --outdir mlst_tree_plot \
    --prefix mlst_on_tree