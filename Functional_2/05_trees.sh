#!/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh
THREADS=60

# Variables
base_dir='/home/senekowitsch/Thesis/Functional_2/05_trees'
output_base="${base_dir}/output"

# Outputs
output_trees="${output_base}/trees"
output_snippy="${output_base}/snippy"

# Inputs
core_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/snippy_alignment/core.full.aln"
recomb_removed_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/core_alignment_noambiguous.aln"
pesi_fastas="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/consensus/pESI"
sweeps_infantis="/home/senekowitsch/Thesis/Functional_2/sweep7_genomes.txt"

# ------------------------------------------
# Create directories
# ------------------------------------------
mkdir -p "${base_dir}"
mkdir -p "${output_base}"
mkdir -p "${output_trees}"
mkdir -p "${output_snippy}"

cd "${base_dir}"

# ------------------------------------------
# Download reference for pESI plasmid & create abricate database
# ------------------------------------------
cd "${base_dir}"
conda activate edirect
efetch -db nucleotide \
  -id NZ_CP016411.1 \
  -format fasta_cds_na \
  > NZ_CP016411.1_cds.fasta

sed 's/>lcl|/>/' NZ_CP016411.1_cds.fasta > NZ_CP016411.1_cds_clean.fa

# ------------------------------------------
# Run snippy on pESI plasmids to create a core alignment
# ------------------------------------------
echo ""
echo "--- Running snippy (per-sample) ---"
conda activate snippy_env

cd "${output_snippy}"
for genome in "${pesi_fastas}"/*.fasta; do
    [[ -e "$genome" ]] || { echo "ERROR: No .fasta files found in ${FOLDER}" >&2; exit 1; }
    sample_name=$(basename "$genome" .fasta)
    echo "  Processing: $sample_name"
    snippy \
        --cpus  "$THREADS" \
        --outdir "${output_snippy}/${sample_name}" \
        --ref    /home/senekowitsch/Thesis/Functional_2/05_trees/NZ_CP016411.1_cds_clean.fa \
        --ctgs   "$genome"
done

# --- snippy-core: merge all per-sample VCFs into a multi-sample alignment ---
# We explicitly cd to output_snippy so that core.* files land there.
# core.full.aln (all reference positions, not just SNPs) is required by Gubbins.
echo ""
echo "--- Running snippy-core ---"
cd "${output_snippy}"
snippy-core --prefix core_pesi --ref "/home/senekowitsch/Thesis/Functional_2/05_trees/NZ_CP016411.1_cds_clean.fa" "${output_snippy}"/*

echo "Snippy complete. Core alignment: ${output_snippy}/core_pesi.full.aln"

# ------------------------------------------
# Build trees from alignments using IQ-TREE
# ------------------------------------------
# Recombination removed tree (core_alignment_noambiguous.aln)
cd "${output_trees}"
echo ""
echo "-- Building phylogenetic tree with IQ-TREE ---"
conda activate phylo_pipeline
iqtree \
    -s "${recomb_removed_alignment}" \
    --prefix "recomb_removed_tree" \
    -m GTR+G+I \
    -T "${THREADS}" \
    -B 1000 \
    --redo

# Recombination raw tree (core.full.aln)
cd "${output_trees}"
echo ""
echo "-- Building phylogenetic tree with IQ-TREE ---"
conda activate phylo_pipeline
iqtree \
    -s "${core_alignment}" \
    --prefix "raw_tree" \
    -m GTR+G+I \
    -T "${THREADS}" \
    -B 1000 \
    --redo

# pESI tree (core_pesi.full.aln)
cd "${output_trees}"
echo ""
echo "-- Building phylogenetic tree with IQ-TREE ---"
conda activate phylo_pipeline
iqtree \
    -s "${output_snippy}/core_pesi.full.aln" \
    --prefix "pESI_tree" \
    -m GTR+G+I \
    -T "${THREADS}" \
    -B 1000 \
    --redo

# ------------------------------------------
# Plot trees
# ------------------------------------------
cd "${base_dir}"
# Raw tree (core.full.aln)
python3 plot_tree.py \
    --tree "${output_trees}/raw_tree.treefile" \
    --highlight "${sweeps_infantis}" "#1B7837" "Sweep 7" \
    --out "${output_base}/raw_tree.png"

# Recombination removed tree (core_alignment_noambiguous.aln)
python3 plot_tree.py \
    --tree "${output_trees}/recomb_removed_tree.treefile" \
    --highlight "${sweeps_infantis}" "#1B7837" "Sweep 7" \
    --out "${output_base}/recomb_removed_tree.png"

# pESI tree (core_pesi.full.aln)
python3 plot_tree.py \
    --tree "${output_trees}/pESI_tree.treefile" \
    --highlight "${sweeps_infantis}" "#1B7837" "Sweep 7" \
    --out "${output_base}/pESI_tree.png"

# Tanglegram
python3 tanglegram.py \
    --left  "${output_trees}/raw_tree.treefile" "Raw" \
    --mid   "${output_trees}/recomb_removed_tree.treefile" "Recombination Removed" \
    --right "${output_trees}/pESI_tree.treefile" "pESI" \
    --highlight "${sweeps_infantis}" "#1B7837" "Sweep 7" \
    --labels \
    --out tanglegram.png