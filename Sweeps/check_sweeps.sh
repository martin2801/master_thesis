#!/bin/bash

# Variables
base_dir='/home/senekowitsch/Thesis/Sweeps'
output_base='/home/senekowitsch/Thesis/Sweeps_output'
output_distances="${output_base}/distances"

treefile="/home/senekowitsch/Time/10_place_on_tree/output/validate_sweeps/full_tree.treefile"

SNP_LIST="${output_distances}/snp_list.tab"
SCRIPT="${base_dir}/check_5x_rule.py"

threads=30


raw_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/snippy_alignment/core.full.aln"
cleaned_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/core_alignment_noambiguous.aln"

raw_tree="/home/senekowitsch/Thesis/Sweeps_output/raw_tree.treefile"
clean_tree="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/validate_sweeps/full_tree.treefile"

sweeps="/home/senekowitsch/Thesis/Sweeps/05_check_distance/output/sweeps_bottomup_clonal_5x.txt"
draw_tree="/home/senekowitsch/Thesis/Sweeps/05_check_distance/plot_tree_with7.py"

conda activate phylo_pipeline
cd "${base_dir}"
mkdir -p "${output_base}"
cd "${output_base}"


iqtree \
    -s "${raw_alignment}" \
    --prefix raw_tree \
    -m GTR+G+I \
    -T "${threads}" \
    -B 1000 \
    --redo

conda activate phylo_pipeline
python3 "${draw_tree}" \
    --tree "${raw_tree}" \
    --sweeps "${sweeps}" \
    --out "${output_base}/raw_tree_sweeps_bottomup_clonal_5x_pres_with7.pdf" \
    --labels \
    --dpi 400

python3 "${draw_tree}" \
    --tree "${clean_tree}" \
    --sweeps "${sweeps}" \
    --out "${output_base}/clean_tree_sweeps_bottomup_clonal_5x_pres_with7.pdf" \
    --labels \
    --dpi 400
