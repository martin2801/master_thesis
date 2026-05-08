#!/bin/bash

# Variables
base_dir='/home/senekowitsch/Time/11_check_distance'
output_base='/home/senekowitsch/Time/11_check_distance/output'
output_distances="${output_base}/distances"
unambiguous_core_alignment="/home/senekowitsch/Time/10_place_on_tree/output/snippy_alignment/core.full.aln"
treefile="/home/senekowitsch/Time/10_place_on_tree/output/validate_sweeps/full_tree.treefile"

SNP_LIST="${output_distances}/snp_list.tab"
SCRIPT="${base_dir}/check_5x_rule.py"

threads=30

# Define sweeps
sweep1_genomes=(
    "out_11452375"
    "out_11452785"
    "out_11453035"
    "out_8046645"
    "out_8091485"
    "out_11452255"
    "out_11452675"
)
sweep1_sister=("out_8790705")

sweep2_genomes=(
    "out_14966075"
    "out_8963375"
    "out_31832435"
    "out_8067535"
    "out_8476865"
    "out_15254495"
    "out_32036575"
    "out_32037335"
    "out_8070535"
    "out_8514935"
    "out_8777535"
    "out_9027665"
    "out_9523435"
)
sweep2_sister=("out_11452515")

sweep3_genomes=(
    "out_20671165"
    "out_20671355"
    "out_20671465"
    "out_20671385"
)
sweep3_sister=("out_19267735")

sweep4_genomes=(
    "out_45809245"
    "out_45808905"
    "out_45808925"
    "out_45808965"
    "out_45808885"
)
sweep4_sister=("out_9523415")

# Create output directories if they don't exist
cd "${base_dir}"
mkdir -p "${output_sweep_alignments}"
mkdir -p "${output_distances}"

# Calculate pairwise SNP distances using snp-dists
conda activate phylo_pipeline
snp-dists -j "${threads}" -m "${unambiguous_core_alignment}" > "${SNP_LIST}"

# Check the 5x rule for each sweep
python3 $SCRIPT --name "Sweep 1" --sister "$sweep1_sister" --snps "$SNP_LIST" --genomes "${sweep1_genomes[@]}"
python3 $SCRIPT --name "Sweep 2" --sister "$sweep2_sister" --snps "$SNP_LIST" --genomes "${sweep2_genomes[@]}"
python3 $SCRIPT --name "Sweep 3" --sister "$sweep3_sister" --snps "$SNP_LIST" --genomes "${sweep3_genomes[@]}"
python3 $SCRIPT --name "Sweep 4" --sister "$sweep4_sister" --snps "$SNP_LIST" --genomes "${sweep4_genomes[@]}"

# Detect sweeps
python3 detect_sweeps.py --tree "${treefile}" --snps "${SNP_LIST}" --threshold 5 --min-tips 3 --out "${output_base}/sweeps.txt"

python3 detect_sweeps_bottomup.py --tree "${treefile}" --snps "${SNP_LIST}" --threshold 2 --min-tips 3 --out "${output_base}/sweeps_bottomup_2x.txt"
python3 detect_sweeps_bottomup.py --tree "${treefile}" --snps "${SNP_LIST}" --threshold 3 --min-tips 3 --out "${output_base}/sweeps_bottomup_3x.txt"
python3 detect_sweeps_bottomup.py --tree "${treefile}" --snps "${SNP_LIST}" --threshold 4 --min-tips 3 --out "${output_base}/sweeps_bottomup_4x.txt"
python3 detect_sweeps_bottomup.py --tree "${treefile}" --snps "${SNP_LIST}" --threshold 5 --min-tips 3 --out "${output_base}/sweeps_bottomup_5x.txt"

python3 detect_sweeps_bottomup_clonal.py --tree "${treefile}" --snps "${SNP_LIST}" --threshold 5 --min-tips 3 --out "${output_base}/sweeps_bottomup_clonal_5x.txt"

# Create image of tree with sweeps highlighted
python3 plot_tree.py --tree "${treefile}" --sweeps "${output_base}/sweeps.txt" --out "${output_base}/tree_sweeps.png" --labels

python3 plot_tree.py --tree "${treefile}" --sweeps "${output_base}/sweeps_bottomup_2x.txt" --out "${output_base}/tree_sweeps_bottomup_2x.png" --labels
python3 plot_tree.py --tree "${treefile}" --sweeps "${output_base}/sweeps_bottomup_3x.txt" --out "${output_base}/tree_sweeps_bottomup_3x.png" --labels
python3 plot_tree.py --tree "${treefile}" --sweeps "${output_base}/sweeps_bottomup_4x.txt" --out "${output_base}/tree_sweeps_bottomup_4x.png" --labels
python3 plot_tree.py --tree "${treefile}" --sweeps "${output_base}/sweeps_bottomup_5x.txt" --out "${output_base}/tree_sweeps_bottomup_5x.png" --labels

python3 plot_tree.py --tree "${treefile}" --sweeps "${output_base}/sweeps_bottomup_clonal_5x.txt" --out "${output_base}/tree_sweeps_bottomup_clonal_5x.png" --labels