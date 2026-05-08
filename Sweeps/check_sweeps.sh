#!/bin/bash

# Variables
base_dir='/home/senekowitsch/Thesis/Sweeps'
output_base='/home/senekowitsch/Thesis/Sweeps_output'
output_distances="${output_base}/distances"

treefile="/home/senekowitsch/Time/10_place_on_tree/output/validate_sweeps/full_tree.treefile"

SNP_LIST="${output_distances}/snp_list.tab"
SCRIPT="${base_dir}/check_5x_rule.py"

threads=30