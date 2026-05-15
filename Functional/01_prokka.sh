#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Functional/01_prokka'
output_base="${base_dir}/output"


mkdir -p "${base_dir}"
mkdir -p "${output_base}"
cd "${base_dir}"
conda activate prokka
python3 run_prokka.py

