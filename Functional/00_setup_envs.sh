#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Create a new environment for prokka if not exists
if ! conda env list | grep -q "prokka"; then
    conda create -n prokka bioconda::prokka -y
    conda activate prokka
    prokka --version
else
    echo "prokka environment already exists."
fi

