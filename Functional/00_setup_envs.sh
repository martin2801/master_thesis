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

# Create a new environment for eggnog-mapper if not exists
if ! conda env list | grep -q "eggnog"; then
    conda create -n eggnog -y
    conda activate eggnog
    conda install -c bioconda -c conda-forge eggnog-mapper -y
    emapper.py --version
else
    echo "eggnog environment already exists."
fi

# Create a new environment for blast if not exists
if ! conda env list | grep -q "blast"; then
    conda create -n blast -y
    conda activate blast
    conda install -c bioconda blast -y
    blastp -version
else
    echo "blast environment already exists."
fi