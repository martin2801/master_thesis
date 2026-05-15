#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Create a new environment for snippy if not exists
if ! conda env list | grep -q "snippy_env"; then
    conda create -n snippy_env -c conda-forge -c bioconda -c defaults snippy -y
    conda activate snippy_env
    snippy --check
else
    echo "snippy_env environment already exists."
fi

# Create a new environment for gubbins if not exists
conda create -n gubbins
conda activate gubbins
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install gubbins
echo "Created gubbins_env environment and installed gubbins."

# Create a new environment for phylo_time if not exists
if ! conda env list | grep -q "phylo_pipeline"; then
    conda create -n phylo_pipeline python=3.10 -y
    conda activate phylo_pipeline
    conda install -c bioconda -c conda-forge iqtree treetime -y
    conda install -c conda-forge -c bioconda phylip -y
    conda install -c conda-forge -c bioconda seqmagick -y
    conda install -c conda-forge -c bioconda trimal -y
    conda install -c conda-forge -c bioconda snp-sites -y
    conda install -c conda-forge -c bioconda snp-dists -y
    conda install -c conda-forge -c bioconda bedtools -y
    conda install -c conda-forge -c bioconda scikit-learn -y
    # Verify installs
    iqtree --version
    treetime --version
    trimal --version
else
    echo "phylo_time environment already exists."
fi

# create new environment for parsnp if not exists
if ! conda env list | grep -q "parsnp"; then
    conda create -n parsnp python=3.10 -y
    conda activate parsnp
    conda install -c bioconda -c conda-forge parsnp -y
else
    echo "parsnp environment already exists."
fi

# create new environment for vcftools if not exists
if ! conda env list | grep -q "vcftools"; then
    conda create -n vcftools -c conda-forge -c bioconda vcftools -y
    conda activate vcftools
else
    echo "vcftools environment already exists."
fi

# create new environment for clonalframeML if not exists
if ! conda env list | grep -q "clonalframeML"; then
    conda create -n clonalframeML -c conda-forge -c bioconda clonalframeml -y
    conda activate clonalframeML
    conda install -c conda-forge -c bioconda phyx -y
else
    echo "clonalframeML environment already exists."
fi

# create new environment for bedtools if not exists
if ! conda env list | grep -q "bedtools"; then
    conda create -n bedtools -c conda-forge -c bioconda bedtools -y
    conda activate bedtools
else
    echo "bedtools environment already exists."
fi

# create new environment for check_sweeps if not exists
if ! conda env list | grep -q "check_sweeps"; then
    conda create -n check_sweeps -y
    conda activate check_sweeps
    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install bioconda::samtools -y
    conda install -c bioconda emboss -y
else
    echo "check_sweeps environment already exists."
fi

# create new environment for epa-ng if not exists
if ! conda env list | grep -q "epa-ng"; then
    conda create -n epa-ng -c conda-forge -c bioconda epa-ng -y
    conda activate epa-ng
    conda install -c conda-forge -c bioconda gappa -y
else
    echo "epa-ng environment already exists."
fi

# create new environment for chewBBACA if not exists
if ! conda env list | grep -q "chewBBACA"; then
    conda create -n chewBBACA -c conda-forge -c bioconda chewbbaca -y
    conda activate chewBBACA
else
    echo "chewBBACA environment already exists."
fi

# move to Time folder
cd /home/senekowitsch/Thesis/Sweeps

# create folder structure
mkdir -p 01_recombination_checks
mkdir -p 02_remove_recombination
mkdir -p 03_check_sweeps
mkdir -p 04_place_on_tree
mkdir -p 05_check_distance
