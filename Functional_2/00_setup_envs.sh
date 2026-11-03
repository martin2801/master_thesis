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

# Create a new environment for pESI if not exists
if ! conda env list | grep -q "pESI"; then
    conda create -n pESI -y
    conda activate pESI
    conda install -c bioconda mob_suite -y
    conda install -c conda-forge -c bioconda abricate -y
    conda install -c conda-forge -c bioconda efetch -y
    mob_recon --version
    abricate --version
    efetch --help
    conda install -c bioconda mlst -y
    mlst --help
else
    echo "pESI environment already exists."
fi

# Create a new environment for pMLST if not exists
if ! conda env list | grep -q "pmlst"; then
    conda create -y -n pmlst
    conda activate pmlst
    conda install -c conda-forge -c bioconda pmlst_ssi -y
    pmlst --help
else
    echo "pmlst environment already exists."
fi

# Create a new environment for hyphy if not exists
if ! conda env list | grep -q "hyphy"; then
    conda create -y -n hyphy
    conda activate hyphy
    conda install -c bioconda hyphy -y
    conda install -c etetoolkit ete3 -y
    conda install -c conda-forge parallel -y
    hyphy --version
else
    echo "hyphy environment already exists."
fi

# Create a new environment for bakta if not exists
if ! conda env list | grep -q "bakta"; then
    conda create -y -n bakta
    conda activate bakta
    conda install -c conda-forge -c bioconda bakta -y
    bakta --version
else
    echo "bakta environment already exists."
fi

# Create a new environment for panaroo if not exists
if ! conda env list | grep -q "panaroo"; then
    conda create -y -n panaroo
    conda activate panaroo
    conda install bioconda::panaroo -y
    conda install -c bioconda dendropy -y
    conda install pandas -y
    panaroo --version
else
    echo "panaroo environment already exists."
fi

# Create a new environment for pyseer if not exists
if ! conda env list | grep -q "pyseer"; then
    conda create -y -n pyseer
    conda activate pyseer
    conda install -c bioconda pyseer -y
    pyseer --version
else
    echo "pyseer environment already exists."
fi

# Remove conda environments if needed
#conda deactivate
#conda env remove -n pmlst

# Create a new environment for minimap2 if not exists
if ! conda env list | grep -q "minimap2"; then
    conda create -y -n minimap2
    conda activate minimap2
    conda install -c bioconda minimap2 -y
    minimap2 --version
else
    echo "minimap2 environment already exists."
fi

# Create a new environment for pESI if not exists
if ! conda env list | grep -q "pESI"; then
    conda create -n pESI -y
    conda activate pESI
    conda install -c bioconda mob_suite -y
    conda install -c conda-forge -c bioconda abricate -y
    conda install -c conda-forge -c bioconda efetch -y
    mob_recon --version
    abricate --version
    efetch --help
    conda install -c bioconda mlst -y
    mlst --help
else
    echo "pESI environment already exists."
fi

# Create a new environment for platon if not exists
if ! conda env list | grep -q "platon"; then
    conda create -y -n platon
    conda activate platon
    conda install -c conda-forge -c bioconda -c defaults platon -y
    platon --version
else
    echo "platon environment already exists."
fi

# Create a new environment for genomad if not exists
if ! conda env list | grep -q "genomad"; then
    conda create -y -n genomad
    conda activate genomad
    conda install -c bioconda -c conda-forge genomad -y
    conda install -c bioconda -c conda-forge blast -y
    genomad --version
else
    echo "genomad environment already exists."
fi

# Create a new environment for edirect if not exists
if ! conda env list | grep -q "edirect"; then
    conda create -y -n edirect
    conda activate edirect
    conda install -c bioconda -c conda-forge entrez-direct -y
    efetch --help
else
    echo "edirect environment already exists."
fi


# Create a new environment for amrfinder if not exists
if ! conda env list | grep -q "amrfinder"; then
    conda create -y -n amrfinder
    conda activate amrfinder
    conda install bioconda::ncbi-amrfinderplus -y
    amrfinder --version
else
    echo "amrfinder environment already exists."
fi

# Create a new environment for snippy if not exists
if ! conda env list | grep -q "snippy_env"; then
    conda create -n snippy_env -c conda-forge -c bioconda -c defaults snippy -y
    conda activate snippy_env
    snippy --check
else
    echo "snippy_env environment already exists."
fi