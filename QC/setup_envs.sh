# !/bin/bash

# go to QC folder
cd /home/senekowitsch/Thesis/QC

# setup conda envs for the different steps of the QC pipeline

# ncbi-datasets
if conda info --envs | grep -q 'ncbi-datasets'; then
    echo "Conda environment 'ncbi-datasets' already exists."
else
    echo "Creating conda environment 'ncbi-datasets'."
    conda create -n ncbi-datasets -c conda-forge ncbi-datasets-cli
fi

# fastANI
if conda info --envs | grep -q 'fastANI'; then
    echo "Conda environment 'fastANI' already exists."
else
    echo "Creating conda environment 'fastANI'."
    conda create -n fastANI -c conda-forge -c bioconda fastani=1.34 -y
    conda install -c conda-forge -c bioconda biopython -y
fi

# ani_analysis
if conda info --envs | grep -q 'ani_analysis'; then
    echo "Conda environment 'ani_analysis' already exists."
else
    echo "Creating conda environment 'ani_analysis'."
    conda create -n ani_analysis -c conda-forge pandas scipy matplotlib -y
fi


# seqsero2
if conda info --envs | grep -q 'seqsero'; then
    echo "Conda environment 'seqsero' already exists."
else
    echo "Creating conda environment 'seqsero'."
    conda create -n seqsero -c conda-forge -c bioconda seqsero2=1.3.2 -y
fi

# checkM
if conda info --envs | grep -q 'checkm2'; then
    echo "Conda environment 'checkm2' already exists."
else
    echo "Creating conda environment 'checkm2'."
    conda create -n checkm2 -c bioconda -c conda-forge checkm2 -y
fi

# QUAST
if conda info --envs | grep -q 'quast'; then
    echo "Conda environment 'quast' already exists."
else
    echo "Creating conda environment 'quast'."
    conda create -n quast -c bioconda quast -y
fi

# ani_heatmap
if conda info --envs | grep -q 'ani_heatmap'; then
    echo "Conda environment 'ani_heatmap' already exists."
else
    echo "Creating conda environment 'ani_heatmap'."
    conda create -n ani_heatmap -c conda-forge pandas scipy matplotlib seaborn -y
    conda install networkx -n ani_heatmap -c conda-forge -y
    conda install -c anaconda scikit-learn -y
fi  
