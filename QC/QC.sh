# !/bin/bash

# go to QC folder
cd /home/senekowitsch/Thesis/QC

# create folder structure if not exists
if [ ! -d "00_data" ]; then
    mkdir -p 00_data
    else
    echo "Directory 00_data already exists."
fi

if [ ! -d "01_ANI/results" ]; then
    mkdir -p 01_ANI/results
    else
    echo "Directory 01_ANI/results already exists."
fi

if [ ! -d "02_seqsero/results" ]; then
    mkdir -p 02_seqsero/results
    else
    echo "Directory 02_seqsero/results already exists."
fi

if [ ! -d "03_checkM/results" ]; then
    mkdir -p 03_checkM/results
    else
    echo "Directory 03_checkM/results already exists."
fi

if [ ! -d "04_QUAST/results" ]; then
    mkdir -p 04_QUAST/results
    else
    echo "Directory 04_QUAST/results already exists."
fi

if [ ! -d "05_filter/results" ]; then
    mkdir -p 05_filter/results
    else
    echo "Directory 05_filter/results already exists."
fi

if [ ! -d "06_pairwise_ANI/results" ]; then
    mkdir -p 06_pairwise_ANI/results
    else
    echo "Directory 06_pairwise_ANI/results already exists."
fi

if [ ! -d "07_filtered_by_clusters/results" ]; then
    mkdir -p 07_filtered_by_clusters/results
    else
    echo "Directory 07_filtered_by_clusters/results already exists."
fi

if [ ! -d "08_filtered_data" ]; then
    mkdir -p 08_filtered_data
    else
    echo "Directory 08_filtered_data already exists."
fi

if [ ! -d "09_check_meta/results" ]; then
    mkdir -p 09_check_meta/results
    else
    echo "Directory 09_check_meta/results already exists."
fi

# make a link to the data in foler 00_data
if [ ! -d "00_data/*.fna" ]; then
    ln -s /home/senekowitsch/raw_data/*.fna 00_data
    else
    echo "Symbolic link 00_data/data already exists."
fi

# Run QC
bash setup_envs.sh
bash 01_ANI.sh
bash 02_seqsero.sh
bash 03_checkM.sh
bash 04_QUAST.sh
bash 05_filter.sh
bash 06_pairwise_ANI.sh
bash 07_filtered_by_clusters.sh
bash 08_parse_metadata.sh
