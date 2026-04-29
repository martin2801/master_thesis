# !/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh

# go to QC folder
cd /home/senekowitsch/Thesis/QC/

# activate conda env checkm2
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ani_analysis

# move to the right folder
cd /home/senekowitsch/Thesis/QC/05_filter

# make folder for filtered results if it does not exist
if [ ! -d "filtered_data" ]; then
    mkdir -p filtered_data
    else
    echo "Directory filtered_data already exists."
fi

# remove everything in that folder
rm -r /home/senekowitsch/Thesis/QC/05_filter/filtered_data/*

# run the filtering script
conda activate ani_analysis
python3 /home/senekowitsch/Thesis/QC/05_filter/01_filter.py


