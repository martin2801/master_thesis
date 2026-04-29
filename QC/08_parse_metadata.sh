# !/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh

# go to QC folder
cd /home/senekowitsch/Thesis/QC/

# move to the right folder
cd /home/senekowitsch/Thesis/QC/09_check_meta

# create a list of all files in 00_data and save as txt file
cd /home/senekowitsch/Thesis/QC/00_data
ls *.fna > genome_list.txt
cd /home/senekowitsch/Thesis/QC/09_check_meta

# set up env for the filtering, clustering and heatmap and run it
conda activate ani_heatmap
python3 01_parse_metadata.py
conda deactivate