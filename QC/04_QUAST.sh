# !/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh

# go to QC folder
cd /home/senekowitsch/Thesis/QC/

# activate conda env quast
source $(conda info --base)/etc/profile.d/conda.sh
conda activate quast

# move to the right folder
cd /home/senekowitsch/Thesis/QC/04_QUAST

# set path to the input folder
INPUT="/home/senekowitsch/Thesis/QC/03_checkM/inputdata"
TEST="/home/senekowitsch/Thesis/QC/04_QUAST/test_input"

# set path to the output folder
OUTPUT="/home/senekowitsch/Thesis/QC/04_QUAST/results"

# run quast
rm -r /home/senekowitsch/Thesis/QC/04_QUAST/results/*
quast.py -o $OUTPUT -t 10 $INPUT/*.fna
