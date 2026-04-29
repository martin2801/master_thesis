# !/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh

# go to QC folder
cd /home/senekowitsch/Thesis/QC/

# activate conda env fastANI
source $(conda info --base)/etc/profile.d/conda.sh
conda activate fastANI

# move to the right folder
cd /home/senekowitsch/Thesis/QC/06_pairwise_ANI

# make an input file with all the paths to the filtered assemblies
cut -f1 ../05_filter/results/filtered_ANI_results.txt > pairwise_ANI_input_assemblies.txt

# make a test input file with only 10 assemblies
head -n 10 pairwise_ANI_input_assemblies.txt > pairwise_ANI_input_assemblies_test.txt

# run fastANI on the real files
fastANI --ql pairwise_ANI_input_assemblies.txt \ 
    --rl pairwise_ANI_input_assemblies.txt \
    --output results/fastANI_all_vs_all_output.txt \
    --threads 30 \
    --matrix

# set up env for the clustering and heatmap and run it
conda activate ani_heatmap
python3 01_cluster_heatmap.py
conda deactivate

cd /home/senekowitsch/Thesis/QC

