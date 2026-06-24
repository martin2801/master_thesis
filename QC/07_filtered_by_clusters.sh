# !/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh

# go to QC folder
cd /home/senekowitsch/Thesis/QC/

# move to the right folder
cd /home/senekowitsch/Thesis/QC/07_filtered_by_clusters

# set up env for the filtering, clustering and heatmap and run it
conda activate ani_heatmap
python3 01_define_threshold.py
#python3 02_filter_cluster_heatmap.py
python3 02_filter_cluster_heatmap_ANI99.97.py
python3 02_filter_cluster_heatmap_ANI99.97_pres.py
python3 03_heatmap_filtered.py
conda deactivate

cd /home/senekowitsch/Thesis/QC/08_filtered_data

# print how many genomes are in the filtered set
num_genomes=$(ls -l | wc -l)
echo "There are $((num_genomes - 1)) genomes in the filtered set." # subtract 1 for the total line

cd /home/senekowitsch/Thesis/QC/07_filtered_by_clusters/results

# create a file with the names of the filtered genomes
ls /home/senekowitsch/Thesis/QC/08_filtered_data > filtered_genomes.txt
awk '/Cluster/ {cluster=$2; sub(/:/,"",cluster)} /\.fna/ {print $1 "\t" cluster}' cluster_membership.txt > clusters_table.tsv

conda deactivate


