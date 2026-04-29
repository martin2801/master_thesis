# !/bin/bash

# go to PopCoGenomeS folder
cd /home/senekowitsch/Thesis/PopCoGenomeS

# create folder structure if not exists
if [ ! -d "00_data_alt" ]; then
    mkdir -p 00_data_alt
    else
    echo "Directory 00_data_alt already exists."
fi

# create folder for output
if [ ! -d "/home/senekowitsch/Thesis/PopCoGenomeS/output_alt" ]; then
    mkdir -p /home/senekowitsch/Thesis/PopCoGenomeS/output_alt
    echo "Created output directory."
else
    echo "Directory /home/senekowitsch/Thesis/PopCoGenomeS/output_alt already exists."
fi


# set up variables
RAW_DATA_DIR="/home/senekowitsch/raw_data"
FILTERED_DATA_DIR="/home/senekowitsch/Thesis/PopCoGenomeS/00_data_alt"
GENOMES_FILE="/home/senekowitsch/QC/07_filtered_by_clusters/results/filtered_genomes.txt"
OUTPUT_DIR="/home/senekowitsch/Thesis/PopCoGenomeS/output_alt"

# Ensure the destination directory exists
mkdir -p "$FILTERED_DATA_DIR"

# clear out folder 00_data if it already exists
rm -r /home/senekowitsch/Thesis/PopCoGenomeS/00_data_alt/*

# filter RAW_DATA_DIR based on GENOMES_FILE and COPY files to FILTERED_DATA_DIR
while IFS= read -r genome; do
    # Skip empty lines in the genome file
    [[ -z "$genome" ]] && continue

    if [ -f "$RAW_DATA_DIR/$genome" ]; then
        # Changed 'ln -s' to 'cp'
        # -p preserves timestamps and permissions
        cp -p "$RAW_DATA_DIR/$genome" "$FILTERED_DATA_DIR/$genome"
        echo "Copied $genome to $FILTERED_DATA_DIR"
    else
        echo "Warning: $genome not found in $RAW_DATA_DIR"
    fi
done < "$GENOMES_FILE"

# OPTIONAL!!
# If file names contain "_", ".", remove them and rename files to [ID].fna
cd /home/senekowitsch/Thesis/PopCoGenomeS/00_data_alt

# Rename files in 00_data_alt to [ID].fna
for file in GCA_*_genomic.fna; do
    [[ -e "$file" ]] || continue

    # 1. Extract the second segment (007618695.1)
    full_id=$(echo "$file" | cut -d'_' -f2)

    # 2. Remove everything from the dot onwards
    # ${variable%.*} is a Bash trick to delete the shortest match of '.' from the end
    short_id=${full_id%.*}

    mv "$file" "${short_id}.fna"

    # 3. Remove leading zeros from the ID
    short_id_no_zeros=$(echo "$short_id" | sed 's/^0*//')
    mv "${short_id}.fna" "${short_id_no_zeros}.fna"

    echo "Renamed: $file -> ${short_id}.fna"
done


# -----------------------------------------
# run PopCoGenomeS part 1
# -----------------------------------------
conda activate popcogenomes
# export PATH="$PATH:/home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS/src/PopCoGenomeS_part_1/MUMmer3.20"
cd /home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS/src/PopCoGenomeS_part_1
# check if the config file was edited with the user and ask for confirmation to run the script
echo "Please confirm that you have edited the config file with your own settings."
read -p "Do you want to proceed with running PopCoGenomeS? (y/n): " confirm
if [[ $confirm != "y" ]]; then
    echo "Aborted."
    exit 1
fi

rm *.log
rm -r proc/
#rm -r output/
#rm -r /home/senekowitsch/Thesis/PopCoGenomeS/output/*
rm /home/senekowitsch/Thesis/PopCoGenomeS/00_exp_data/*.renamed.mugsy
rm /home/senekowitsch/Thesis/PopCoGenomeS/00_data/*.renamed.mugsy
rm /home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS/example_genomes/*.renamed.mugsy


bash PopCOGenomeS.sh

# -----------------------------------------
# run PopCoGenomeS part 2
# -----------------------------------------
cd /home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS/src/PopCoGenomeS_part_2
source phybreak_config.sh

# remove old output files if they exist
rm -rf /home/senekowitsch/Thesis/PopCoGenomeS/${basename}_cf_size_3.list
rm -rf /home/senekowitsch/Thesis/PopCoGenomeS/example_cf_001
#rm -rf /data/Unit_LMM/selberherr-group/senekowitsch/Thesis/PopCoGenomeS/temp_2/*


while read line; do
 echo "Processing cluster: ${line}"
 bash align_and_construct_trees.sh ${line}
done < ${output_dir}/${basename}_cf_size_3.list

echo "Trees constructed, now running 5x rule analyses"


conda activate popcogenomes_r
source phybreak_config.sh
cd /home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS/src/PopCoGenomeS_part_2
Rscript find_sweeps.R ${basename} ${output_dir}

