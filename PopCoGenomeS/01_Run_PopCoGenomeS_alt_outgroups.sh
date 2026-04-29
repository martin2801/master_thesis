# !/bin/bash

# go to PopCoGenomeS folder
cd /home/senekowitsch/Thesis/PopCoGenomeS

# create folder structure if not exists
if [ ! -d "00_data_alt_outgroups" ]; then
    mkdir -p 00_data_alt_outgroups
    else
    echo "Directory 00_data_alt_outgroups already exists."
fi

# create folder for output
if [ ! -d "/home/senekowitsch/Thesis/PopCoGenomeS/output_alt_outgroups" ]; then
    mkdir -p /home/senekowitsch/Thesis/PopCoGenomeS/output_alt_outgroups
    echo "Created output directory."
else
    echo "Directory /home/senekowitsch/Thesis/PopCoGenomeS/output_alt_outgroups already exists."
fi


# set up variables
RAW_DATA_DIR="/home/senekowitsch/raw_data"
FILTERED_DATA_DIR="/home/senekowitsch/Thesis/PopCoGenomeS/00_data_alt_outgroups"
OUTGROUPS_DIR="/home/senekowitsch/Thesis/PopCoGenomeS/00_outgroups"
GENOMES_FILE="/home/senekowitsch/QC/07_filtered_by_clusters/results/filtered_genomes.txt"
OUTPUT_DIR="/home/senekowitsch/Thesis/PopCoGenomeS/output_alt_outgroups"

# Ensure the destination directory exists
mkdir -p "$FILTERED_DATA_DIR"
mkdir -p "$OUTGROUPS_DIR"

# clear out folder 00_data_alt_outgroups if it already exists
rm -r /home/senekowitsch/Thesis/PopCoGenomeS/00_data_alt_outgroups/*

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

# Download outgroup genomes from NCBI and place them in 00_data_alt_outgroups
# activate conda environment with ncbi-genome-download installed
conda activate ncbi-datasets
# download outgroup genomes from NCBI and place them in 00_data_alt_outgroups
# Salmonella Virchow
datasets download genome accession GCA_016852765.1 --filename $OUTGROUPS_DIR/GCA_016852765.1.zip
datasets download genome accession GCF_000516855.1 --filename $OUTGROUPS_DIR/GCF_000516855.1.zip
datasets download genome accession GCF_000171535.2 --filename $OUTGROUPS_DIR/GCF_000171535.2.zip
# unzip the downloaded files, if not already unzipped
if [ ! -d "$OUTGROUPS_DIR/GCA_016852765.1" ]; then
    unzip $OUTGROUPS_DIR/GCA_016852765.1.zip -d $OUTGROUPS_DIR/GCA_016852765.1
else
    echo "Directory $OUTGROUPS_DIR/GCA_016852765.1 already exists."
fi

if [ ! -d "$OUTGROUPS_DIR/GCF_000516855.1" ]; then
    unzip $OUTGROUPS_DIR/GCF_000516855.1.zip -d $OUTGROUPS_DIR/GCF_000516855.1
else
    echo "Directory $OUTGROUPS_DIR/GCF_000516855.1 already exists."
fi

if [ ! -d "$OUTGROUPS_DIR/GCF_000171535.2" ]; then
    unzip $OUTGROUPS_DIR/GCF_000171535.2.zip -d $OUTGROUPS_DIR/GCF_000171535.2
else
    echo "Directory $OUTGROUPS_DIR/GCF_000171535.2 already exists."
fi

# copy the outgroup genomes to 00_data_alt_outgroups
cp $OUTGROUPS_DIR/GCA_016852765.1/ncbi_dataset/data/GCA_016852765.1/GCA_016852765.1_PDT000935349.1_genomic.fna $FILTERED_DATA_DIR/
cp $OUTGROUPS_DIR/GCF_000171535.2/ncbi_dataset/data/GCF_000171535.2/GCF_000171535.2_ASM17153v2_genomic.fna $FILTERED_DATA_DIR/
cp $OUTGROUPS_DIR/GCF_000516855.1/ncbi_dataset/data/GCF_000516855.1/GCF_000516855.1_SvirSVQ1_v1.0_genomic.fna $FILTERED_DATA_DIR/

# OPTIONAL!!
# If file names contain "_", ".", remove them and rename files to [ID].fna
cd /home/senekowitsch/Thesis/PopCoGenomeS/00_data_alt_outgroups

# Rename files in 00_data_alt_outgroups to [ID].fna
for file in GC*_*_genomic.fna; do
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

    echo "Renamed: $file -> ${short_id_no_zeros}.fna"
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
rm -r /home/senekowitsch/Thesis/PopCoGenomeS/output_alt_outgroups
rm /home/senekowitsch/Thesis/PopCoGenomeS/00_data_alt_outgroups/*.renamed.mugsy
#rm /home/senekowitsch/Thesis/PopCoGenomeS/00_exp_data/*.renamed.mugsy
#rm /home/senekowitsch/Thesis/PopCoGenomeS/00_data/*.renamed.mugsy
#rm /home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS/example_genomes/*.renamed.mugsy

bash PopCOGenomeS.sh

# -----------------------------------------
# run PopCoGenomeS part 2
# -----------------------------------------
cd /home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS/src/PopCoGenomeS_part_2
source phybreak_config.sh

# remove old output files if they exist
#rm -rf /data/Unit_LMM/selberherr-group/senekowitsch/Thesis/PopCoGenomeS/temp_3/*

conda activate popcogenomes
while read line; do
 echo "Processing cluster: ${line}"
 bash align_and_construct_trees.sh ${line}
done < ${output_dir}/${basename}_cf_size_3.list

echo "Trees constructed, now running 5x rule analyses"

conda activate popcogenomes_r
source phybreak_config.sh
cd /home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS/src/PopCoGenomeS_part_2
Rscript find_sweeps.R ${basename} ${output_dir}

