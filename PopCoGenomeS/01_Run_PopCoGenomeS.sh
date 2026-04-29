# !/bin/bash

# go to PopCoGenomeS folder
cd /home/senekowitsch/Thesis/PopCoGenomeS

# create folder structure if not exists
if [ ! -d "00_data" ]; then
    mkdir -p 00_data
    else
    echo "Directory 00_data already exists."
fi

# create folder for output
if [ ! -d "/home/senekowitsch/Thesis/PopCoGenomeS/output" ]; then
    mkdir -p /home/senekowitsch/Thesis/PopCoGenomeS/output
    echo "Created output directory."
else
    echo "Directory /home/senekowitsch/Thesis/PopCoGenomeS/output already exists."
fi


# set up variables
RAW_DATA_DIR="/home/senekowitsch/raw_data"
FILTERED_DATA_DIR="/home/senekowitsch/Thesis/PopCoGenomeS/00_data"
GENOMES_FILE="/home/senekowitsch/QC/07_filtered_by_clusters/results/filtered_genomes.txt"

# Ensure the destination directory exists
mkdir -p "$FILTERED_DATA_DIR"

# clear out folder 00_data if it already exists
rm -rf "$FILTERED_DATA_DIR"/*

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

conda activate popcogenomes

# -----------------------------------------
# run PopCoGenomeS part 1
# -----------------------------------------
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
rm -rf proc/
#rm -rf output/
#rm -rf /home/senekowitsch/Thesis/PopCoGenomeS/output/*
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
#rm -rf /data/Unit_LMM/selberherr-group/senekowitsch/Thesis/PopCoGenomeS/temp/*


while read line; do
 echo "Processing cluster: ${line}"
 bash align_and_construct_trees.sh ${line}
done < ${output_dir}/${basename}_cf_size_3.list

echo "If it fails (00_prepare_phybreak.sh: 20: [[: not found) make sure to change line 33 (sh 00_prepare_phybreak.sh) in align_and_construct_trees.sh to bash 00_prepare_phybreak.sh"


'
Failed with:
Done reading maf file. Starting alignment filtering.
Done filtering, finding gap-columns and SNP locations in alignment.
Done finding gap-columns and counting SNPs
Done removing gap-columns
Writing information about alignment blocks to file.
Traceback (most recent call last):
File "phybreak2.maf_to_fasta.py", line 327, in ‹module> corefile write(">"+iso +"\n"+ full segdict [iso] +"\n")
KeyError:
GCA_0076186951_PDT0005387451_genomic
Loading required package: ape
Error in seq-start:seq-end : NA/NaN argument Calls: read.fasta →> gsub →> is. factor -> paste
Execution halted
'

cd /data/Unit_LMM/selberherr-group/senekowitsch/Thesis/PopCoGenomeS/temp/example_cf_001/align_and_construct_trees
ls phybreak_parameters.txt  # confirm it's there
python /home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS/src/PopCoGenomeS_part_2/align_and_construct_trees/phybreak2.maf_to_fasta.py