# !/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh

# go to QC folder
cd /home/senekowitsch/Thesis/QC

# activate conda env seqsero
source $(conda info --base)/etc/profile.d/conda.sh
conda activate seqsero

# move to the right folder
cd /home/senekowitsch/Thesis/QC

# set path to the output folder
OUTPUT="/home/senekowitsch/Thesis/QC/02_seqsero/results"

# run seqsero2 on files mentioned in file 01_ANI/results/fastANI_all_vs_infantis_reference_output_below_99.txt
while IFS= read -r line; do
    file=$(echo $line | awk '{print $1}')
    echo "Processing file: $file"
    basename=$(basename $file .fna)
    mkdir -p $OUTPUT/$basename
    SeqSero2_package.py -m k -t 4 -i $file -d $OUTPUT/$basename
done < /home/senekowitsch/Thesis/QC/01_ANI/results/fastANI_all_vs_infantis_reference_output_below_99.txt

# concatenate the results into one file, first line is header, so we need to skip it for all files except the first one
first_file=true
for file in $OUTPUT/*/SeqSero_result.tsv; do
    if [ "$first_file" = true ]; then
        cat "$file" > "$OUTPUT/SeqSero_result_all.tsv"
        first_file=false
    else
        tail -n +2 "$file" >> "$OUTPUT/SeqSero_result_all.tsv"
    fi
done

