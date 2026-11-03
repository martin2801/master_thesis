

# set up variables
RAW_DATA_DIR="/data/Unit_LMM/selberherr-group/senekowitsch/Thesis/raw_data"
FILTERED_DATA_DIR="/home/senekowitsch/Thesis/Functional_2/00_data"
GENOMES_FILE="/home/senekowitsch/Thesis/QC/05_filter/results/filtered_samples.txt"


mkdir -p "$FILTERED_DATA_DIR"
# clear out folder 00_data if it already exists
rm -r ${FILTERED_DATA_DIR}/*

# turn any name form into the short numeric ID
short_id() {
    local x=${1##*/}       # drop any path
    x=${x#GCA_}            # drop GCA_ prefix
    x=${x%%[._]*}          # keep up to the first . or _
    echo $((10#$x))        # drop leading zeros (10# avoids octal parsing)
}

# map short ID -> raw file
declare -A raw
for f in "$RAW_DATA_DIR"/GCA_*_genomic.fna; do
    raw[$(short_id "$f")]=$f
done

# link every genome listed in GENOMES_FILE
n_ok=0; n_missing=0
while IFS= read -r line || [[ -n $line ]]; do
    line=${line//$'\r'/}; line=${line//[[:space:]]/}   # strip CR/whitespace
    [[ -z $line || $line == \#* ]] && continue
    id=$(short_id "$line")
    if [[ -n ${raw[$id]} ]]; then
        ln -sn "${raw[$id]}" "$FILTERED_DATA_DIR/$id.fna" && ((n_ok++))
    else
        echo "not in raw_data: $line" >&2; ((n_missing++))
    fi
done < "$GENOMES_FILE"

echo "linked: $n_ok   missing: $n_missing"
















# Ensure the destination directory exists
mkdir -p "$FILTERED_DATA_DIR"



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
