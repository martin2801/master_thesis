#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Functional/03_blast'
output_base="${base_dir}/output"
PROKKA_OUT='/home/senekowitsch/Thesis/Functional/01_prokka/output'
GENOMES='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/genomes'
THREADS=30

mkdir -p "${base_dir}"
mkdir -p "${output_base}"
cd "${base_dir}"

# Step 1: Build BLAST database with genome IDs in headers
echo "=== Building BLAST database ==="
#cat "${GENOMES}"/*.fna > all_genomes.fna

for fna in "${GENOMES}"/*.fna; do
    genome=$(basename "${fna}" .fna)
    sed "s/^>/>${genome}|/" "${fna}"
done > all_genomes.fna

echo "=== Building BLAST database ==="
conda activate blast
makeblastdb \
    -in all_genomes.fna \
    -dbtype nucl \
    -out all_genomes_db \
    -title "all_genomes"
echo "Done."

# Step 2: Loop over each representative genome's .ffn and BLAST separately
echo "=== Running blastn ==="
for ffn in "${PROKKA_OUT}"/*/*/*.ffn; do
    genome=$(basename "${ffn}" .ffn)
    label=$(basename $(dirname $(dirname "${ffn}")))
    outfile="${base_dir}/${label}_${genome}_blast.tsv"

    if [[ -f "${outfile}" ]]; then
        echo "Skipping ${genome} (${label}) — already done."
        continue
    fi

    echo "BLASTing ${genome} (${label})..."
    blastn \
        -query "${ffn}" \
        -db all_genomes_db \
        -out "${outfile}" \
        -outfmt "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore" \
        -max_target_seqs 500 \
        -max_hsps 1 \
        -perc_identity 90 \
        -qcov_hsp_perc 80 \
        -num_threads "${THREADS}" \
        -evalue 1e-10
    echo "Done: ${genome}"
done 

mv *.tsv "${output_base}/"
echo "All BLAST searches completed. Results moved to ${output_base}."