#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Functional/02_eggnog'
output_base="${base_dir}/output"
eggnogg_db_dir='/data/Unit_LMM/selberherr-group/senekowitsch/Thesis/eggnog_db'

PROKKA_OUT='/home/senekowitsch/Thesis/Functional/01_prokka/output'
THREADS=30

mkdir -p "${base_dir}"
cd "${base_dir}"
mkdir -p "${output_base}"
mkdir -p "${eggnogg_db_dir}"

conda activate eggnog

# Update environment variable
export EGGNOG_DATA_DIR=/data/Unit_LMM/selberherr-group/senekowitsch/Thesis/eggnog_db


# Download eggnog-mapper database
cd "${eggnogg_db_dir}"
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz && gunzip eggnog.db.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz && tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz && gunzip eggnog_proteins.dmnd.gz

cd "${base_dir}"

# Run eggnog-mapper
for faa in "${PROKKA_OUT}"/*/*/*.faa; do
    # Extract genome ID from filename
    genome=$(basename "${faa}" .faa)
    # Extract sweep label from directory structure (e.g. output/sweep_1/32036575/32036575.faa)
    label=$(basename $(dirname $(dirname "${faa}")))
    outdir="${output_base}/${label}/${genome}"

    if [[ -d "${outdir}" ]]; then
        echo "Skipping ${genome} (${label}) — already done."
        continue
    fi

    mkdir -p "${outdir}"
    echo "Running eggNOG-mapper on ${genome} (${label})..."

    emapper.py \
        -i "${faa}" \
        --itype proteins \
        -o "${genome}" \
        --output_dir "${outdir}" \
        --cpu "${THREADS}" \
        --tax_scope bacteria \
        --override

    echo "Done: ${genome}"
done

echo "=== All done! ==="‚