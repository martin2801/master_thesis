#!/bin/bash
#set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Functional/01_prokka'
output_base="${base_dir}/output"

GENOMES='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/genomes'
SWEEPS='/home/senekowitsch/Thesis/Sweeps/05_check_distance/output/sweeps_bottomup_clonal_5x.txt'
extra_sweeps='/home/senekowitsch/Thesis/Functional/01_prokka/sweep_7.txt'
THREADS=30

# Create output directory
mkdir -p "${base_dir}"
mkdir -p "${output_base}"
cd "${base_dir}"

# Output file
OUTPUT="${base_dir}/genome_sweep_labels.txt"
mkdir -p "${base_dir}"
 
# ---------------------------------------------------------------------
# Create a mapping of genome IDs to sweep labels based on the sweeps file, then label each genome accordingly.
# ---------------------------------------------------------------------
# Build a lookup: genome_id -> sweep_label from the sweeps file
# Skip header, strip "out_" prefix, format as sweep_N
declare -A sweep_map
while IFS=$'\t' read -r sweep_id genome rest; do
    [[ "${sweep_id}" == "sweep_id" ]] && continue   # skip header
    genome_clean="${genome#out_}"                    # strip "out_" prefix
    sweep_map["${genome_clean}"]="sweep_${sweep_id}"
done < "${SWEEPS}"
 
# Write header
echo -e "genome\tsweep" > "${OUTPUT}"
 
# Iterate over all .fna files in GENOMES directory
for fna in "${GENOMES}"/*.fna; do
    fname=$(basename "${fna}" .fna)   # e.g. 32036575
    if [[ -v sweep_map["${fname}"] ]]; then
        label="${sweep_map[${fname}]}"
    else
        label="no_sweep"
    fi
    echo -e "${fname}\t${label}"
done >> "${OUTPUT}"
 
echo "Done. Output written to ${OUTPUT}"
 
# Quick summary
echo ""
echo "=== Label counts ==="
tail -n +2 "${OUTPUT}" | awk '{print $2}' | sort | uniq -c | sort -k2

# Add sweep_7 labels to the existing genome_sweep_labels.txt file, overwriting any existing label for those genomes.
python3 << 'EOF'
sweep7 = set()
with open('/home/senekowitsch/Thesis/Functional/01_prokka/sweep_7.txt') as f:
    for line in f:
        g = line.strip()
        if g:
            sweep7.add(g)

# Read existing labels
rows = []
with open('/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt') as f:
    header = f.readline()
    for line in f:
        line = line.strip()
        if not line or '\t' not in line:
            continue
        genome, label = line.split('\t')
        rows.append([genome, label])

# Apply sweep_7 labels
updated = 0
for row in rows:
    if row[0] in sweep7:
        row[1] = 'sweep_7'
        updated += 1

# Write back
with open('/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt', 'w') as f:
    f.write('genome\tsweep\n')
    for row in rows:
        f.write(f'{row[0]}\t{row[1]}\n')

print(f"Updated {updated} genomes to sweep_7")
print(f"Total genomes written: {len(rows)}")

from collections import Counter
counts = Counter(r[1] for r in rows)
for label, count in sorted(counts.items()):
    print(f"  {label}: {count}")
EOF

# ---------------------------------------------------------------------
# Create a mapping of genome IDs to sweep labels based on the sweeps file, then label each genome accordingly.
# ---------------------------------------------------------------------

LABELS="${base_dir}/genome_sweep_labels.txt"
REPRESENTATIVES="${base_dir}/prokka_representatives.txt"

python3 << 'EOF'
import random
import sys
 
random.seed(42)  # reproducible selection
 
labels_file = '/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt'
output_file = '/home/senekowitsch/Thesis/Functional/01_prokka/prokka_representatives.txt'
 
# Group genomes by sweep label
from collections import defaultdict
groups = defaultdict(list)
with open(labels_file) as f:
    next(f)  # skip header
    for line in f:
        genome, label = line.strip().split('\t')
        groups[label].append(genome)
 
selected = []
for label, genomes in sorted(groups.items()):
    n = min(3, len(genomes))
    chosen = random.sample(genomes, n)
    for g in chosen:
        selected.append((g, label))
    print(f"{label}: {chosen}")
 
with open(output_file, 'w') as f:
    f.write('genome\tsweep\n')
    for genome, label in selected:
        f.write(f'{genome}\t{label}\n')
 
print(f"\nTotal representatives selected: {len(selected)}")
EOF


echo ""
echo "=== Running Prokka on representatives ==="
 
# Step 2: Run Prokka on each representative
while IFS=$'\t' read -r genome label; do
    [[ "${genome}" == "genome" ]] && continue  # skip header
 
    fna="${GENOMES}/${genome}.fna"
    outdir="${output_base}/${label}/${genome}"
 
    if [[ ! -f "${fna}" ]]; then
        echo "WARNING: ${fna} not found, skipping."
        continue
    fi
 
    if [[ -d "${outdir}" ]]; then
        echo "Skipping ${genome} (${label}) — output already exists."
        continue
    fi
 
    echo "Annotating ${genome} (${label})..."
    mkdir -p "${outdir}"
 
    prokka \
        --outdir "${outdir}" \
        --prefix "${genome}" \
        --genus Salmonella \
        --species Infantis \
        --kingdom Bacteria \
        --cpus "${THREADS}" \
        --force \
        "${fna}"
 
    echo "Done: ${genome}"
done < "${REPRESENTATIVES}"
 
echo ""
echo "=== All done! ==="
echo "Output in: ${output_base}"
echo ""
echo "Key files per genome:"
echo "  *.faa  -> proteins for eggNOG-mapper"
echo "  *.ffn  -> nucleotide gene sequences for BLAST"
echo "  *.gff  -> annotation"
echo "  *.gbk  -> GenBank format"