#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Functional/07_pESI'
output_base="${base_dir}/output"
GENOMES='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/genomes'
abricate_results_dir="${output_base}/abricate_results"
abricate_summary_file="${output_base}/abricate_linecounts.txt"
sweep_labels='/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt'
PMLST_DB='/home/senekowitsch/miniconda3/envs/pmlst/share/pmlst/db'
THREADS=30
gene_count_cutoff=220
# for tree visualization
treefile="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/validate_sweeps/full_tree.treefile"
SWEEPS="/home/senekowitsch/Thesis/Sweeps/05_check_distance/output/sweeps_bottomup_clonal_5x.txt"

mkdir -p "${base_dir}"
mkdir -p "${output_base}"

cd "${base_dir}"

# ------------------------------------------
# Identify pESI plasmid in genome assemblies
# ------------------------------------------

# Identify plasmid sequences from genomes
conda activate pESI

# Database stored here:
ls -lh /home/senekowitsch/miniconda3/envs/pESI/lib/python3.11/site-packages/mob_suite/databases

# Run mob_recon for each genome
for fna in "${GENOMES}"/*.fna; do
    genome=$(basename "${fna}" .fna)
    echo "Running mob_recon for ${genome}..."
    mob_recon \
        -i "${fna}" \
        -o "${output_base}/${genome}_mob_recon_output" \
        -n "${THREADS}" \
        -f
    echo "Done: ${genome}"
done

mkdir -p "${output_base}/mob_recon_outputs"
mv "${output_base}"/*_mob_recon_output "${output_base}/mob_recon_outputs/"

# Download reference for pESI plasmid
cd "${base_dir}"
efetch -db nucleotide \
  -id NZ_CP016411.1 \
  -format fasta_cds_na \
  > NZ_CP016411.1_cds.fasta

sed 's/>lcl|/>/' NZ_CP016411.1_cds.fasta > NZ_CP016411.1_cds_clean.fa

# Create db for abricate from pESI reference
cd /home/senekowitsch/miniconda3/envs/pESI/db
mkdir -p pESI_db
cd pESI_db
cp /home/senekowitsch/Thesis/Functional/07_pESI/NZ_CP016411.1_cds_clean.fa sequences
ls
abricate --setupdb
abricate --list
cd "${base_dir}"

# Run abricate on all identified plasmids with pESI database
mkdir -p "${output_base}/abricate_results"

for plasmid in "${output_base}/mob_recon_outputs"/*_mob_recon_output/plasmid_*.fasta; do
    genome=$(basename $(dirname "${plasmid}") _mob_recon_output)
    plasmid_id=$(basename "${plasmid}" .fasta)
    echo "Running abricate for ${genome} ${plasmid_id}..."
    abricate \
        --db pESI_db \
        "${plasmid}" \
        > "${output_base}/abricate_results/${genome}_${plasmid_id}_abricate.txt"
done

abricate --summary "${output_base}/abricate_results/"*_abricate.txt \
    > "${output_base}/abricate_summary.txt"

ls "${output_base}/mob_recon_outputs"/*_mob_recon_output/plasmid_*.fasta | wc -l


echo -e "genome\hits" > "${abricate_summary_file}"

for f in "${abricate_results_dir}"/*_abricate.txt; do
    genome=$(basename "${f}" _abricate.txt)
    lines=$(( $(wc -l < "${f}") - 1 ))
    echo -e "${genome}\t${lines}"
done >> "${abricate_summary_file}"

cat "${abricate_summary_file}"

# Keep only lines where the cut off is met
awk -F'\t' -v cutoff="${gene_count_cutoff}" \
    'NR==1 || $2 >= cutoff' \
    "${abricate_summary_file}" \
    > "${output_base}/abricate_pESI_present.txt"

cat "${output_base}/abricate_pESI_present.txt"

# Create a list of genomes with pESI present
awk -F'\t' 'NR>1 {split($1, a, "_"); print a[1]}' \
    "${output_base}/abricate_pESI_present.txt" \
    > "${output_base}/abricate_pESI_genomes.txt"

cat "${output_base}/abricate_pESI_genomes.txt"

# Compare to sweep labels and add labels to it
awk '
  NR==FNR {
    if (FNR>1) { id=$1; sub(/_plasmid_.*/, "", id); present[id]=1 }
    next
  }
  FNR==1 { print $0 "\tpESI"; next }
  { print $0 "\t" ($1 in present ? "YES" : "NO") }
' "${output_base}/abricate_pESI_present.txt" "${sweep_labels}" \
  > "${output_base}/sweep_labels_pESI.txt"

cat "${output_base}/sweep_labels_pESI.txt"

# Create a tree that highlights the genomes with pESI present
pESI_present_file="${output_base}/sweep_labels_pESI.txt"

awk -F'\t' '$NF=="YES" {print "out_" $1}' "${output_base}/sweep_labels_pESI.txt" \
    > "${output_base}/pesi_genomes.txt"

cd "${base_dir}"
conda activate phylo_pipeline
python3 plot_tree_pesi.py \
    --tree "${treefile}" \
    --sweeps "${SWEEPS}" \
    --pesi "${output_base}/pesi_genomes.txt" \
    --out "${output_base}/tree_pESI.png" \
    --labels \
    --dpi 400



# ------------------------------------------
# Identify other plasmids present in genomes
# ------------------------------------------

# Create a combined table of all mobtyper results
out="${output_base}/mob_recon_outputs"
{ head -n1 "$(ls ${out}/*_mob_recon_output/mobtyper_results.txt | head -1)"
  for f in ${out}/*_mob_recon_output/mobtyper_results.txt; do
      tail -n +2 "$f"
  done
} > "${output_base}/all_mobtyper_results.tsv"

# Run abricate against plasmidfinder database for all identified plasmids
abricate --db plasmidfinder "${out}"/*_mob_recon_output/plasmid_*.fasta \
    > "${output_base}/plasmidfinder_results.tsv"

abricate --summary "${output_base}/plasmidfinder_results.tsv" \
    > "${output_base}/plasmidfinder_summary.tsv"


# One line per plasmid with its replicon(s) listed
awk -F'\t' '
  NR==1 {
    for (i=3; i<=NF; i++) col[i]=$i          # replicon names from header
    print "genome\tplasmid\treplicons"
    next
  }
  {
    n=split($1, parts, "/")
    genome="NA"
    for (i=1; i<=n; i++)
      if (parts[i] ~ /_mob_recon_output$/) {
        genome=parts[i]; sub(/_mob_recon_output$/, "", genome)
      }
    plasmid=parts[n]; sub(/\.fasta$/, "", plasmid)

    reps=""
    for (i=3; i<=NF; i++)
      if ($i != "." && $i != "")
        reps = (reps=="" ? col[i] : reps ";" col[i])
    if (reps=="") reps="none"

    print genome "\t" plasmid "\t" reps
  }
' "${output_base}/plasmidfinder_summary.tsv" \
  > "${output_base}/plasmid_inventory.tsv"

# Count how many plasmids carry each replicon type
awk -F'\t' 'NR>1 {
    n=split($3, a, ";")
    for (i=1; i<=n; i++) count[a[i]]++
  }
  END { for (r in count) print count[r] "\t" r }' \
  "${output_base}/plasmid_inventory.tsv" \
  | sort -rn > "${output_base}/replicon_counts.tsv"

# Drop the pESI plasmids (handled separately above)
awk -F'\t' 'NR==1 || $3 !~ /IncFIB\(pN55391\)/' \
    "${output_base}/plasmid_inventory.tsv" \
    > "${output_base}/plasmid_inventory_no_pESI.tsv"


# Lets look at IncI1_1_Alpha replicon, as it is in all 3 sweep 1 genomes
# Pull plasmids out that have this replicon
awk -F'\t' 'NR==1 || $3 ~ /IncI1_1_Alpha/' \
    "${output_base}/plasmid_inventory.tsv" \
    > "${output_base}/incI1_plasmids.tsv"

cat "${output_base}/incI1_plasmids.tsv"

# Collect the coresponding fasta files
mkdir -p "${output_base}/incI1_fastas"
awk -F'\t' 'NR>1 {print $1"\t"$2}' "${output_base}/incI1_plasmids.tsv" \
  | while IFS=$'\t' read genome plasmid; do
      src="${output_base}/mob_recon_outputs/${genome}_mob_recon_output/${plasmid}.fasta"
      cp "$src" "${output_base}/incI1_fastas/${genome}_${plasmid}.fasta"
  done

conda activate pmlst
pmlst-download-db
ls /home/senekowitsch/miniconda3/envs/pmlst/share/pmlst/db/

mkdir -p "${output_base}/incI1_pmlst"

for fa in "${output_base}/incI1_fastas/"*.fasta; do
    name=$(basename "$fa" .fasta)
    echo "pMLST: ${name}"
    mkdir -p "${output_base}/incI1_pmlst/${name}"
    pmlst -i "$fa" \
          -o "${output_base}/incI1_pmlst/${name}" \
          -s inci1 \
          -p "${PMLST_DB}" \
          -mp $(which blastn)
done

python3 << 'EOF'
import json, glob, os, csv

output_base = '/home/senekowitsch/Thesis/Functional/07_pESI/output'
pmlst_dir   = f'{output_base}/incI1_pmlst'
loci        = ['repI1', 'ardA', 'sogS', 'pilL', 'trbA']

rows = []
for jf in sorted(glob.glob(f'{pmlst_dir}/*/data.json')):
    name = os.path.basename(os.path.dirname(jf))   # folder name = genome_plasmid
    with open(jf) as fh:
        d = json.load(fh)['pmlst']['results']
    alleles = d.get('allele_profile', {})
    # split "10177815_plasmid_AA474" -> genome, plasmid
    genome  = name.split('_plasmid_')[0]
    plasmid = 'plasmid_' + name.split('_plasmid_')[1] if '_plasmid_' in name else name
    row = {
        'genome'        : genome,
        'plasmid'       : plasmid,
        'sequence_type' : d.get('sequence_type', 'NA'),
        'nearest_sts'   : d.get('nearest_sts', 'NA'),
        'clonal_complex': d.get('clonal_complex', 'NA'),
    }
    for locus in loci:
        row[locus] = alleles.get(locus, {}).get('allele', '-')
    rows.append(row)

cols = ['genome','plasmid','sequence_type','nearest_sts','clonal_complex'] + loci
out = f'{output_base}/incI1_pmlst_summary.tsv'
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
    w.writeheader()
    w.writerows(rows)

print(f'Wrote {len(rows)} plasmids -> {out}\n')
# quick console view
widths = {c: max(len(c), max((len(str(r[c])) for r in rows), default=0)) for c in cols}
print('  '.join(c.ljust(widths[c]) for c in cols))
for r in rows:
    print('  '.join(str(r[c]).ljust(widths[c]) for c in cols))

# tally of STs / clonal complexes
from collections import Counter
print('\nClonal complex tally:')
for cc, n in Counter(r['clonal_complex'] for r in rows).most_common():
    print(f'  {cc or "none":<12} {n}')
EOF

