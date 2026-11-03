#!/bin/bash
# ------------------------------------------------------------------
# Plasmid identification on FASTA assemblies (Infantis)
#   1. minimap2 against a closed pESI reference  -> pESI contigs + pESI coverage per genome
#   2. Platon                                    -> plasmid vote
#   3. geNomad                                   -> plasmid vote (+ prophage annotation for later)
#   4. MOB-suite (existing mob_recon outputs)    -> plasmid vote
#   5. plasmid_consensus.py                      -> contig labels, genome summary, cleaned chromosomes
# ------------------------------------------------------------------
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh
THREADS=30

# ---- inputs (edit these) ----
GENOMES="/home/senekowitsch/Thesis/Functional_2/00_data"   # *.fna
EXT=".fna"
PESI_REF="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/pESI_CP016409.fasta"     # a closed, complete pESI plasmid sequence

PLATON_DB="/home/senekowitsch/miniconda3/envs/platon/DB/db"                      # platon db (download once, see Platon README)
GENOMAD_DB="/home/senekowitsch/miniconda3/envs/genomad/DB/genomad_db"                    # genomad download-database <dir>

# ---- outputs ----
base_dir='/home/senekowitsch/Thesis/Functional_2/01_plasmid_id'
output_base="${base_dir}/output"

output_mob_recon="${output_base}/mob_recon"
output_minimap2="${output_base}/minimap2"
output_platon="${output_base}/platon"
output_genomad="${output_base}/genomad"
output_consensus="${output_base}/consensus"

# ---- create folders ----
mkdir -p "${output_base}"
mkdir -p "${output_mob_recon}"
mkdir -p "${output_minimap2}"
mkdir -p "${output_platon}"
mkdir -p "${output_genomad}"
mkdir -p "${output_consensus}"

# ------------------------------------------------------------------
# Download the pESI reference sequence from NCBI (if not already done)
# ------------------------------------------------------------------
cd ${output_base}
conda activate edirect
[ -s "${PESI_REF}" ] || efetch -db nuccore -id CP016409.1 -format fasta > pESI_CP016409.fasta
efetch -db nuccore -id CP016409.1 -format fasta > pESI_CP016409.fasta
# or: datasets download genome / NCBI web "Send to → FASTA"

# ------------------------------------------------------------------
# 1. minimap2 vs pESI reference
# ------------------------------------------------------------------
conda activate minimap2          # env with minimap2
cd ${output_base}
for fna in "${GENOMES}"/*${EXT}; do
    name=$(basename "${fna}" "${EXT}")
    paf="${output_minimap2}/${name}.paf"
    [ -s "${paf}" ] && continue
    minimap2 -x asm10 -c --secondary=no -t "${THREADS}" "${PESI_REF}" "${fna}" > "${paf}"
done

# ------------------------------------------------------------------
# 2. Platon
# ------------------------------------------------------------------
conda activate platon

# download the Platon database once (see Platon README)
cd /home/senekowitsch/miniconda3/envs/platon
mkdir -p DB
cd DB
wget https://zenodo.org/record/4066768/files/db.tar.gz
tar -xzf db.tar.gz
rm db.tar.gz


for fna in "${GENOMES}"/*${EXT}; do
    name=$(basename "${fna}" "${EXT}")
    [ -s "${output_platon}/${name}/${name}.tsv" ] && continue
    platon --db "${PLATON_DB}" \
        --output "${output_platon}/${name}" \
        --prefix "${name}" \
        --threads "${THREADS}" \
        "${fna}"
done

# ------------------------------------------------------------------
# 3. geNomad (also gives provirus calls -> reuse for prophage analysis)
# ------------------------------------------------------------------
conda activate genomad

cd /home/senekowitsch/miniconda3/envs/genomad
mkdir DB
cd DB
genomad download-database .

cd ${output_base}
for fna in "${GENOMES}"/*${EXT}; do
    name=$(basename "${fna}" "${EXT}")
    [ -d "${output_genomad}/${name}/${name}_summary" ] && continue
    genomad end-to-end --cleanup --threads "${THREADS}" \
        "${fna}" "${output_genomad}/${name}" "${GENOMAD_DB}"
done

# ------------------------------------------------------------------
# 4. MOB-suite (mob_recon)
# ------------------------------------------------------------------
conda activate pESI
for fna in "${GENOMES}"/*${EXT}; do
    genome=$(basename "${fna}" .fna)
    echo "Running mob_recon for ${genome}..."
    mob_recon \
        -i "${fna}" \
        -o "${output_mob_recon}/${genome}_mob_recon_output" \
        -n "${THREADS}" \
        -f
    echo "Done: ${genome}"
done

# ------------------------------------------------------------------
# 4. Consensus
# ------------------------------------------------------------------
conda activate phylo_pipeline    # any env with python3 >= 3.8 (standard library only)
'''
python3 /home/senekowitsch/Salmonella/AGES/Functional/07_plasmid_id/plasmid_consensus.py \
    --genomes       "${GENOMES}" \
    --ext           "${EXT}" \
    --paf-dir       "${OUT}/minimap2" \
    --platon-dir    "${OUT}/platon" \
    --genomad-dir   "${OUT}/genomad" \
    --mob-dir       "${MOB_DIR}" \
    --abricate-list "${ABRICATE_LIST}" \
    --pos-cov       0.8 \
    --partial-cov   0.5 \
    --out           "${OUT}/consensus"
'''

cd ${base_dir}
python3 plasmid_consensus_2.py \
    --genomes       "${GENOMES}" \
    --ext           "${EXT}" \
    --paf-dir       "${output_minimap2}" \
    --platon-dir    "${output_platon}" \
    --genomad-dir   "${output_genomad}" \
    --mob-dir       "${output_mob_recon}" \
    --pos-cov       0.8 \
    --partial-cov   0.5 \
    --pesi-include-short \
    --out           "${output_consensus}"

# Main outputs:
#   consensus/contig_labels.tsv      one row per contig: label + evidence
#   consensus/genome_summary.tsv     pESI coverage + status per genome, comparison to abricate
#   consensus/chromosomes/<g>.fasta  chromosome contigs only (pESI + other plasmids removed)
#   -> use chromosomes/*.fasta as Bakta input for the new pangenome

# ------------------------------------------------------------------
# 5. List of pESI+ genomes (minimap2 / pESI reference coverage)
# ------------------------------------------------------------------
SUMMARY="${output_consensus}/genome_summary.tsv"

awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="status") c=i; next}
            $c=="pESI_pos" {print $1}' "${SUMMARY}" \
    | sort > "${output_base}/pESI_genomes_infantis_minimap2.txt"

# keep the uncertain ones separately, for checking
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="status") c=i; next}
            $c=="pESI_partial" {print $1}' "${SUMMARY}" \
    | sort > "${output_base}/pESI_partial_genomes_infantis_minimap2.txt"

wc -l "${output_base}"/pESI_*genomes_infantis_minimap2.txt

# compare with the old abricate list
echo "only abricate: $(comm -23 <(sort "${ABRICATE_LIST}") "${output_base}/pESI_genomes_infantis_minimap2.txt" | wc -l)"
echo "only minimap2: $(comm -13 <(sort "${ABRICATE_LIST}") "${output_base}/pESI_genomes_infantis_minimap2.txt" | wc -l)"

# ------------------------------------------------------------------
# visualization of the tree with highlighted AGES and sweep7 clades
# ------------------------------------------------------------------
infantis_tree="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/validate_sweeps/full_tree.treefile"
infantis_tree_renamed="${output_base}/tree/full_tree_renamed.treefile"
infantis_tree_pesi="${output_base}/tree/infantis_pESI_tree.treefile"
pesi_genomes="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/pESI_genomes_infantis_minimap2.txt"
partial_pesi_genomes="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/pESI_partial_genomes_infantis_minimap2.txt"

# copy the tree and strip the "out_" prefix from leaf names
# (only at the start of a leaf name, i.e. after "(" or ",", optionally quoted)
mkdir -p "${output_base}/tree"
sed -E "s/([(,]'?)out_/\1/g" "${infantis_tree}" > "${infantis_tree_renamed}"

# sanity check: no "out_" left, and how many pESI genomes are found as leaves
echo "remaining 'out_' labels: $(grep -o 'out_' "${infantis_tree_renamed}" | wc -l)"
grep -oE "[(,]'?[^(),:;']+" "${infantis_tree_renamed}" | sed -E "s/^[(,]'?//" | sort -u > "${output_base}/tree/leaves.txt"
echo "pESI+ genomes in tree:     $(comm -12 "${output_base}/tree/leaves.txt" <(sort -u "${pesi_genomes}") | wc -l) / $(wc -l < "${pesi_genomes}")"
echo "partial pESI in tree:      $(comm -12 "${output_base}/tree/leaves.txt" <(sort -u "${partial_pesi_genomes}") | wc -l) / $(wc -l < "${partial_pesi_genomes}")"


cd ${base_dir}
conda activate phylo_pipeline
python3 plot_tree.py \
    --tree "${infantis_tree_renamed}" \
    --highlight "${pesi_genomes}" "#00BCD4" "pESI+" \
    --highlight "${partial_pesi_genomes}" "#1B7837" "partial pESI+" \
    --labels-all \
    --out "${output_base}"/infantis_norecomb_tree_with_pesi.png

python3 /home/senekowitsch/Salmonella/AGES/Functional/plot_tree.py \
    --tree "${infantis_tree_pesi}" \
    --highlight "${pesi_genomes}" "#00BCD4" "pESI+" \
    --highlight "${partial_pesi_genomes}" "#1B7837" "partial pESI+" \
    --labels-all \
    --out "${output_base}"/infantis_pesi_tree_with_pesi.png