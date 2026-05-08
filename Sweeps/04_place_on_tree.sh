#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree'
output_base='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output'

input_tree='/home/senekowitsch/Thesis/Sweeps/02_remove_recombination/output/trees/clean_tree_noambiguous.treefile'
input_alignment='/home/senekowitsch/Thesis/Sweeps/02_remove_recombination/output/core_alignment_noambiguous.aln'
input_bed="/home/senekowitsch/Thesis/Sweeps/01_recombination_checks/output/gubbins_fixed.bed"
input_cfml="/home/senekowitsch/Thesis/Sweeps/01_recombination_checks/output/cfml_fixed.bed"

output_full_alignment="${output_base}/snippy_alignment"
output_clonal_alignment="${output_base}/masked.aln"
output_alignment_noambiguous="${output_base}/core_alignment_noambiguous.aln"
output_trees="${output_base}/trees"

RAW_DATA="/home/senekowitsch/raw_data"
GENOMES="${base_dir}/genomes"
genomes_filter="/home/senekowitsch/Thesis/QC/05_filter/results/filtered_samples.txt"
reference_genome='/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/infantis/GCF_000506925.1/ncbi_dataset/data/GCF_000506925.1/GCF_000506925.1_SI119944_genomic.fna'

full_alignment="${output_full_alignment}/core.full.aln"
REF='/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/infantis/GCF_000506925.1/ncbi_dataset/data/GCF_000506925.1/GCF_000506925.1_SI119944_genomic.fna'
threads=40

# =============================================================================
# PREPARATION
# =============================================================================

# Define sweep genomes
sweep1_genomes=(
    "out_11452375"
    "out_11452785"
    "out_11453035"
    "out_8046645"
    "out_8091485"
)
sweep2_genomes=(
    "out_14966075"
    "out_8963375"
    "out_31832435"
    "out_8067535"
    "out_8476865"
)

# Create a list of all sweep genomes for easy checking later
cd "${base_dir}"
for genome in "${sweep1_genomes[@]}"; do
    echo -e "${genome}\tsweep1"
done > sweeps.txt

for genome in "${sweep2_genomes[@]}"; do
    echo -e "${genome}\tsweep2"
done >> sweeps.txt

# create output directories
cd "${base_dir}"
mkdir -p "${GENOMES}"
mkdir -p "${output_base}"
mkdir -p "${output_full_alignment}"
mkdir -p "${output_trees}"

# copy genomes from genomes_filtered to GENOMES directory
while IFS= read -r sample; do
    src="${RAW_DATA}/${sample}.fna"
    if [[ -f "$src" ]]; then
        cp "$src" "${GENOMES}/"
        echo "Copied: ${sample}.fna"
    else
        echo "WARNING: Not found: ${src}"
    fi
done < "${genomes_filter}"

echo "Done. Genomes copied: $(ls ${GENOMES}/*.fna | wc -l)"


# Renameing according to popcogenomes workflow!!!!!
cd "${GENOMES}"
for file in GC*_*_genomic.fna; do
    [[ -e "$file" ]] || continue

    # 1. Extract the second segment (007618695.1)
    full_id=$(echo "$file" | cut -d'_' -f2)

    # 2. Remove everything from the dot onwards
    short_id=${full_id%.*}

    # 3. Remove leading zeros from the ID
    short_id_no_zeros=$(echo "$short_id" | sed 's/^0*//')

    # 4. Rename directly to final format
    mv "$file" "${short_id_no_zeros}.fna"

    echo "Renamed: $file -> ${short_id_no_zeros}.fna"
done
echo "Done. Genomes renamed: $(ls ${GENOMES}/*.fna | wc -l)"
cd "${base_dir}"

# =============================================================================
# RUN SNIPPY
# =============================================================================
# Each genome is processed independently against the named reference.
# Output is written to a sample-specific subdirectory (out_<sample_name>).
conda activate snippy_env
cd "${output_full_alignment}"
for genome in "$GENOMES"/*.fna; do
    # Extract the filename without the path and extension for the folder name
    sample_name=$(basename "$genome" .fna)
    
    echo "Processing $sample_name..."

    # Key flags:
    #   --cpus   number of CPU threads
    #   --outdir per-sample output directory
    #   --ref    reference genome (FASTA or GenBank)
    #   --ctgs   input assembled contigs (instead of reads)
    snippy --cpus "$threads" \
           --outdir "$output_full_alignment/out_$sample_name" \
           --ref "$reference_genome" \
           --ctgs "$genome"
done


# --- Core alignment ---
# snippy-core merges all per-sample VCFs and produces a multi-sample alignment.
# core.full.aln includes ALL reference positions (not just polymorphic sites),
# which is required by Gubbins (it needs invariant sites to model substitution
# rates correctly).
snippy-core --prefix core --ref "$reference_genome" $output_full_alignment/*

# set variable to the core alignment for the next step
full_alignment="${output_full_alignment}/core.full.aln"

# remove the temporary per-sample output directories to save space
rm -rf "${output_full_alignment}/out_"*

echo "Snippy complete. Core alignment generated at ${full_alignment}"

# =============================================================================
# RECOMBINATION MASKING
# =============================================================================
cd "${base_dir}"
conda activate phylo_pipeline
python3 mask_recombination.py \
    -i "${full_alignment}" \
    -b "${input_bed}" \
    -o "${output_base}/masked.aln"

# set variable to the masked alignment for the next step
output_clonal_alignment="${output_base}/masked.aln"
echo "Recombination masking complete. Masked alignment generated at ${output_clonal_alignment}"

# =============================================================================
# VERIFY MASKING
# =============================================================================
conda activate phylo_pipeline
cd "${base_dir}"
python3 snp_density_plot.py \
    -r "${full_alignment}" \
    -c "${output_clonal_alignment}" \
    -g "${input_bed}" \
    -f "${input_cfml}" \
    -o "${output_base}/snp_density.png"

# =============================================================================
# REMOVE AMBIGUOUS SITES
# =============================================================================
cd "${base_dir}"
conda activate phylo_pipeline

# Replace N with - so trimAl sees them as gaps
sed '/^>/! s/N/-/g' "${output_clonal_alignment}" > "${output_base}/core_alignment_masked.aln"

# Then remove all gap columns
trimal -in "${output_base}/core_alignment_masked.aln" -out "${output_base}/core_alignment_noambiguous.aln" -nogaps
output_alignment_noambiguous="${output_base}/core_alignment_noambiguous.aln"
echo "Ambiguous site removal complete. Final alignment generated at ${output_alignment_noambiguous}"

# Sanity checks:
# Check alignment length before and after
echo "Alignment length with ambiguous sites:"
awk '/^>/{if(NR>1) exit} !/^>/{count+=length($0)} END{print count}' "${output_clonal_alignment}"
echo "Alignment length with Ns masked:"
awk '/^>/{if(NR>1) exit} !/^>/{count+=length($0)} END{print count}' "${output_base}/core_alignment_masked.aln"
echo "Alignment length without ambiguous sites:"
awk '/^>/{if(NR>1) exit} !/^>/{count+=length($0)} END{print count}' "${output_alignment_noambiguous}"

# Check N count before, during and after
echo "Ns in clonal alignment: $(grep -v ">" "${output_clonal_alignment}" | tr -cd 'Nn' | wc -c)"
echo "Ns in masked alignment: $(grep -v ">" "${output_base}/core_alignment_masked.aln" | tr -cd 'Nn' | wc -c)"
echo "Ns in no-ambiguous alignment: $(grep -v ">" "${output_alignment_noambiguous}" | tr -cd 'Nn' | wc -c)"

# Check how many letters there are, that are not A, T, C, G or N
echo "ATCGN characters in clonal alignment: $(grep -v ">" "${output_clonal_alignment}" | tr -cd 'ATCGNatcgn' | wc -c)"
echo "ATCGN characters in masked alignment: $(grep -v ">" "${output_base}/core_alignment_masked.aln" | tr -cd 'ATCGNatcgn' | wc -c)"
echo "ATCGN characters in no-ambiguous alignment: $(grep -v ">" "${output_alignment_noambiguous}" | tr -cd 'ATCGNatcgn' | wc -c)"

# Show count and identity of unexpected characters
# Check for unexpected characters (anything other than ATCG, N, gaps, and newlines)
count=$(grep -v "^>" "${output_clonal_alignment}" | tr -d 'ATCGNatcgn-' | tr -d '\n' | wc -c)
echo "Unexpected characters in clonal alignment: ${count}"
count=$(grep -v "^>" "${output_base}/core_alignment_masked.aln" | tr -d 'ATCGNatcgn-' | tr -d '\n' | wc -c)
echo "Unexpected characters in masked alignment: ${count}"
count=$(grep -v "^>" "${output_alignment_noambiguous}" | tr -d 'ATCGNatcgn-' | tr -d '\n' | wc -c)
echo "Unexpected characters in no-ambiguous alignment: ${count}"

rm "${output_base}/core_alignment_masked.aln"
rm "${output_base}/masked.aln"

# =============================================================================
# SPLIT ALIGNMENT INTO REFERENCE (117+ref) AND QUERY SEQUENCES
# =============================================================================
conda activate phylo_pipeline

# Extract sequence names from the original 117+ref alignment
grep ">" "${input_alignment}" | sed 's/>//' > "${output_base}/ref_names.txt"
echo "Reference sequences (117 + ref): $(wc -l < ${output_base}/ref_names.txt)"

# Split into reference and query alignments using python
python3 - <<EOF
from Bio import SeqIO

ref_names = set(line.strip() for line in open("${output_base}/ref_names.txt"))

all_seqs = list(SeqIO.parse("${output_alignment_noambiguous}", "fasta"))
print(f"Total sequences in full alignment: {len(all_seqs)}")

ref_seqs   = [rec for rec in all_seqs if rec.id in ref_names]
query_seqs = [rec for rec in all_seqs if rec.id not in ref_names]

# Check for any reference sequences missing from the new alignment
found_refs = set(rec.id for rec in ref_seqs)
missing    = ref_names - found_refs
if missing:
    print(f"WARNING: {len(missing)} reference sequences not found in new alignment:")
    for name in sorted(missing):
        print(f"  {name}")

SeqIO.write(ref_seqs,   "${output_base}/ref_alignment.aln",   "fasta")
SeqIO.write(query_seqs, "${output_base}/query_alignment.aln", "fasta")

print(f"Reference sequences written: {len(ref_seqs)}")
print(f"Query sequences written:     {len(query_seqs)}")
EOF

echo "Split complete."
echo "  Reference alignment: ${output_base}/ref_alignment.aln"
echo "  Query alignment:     ${output_base}/query_alignment.aln"

# =============================================================================
# RECREATE TREE WITHOUT AMBIGUOUS SITES
# =============================================================================
conda activate phylo_pipeline
cd "${output_trees}"
# Tree 3 — raw snippy alignment (pre-recombination removal)
iqtree \
    -s "${output_base}/ref_alignment.aln" \
    --prefix ref_clean_tree_noambiguous \
    -m GTR+G+I \
    -T ${threads} \
    -B 1000 \
    --redo

# Compare the trees using Robinson-Foulds distance
iqtree -rf ${input_tree} \
           ${output_trees}/ref_clean_tree_noambiguous.treefile 

python3 "${base_dir}/tanglegram.py" \
        ${input_tree} \
        ${output_trees}/ref_clean_tree_noambiguous.treefile \
        tanglegram_rawold_to_rawnew_noambiguous.png

cd "${base_dir}"

# =============================================================================
# EPA-NG — PLACE QUERY SEQUENCES ON REFERENCE TREE
# =============================================================================
conda activate epa-ng
mkdir -p "${output_base}/epa_output/"

epa-ng \
    --tree "${output_trees}/ref_clean_tree_noambiguous.treefile" \
    --ref-msa "${output_base}/ref_alignment.aln" \
    --query "${output_base}/query_alignment.aln" \
    --model "${output_trees}/ref_clean_tree_noambiguous.iqtree" \
    --outdir "${output_base}/epa_output/" \
    --redo \
    --threads "${threads}"


echo "EPA-ng complete. Results in ${output_base}/epa_output/"

gappa examine graft \
    --jplace-path "${output_base}/epa_output/epa_result.jplace" \
    --out-dir "${output_base}/gappa_output/"

gappa edit extract \
    --jplace-path "${output_base}/epa_output/epa_result.jplace" \
    --clade-list-file "${base_dir}/sweeps.txt" \
    --fasta-path "${output_base}/query_alignment.aln" \
    --threshold 0.95 \
    --color-tree-file "${output_base}/gappa_output/clade_tree.svg" \
    --samples-out-dir "${output_base}/gappa_output/samples/" \
    --sequences-out-dir "${output_base}/gappa_output/sequences/" \
    --allow-file-overwriting \
    --threads "${threads}"

# Check what sequences are in each clade fasta file
grep ">" "${output_base}/gappa_output/sequences/sweep1.fasta" >> "${output_base}/gappa_output/sweep1_samples.txt"
grep ">" "${output_base}/gappa_output/sequences/sweep2.fasta" >> "${output_base}/gappa_output/sweep2_samples.txt"

# =============================================================================
# VALIDATE ADDITIONS TO SWEEPS
# =============================================================================
conda activate phylo_pipeline
mkdir -p "${output_base}/validate_sweeps/"
mkdir -p "${output_base}/validate_sweeps/genomes_sweep1/"
mkdir -p "${output_base}/validate_sweeps/genomes_sweep2/"

# Copy the newly identified sweep genomes to the validation folders
# Copy the newly identified sweep genomes to the validation folders
while IFS= read -r sample; do
    # Strip ">out_" from the beginning of the string
    clean_sample="${sample#>out_}"
    
    src="${GENOMES}/${clean_sample}.fna"
    
    if [[ -f "$src" ]]; then
        cp "$src" "${output_base}/validate_sweeps/genomes_sweep1/"
        echo "Copied to sweep1 validation: ${clean_sample}.fna"
    else
        echo "WARNING: Not found for sweep1 validation: ${src}"
    fi
done < "${output_base}/gappa_output/sweep1_samples.txt"

while IFS= read -r sample; do
    # Strip ">out_" from the beginning of the string
    clean_sample="${sample#>out_}"
    
    src="${GENOMES}/${clean_sample}.fna"
    
    if [[ -f "$src" ]]; then
        cp "$src" "${output_base}/validate_sweeps/genomes_sweep2/"
        echo "Copied to sweep2 validation: ${clean_sample}.fna"
    else
        echo "WARNING: Not found for sweep2 validation: ${src}"
    fi
done < "${output_base}/gappa_output/sweep2_samples.txt"

# Copy original sweep genomes to the validation folders
for sample in "${sweep1_genomes[@]}"; do
    # This removes the prefix "out_" from the string
    clean_sample="${sample#out_}"
    
    src="${GENOMES}/${clean_sample}.fna"
    if [[ -f "$src" ]]; then
        cp "$src" "${output_base}/validate_sweeps/genomes_sweep1/"
        echo "Copied original sweep1 genome: ${clean_sample}.fna"
    else
        echo "WARNING: Not found for original sweep1 genome: ${src}"
    fi
done

for sample in "${sweep2_genomes[@]}"; do
    # This removes the prefix "out_" from the string
    clean_sample="${sample#out_}"
    
    src="${GENOMES}/${clean_sample}.fna"
    if [[ -f "$src" ]]; then
        cp "$src" "${output_base}/validate_sweeps/genomes_sweep2/"
        echo "Copied original sweep2 genome: ${clean_sample}.fna"
    else
        echo "WARNING: Not found for original sweep2 genome: ${src}"
    fi
done

# Add genomes up to 25 total in each sweep folder by randomly sampling from the remaining genomes
# Target total sample size
TARGET_SIZE=25

for sweep in sweep1 sweep2; do
    echo "Filling ${sweep} validation folder to ${TARGET_SIZE} genomes..."
    
    # 1. Count how many genomes we already have
    current_count=$(ls -1 "${output_base}/validate_sweeps/genomes_${sweep}/"*.fna 2>/dev/null | wc -l)
    needed=$((TARGET_SIZE - current_count))
    
    if [ "$needed" -le 0 ]; then
        echo "Folder ${sweep} already has ${current_count} genomes. No additions needed."
        continue
    fi
    
    echo "Current count: ${current_count}. Adding ${needed} random genomes."

    # 2. Get a list of all genomes, exclude those already in the folder, 
    # shuffle them, and take the top N needed.
    # Note: We strip the path to compare just the filenames.
    ls -1 "${GENOMES}/"*.fna | while read -r full_path; do
        filename=$(basename "$full_path")
        if [[ ! -f "${output_base}/validate_sweeps/genomes_${sweep}/${filename}" ]]; then
            echo "$full_path"
        fi
    done | shuf -n "$needed" | while read -r random_src; do
        cp "$random_src" "${output_base}/validate_sweeps/genomes_${sweep}/"
        echo "Added random genome to ${sweep}: $(basename "$random_src")"
    done
done

# Run tree pipeline for the sweep genomes to validate their placement on the tree
# Run for all folders in the validate_sweeps directory
cd "${output_base}/validate_sweeps/"
for folder in "${output_base}/validate_sweeps/"*; do
    bash /home/senekowitsch/Thesis/Sweeps/03_check_sweeps/run_tree_pipeline.sh \
        -f "$folder" \
        -r "${REF}" \
        -n NO \
        -t "${threads}"
done

# Combine the genomes and run the tree again
mkdir -p "${output_base}/validate_sweeps/combined/"
cp "${output_base}/validate_sweeps/genomes_sweep1/"*.fna "${output_base}/validate_sweeps/combined/"
cp "${output_base}/validate_sweeps/genomes_sweep2/"*.fna "${output_base}/validate_sweeps/combined/"

# Check number of genomes in combined folder
echo "Total genomes in combined folder: $(ls -1 "${output_base}/validate_sweeps/combined/"*.fna | wc -l)"

bash /home/senekowitsch/Thesis/Sweeps/03_check_sweeps/run_tree_pipeline.sh \
    -f "${output_base}/validate_sweeps/combined/" \
    -r "${REF}" \
    -n NO \
    -t "${threads}"

# Remove the snippy alignments to save space
rm -rf "${output_base}"/validate_sweeps/genomes_sweep*/snippy/out_*
rm -rf "${output_base}"/validate_sweeps/combined/snippy/out_*

# Remove genomes from the combined folder to save space
rm -rf "${output_base}/validate_sweeps/combined/"*.fna

# Remove genomes from the individual sweep folders to save space
rm -rf "${output_base}"/validate_sweeps/genomes_sweep*/*.fna


# =============================================================================
# RUN IQ-TREE ON FULL ALIGNMENT
# =============================================================================

conda activate phylo_pipeline

iqtree \
    -s "${output_base}/core_alignment_noambiguous.aln" \
    --prefix full_tree \
    -m GTR+G+I \
    -T "${threads}" \
    -B 1000 \
    --redo

