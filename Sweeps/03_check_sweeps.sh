#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Sweeps/03_check_sweeps'
output_base="${base_dir}/output"
output_sweep_alignments="${output_base}/sweep_alignments"
output_fastANI="${output_base}/fastANI"
output_genomes="${output_base}/genomes"
output_trees="${output_base}/trees"
unambiguous_core_alignment="/home/senekowitsch/Thesis/Sweeps/02_remove_recombination/output/core_alignment_noambiguous.aln"
GENOMES="/home/senekowitsch/Thesis/QC/00_data/"
threads=30

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

# Create output directory
mkdir -p "${output_base}"
mkdir -p "${output_sweep_alignments}"
cd "${output_base}"


# =============================================================================
# EXTRACT SWEEP ALIGNMENTS USING SAMTOOLS
# =============================================================================
conda activate check_sweeps

# Index the alignment
samtools faidx "${unambiguous_core_alignment}"

# Extract sweep alignments
samtools faidx "${unambiguous_core_alignment}" "${sweep1_genomes[@]}" > "${output_sweep_alignments}/sweep1.aln"
samtools faidx "${unambiguous_core_alignment}" "${sweep2_genomes[@]}" > "${output_sweep_alignments}/sweep2.aln"


# =============================================================================
# SANITY CHECKS - ARE ALL SEQUENCES THE SAME LENGTH?
# =============================================================================
echo "=== Sweep Alignment Length Checks ==="

for sweep_aln in "${output_sweep_alignments}/sweep1.aln" "${output_sweep_alignments}/sweep2.aln"; do
    echo ""
    echo "--- ${sweep_aln} ---"

    # Get length of every sequence
    lengths=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{if(seq) print length(seq)}' "${sweep_aln}")

    # Count unique lengths
    unique_lengths=$(echo "${lengths}" | sort -u)
    n_unique=$(echo "${unique_lengths}" | wc -l)

    if [ "${n_unique}" -eq 1 ]; then
        echo "PASS: All sequences are the same length ($(echo "${unique_lengths}") bp)"
    else
        echo "FAIL: Sequences have different lengths!"
        # Show which sequence has which length
        awk '/^>/{if(seq) print name"\t"length(seq); name=$0; seq=""} !/^>/{seq=seq$0} END{if(seq) print name"\t"length(seq)}' "${sweep_aln}"
    fi
done


# =============================================================================
# GET CONSENSUS SEQUENCE FOR EACH SWEEP
# =============================================================================
cons -sequence "${output_sweep_alignments}/sweep1.aln" -outseq "${output_sweep_alignments}/sweep1_consensus.fasta" -plurality 0.5 -name "sweep1_consensus"
cons -sequence "${output_sweep_alignments}/sweep2.aln" -outseq "${output_sweep_alignments}/sweep2_consensus.fasta" -plurality 0.5 -name "sweep2_consensus"


# =============================================================================
# RUN FASTANI BETWEEN THE CONSENSUS SEQUENCES AND ALL OTHER GENOMES TO CHECK FOR SIMILARITY
# =============================================================================
conda activate fastANI
echo "Starting ANI comparisons with sweep1 consensus sequence..."
mkdir -p "${output_fastANI}"

# Build genome list
find "${GENOMES}" -name "*.fasta" -o -name "*.fa" -o -name "*.fna" > "${output_fastANI}/genome_list.txt"
echo "Found $(wc -l < "${output_fastANI}/genome_list.txt") genomes"

# Run fastANI for sweep1 consensus vs all genomes
echo "Starting ANI comparison with sweep1 consensus..."
fastANI \
    -q "${output_sweep_alignments}/sweep1_consensus.fasta" \
    --rl "${output_fastANI}/genome_list.txt" \
    -o "${output_fastANI}/sweep1_vs_all.tsv" \
    -t ${threads}

# Run fastANI for sweep2 consensus vs all genomes
echo "Starting ANI comparison with sweep2 consensus..."
fastANI \
    -q "${output_sweep_alignments}/sweep2_consensus.fasta" \
    --rl "${output_fastANI}/genome_list.txt" \
    -o "${output_fastANI}/sweep2_vs_all.tsv" \
    -t ${threads}

echo "Done. Results in ${output_fastANI}"

# =============================================================================
# EXTRACT HIGHEST AND LOWEST ANI SCORING GENOMES FOR EACH SWEEP
# =============================================================================
# check the ANI values for each sweep and find the lowest and highest similarity to other genomes
echo "=== ANI Similarity Checks ==="
for sweep in sweep1 sweep2; do
    echo ""
    echo "--- ${sweep} ---"
    tsv_file="${output_fastANI}/${sweep}_vs_all.tsv"

    if [ -f "${tsv_file}" ]; then
        # Extract ANI values (3rd column) and find min/max
        awk '{print $3}' "${tsv_file}" | sort -n | awk 'NR==1{print "Min ANI: "$0} END{print "Max ANI: "$0}'
    else
        echo "Error: ${tsv_file} not found!"
    fi
done

# write highest and lowest ANI scoring genomes to separate files
echo "=== Top and Bottom ANI Scoring Genomes ==="
for sweep in sweep1 sweep2; do
    tsv_file="${output_fastANI}/${sweep}_vs_all.tsv"

    if [ ! -f "${tsv_file}" ]; then
        echo "Error: ${tsv_file} not found!"
        continue
    fi

    # Sort by ANI (3rd column) descending, write top 20 genome names
    sort -t$'\t' -k3 -rn "${tsv_file}" | head -20 | awk '{print $2}' \
        > "${output_fastANI}/${sweep}_top20_genomes.txt"

    # Sort by ANI (3rd column) ascending, write bottom 20 genome names
    sort -t$'\t' -k3 -n "${tsv_file}" | head -20 | awk '{print $2}' \
        > "${output_fastANI}/${sweep}_bottom20_genomes.txt"

    echo "${sweep}: top and bottom 20 genomes written"
    echo "  Top 20 ANI scores:    $(sort -t$'\t' -k3 -rn "${tsv_file}" | head -20 | awk '{print $3}' | tr '\n' ' ')"
    echo "  Bottom 20 ANI scores: $(sort -t$'\t' -k3 -n "${tsv_file}" | head -20 | awk '{print $3}' | tr '\n' ' ')"
done


# =============================================================================
# MAKE FOLDERS FOR THE HIGHEST AND LOWEST ANI SCORING GENOMES AND COPY THE GENOMES INTO THEM
# =============================================================================
mkdir -p "${output_genomes}/sweep1_top20" "${output_genomes}/sweep1_bottom20" "${output_genomes}/sweep2_top20" "${output_genomes}/sweep2_bottom20"
for sweep in sweep1 sweep2; do
    # Copy top 20 genomes
    while read -r genome_path; do
        cp "${genome_path}" "${output_genomes}/${sweep}_top20/"
    done < "${output_fastANI}/${sweep}_top20_genomes.txt"

    # Copy bottom 20 genomes
    while read -r genome_path; do
        cp "${genome_path}" "${output_genomes}/${sweep}_bottom20/"
    done < "${output_fastANI}/${sweep}_bottom20_genomes.txt"
done

# copy the sweep genomes into the same folders
for sweep in sweep1 sweep2; do
    # Get the sweep genomes array dynamically
    sweep_array="${sweep}_genomes[@]"

    for genome_id in "${!sweep_array}"; do
        # Strip the "out_" prefix to get just the number
        number="${genome_id#out_}"

        # Find matching file in GENOMES folder
        match=$(find "${GENOMES}" -name "*${number}*" \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \))

        if [ -z "${match}" ]; then
            echo "WARN: No match found for ${genome_id} (searching for ${number})"
        else
            echo "Copying ${match} → ${output_genomes}/${sweep}_top20/"
            cp "${match}" "${output_genomes}/${sweep}_top20/"
            cp "${match}" "${output_genomes}/${sweep}_bottom20/"
        fi
    done
done


# =============================================================================
# RUN TREE PIPELINE FOR TOP AND BOTTOM ANI SCORING GENOMES
# =============================================================================

REF='/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/infantis/GCF_000506925.1/ncbi_dataset/data/GCF_000506925.1/GCF_000506925.1_SI119944_genomic.fna'

# Run for all folders in the output_genomes directory
for folder in "${output_genomes}"/*; do
    bash ../run_tree_pipeline.sh \
        -f "$folder" \
        -r "${REF}" \
        -n YES \
        -t 20
done

# Remove temporary genome folders to save space
rm -rf "${output_genomes}/sweep1_top20" "${output_genomes}/sweep1_bottom20" "${output_genomes}/sweep2_top20" "${output_genomes}/sweep2_bottom20"

# Remove snippy output folders to save space
rm -rf "${output_base}"/sweep*/snippy/out_*
