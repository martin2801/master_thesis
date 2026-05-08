#!/usr/bin/env bash
# =============================================================================
# run_tree_pipeline.sh
#
# Phylogenomics pipeline: snippy → gubbins → masking → trimAl → IQ-TREE → snp-dists
#
# Usage:
#   bash run_tree_pipeline.sh \
#       -f /path/to/genome_folder \
#       -r /path/to/reference.fna \
#       -n YES   # or NO — whether genomes need NCBI-style renaming
#
# Options:
#   -f | --folder        Folder containing input genomes (.fna files)
#   -r | --reference     Path to reference genome (FASTA or GenBank)
#   -n | --rename        Rename NCBI-style filenames? YES or NO
#   -t | --threads       Number of CPU threads (default: 20)
#   -h | --help          Show this help message
#
# Dependencies (must be available via conda environments):
#   snippy_env    : snippy, snippy-core
#   gubbins       : run_gubbins.py
#   bedtools      : bedtools
#   phylo_pipeline: mask_recombination.py, trimal, snp-sites, snp-dists, iqtree
# =============================================================================

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Resolve the directory where this script lives, regardless of where it is
# called from. Must be done before any cd commands change the working directory.
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

# =============================================================================
# DEFAULTS
# =============================================================================
THREADS=20
FOLDER=""
REFERENCE=""
RENAME=""

# =============================================================================
# USAGE
# =============================================================================
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 0
}

# =============================================================================
# ARGUMENT PARSING
# =============================================================================
if [[ $# -eq 0 ]]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--folder)
            [[ -n "${2-}" ]] || { echo "ERROR: -f requires an argument" >&2; exit 1; }
            FOLDER="$2"; shift 2 ;;
        -r|--reference)
            [[ -n "${2-}" ]] || { echo "ERROR: -r requires an argument" >&2; exit 1; }
            REFERENCE="$2"; shift 2 ;;
        -n|--rename)
            [[ -n "${2-}" ]] || { echo "ERROR: -n requires YES or NO" >&2; exit 1; }
            RENAME="$2"; shift 2 ;;
        -t|--threads)
            [[ -n "${2-}" ]] || { echo "ERROR: -t requires a number" >&2; exit 1; }
            THREADS="$2"; shift 2 ;;
        -h|--help)
            usage ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            echo "Run with -h for help." >&2
            exit 1 ;;
    esac
done

# =============================================================================
# VALIDATE INPUTS
# =============================================================================
errors=0

if [[ -z "$FOLDER" ]]; then
    echo "ERROR: -f (genome folder) is required." >&2
    errors=$((errors + 1))
fi

if [[ -z "$REFERENCE" ]]; then
    echo "ERROR: -r (reference genome) is required." >&2
    errors=$((errors + 1))
fi

if [[ -z "$RENAME" ]]; then
    echo "ERROR: -n (rename YES/NO) is required." >&2
    errors=$((errors + 1))
fi

if [[ -n "$FOLDER" && ! -d "$FOLDER" ]]; then
    echo "ERROR: Genome folder does not exist: $FOLDER" >&2
    errors=$((errors + 1))
fi

if [[ -n "$REFERENCE" && ! -f "$REFERENCE" ]]; then
    echo "ERROR: Reference genome file does not exist: $REFERENCE" >&2
    errors=$((errors + 1))
fi

if [[ -n "$RENAME" && "$RENAME" != "YES" && "$RENAME" != "NO" ]]; then
    echo "ERROR: -n must be YES or NO, got: $RENAME" >&2
    errors=$((errors + 1))
fi

if [[ $errors -gt 0 ]]; then
    echo "Run with -h for help." >&2
    exit 1
fi

# Resolve to absolute paths so cd calls don't break relative references
FOLDER=$(realpath "$FOLDER")
REFERENCE=$(realpath "$REFERENCE")

# Strip any trailing slash, then derive prefix from the last path component
FOLDER="${FOLDER%/}"
PREFIX=$(basename "$FOLDER")

# =============================================================================
# CONDA HELPER
# =============================================================================
# conda activate is not available in non-interactive shells without sourcing
# the conda initialisation script first. We do that here once.
_conda_init_done=0
activate_conda() {
    local env_name="$1"
    if [[ $_conda_init_done -eq 0 ]]; then
        local conda_base
        conda_base=$(conda info --base 2>/dev/null) \
            || { echo "ERROR: conda not found in PATH." >&2; exit 1; }
        # shellcheck source=/dev/null
        source "${conda_base}/etc/profile.d/conda.sh"
        _conda_init_done=1
    fi
    conda activate "$env_name"
}

# =============================================================================
# SET UP OUTPUT DIRECTORIES
# =============================================================================
mkdir -p "${PREFIX}"
cd "${PREFIX}"
output_base=$(pwd)
echo "Output base directory: ${output_base}"

output_snippy="${output_base}/snippy"
output_gubbins="${output_base}/gubbins"
output_iqtree="${output_base}/iqtree"
output_clean_alignments="${output_base}/clean_alignments"
output_stats="${output_base}/stats"

mkdir -p \
    "${output_snippy}" \
    "${output_gubbins}" \
    "${output_iqtree}" \
    "${output_clean_alignments}" \
    "${output_stats}"

# =============================================================================
# LOGGING
# Redirect all stdout and stderr — from both echo and external programs — to
# both the terminal and a log file simultaneously, for the rest of the script.
# =============================================================================
LOG_FILE="${output_base}/${PREFIX}_pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1
echo "Logging to: ${LOG_FILE}"
echo "Pipeline started: $(date)"
echo ""

# =============================================================================
# OPTIONAL RENAMING
# Strips NCBI-style GCF_XXXXXXXXX.X_... prefix down to a short numeric ID.
# Example: GCF_000618695.1_genomic.fna -> 618695.fna
# =============================================================================
if [[ "$RENAME" == "YES" ]]; then
    echo "--- Renaming NCBI-style filenames in ${FOLDER} ---"
    for file in "${FOLDER}"/GC*_*_genomic.fna; do
        [[ -e "$file" ]] || { echo "WARNING: No GC*_*_genomic.fna files found — skipping rename."; break; }

        filename=$(basename "$file")

        # Extract the second underscore-delimited segment (e.g. 000618695.1)
        full_id=$(echo "$filename" | cut -d'_' -f2)

        # Drop the dot-suffix (.1)
        short_id="${full_id%.*}"

        # Strip leading zeros
        short_id_no_zeros=$(echo "$short_id" | sed 's/^0*//')

        target="${FOLDER}/${short_id_no_zeros}.fna"
        mv "$file" "$target"
        echo "  Renamed: $filename -> ${short_id_no_zeros}.fna"
    done
else
    echo "--- Skipping rename (--rename NO) ---"
fi

# =============================================================================
# SNIPPY — per-sample runs
# =============================================================================
echo ""
echo "--- Running snippy (per-sample) ---"
activate_conda snippy_env

cd "${output_snippy}"
for genome in "${FOLDER}"/*.fna; do
    [[ -e "$genome" ]] || { echo "ERROR: No .fna files found in ${FOLDER}" >&2; exit 1; }
    sample_name=$(basename "$genome" .fna)
    echo "  Processing: $sample_name"
    snippy \
        --cpus  "$THREADS" \
        --outdir "${output_snippy}/out_${sample_name}" \
        --ref    "$REFERENCE" \
        --ctgs   "$genome"
done

# --- snippy-core: merge all per-sample VCFs into a multi-sample alignment ---
# We explicitly cd to output_snippy so that core.* files land there.
# core.full.aln (all reference positions, not just SNPs) is required by Gubbins.
echo ""
echo "--- Running snippy-core ---"
cd "${output_snippy}"
snippy-core --prefix core --ref "$REFERENCE" "${output_snippy}"/out_*

echo "Snippy complete. Core alignment: ${output_snippy}/core.full.aln"

# =============================================================================
# GUBBINS — recombination detection
# =============================================================================
# Note: Gubbins is most reliable with ≥20 samples. With fewer samples,
# treat recombination predictions with caution.
echo ""
echo "--- Running Gubbins ---"
activate_conda gubbins
cd "${output_gubbins}"

run_gubbins.py \
    --prefix  "${PREFIX}" \
    --threads "$THREADS" \
    --min-snps 3 \
    "${output_snippy}/core.full.aln"

echo "Gubbins complete. Filtered alignment: ${output_gubbins}/${PREFIX}.filtered_polymorphic_sites.aln"

# Convert Gubbins GFF to BED, then fix chromosome label to match snippy-core's
# default label ("core"), which is required for bedtools intersection.
grep -v "^#" "${output_gubbins}/${PREFIX}.recombination_predictions.gff" \
    | awk '{print $1"\t"$4"\t"$5"\t"$9}' \
    | sort -k2,2n \
    > "${output_gubbins}/${PREFIX}_gubbins_regions.bed"

awk '{$1="core"; print}' OFS="\t" \
    "${output_gubbins}/${PREFIX}_gubbins_regions.bed" \
    > "${output_gubbins}/${PREFIX}_gubbins_fixed.bed"

echo "Recombinant regions found by Gubbins: $(wc -l < "${output_gubbins}/${PREFIX}_gubbins_fixed.bed")"
awk '{sum+=$3-$2} END{print "Gubbins total bp (raw, may overlap):", sum}' \
    "${output_gubbins}/${PREFIX}_gubbins_fixed.bed"

# =============================================================================
# MASK RECOMBINATION
# =============================================================================
echo ""
echo "--- Masking recombinant regions ---"
activate_conda phylo_pipeline
cd "${output_base}"

# mask_recombination.py is expected in the parent directory of the output base.
# Adjust the path below if it lives elsewhere.
python3 "${SCRIPT_DIR}/mask_recombination.py" \
    -i "${output_snippy}/core.full.aln" \
    -b "${output_gubbins}/${PREFIX}_gubbins_fixed.bed" \
    -o "${output_clean_alignments}/${PREFIX}_core_alignment.aln"

# --- Verify masking ---
echo ""
echo "--- Verifying masking ---"
activate_conda bedtools

input_core_alignment="${output_snippy}/core.full.aln"
output_clonal_alignment="${output_clean_alignments}/${PREFIX}_core_alignment.aln"
gubbins_bed="${output_gubbins}/${PREFIX}_gubbins_fixed.bed"

aln_len=$(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' \
    "${input_core_alignment}")
n_seqs=$(grep -c '>' "${input_core_alignment}")

echo "=== Alignment Info ==="
echo "Alignment length:     ${aln_len} bp"
echo "Number of sequences:  ${n_seqs}"

raw_n=$(grep -v ">" "${input_core_alignment}"  | tr -cd 'Nn' | wc -c)
clean_n=$(grep -v ">" "${output_clonal_alignment}" | tr -cd 'Nn' | wc -c)
new_ns=$((clean_n - raw_n))

echo ""
echo "=== N Count Comparison ==="
echo "Ns in raw alignment:   ${raw_n}"
echo "Ns in clean alignment: ${clean_n}"
echo "New Ns added:          ${new_ns}"

merged_bases=$(sort -k1,1 -k2,2n "${gubbins_bed}" \
    | bedtools merge -i stdin \
    | awk '{sum += $3 - $2} END{print sum}')
expected_new_ns=$((merged_bases * n_seqs))

echo ""
echo "=== BED File Verification ==="
echo "Max BED coordinate:           $(awk '{if($3>max) max=$3} END{print max}' "${gubbins_bed}") bp"
echo "Merged BED covered bases:     ${merged_bases} bp  (overlaps removed)"
echo "Expected new Ns (${merged_bases} x ${n_seqs} seqs): ${expected_new_ns}"

echo ""
echo "=== Sanity Check ==="
diff=$(( new_ns - expected_new_ns ))
diff=${diff#-}  # absolute value
if [[ $expected_new_ns -gt 0 ]]; then
    pct=$(( 100 * diff / expected_new_ns ))
else
    pct=0
fi

if [[ "${pct}" -lt 5 ]]; then
    echo "PASS: New Ns (${new_ns}) are within 5% of expected (${expected_new_ns})."
    echo "      Difference of ${diff} bp is explained by positions already N in the raw alignment."
else
    echo "WARN: New Ns (${new_ns}) differ from expected (${expected_new_ns}) by ${pct}%."
    echo "      Difference: ${diff} bp — investigate further."
fi

# =============================================================================
# REMOVE AMBIGUOUS SITES
# =============================================================================
echo ""
echo "--- Removing ambiguous sites with trimAl ---"
activate_conda phylo_pipeline

# Replace N with - so trimAl treats masked positions as gaps, then strip all
# gap-only columns. The intermediate masked file is removed afterwards.
sed '/^>/! s/N/-/g' \
    "${output_clean_alignments}/${PREFIX}_core_alignment.aln" \
    > "${output_clean_alignments}/${PREFIX}_core_alignment_masked.aln"

trimal \
    -in  "${output_clean_alignments}/${PREFIX}_core_alignment_masked.aln" \
    -out "${output_clean_alignments}/${PREFIX}_core_alignment_noambiguous.aln" \
    -nogaps

# --- Sanity checks ---
echo ""
echo "=== Alignment Length ==="
echo "  With ambiguous sites: $(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' \
    "${output_clean_alignments}/${PREFIX}_core_alignment.aln") bp"
echo "  With Ns masked:       $(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' \
    "${output_clean_alignments}/${PREFIX}_core_alignment_masked.aln") bp"
echo "  Without ambiguous:    $(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' \
    "${output_clean_alignments}/${PREFIX}_core_alignment_noambiguous.aln") bp"

echo ""
echo "=== N Counts ==="
echo "  Clonal alignment:     $(grep -v ">" "${output_clean_alignments}/${PREFIX}_core_alignment.aln" \
    | tr -cd 'Nn' | wc -c)"
echo "  Masked alignment:     $(grep -v ">" "${output_clean_alignments}/${PREFIX}_core_alignment_masked.aln" \
    | tr -cd 'Nn' | wc -c)"
echo "  No-ambiguous:         $(grep -v ">" "${output_clean_alignments}/${PREFIX}_core_alignment_noambiguous.aln" \
    | tr -cd 'Nn' | wc -c)"

echo ""
echo "=== Gap Counts ==="
echo "  Clonal alignment:     $(grep -v ">" "${output_clean_alignments}/${PREFIX}_core_alignment.aln" \
    | tr -cd '-' | wc -c)"
echo "  Masked alignment:     $(grep -v ">" "${output_clean_alignments}/${PREFIX}_core_alignment_masked.aln" \
    | tr -cd '-' | wc -c)"
echo "  No-ambiguous:         $(grep -v ">" "${output_clean_alignments}/${PREFIX}_core_alignment_noambiguous.aln" \
    | tr -cd '-' | wc -c)"

echo ""
echo "=== Non-ATCGN Characters ==="
for label_file in \
    "Clonal:${output_clean_alignments}/${PREFIX}_core_alignment.aln" \
    "Masked:${output_clean_alignments}/${PREFIX}_core_alignment_masked.aln" \
    "No-ambiguous:${output_clean_alignments}/${PREFIX}_core_alignment_noambiguous.aln"; do
    label="${label_file%%:*}"
    fpath="${label_file##*:}"
    count=$(grep -v ">" "$fpath" | tr -cd 'ATCGNatcgn-' | wc -c)
    total=$(grep -v ">" "$fpath" | tr -cd '[:alpha:]-' | wc -c)
    other=$((total - count))
    echo "  ${label}: ${other} non-ATCGN/gap characters"
done

# Remove intermediate masked alignment to save space
rm "${output_clean_alignments}/${PREFIX}_core_alignment_masked.aln"

# =============================================================================
# IQ-TREE — phylogenetic tree
# =============================================================================
echo ""
echo "--- Building phylogenetic tree with IQ-TREE ---"
activate_conda phylo_pipeline
cd "${output_iqtree}"

iqtree \
    -s "${output_clean_alignments}/${PREFIX}_core_alignment_noambiguous.aln" \
    --prefix "${PREFIX}_clean_tree_noambiguous" \
    -m GTR+G+I \
    -T "${THREADS}" \
    -B 1000 \
    --redo

echo "IQ-TREE complete. Tree: ${output_iqtree}/${PREFIX}_clean_tree_noambiguous.treefile"

# =============================================================================
# SNP DISTANCES
# =============================================================================
echo ""
echo "--- Computing pairwise SNP distances ---"
activate_conda phylo_pipeline

# Extract SNP-only sites (removes invariant columns)
snp-sites \
    -o "${output_stats}/${PREFIX}_core_alignment_noambiguous_snps.aln" \
    "${output_clean_alignments}/${PREFIX}_core_alignment_noambiguous.aln"

echo "SNP-only alignment length: $(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' \
    "${output_stats}/${PREFIX}_core_alignment_noambiguous_snps.aln") bp"

# Pairwise distance matrix
snp-dists \
    "${output_stats}/${PREFIX}_core_alignment_noambiguous_snps.aln" \
    > "${output_stats}/${PREFIX}_snp_matrix.tab"

awk 'NR>1 {
    for(i=2; i<=NF; i++) {
        if ($i+0 != 0) {
            if (min=="" || $i+0 < min+0) min=$i
            if (max=="" || $i+0 > max+0) max=$i
        }
    }
} END {
    print "Min SNP distance: " min
    print "Max SNP distance: " max
}' "${output_stats}/${PREFIX}_snp_matrix.tab"

echo ""
echo "=== Pipeline complete ==="
echo "Results in: ${output_base}"
echo "Pipeline finished: $(date)"
echo "Log saved to: ${LOG_FILE}"