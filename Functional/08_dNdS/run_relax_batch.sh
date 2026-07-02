#!/usr/bin/env bash
# Build a job list from relax_inputs/_manifest.tsv and run hyphy relax
# in parallel for every "ok" (gene, sweep) pair.
#
# Assumes extract_matching_alignments.py has already been run, so each
# relax_inputs/<sweep>/<gene>.fasta has a matching .nwk with identical
# tip sets (RELAX requires this exact match).
#
# Usage:
#   ./run_relax_batch.sh <manifest.tsv> <relax_inputs_dir> <output_dir> <n_jobs>
#
# Example:
#   ./run_relax_batch.sh relax_inputs/_manifest.tsv relax_inputs relax_results 8

set -euo pipefail

MANIFEST="${1:?manifest tsv required}"
TREE_DIR="${2:?relax_inputs dir required (also contains matching .fasta files)}"
OUT_DIR="${3:?output dir required}"
N_JOBS="${4:-4}"

mkdir -p "$OUT_DIR"
JOBLIST="$(mktemp)"

# Build one shell command per "ok" row, skip header
awk -F'\t' 'NR>1 && $3=="ok"{print $1"\t"$2}' "$MANIFEST" | while IFS=$'\t' read -r SWEEP GENE; do
    mkdir -p "${OUT_DIR}/${SWEEP}"
    OUT_JSON="${OUT_DIR}/${SWEEP}/${GENE}.RELAX.json"
    # skip if already done (lets you resume an interrupted run)
    if [[ -s "$OUT_JSON" ]]; then
        continue
    fi
    ALN="${TREE_DIR}/${SWEEP}/${GENE}.fasta"
    TREE="${TREE_DIR}/${SWEEP}/${GENE}.nwk"
    echo "hyphy relax --alignment '${ALN}' --tree '${TREE}' --test Foreground --code Universal --output '${OUT_JSON}' > '${OUT_DIR}/${SWEEP}/${GENE}.log' 2>&1" >> "$JOBLIST"
done

N_TOTAL=$(wc -l < "$JOBLIST")
echo "Built ${N_TOTAL} jobs (already-completed genes skipped)."
echo "Running with ${N_JOBS} parallel workers..."

parallel -j "$N_JOBS" --bar --joblog "${OUT_DIR}/_parallel_joblog.tsv" < "$JOBLIST"

rm -f "$JOBLIST"

echo "Done. Check ${OUT_DIR}/_parallel_joblog.tsv for per-job exit codes/timing."
echo "Non-zero exit codes indicate failed hyphy runs - check the matching .log file."