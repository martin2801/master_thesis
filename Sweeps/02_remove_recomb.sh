#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Sweeps/02_remove_recombination'
output_base="${base_dir}/output"
output_trees="${output_base}/trees"
output_clonal_alignment="${output_base}/core_alignment.aln"
output_alignment_noambiguous="${output_base}/core_alignment_noambiguous.aln"

input_core_alignment="/home/senekowitsch/Thesis/Sweeps/01_recombination_checks/output/snippy/core.full.aln"
gubbins_bed="/home/senekowitsch/Thesis/Sweeps/01_recombination_checks/output/gubbins_fixed.bed"
cfml_bed="/home/senekowitsch/Thesis/Sweeps/01_recombination_checks/output/cfml_fixed.bed"

tree_prefix="clean_tree"
threads=30

# Setup output directory
cd "${base_dir}"
mkdir -p "${output_base}"

# =============================================================================
# RECOMBINATION MASKING
# =============================================================================
conda activate phylo_pipeline
python3 mask_recombination.py \
  -i "${input_core_alignment}" \
  -b "${gubbins_bed}" \
  -o "${output_clonal_alignment}"

grep -v ">" "${output_clonal_alignment}" | tr -cd 'Nn' | wc -c
grep -v ">" "${input_core_alignment}" | tr -cd 'Nn' | wc -c
awk '/^>out_10177815/{p=1} /^>/ && !/^>out_10177815/{p=0} p && !/^>/{printf $0}' \
    "${output_clonal_alignment}" | tr -cd 'Nn' | wc -c

# =============================================================================
# VERIFY MASKING
# =============================================================================
 
conda activate bedtools
# Get true alignment length (handles multi-line wrapped FASTA)
aln_len=$(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' "${input_core_alignment}")
n_seqs=$(grep -c '>' "${input_core_alignment}")

echo "=== Alignment Info ==="
echo "Alignment length:     ${aln_len} bp"
echo "Number of sequences:  ${n_seqs}"

echo ""
echo "=== N Count Comparison ==="
raw_n=$(grep -v ">" "${input_core_alignment}" | tr -cd 'Nn' | wc -c)
clean_n=$(grep -v ">" "${output_clonal_alignment}" | tr -cd 'Nn' | wc -c)
new_ns=$((clean_n - raw_n))
echo "Ns in raw alignment:   ${raw_n}"
echo "Ns in clean alignment: ${clean_n}"
echo "New Ns added:          ${new_ns}"

echo ""
echo "=== BED File Verification ==="
echo "Max BED coordinate:    $(awk '{if($3>max) max=$3} END{print max}' ${gubbins_bed}) bp"
echo "Raw BED covered bases: $(awk '{sum += $3 - $2} END{print sum}' ${gubbins_bed}) bp  (overlapping regions, not meaningful)"

merged_bases=$(bedtools merge -i "${gubbins_bed}" | awk '{sum += $3 - $2} END{print sum}')
expected_new_ns=$((merged_bases * n_seqs))
echo "Merged BED covered bases:  ${merged_bases} bp  (overlaps removed)"
echo "Expected new Ns (${merged_bases} x ${n_seqs} sequences): ${expected_new_ns}"

echo ""
echo "=== Sanity Check ==="
diff=$(( new_ns - expected_new_ns ))
diff=${diff#-}  # absolute value
pct=$(( 100 * diff / expected_new_ns ))
if [ "${pct}" -lt 5 ]; then
    echo "PASS: New Ns (${new_ns}) are within 5% of expected (${expected_new_ns})."
    echo "      Difference of ${diff} bp is explained by positions already N in the raw alignment."
else
    echo "WARN: New Ns (${new_ns}) differ from expected (${expected_new_ns}) by ${pct}%."
    echo "      Difference: ${diff} bp — investigate further."
fi

# =============================================================================
# VISUALIZE SNP DISTRIBUTION BEFORE AND AFTER RECOMBINATION REMOVAL
# =============================================================================

conda activate phylo_pipeline
python3 snp_density_plot.py \
    -r "${input_core_alignment}" \
    -c "${output_clonal_alignment}" \
    -g "${gubbins_bed}" \
    -f "${cfml_bed}" \
    -o "${output_base}/snp_density.png"

python3 snp_density_plot_pres.py \
    -r "${input_core_alignment}" \
    -c "${output_clonal_alignment}" \
    -g "${gubbins_bed}" \
    -f "${cfml_bed}" \
    -o "${output_base}/snp_density_pres.png"

# =============================================================================
# TREE COMPARISON BEFORE AND AFTER RECOMBINATION REMOVAL
# =============================================================================

conda activate phylo_pipeline

# Tree 1 — raw snippy alignment (pre-recombination removal)
iqtree \
    -s "${input_core_alignment}" \
    --prefix raw_tree \
    -m GTR+G+I \
    -T ${threads} \
    -B 1000 \
    --redo

# Tree 2 — recombination-masked alignment (post-removal)
iqtree \
    -s "${output_clonal_alignment}" \
    --prefix clean_tree \
    -m GTR+G+I \
    -T ${threads} \
    -B 1000 \
    --redo

# Compare the trees using Robinson-Foulds distance
iqtree -rf raw_tree.treefile clean_tree.treefile
# RF 152: 152 splits differ between the two trees, out of a total of 2x117=234 possible splits in the larger tree. 
# An RF distance of 152 / 234 (~0.65) means the trees differ in over half of their splits,
# indicating a substantial topological change, very likely driven by recombination or filtering.

# Visualize the trees side by side using FigTree or iTOL to see how the topology has changed after recombination removal.
# Visualize the trees as tanglegram
python3 tanglegram.py raw_tree.treefile clean_tree.treefile tanglegram_raw_to_clean.png


# =============================================================================
# REMOVE AMBIGUOUS SITES
# =============================================================================

conda activate phylo_pipeline

# Replace N with - so trimAl sees them as gaps
sed '/^>/! s/N/-/g' "${output_clonal_alignment}" > "${output_base}/core_alignment_masked.aln"

# Then remove all gap columns
trimal -in "${output_base}/core_alignment_masked.aln" -out "${output_base}/core_alignment_noambiguous.aln" -nogaps

# Sanity checks:
# Check alignment length before and after
echo "Alignment length with ambiguous sites: $(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' "${output_clonal_alignment}") bp"
echo "Alignment length with Ns masked: $(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' "${output_base}/core_alignment_masked.aln") bp"
echo "Alignment length without ambiguous sites: $(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' "${output_alignment_noambiguous}") bp"

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
unexpected=$(grep -v "^>" "${output_clonal_alignment}" | tr -cd 'ATCGNatcgn-\n' | wc -c)
echo "Unexpected characters in clonal alignment: $(grep -v "^>" "${output_clonal_alignment}" | tr -d 'ATCGNatcgn-\n' | wc -c)"
unexpected_masked=$(grep -v "^>" "${output_base}/core_alignment_masked.aln" | tr -cd 'ATCGNatcgn-\n' | wc -c)
echo "Unexpected characters in masked alignment: $(grep -v "^>" "${output_base}/core_alignment_masked.aln}" | tr -d 'ATCGNatcgn-\n' | wc -c)"
unexpected_noambiguous=$(grep -v "^>" "${output_alignment_noambiguous}" | tr -cd 'ATCGNatcgn-\n' | wc -c)
echo "Unexpected characters in no-ambiguous alignment: $(grep -v "^>" "${output_alignment_noambiguous}" | tr -d 'ATCGNatcgn-\n' | wc -c)"

rm "${output_base}/core_alignment_masked.aln"

# =============================================================================
# RECREATE TREE WITHOUT AMBIGUOUS SITES
# =============================================================================

conda activate phylo_pipeline

# Tree 3 — raw snippy alignment (pre-recombination removal)
iqtree \
    -s "${output_base}/core_alignment_noambiguous.aln" \
    --prefix clean_tree_noambiguous \
    -m GTR+G+I \
    -T ${threads} \
    -B 1000 \
    --redo

# Compare the trees using Robinson-Foulds distance
iqtree -rf raw_tree.treefile clean_tree_noambiguous.treefile
iqtree -rf clean_tree.treefile clean_tree_noambiguous.treefile

# Visualize the trees side by side using FigTree or iTOL to see how the topology has changed after removing ambiguous sites.
python3 tanglegram.py raw_tree.treefile clean_tree_noambiguous.treefile tanglegram_raw_to_noambiguous.png
python3 tanglegram.py clean_tree.treefile clean_tree_noambiguous.treefile tanglegram_clean_to_noambiguous.png

# =============================================================================
# HOUSEKEEPING
# =============================================================================

# move all tree files to a dedicated folder
mkdir -p "${output_trees}"
mv raw_tree* "${output_trees}/"
mv clean_tree* "${output_trees}/"
mv clean_tree_noambiguous* "${output_trees}/"

mv tanglegram* "${output_base}/"

# =============================================================================
# EXTRACT ONLY SNPs FROM THE ALIGNMENT
# =============================================================================

conda activate phylo_pipeline
snp-sites -o "${output_base}/core_alignment_noambiguous_snps.aln" "${output_alignment_noambiguous}"

# check the length of the SNP-only alignment
echo "Length of SNP-only alignment: $(awk '/^>/{if(seq) {print length(seq); exit} seq=""} !/^>/{seq=seq$0}' "${output_base}/core_alignment_noambiguous_snps.aln") bp"

# =============================================================================
# DO PAIRWISE SNP DISTANCE COMPARISON of the SNPs ONLY ALIGNMENT
# =============================================================================

conda activate phylo_pipeline
snp-dists "${output_base}/core_alignment_noambiguous_snps.aln" > "${output_base}/snp_matrix.tab"

# find lowest and highest value in the matrix (excluding self-comparisons)
awk 'NR>1 {for(i=2;i<=NF;i++) if($i != 0) {if(min=="") min=$i; if($i<min) min=$i; if(max=="") max=$i; if($i>max) max=$i}} END{print "Min SNP distance: " min "\nMax SNP distance: " max}' "${output_base}/snp_matrix.tab"

