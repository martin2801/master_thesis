#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Functional/09_Fst'
output_base="${base_dir}/output"
GENOMES='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/genomes'
sweep_labels='/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt'
no_rec_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/core_alignment_noambiguous.aln"
raw_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/snippy_alignment/core.full.aln"
masked_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/masked.aln"
THREADS=30
reference_genome='/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/infantis/GCF_000506925.1/ncbi_dataset/data/GCF_000506925.1/GCF_000506925.1_SI119944_genomic.fna'

SWEEPS="/home/senekowitsch/Thesis/Sweeps/05_check_distance/output/sweeps_bottomup_clonal_5x.txt"

mkdir -p "${base_dir}"
mkdir -p "${output_base}"

cd "${base_dir}"

# ------------------------------------------
# Prepare vcf file for Fst calculation
# ------------------------------------------
conda activate vcftools
snp-sites -v -o "${output_base}/core_alignment_noambiguous.vcf" "${no_rec_alignment}"
snp-sites -v -o "${output_base}/core_alignment_full.vcf" "${raw_alignment}"
snp-sites -v -o "${output_base}/core_alignment_masked.vcf" "${masked_alignment}"

conda activate bcftools
bcftools +fixploidy "${output_base}/core_alignment_noambiguous.vcf" -- -f 2 > "${output_base}/pseudo_diploid.vcf"
bcftools +fixploidy "${output_base}/core_alignment_full.vcf" -- -f 2 > "${output_base}/pseudo_diploid_full.vcf"
bcftools +fixploidy "${output_base}/core_alignment_masked.vcf" -- -f 2 > "${output_base}/pseudo_diploid_masked.vcf"

# ------------------------------------------
# Prepare grouping files
# ------------------------------------------
cd "${output_base}"
sweep_no="7"
awk '{if ($2=="sweep_'${sweep_no}'") print "out_"$1 > "sweep_'${sweep_no}'.txt"; else print "out_"$1 > "other_sweep_'${sweep_no}'.txt"}' "${sweep_labels}"
wc -l sweep_'${sweep_no}'.txt other_sweep_${sweep_no}.txt
ls "sweep_${sweep_no}.txt"
# Sanity checks
grep -v "^out_" "sweep_${sweep_no}.txt"   # should return nothing if prefixing worked
conda activate bcftools
bcftools query -l pseudo_diploid.vcf > vcf_samples.txt
comm -23 <(sort "sweep_${sweep_no}.txt") <(sort vcf_samples.txt)   # IDs in sweep_${sweep_no} but not in VCF

# ------------------------------------------
# Run Fst calculation
# ------------------------------------------
# No recombination alignment
conda activate vcftools
vcftools --vcf pseudo_diploid.vcf \
  --weir-fst-pop "sweep_${sweep_no}.txt" \
  --weir-fst-pop "other_sweep_${sweep_no}.txt" \
  --out "${output_base}"/fst_results_sweep_${sweep_no}
# visualize results
conda activate ani_heatmap
cd "${output_base}"
python3 ../visualize_fst.py --input fst_results_sweep_${sweep_no}.weir.fst --out fst_plot_sweep_${sweep_no}.png

# Full alignment
conda activate vcftools
vcftools --vcf pseudo_diploid_full.vcf \
  --weir-fst-pop "sweep_${sweep_no}.txt" \
  --weir-fst-pop "other_sweep_${sweep_no}.txt" \
  --out "${output_base}"/fst_results_sweep_${sweep_no}_full
# visualize results
conda activate ani_heatmap
cd "${output_base}"
python3 ../visualize_fst.py --input fst_results_sweep_${sweep_no}_full.weir.fst --out fst_plot_sweep_${sweep_no}_full.png

# Masked alignment
conda activate vcftools
vcftools --vcf pseudo_diploid_masked.vcf \
  --weir-fst-pop "sweep_${sweep_no}.txt" \
  --weir-fst-pop "other_sweep_${sweep_no}.txt" \
  --out "${output_base}"/fst_results_sweep_${sweep_no}_masked
# visualize results
conda activate ani_heatmap
cd "${output_base}"
python3 ../visualize_fst.py --input fst_results_sweep_${sweep_no}_masked.weir.fst --out fst_plot_sweep_${sweep_no}_masked.png

# ------------------------------------------
# Run Bakta annotation against reference
# ------------------------------------------
conda activate bakta
bakta_db list
#mkdir /data/Unit_LMM/selberherr-group/senekowitsch/db/bakta
#bakta_db download --output /data/Unit_LMM/selberherr-group/senekowitsch/db/bakta --type full

# Extract reference genome from the masked alignment
awk '/^>Reference/{print ">1"; next} /^>/{exit} {gsub(/-/,""); print}' "${raw_alignment}" > reference_extracted.fasta

echo "Running bakta for Reference genome"
bakta --db /data/Unit_LMM/selberherr-group/senekowitsch/db/bakta/db \
    --output "${output_base}/bakta" \
    --prefix reference \
    --threads 60 \
    --force \
    reference_extracted.fasta

# ------------------------------------------
# Overlap gene annotations with Fst results
# ------------------------------------------
conda activate ani_heatmap
python3 annotate_fst_genes.py \
  --fst "${output_base}/fst_results_sweep_${sweep_no}_masked.weir.fst" \
  --gff "${output_base}/bakta/reference.gff3" \
  --out-snps "${output_base}/sweep_${sweep_no}_fst_annotated_snps.tsv" \
  --out-genes "${output_base}/sweep_${sweep_no}_fst_gene_summary.tsv" \
  --fst-threshold 0.99 \
  --floor-negative


# check headers of the alignment file (lines with >)
grep ">" "${no_rec_alignment}" | head -n 10
cat "${sweep_labels}" | head -n 10