#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Functional/08_dNdS'
output_base="${base_dir}/output"
ALN_full='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/snippy_alignment/core.full.aln'
ALN_masked='/home/senekowitsch/Thesis/Sweeps/08_dNdS/output/masked.aln'
input_bed='/home/senekowitsch/Thesis/Sweeps/01_recombination_checks/output/gubbins_fixed.bed'
input_cfml='/home/senekowitsch/Thesis/Sweeps/01_recombination_checks/output/cfml_fixed.bed'
reference_genome='/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/infantis/GCF_000506925.1/ncbi_dataset/data/GCF_000506925.1/GCF_000506925.1_SI119944_genomic.fna'
sweep_labels='/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt'
treefile="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/validate_sweeps/full_tree.treefile"

MASKING='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/mask_recombination.py'
SNP_DENSITY='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/snp_density_plot.py'

THREADS=30

mkdir -p "${output_base}"

# =============================================================================
# RECOMBINATION MASKING
# =============================================================================
cd "${base_dir}"
conda activate phylo_pipeline
python3 "${MASKING}" \
    -i "${ALN_full}" \
    -b "${input_bed}" \
    -o "${output_base}/masked.aln"

# set variable to the masked alignment for the next step
echo "Recombination masking complete. Masked alignment generated at ${ALN_masked}"

# =============================================================================
# VERIFY MASKING
# =============================================================================
conda activate phylo_pipeline
cd "${base_dir}"
python3 "${SNP_DENSITY}" \
    -r "${ALN_full}" \
    -c "${output_base}/masked.aln" \
    -g "${input_bed}" \
    -f "${input_cfml}" \
    -o "${output_base}/snp_density.png"

# =============================================================================
# Detect genes
# =============================================================================
conda activate prokka
prokka \
  --outdir prokka_reference \
  --prefix ref \
  --genus Salmonella \
  --species enterica \
  --kingdom Bacteria \
  --gcode 11 \
  --cpus ${THREADS} \
  "$reference_genome"

mv prokka_reference "${output_base}"

samtools faidx "$reference_genome"
cat "${reference_genome}.fai"

conda activate phylo_pipeline
python3 slice_genes.py \
  --aln "${output_base}/masked.aln" \
  --gff "${output_base}/prokka_reference/ref.gff" \
  --fai "${reference_genome}.fai" \
  --outdir per_gene_alignments

head "$sweep_labels"
head -c 500 "$treefile"
cut -f2 "$sweep_labels" | tail -n +2 | sort | uniq -c

# =============================================================================
# Prepare RELAX inputs
# =============================================================================
cd "${base_dir}"
conda activate phylo_pipeline
python3 prepare_relax_trees.py \
  --tree "$treefile" \
  --labels "$sweep_labels" \
  --gene-dir per_gene_alignments \
  --outdir relax_inputs \
  --min-test-tips 3

python3 strip_stop_codons.py \
  --in-dir per_gene_alignments \
  --out-dir per_gene_alignments_nostop

python3 extract_matching_alignments.py \
  --manifest relax_inputs/_manifest.tsv \
  --gene-dir per_gene_alignments_nostop \
  --tree-dir relax_inputs \
  --outdir relax_inputs


# =============================================================================
# Run RELAX
# =============================================================================
cd "${base_dir}"
conda activate hyphy

mkdir -p relax_results/sweep_7

hyphy relax \
  --alignment relax_inputs/sweep_7/OHFJIJCL_00001.fasta \
  --tree relax_inputs/sweep_7/OHFJIJCL_00001.nwk \
  --test Foreground \
  --code Universal \
  --output relax_results/sweep_7/OHFJIJCL_00001.RELAX.json

cd ~/Thesis/Functional/08_dNdS
chmod +x run_relax_batch.sh
./run_relax_batch.sh relax_inputs/_manifest.tsv relax_inputs relax_results 30

wc -l relax_results/_parallel_joblog.tsv   # completed jobs so far (minus 1 for header)

# =============================================================================
# Parse results and FDR correction
# =============================================================================
cd "${base_dir}"
conda activate phylo_pipeline

python3 aggregate_relax.py --results-dir relax_results --out relax_summary.csv
# -> no significant results, lets look at dN/dS whole genome wide

# =============================================================================
# dN/dS whole genome wide
# =============================================================================
cd "${base_dir}"
conda activate phylo_pipeline

python3 concatenate_genes.py \
  --gene-dir per_gene_alignments_nostop \
  --out-fasta concatenated_core.fasta \
  --out-manifest concatenated_core.manifest.tsv

mkdir -p relax_inputs_concat/gene_dir
mv concatenated_core.fasta relax_inputs_concat/gene_dir/

python3 prepare_relax_trees.py \
  --tree "$treefile" \
  --labels "$sweep_labels" \
  --gene-dir relax_inputs_concat/gene_dir \
  --outdir relax_inputs_concat \
  --min-test-tips 3

for sweep in relax_inputs_concat/sweep_*/; do
  name=$(basename "$sweep")
  python3 match_alignment_to_tree.py \
    --alignment relax_inputs_concat/gene_dir/concatenated_core.fasta \
    --tree "${sweep}concatenated_core.nwk" \
    --out "${sweep}concatenated_core.fasta"
done

conda activate hyphy
for sweep in relax_inputs_concat/sweep_*/; do
  name=$(basename "$sweep")
  mkdir -p "relax_results_concat/${name}"
  hyphy relax \
    --alignment "${sweep}concatenated_core.fasta" \
    --tree "${sweep}concatenated_core.nwk" \
    --test Foreground \
    --code Universal \
    --output "relax_results_concat/${name}/concatenated_core.RELAX.json"
done

# Diagnose RELAX results
cd "${base_dir}"
conda activate phylo_pipeline
python3 diagnose_relax_power.py --results-dir relax_results_concat --out power_diagnostics.csv

python3 diagnose_relax_power.py --file relax_results_concat/sweep_1/concatenated_core.RELAX.json --out sweep_1_diag.csv
python3 diagnose_relax_power.py --file relax_results_concat/sweep_2/concatenated_core.RELAX.json --out sweep_2_diag.csv
python3 diagnose_relax_power.py --file relax_results_concat/sweep_3/concatenated_core.RELAX.json --out sweep_3_diag.csv
python3 diagnose_relax_power.py --file relax_results_concat/sweep_4/concatenated_core.RELAX.json --out sweep_4_diag.csv
python3 diagnose_relax_power.py --file relax_results_concat/sweep_5/concatenated_core.RELAX.json --out sweep_5_diag.csv
python3 diagnose_relax_power.py --file relax_results_concat/sweep_6/concatenated_core.RELAX.json --out sweep_6_diag.csv
python3 diagnose_relax_power.py --file relax_results_concat/sweep_7/concatenated_core.RELAX.json --out sweep_7_diag.csv

# Parse results and FDR correction
cd "${base_dir}"
conda activate phylo_pipeline
python3 aggregate_relax.py --results-dir relax_results_concat --out relax_summary_concat.csv



