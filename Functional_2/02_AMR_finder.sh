#!/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh
THREADS=60

# Variables
GENOMES="/home/senekowitsch/Thesis/Functional_2/00_data"   # *.fna
EXT=".fna"

base_dir='/home/senekowitsch/Thesis/Functional_2/02_AMR_finder'
output_base="${base_dir}/output"
output_AMR_full="${output_base}/AMRFinder_full"
output_AMR_chromosomes="${output_base}/AMRFinder_chromosomes"
output_AMR_pESI="${output_base}/AMRFinder_pESI"
output_AMR_plasmids="${output_base}/AMRFinder_plasmids"
infantis_tree="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/tree/full_tree_renamed.treefile"

# Create output directory if it doesn't exist
mkdir -p "${output_AMR_full}"
mkdir -p "${output_AMR_chromosomes}"
mkdir -p "${output_AMR_plasmids}"
mkdir -p "${output_AMR_pESI}"


# ------------------------------------------------------------------
# 1. Activate AMRFinder and update the database and run AMRFinder on all genomes
# ------------------------------------------------------------------
conda activate amrfinder
# update the AMRFinder database
amrfinder -u

cd "${output_AMR_full}"
mkdir -p amrfinder
cd amrfinder
for genome in "${GENOMES}"/*${EXT}; do
    sample_name=$(basename "$genome" .fna)
    echo "Processing: $sample_name"
    amrfinder -n "$genome" -O Salmonella --plus -o "${sample_name}_amrfinder.txt"
done

# ------------------------------------------------------------------
# 2. Create metadata file for AMRFinder results
# ------------------------------------------------------------------
cd "${base_dir}"
python3 make_metadata.py \
    --genome-dir /home/senekowitsch/Thesis/Functional_2/00_data \
    --pesi-full /home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/pESI_genomes_infantis_minimap2.txt \
    --pesi-partial /home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/pESI_partial_genomes_infantis_minimap2.txt \
    --sweep7 /home/senekowitsch/Thesis/Functional_2/sweep7_genomes.txt \
    --out metadata.tsv

# ------------------------------------------------------------------
# 3. Run summary statistics on the AMRFinder results
# ------------------------------------------------------------------
cd "${base_dir}"
conda activate phylo_pipeline

python3 amrfinder_summary.py \
    --indir     "${output_AMR_full}/amrfinder/" \
    --tree      "${infantis_tree}" \
    --midpoint \
    --metadata  metadata.tsv \
    --color-by  pESI pESI_detail group   \
    --outdir    "${output_AMR_full}"/amrfinder_summary

# ------------------------------------------------------------------
# 4. Run AMRFinder on all plasmids and chromosomes seperately
# ------------------------------------------------------------------
conda activate amrfinder
# update the AMRFinder database
amrfinder -u

cd "${output_AMR_chromosomes}"
mkdir -p amrfinder
cd amrfinder
for genome in /home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/consensus/chromosomes/*.fasta; do
    sample_name=$(basename "$genome" .fasta)
    echo "Processing: $sample_name"
    amrfinder -n "$genome" -O Salmonella --plus -o "${sample_name}_chromosome_amrfinder.txt"
done

cd "${output_AMR_pESI}"
mkdir -p amrfinder
cd amrfinder
for genome in /home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/consensus/pESI/*.fasta; do
    sample_name=$(basename "$genome" .fasta)
    echo "Processing: $sample_name"
    amrfinder -n "$genome" -O Salmonella --plus -o "${sample_name}_pESI_amrfinder.txt"
done

cd "${output_AMR_plasmids}"
mkdir -p amrfinder
cd amrfinder
for genome in /home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/consensus/plasmid_other/*.fasta; do
    sample_name=$(basename "$genome" .fasta)
    echo "Processing: $sample_name"
    amrfinder -n "$genome" -O Salmonella --plus -o "${sample_name}_plasmid_amrfinder.txt"
done

# ------------------------------------------------------------------
# 5. Run summary statistics on the AMRFinder results for plasmids and chromosomes seperately
# ------------------------------------------------------------------
cd "${base_dir}"
conda activate phylo_pipeline

python3 amrfinder_summary.py \
    --indir     "${output_AMR_chromosomes}/amrfinder/" \
    --suffix  _chromosome_amrfinder \
    --tree      "${infantis_tree}" \
    --midpoint \
    --metadata  metadata.tsv \
    --color-by  pESI pESI_detail group   \
    --outdir    "${output_AMR_chromosomes}"/amrfinder_summary

python3 amrfinder_summary.py \
    --indir     "${output_AMR_pESI}/amrfinder/" \
    --suffix  _pESI_amrfinder \
    --tree      "${infantis_tree}" \
    --midpoint \
    --metadata  metadata.tsv \
    --color-by  pESI pESI_detail group   \
    --outdir    "${output_AMR_pESI}"/amrfinder_summary

python3 amrfinder_summary.py \
    --indir     "${output_AMR_plasmids}/amrfinder/" \
    --suffix  _plasmid_amrfinder \
    --tree      "${infantis_tree}" \
    --midpoint \
    --metadata  metadata.tsv \
    --color-by  pESI pESI_detail group   \
    --outdir    "${output_AMR_plasmids}"/amrfinder_summary


