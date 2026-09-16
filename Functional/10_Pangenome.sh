#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

# Variables
base_dir='/home/senekowitsch/Thesis/Functional/10_PanGenome'
output_base="${base_dir}/output"
output_bakta="/data/Unit_LMM/selberherr-group/senekowitsch/Thesis/bakta_temp"
output_panaro="/data/Unit_LMM/selberherr-group/senekowitsch/Thesis/panaroo_temp"

GENOMES='/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/genomes'
sweep_labels='/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt'
no_rec_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/core_alignment_noambiguous.aln"
raw_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/snippy_alignment/core.full.aln"
masked_alignment="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/masked.aln"
THREADS=30
reference_genome='/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/infantis/GCF_000506925.1/ncbi_dataset/data/GCF_000506925.1/GCF_000506925.1_SI119944_genomic.fna'
bakta_db='/data/Unit_LMM/selberherr-group/senekowitsch/db/bakta/db'
treefile="/home/senekowitsch/Thesis/Sweeps/04_place_on_tree/output/validate_sweeps/full_tree.treefile"

SWEEPS="/home/senekowitsch/Thesis/Sweeps/05_check_distance/output/sweeps_bottomup_clonal_5x.txt"

# Make output directories
mkdir -p "${base_dir}"
mkdir -p "${output_base}"
mkdir -p "${output_bakta}"
mkdir -p "${output_panaro}"

cd "${base_dir}"

# ------------------------------------------
# Run Bakta for functional annotation
# ------------------------------------------
conda activate bakta
bakta_db list
# Bakta, one genome at a time (or GNU parallel / array job across your cluster)
for fna in "${GENOMES}"/*.fna; do
    genome=$(basename "${fna}" .fna)
    bakta --db "${bakta_db}" \
          --prefix "${genome}" \
          --output "${output_bakta}/${genome}" \
          --threads ${THREADS} \
          --force \
          "${fna}"
done

# Remove all unnecessary files from Bakta output to save space
find . -type f ! \( -name "*.log" -o -name "*.gff3" \) -delete

# ------------------------------------------
# Run Panaroo for pangenome analysis
# ------------------------------------------
panaroo -i "${output_bakta}"/*.gff3 \
        -o "${output_panaro}" \
        --clean-mode strict \
        --remove-invalid-genes \
        -a core \
        -t ${THREADS}

# ------------------------------------------
# build a kinship/similarity matrix from tree
# ------------------------------------------
conda activate pyseer
python phylogeny_distance.py --lmm "${treefile}" > phylogeny_similarity.tsv

