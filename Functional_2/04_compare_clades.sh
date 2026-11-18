#!/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh
THREADS=60

# Variables
base_dir='/home/senekowitsch/Thesis/Functional_2/04_compare_clades'
output_base="${base_dir}/output"
THREADS=60

# Outputs
# Bakta
output_bakta="${output_base}/bakta"
# Panaroo
output_panaroo="${output_base}/panaroo"
# Pyseer
output_pyseer="${output_base}/pyseer"

# Inputs
chromosomes="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/consensus/chromosomes"
pesi_fastas="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/consensus/pESI"
genome_summary="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/consensus/genome_summary.tsv"

# for tree visualization
infantis_tree="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/tree/full_tree_renamed.treefile"
sweeps_infantis="/home/senekowitsch/Thesis/Functional_2/sweep7_genomes.txt"


# ------------------------------------------
# Create directories
# ------------------------------------------
mkdir -p "${base_dir}"
mkdir -p "${output_base}"
mkdir -p "${output_bakta}"
mkdir -p "${output_panaroo}"
mkdir -p "${output_pyseer}"
cd "${base_dir}"

# ------------------------------------------
# Run bakta on chromosomes only
# ------------------------------------------
conda activate bakta

# download the bakta db
#bakta_db list
#mkdir /data/Unit_LMM/selberherr-group/senekowitsch/db/bakta
#bakta_db download --output /data/Unit_LMM/selberherr-group/senekowitsch/db/bakta --type full

for chrom in "${chromosomes}"/*.fasta; do
    genome=$(basename "${chrom}" .fasta)
    bakta --db /data/Unit_LMM/selberherr-group/senekowitsch/db/bakta/db \
        --output "${output_bakta}/bakta_chromosome/${genome}" \
        --prefix "${genome}" \
        --threads 20 \
        --force \
        "${chrom}"
done

#rm -rf "/home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/bakta_chromosome/*"

# ------------------------------------------
# Run bakta on pESI only
# ------------------------------------------
conda activate bakta

# download the bakta db
#bakta_db list
#mkdir /data/Unit_LMM/selberherr-group/senekowitsch/db/bakta
#bakta_db download --output /data/Unit_LMM/selberherr-group/senekowitsch/db/bakta --type full

cd "${output_bakta}"
mkdir -p "${output_bakta}/bakta_pESI"
cd "${output_bakta}/bakta_pESI"

for plasmid in "${pesi_fastas}"/*.fasta; do
    genome=$(basename "${plasmid}" _pESI.fasta)
    bakta --db /data/Unit_LMM/selberherr-group/senekowitsch/db/bakta/db \
        --output "${output_bakta}/bakta_pESI/${genome}" \
        --prefix "${genome}" \
        --threads 20 \
        --force \
        "${plasmid}"
done
#rm -rf "/home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/bakta_pESI/*"

# ------------------------------------------
# Sanity checks
# ------------------------------------------
# how many genomes do you expect total (all infantis genomes, not just pESI+)?
ls /home/senekowitsch/Thesis/Functional_2/00_data/*.fna | wc -l

# how many chromosome GFF3s actually got produced?
find /home/senekowitsch/Thesis/Functional_2/04_compare_clades/output/bakta/bakta_chromosome -name "*.gff3" | wc -l

# which genomes are missing a GFF3 (failed or still running)?
comm -23 \
  <(ls /home/senekowitsch/Thesis/Functional_2/00_data/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort) \
  <(find /home/senekowitsch/Thesis/Functional_2/04_compare_clades/output/bakta/bakta_chromosome -name "*.gff3" | xargs -n1 basename | sed 's/\.gff3$//' | sort)

# ------------------------------------------
# Run panaroo on the infantis chromosomes
# ------------------------------------------
conda activate panaroo
cd "${output_panaroo}"
mkdir temp
panaroo -i /home/senekowitsch/Thesis/Functional_2/04_compare_clades/output/bakta/bakta_chromosome/*/*.gff3 \
    -o "${output_panaroo}" \
    --clean-mode strict \
    -t "${THREADS}" \
    --remove-invalid-genes \
    -a core \
    --core_threshold 0.95 \
    2>&1 | tee panaroo_run.log

# ------------------------------------------
# Prepare pESI phenotype table for GWAS
# ------------------------------------------
cd "${output_base}"

cat /home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/pESI_genomes_infantis_minimap2.txt > pesi_pos.tsv
cat /home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/pESI_partial_genomes_infantis_minimap2.txt >> pesi_pos.tsv

LIST_pESI_ALL="/home/senekowitsch/Thesis/Functional_2/04_compare_clades/output/pesi_pos.tsv"

conda activate phylo_pipeline
python3 /home/senekowitsch/Thesis/Functional_2/01_plasmid_id/plot_tree.py \
    --tree "${infantis_tree}" \
    --highlight "${sweeps_infantis}" "#00BCD4" "Sweep 7"\
    --highlight "${LIST_pESI_ALL}" "#1B7837" "pESI+" \
    --out "${output_base}/pESI_tree_all.pdf" \
    --labels-all 


conda activate phylo_pipeline
cd "${output_base}"
python3 "${base_dir}/make_pesi_phenotype_v2.py" \
    --summary   "${genome_summary}" \
    --bakta-dir "${output_bakta}/bakta_chromosome" \
    --out       "${output_base}/pESI_phenotype.tsv" \
    --pos-list  "${LIST_pESI_ALL}" \
    --partial   pos \
    | tee "${output_base}/pESI_phenotype.log"


cd "${base_dir}"

if [ ! -f phylogeny_distance.py ]; then
    wget https://raw.githubusercontent.com/mgalardini/pyseer/master/scripts/phylogeny_distance.py
fi

python3 phylogeny_distance.py --lmm "${infantis_tree}" > phylogeny_similarity.tsv


# ------------------------------------------
# Run check_sample_overlap.py to check sample overlap between phenotype, gene presence/absence, and similarity matrix
# -----------------------------------------
cd "${base_dir}"
python3 "${base_dir}/check_sample_overlap.py" \
    --phenotype  "${output_base}/pESI_phenotype.tsv" \
    --rtab       "${output_panaroo}/gene_presence_absence.Rtab" \
    --similarity "${base_dir}/phylogeny_similarity.tsv" \
    --show 5

# ------------------------------------------
# Run pyseer for GWAS on chromosome and pESI presence/absence phenotype
# ------------------------------------------
conda activate pyseer
cd "${output_pyseer}"
pyseer --lmm \
    --phenotypes /home/senekowitsch/Thesis/Functional_2/04_compare_clades/output/pESI_phenotype.tsv \
    --pres /home/senekowitsch/Thesis/Functional_2/04_compare_clades/output/panaroo/gene_presence_absence.Rtab \
    --similarity /home/senekowitsch/Thesis/Functional_2/04_compare_clades/phylogeny_similarity.tsv \
    --output-patterns patterns.txt \
    --cpu ${THREADS} \
    > pyseer_results.txt


wc -l pyseer_results.txt
sort -g -k4,4 pyseer_results.txt | head -20

# ------------------------------------------
# Run fishers exact test for clades pESI+ vx pESI-
# ------------------------------------------
conda activate phylo_pipeline

cd "${output_base}"
python3 "${base_dir}/fisher_gene_test.py" \
    --rtab      /home/senekowitsch/Thesis/Functional_2/04_compare_clades/output/panaroo/gene_presence_absence.Rtab \
    --phenotype /home/senekowitsch/Thesis/Functional_2/04_compare_clades/output/pESI_phenotype.tsv \
    --pyseer    "${output_pyseer}/pyseer_results.txt" \
    --out       fisher_results.tsv

# ------------------------------------------
# Merge annotations from Panaroo into fisher results
# ------------------------------------------
conda activate phylo_pipeline

cd "${output_base}"
python3 "${base_dir}/merge_annotations.py" \
    --results     fisher_results.tsv \
    --panaroo-csv /home/senekowitsch/Thesis/Functional_2/04_compare_clades/output/panaroo/gene_presence_absence.csv \
    --out         fisher_results_annotated.tsv

# ------------------------------------------
# Create heatmap
# ------------------------------------------
conda activate phylo_pipeline
cd "${base_dir}"

python3 pangenome_heatmap.py \
    --rtab      "${output_panaroo}/gene_presence_absence.Rtab" \
    --results   "${output_base}/fisher_results_annotated.tsv" \
    --tree      "${infantis_tree}" \
    --midpoint \
    --metadata  /home/senekowitsch/Thesis/Functional_2/02_AMR_finder/metadata.tsv \
    --color-by  pESI pESI_detail group \
    --top-n     500 \
    --outdir    "${output_base}/pangenome_figs"


# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------
# Look at excluded genes becasue they had stop codone
# ------------------------------------------
cd "${output_infantis_panaroo}"
conda activate phylo_pipeline

python3 ./../parse_invalid_genes.py \
    --log panaroo_run.log \
    --out removed_invalid_genes.tsv

python3 ./../lookup_invalid_gene_products.py \
    --in  removed_invalid_genes.tsv \
    --out removed_invalid_genes_annotated.tsv
