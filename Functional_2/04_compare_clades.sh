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
        --threads "${THREADS}" \
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
        --threads 30 \
        --force \
        "${plasmid}"
done
#rm -rf "/home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/bakta_pESI/*"

# ------------------------------------------
# Sanity checks
# ------------------------------------------
# how many genomes do you expect total (all infantis genomes, not just pESI+)?
ls /home/senekowitsch/Salmonella/AGES/Sweep/00_Data/infantis/all_genomes/*.fna | wc -l

# how many chromosome GFF3s actually got produced?
find /home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/bakta_chromosome -name "*.gff3" | wc -l

# which genomes are missing a GFF3 (failed or still running)?
comm -23 \
  <(ls /home/senekowitsch/Salmonella/AGES/Sweep/00_Data/infantis/all_genomes/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort) \
  <(find /home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/bakta_chromosome -name "*.gff3" | xargs -n1 basename | sed 's/\.gff3$//' | sort)

# ------------------------------------------
# Run panaroo on the infantis chromosomes
# ------------------------------------------
conda activate panaroo
cd "${output_infantis_panaroo}"
mkdir temp
panaroo -i /home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/bakta_chromosome/*/*.gff3 \
    -o "${output_infantis_panaroo}" \
    --clean-mode strict \
    -t "${THREADS}" \
    --remove-invalid-genes \
    -a core \
    --core_threshold 0.95 \
    2>&1 | tee panaroo_run.log

# ------------------------------------------
# Prepare pESI phenotype table for GWAS
# ------------------------------------------
cd "${output_infantis}"

LIST_pESI_ALL="/home/senekowitsch/Salmonella/AGES/Functional/01_pESI/output/infantis/abricate_pESI_genomes_infantis.txt"
LIST_pESI_AGES="/home/senekowitsch/Salmonella/AGES/Functional/03_pESI_tree/output/infantis/clade_distances/AGES_tips_filtered.txt"
LIST_pESI_NCBI="/home/senekowitsch/Salmonella/AGES/Functional/03_pESI_tree/output/infantis/clade_distances/sweep7_tips.txt"

python3 /home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/06_plot_tree_with7_AGES.py \
    --tree "/home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/output/infantis/tree/infantis_clean_tree_noambiguous.treefile" \
    --sweeps "/home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/output/infantis_sweeps_bottomup_clonal_5x.txt" \
    --extra-list $LIST_pESI_ALL \
    --extra-label "pESI_all" \
    --extra-color "#00BCD4" \
    --extra-overwrite \
    --out "${output_base}/infantis_pESI_tree_all.pdf" \
    --labels \
    --dpi 400

python3 /home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/06_plot_tree_with7_AGES.py \
    --tree "/home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/output/infantis/tree/infantis_clean_tree_noambiguous.treefile" \
    --sweeps "/home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/output/infantis_sweeps_bottomup_clonal_5x.txt" \
    --extra-list $LIST_pESI_AGES \
    --extra-label "pESI_ages" \
    --extra-color "#00BCD4" \
    --extra-overwrite \
    --out "${output_base}/infantis_pESI_tree_ages.pdf" \
    --labels \
    --dpi 400

python3 /home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/06_plot_tree_with7_AGES.py \
    --tree "/home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/output/infantis/tree/infantis_clean_tree_noambiguous.treefile" \
    --sweeps "/home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/output/infantis_sweeps_bottomup_clonal_5x.txt" \
    --extra-list $LIST_pESI_NCBI \
    --extra-label "pESI_ncbi" \
    --extra-color "#00BCD4" \
    --extra-overwrite \
    --out "${output_base}/infantis_pESI_tree_ncbi.pdf" \
    --labels \
    --dpi 400

python3 make_pesi_phenotype.py

cd "${base_dir}"

if [ ! -f phylogeny_distance.py ]; then
    wget https://raw.githubusercontent.com/mgalardini/pyseer/master/scripts/phylogeny_distance.py
fi

python phylogeny_distance.py --lmm \
    /home/senekowitsch/Salmonella/AGES/Sweep/01_Sweeps/output/infantis/tree/infantis_clean_tree_noambiguous.treefile \
    > phylogeny_similarity.tsv



# ------------------------------------------
# Run check_sample_overlap.py to check sample overlap between phenotype, gene presence/absence, and similarity matrix
# -----------------------------------------
cd "${base_dir}"
python3 check_sample_overlap.py

# ------------------------------------------
# Run pyseer for GWAS on chromosome and pESI presence/absence phenotype
# ------------------------------------------
conda activate pyseer
cd "${output_infantis_pyseer}"
pyseer --lmm \
    --phenotypes /home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/pESI_phenotype.tsv \
    --pres /home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/panaroo/gene_presence_absence.Rtab \
    --similarity /home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/phylogeny_similarity.tsv \
    --output-patterns patterns.txt \
    --cpu ${THREADS} \
    > pyseer_results.txt

wc -l pyseer_results.txt
sort -g -k4,4 pyseer_results.txt | head -20

# ------------------------------------------
# Run fishers exact test for clades pESI+ vx pESI-
# ------------------------------------------
conda activate phylo_pipeline

cd "${output_infantis}"
python3 "${base_dir}/fisher_gene_test.py" \
    --rtab      "${output_infantis_panaroo}/gene_presence_absence.Rtab" \
    --phenotype "${output_infantis}/pESI_phenotype.tsv" \
    --pyseer    "${output_infantis_pyseer}/pyseer_results.txt" \
    --out       fisher_results.tsv

# ------------------------------------------
# Merge annotations from Panaroo into fisher results
# ------------------------------------------
conda activate phylo_pipeline

cd "${output_infantis}"
python3 /home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/merge_annotations.py \
    --results     /home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/fisher_results.tsv \
    --panaroo-csv /home/senekowitsch/Salmonella/AGES/Functional/04_compare_clades/output/infantis/panaroo/gene_presence_absence.csv \
    --out         fisher_results_annotated.tsv


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
