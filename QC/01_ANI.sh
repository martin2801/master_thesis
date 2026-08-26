# !/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh

# Later in the script, when switching:
conda activate fastANI

# go to folder ANI
cd /home/senekowitsch/Thesis/QC/01_ANI

# Create folder structure
# for reference genomes
mkdir -p ./Genomes/Reference/infantis
mkdir -p ./Genomes/Reference/typhimurium
mkdir -p ./Genomes/Reference/enteritidis

# Path to the genomes files
GENOMES="/home/senekowitsch/Thesis/QC/00_data/"

# Path to the reference folder
REF="/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/"

# Path to the output folder
OUTPUT="/home/senekowitsch/Thesis/QC/01_ANI/results"

# create a file with all genome paths
find $GENOMES -name "*.fna" > genome_paths.txt

# check the line count of the genome paths file
echo "Number of genome paths in file genome_paths.txt:"
wc -l genome_paths.txt

# compare to the number of genome files in the genomes folder
echo "Number of .fna files in the genomes folder:"
find $GENOMES -name "*.fna" | wc -l

## download reference genomes
conda init
conda activate ncbi-datasets

# Infantis
# from Cohen E. et al. 2020. Genome Sequence of an Emerging Salmonella enterica 
# Serovar Infantis and Genomic Comparison with Other
# S. Infantis Strains
# Strain 119944
# GenBank accession GCF_000506925.1
/home/senekowitsch/miniconda3/envs/fastANI/bin/datasets download genome accession GCF_000506925.1 --filename $REF/infantis/GCF_000506925.1.zip

# Typhimurium
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006945.2/
/home/senekowitsch/miniconda3/envs/fastANI/bin/datasets download genome accession GCF_000006945.2 --filename $REF/typhimurium/GCF_000006945.2.zip

# Enteritidis
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009505.1/
/home/senekowitsch/miniconda3/envs/fastANI/bin/datasets download genome accession GCF_000009505.1 --filename $REF/enteritidis/GCF_000009505.1.zip

# unzip the downloaded files, if not already unzipped
if [ ! -d "./Genomes/Reference/typhimurium/GCF_000006945.2" ]; then
    unzip $REF/typhimurium/GCF_000006945.2.zip -d ./Genomes/Reference/typhimurium/GCF_000006945.2
else
    echo "Directory ./Genomes/Reference/typhimurium/GCF_000006945.2 already exists."
fi

if [ ! -d "./Genomes/Reference/enteritidis/GCF_000009505.1" ]; then
    unzip $REF/enteritidis/GCF_000009505.1.zip -d ./Genomes/Reference/enteritidis/GCF_000009505.1
else
    echo "Directory ./Genomes/Reference/enteritidis/GCF_000009505.1 already exists."
fi

if [ ! -d "./Genomes/Reference/infantis/GCF_000506925.1" ]; then
    unzip $REF/infantis/GCF_000506925.1.zip -d ./Genomes/Reference/infantis/GCF_000506925.1
else
    echo "Directory ./Genomes/Reference/infantis/GCF_000506925.1 already exists."
fi


# set paths to reference genomes
REF_INFANTIS="/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/infantis/GCF_000506925.1/ncbi_dataset/data/GCF_000506925.1/GCF_000506925.1_SI119944_genomic.fna"
REF_TYPHIMURIUM="/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/typhimurium/GCF_000006945.2/ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna"
REF_ENTERITIDIS="/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/enteritidis/GCF_000009505.1/ncbi_dataset/data/GCF_000009505.1/GCF_000009505.1_ASM950v1_genomic.fna"

# check that the reference genomes files exist
if [ -f "$REF_INFANTIS" ] && [ -f "$REF_TYPHIMURIUM" ] && [ -f "$REF_ENTERITIDIS" ]; then
    echo "All reference genome files found."
else
    echo "One or more reference genome files not found."
    exit 1
fi

echo "Starting ANI comparisons with reference genomes..."
echo "Activating fastANI environment..."
conda activate fastANI

# do fastANI against infantis reference genome
start=$SECONDS
conda run -n fastANI fastANI --ql genome_paths.txt \
    -r $REF_INFANTIS \
    --output $OUTPUT/fastANI_all_vs_infantis_reference_output.txt \
    --threads 20
end=$SECONDS
echo "Elapsed time: $(($end - $start)) seconds"

# do fastANI against typhimurium reference genome
start=$SECONDS
conda run -n fastANI fastANI --ql genome_paths.txt \
    -r $REF_TYPHIMURIUM \
    --output $OUTPUT/fastANI_all_vs_typhimurium_reference_output.txt \
    --threads 20
end=$SECONDS
echo "Elapsed time: $(($end - $start)) seconds"

# do fastANI against enteritidis reference genome
start=$SECONDS
conda run -n fastANI fastANI --ql genome_paths.txt \
    -r $REF_ENTERITIDIS \
    --output $OUTPUT/fastANI_all_vs_enteritidis_reference_output.txt \
    --threads 20
end=$SECONDS
echo "Elapsed time: $(($end - $start)) seconds"

# sort the output file by the 3rd column (ANI value) in descending order
sort -k3,3nr $OUTPUT/fastANI_all_vs_infantis_reference_output.txt > $OUTPUT/fastANI_all_vs_infantis_reference_output_sorted.txt
sort -k3,3nr $OUTPUT/fastANI_all_vs_typhimurium_reference_output.txt > $OUTPUT/fastANI_all_vs_typhimurium_reference_output_sorted.txt
sort -k3,3nr $OUTPUT/fastANI_all_vs_enteritidis_reference_output.txt > $OUTPUT/fastANI_all_vs_enteritidis_reference_output_sorted.txt

conda activate ani_analysis
python3 01_ANI_compare_scores.py
python3 01_ANI_compare_scores_pres.py
python3 01_ANI_compare_scores_vector.py
conda deactivate

cd /home/senekowitsch/Thesis/QC