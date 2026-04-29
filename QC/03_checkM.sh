# !/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh

# go to QC folder
cd /home/senekowitsch/Thesis/QC/

# activate conda env checkm2
source $(conda info --base)/etc/profile.d/conda.sh
conda activate checkm2

# move to the right folder
cd /home/senekowitsch/Thesis/QC/03_checkM

# set path to the output folder
OUTPUT="/home/senekowitsch/Thesis/QC/03_checkM/results"

# set up folder for database
mkdir -p /home/senekowitsch/Thesis/QC/03_checkM/checkm2_db
DB="/home/senekowitsch/Thesis/QC/03_checkM/checkm2_db"

# download checkm2 database if not already downloaded
if [ -d "$DB/CheckM2_database" ]; then
    echo "CheckM2 database already exists."
else
    echo "Downloading CheckM2 database..."
    checkm2 database --download --path $DB
fi

# make a folder for input data
mkdir -p /home/senekowitsch/Thesis/QC/03_checkM/inputdata

# make a link to the data in folder 00_data and remove any that is below 99% ANI
rm /home/senekowitsch/Thesis/QC/03_checkM/inputdata/*.fna
ln -s /home/senekowitsch/Thesis/QC/00_data/* /home/senekowitsch/Thesis/QC/03_checkM/inputdata/
echo "Removing files with ANI below 99%..."
echo "number of files before removal:"
ls /home/senekowitsch/Thesis/QC/03_checkM/inputdata/*.fna | wc -l
cd inputdata
awk '{print $1}' /home/senekowitsch/Thesis/QC/01_ANI/results/fastANI_all_vs_infantis_reference_output_below_99.txt | xargs -n 1 basename | while read f; do
    rm -v "$f"
done
cd ..
echo "number of files after removal:"
ls /home/senekowitsch/Thesis/QC/03_checkM/inputdata/*.fna | wc -l
echo "Files with ANI below 99% have been removed."

# run checkM2
rm -r /home/senekowitsch/Thesis/QC/03_checkM/results/*
checkm2 predict --threads 10 --input /home/senekowitsch/Thesis/QC/03_checkM/inputdata/ --output-directory $OUTPUT --database_path $DB/CheckM2_database/uniref100.KO.1.dmnd

conda activate ani_analysis
python3 /home/senekowitsch/Thesis/QC/03_checkM/01_checkM_check_results.py
