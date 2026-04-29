# !/bin/bash

# remove PopCoGenomeS
#rm -rf /home/senekowitsch/Thesis/PopCoGenomeS/Software/PopCoGenomeS
#conda activate base
#conda remove --name popcogenomes --all
#conda remove --name popcogenomes_r --all
#conda clean --all

# go to PopCoGenomeS folder
cd /home/senekowitsch/Thesis/PopCoGenomeS/Software

# clone the PopCoGenomeS repository if it does not exist
if [ ! -d "PopCoGenomeS" ]; then
    git clone https://github.com/cusoiv/PopCoGenomeS.git
fi

# go to the cloned repository
cd PopCoGenomeS

# set up conda channels and priorities
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible

# setup conda envs for the different steps of the PopCoGenomeS
if conda info --envs | grep -q 'popcogenomes'; then
    echo "Conda environment 'popcogenomes' already exists."
else
    echo "Creating conda environment 'popcogenomes'."
    conda env create -f PopCoGenomeS.yml
fi

if conda info --envs | grep -q 'popcogenomes_r'; then
    echo "Conda environment 'popcogenomes_r' already exists."
else
    echo "Creating conda environment 'popcogenomes_r'."
    conda env create -f PopCoGenomeS_R.yml
fi

# activate both environments to check if they were created successfully
echo "Activating conda environment 'popcogenomes'."
conda activate popcogenomes
echo "Conda environment 'popcogenomes' activated successfully." 

echo "Activating conda environment 'popcogenomes_r'."
conda activate popcogenomes_r
echo "Conda environment 'popcogenomes_r' activated successfully."

conda activate base

echo "Make sure to exchange the "synchain-mugsy" file in miniconda3/envs/popcogenomes/bin/"
