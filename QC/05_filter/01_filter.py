# !/usr/bin/env python3

"""
Multi-Criteria Quality Control and Genomic Filtering.

This script integrates outputs from CheckM, FastANI, and QUAST to filter a 
genomic dataset based on assembly quality and taxonomic relevance. Samples 
must pass all four user-defined thresholds to be included in the final 
filtered set.

QC Criteria:
    - Completeness (CheckM): >= 90.0%
    - Contamination (CheckM): <= 5.0%
    - N50 (QUAST): >= 150,000 bp
    - Taxonomic Identity (FastANI): >= 99.4% ANI against reference.

Workflow:
    1.  Data Import: Loads CheckM (TSV), FastANI (TXT), and QUAST (TSV) reports.
    2.  Stringent Filtering: Independently identifies samples passing each 
        threshold.
    3.  Set Intersection: Computes the intersection of all four passing sets 
        to ensure multi-dimensional quality.
    4.  Path Normalization: Standardizes genomic filenames across different 
        tool outputs (handling extensions and absolute paths).
    5.  Export: Writes the list of passing sample names and a filtered version 
        of the ANI report to the results directory.

Output:
    - results/filtered_samples.txt: Simple list of passing genome IDs.
    - results/filtered_ANI_results.txt: ANI data for passing samples only.

Dependencies:
    - pandas: Dataframe manipulation and filtering.
    - scipy: (Imported but currently unused in this logic).
    - os: File path and extension handling.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import sys

# set thresholds for filtering
completeness_threshold = 90.0
contamination_threshold = 5.0
ANI_threshold = 99.4
N50_threshold = 150000

# set input and output folders
input_folder = "/home/senekowitsch/Thesis/QC/03_checkM/results/"
results_folder = "/home/senekowitsch/Thesis/QC/05_filter/results/"
output_folder = "/home/senekowitsch/Thesis/QC/05_filter/filtered_data/"

# read in the files for filtering (checkM results, ANI results, N50 results)
ANI_file = "/home/senekowitsch/Thesis/QC/01_ANI/results/fastANI_all_vs_infantis_reference_output.txt"
checkm_file = "/home/senekowitsch/Thesis/QC/03_checkM/results/quality_report.tsv"
N50_file = "/home/senekowitsch/Thesis/QC/04_QUAST/results/transposed_report.tsv"

ANI_df = pd.read_csv(ANI_file, sep="\t", header=None, names=["query", "reference", "ANI", "fragments_mapped", "fragments_total"])
checkm_df = pd.read_csv(checkm_file, sep="\t")
N50_df = pd.read_csv(N50_file, sep="\t")

# filter by completeness and contamination
filtered_checkm_df = checkm_df[(checkm_df["Completeness"] >= completeness_threshold) & (checkm_df["Contamination"] <= contamination_threshold)]
# print(filtered_checkm_df)

# filter by ANI
# clean the ANI query column to get only the filename without the path and extension
filtered_ANI_df = ANI_df[ANI_df["ANI"] >= ANI_threshold]
# print(filtered_ANI_df)

# filter by N50
filtered_N50_df = N50_df[N50_df["N50"] >= N50_threshold]
# print(filtered_N50_df)

# get the list of samples that pass all filters
#filtered_samples = set(filtered_checkm_df["Name"]).intersection(set(filtered_ANI_df["query"])).intersection(set(filtered_N50_df["Assembly"]))

# Clean the ANI 'query' column: remove path AND remove '.fna'
ani_basenames = filtered_ANI_df["query"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

filtered_samples = (
    set(filtered_checkm_df["Name"]) & 
    set(ani_basenames) & 
    set(filtered_N50_df["Assembly"])
)

# write the filtered samples to a file
with open(os.path.join(results_folder, "filtered_samples.txt"), "w") as f:
    for sample in filtered_samples:
        f.write(sample + "\n")

print("Filtered samples:", filtered_samples)
print("Number of filtered samples:", len(filtered_samples))

# filter the original ANI file for the samples that are in filtered_samples.txt
new_ANI_df = ANI_df[ANI_df["query"].apply(lambda x: os.path.splitext(os.path.basename(x))[0] in filtered_samples)]
# write the filtered ANI dataframe to a new file
new_ANI_df.to_csv(os.path.join(results_folder, "filtered_ANI_results.txt"), sep="\t", index=False, header=False)

print("Filtered ANI DataFrame:")
print(new_ANI_df)


