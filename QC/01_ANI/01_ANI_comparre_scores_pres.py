# !/usr/bin/env python3

"""
Comparative ANI Profiling and Statistical Significance Testing.

This script evaluates genomic identity across multiple reference lineages 
(Enteritidis, Typhimurium, and Infantis) to determine the closest taxonomic 
affiliation for a dataset of Salmonella isolates. It combines descriptive 
statistics with rigorous non-parametric hypothesis testing.

Statistical Workflow:
    1.  Shapiro-Wilk Test: Evaluates the normality of ANI score distributions 
        for each reference group to guide subsequent test selection.
    2.  Kruskal-Wallis H-test: A non-parametric alternative to ANOVA used to 
        determine if there are statistically significant differences in 
        genomic identity between the reference groups.
    3.  Descriptive Metrics: Calculates mean and median ANI to identify 
        the "best-fit" reference lineage.

Features:
    - Boxplot Visualization: Generates a comparative distribution plot 
      to visualize identity variance and outliers across references.
    - Outlier Filtering: Specifically isolates samples with < 99.0% ANI 
      against the Infantis reference for potential exclusion or re-typing.
    - Integrated Logging: Captures all statistical test outputs and 
      summaries into a timestamped analysis log.

Dependencies:
    - scipy.stats: For Shapiro-Wilk and Kruskal-Wallis implementations.
    - pandas: For parsing multi-column FastANI text outputs.
    - matplotlib: For distribution boxplot generation.
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import sys


ref1 = "/home/senekowitsch/Thesis/QC/01_ANI/results/fastANI_all_vs_enteritidis_reference_output.txt"
ref2 = "/home/senekowitsch/Thesis/QC/01_ANI/results/fastANI_all_vs_typhimurium_reference_output.txt"
ref3 = "/home/senekowitsch/Thesis/QC/01_ANI/results/fastANI_all_vs_infantis_reference_output.txt"

# Update these filenames with your actual file paths
files = [
    ref1,
    ref2, 
    ref3
]

# creat dict that maps filenames to their location (all vs enteritidis, all vs typhimurium, all vs infantis)
comparison_labels = {
    ref1: "All vs Enteritidis",
    ref2: "All vs Typhimurium",
    ref3: "All vs Infantis"
}

# define resutls folder
results_folder = "/home/senekowitsch/Thesis/QC/01_ANI/results/"
output_log_file = results_folder + "ani_comparison_analysis_log.txt"

# --- Redirection Class ---
class Logger(object):
    """Redirects stdout to both the console and a file."""
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # Necessary for compatibility with some python versions
        self.terminal.flush()
        self.log.flush()

# --- Analysis Functions ---
def analyze_ani_files(file_paths):
    data_groups = []
    group_names = []

    # 1. Import the files and extract the 3rd column (index 2)
    for path in file_paths:
        try:
            # Files appear to be tab-separated based on the content structure
            df = pd.read_csv(path, sep='\t', header=None)
            # The 3rd column contains the identity scores (e.g., 98.79)
            scores = df[2].dropna()
            data_groups.append(scores)
            group_names.append(path)
        except Exception as e:
            print(f"Error loading {path}: {e}")

    if len(data_groups) < 3:
        print("Error: Need 3 valid files for comparison.")
        return

    print("--- Normality Test (Shapiro-Wilk) ---")
    # A p-value < 0.05 suggests the data is NOT normally distributed
    all_normal = True
    for name, data in zip(group_names, data_groups):
        stat, p = stats.shapiro(data)
        status = "Normal" if p > 0.05 else "NOT Normal"
        print(f"{name}: p-value = {p:.4f} ({status})")
        if p <= 0.05:
            all_normal = False

    print("\n--- Kruskal-Wallis H-test ---")
    # This is a non-parametric version of ANOVA
    h_stat, p_val = stats.kruskal(*data_groups)
    
    print(f"H-statistic: {h_stat:.4f}")
    print(f"p-value: {p_val:.4e}")

    if p_val < 0.05:
        print("\nResult: There IS a statistically significant difference between the groups.")
    else:
        print("\nResult: There is NO statistically significant difference between the groups.")


# analyze_ani_files(files)


def calculate_stats(file_paths):
    for path in file_paths:
        try:
            # Import the file (assuming tab-separated values)
            df = pd.read_csv(path, sep='\t', header=None)
            
            # The 3rd column (index 2) contains the identity scores
            scores = df[2].dropna()
            
            # Calculate Mean and Median
            mean_val = scores.mean()
            median_val = scores.median()
            
            print(f"File: {path}")
            print(f"  Mean ANI score:   {mean_val:.4f}")
            print(f"  Median ANI score: {median_val:.4f}")
            print("-" * 30)
            
        except Exception as e:
            print(f"Error processing {path}: {e}")


# calculate_stats(files)

# --- Main Execution ---
# This block starts the logging and calls your functions
sys.stdout = Logger(output_log_file)

try:
    analyze_ani_files(files)
    calculate_stats(files)
    print(f"\nAnalysis complete. Results saved to: {output_log_file}")
finally:
    # Always restore original stdout and close the file
    sys.stdout.log.close()
    sys.stdout = sys.stdout.terminal


# --- Plot the scores ---
def plot_ani_boxplot(file_paths, labels_dict):
    data_groups = []
    labels = []
    
    # Load data from each file
    for path in file_paths:
        try:
            df = pd.read_csv(path, sep='\t', header=None)
            scores = df[2].dropna()
            data_groups.append(scores)
            labels.append(labels_dict.get(path, path))
        except Exception as e:
            print(f"Error loading {path}: {e}")
    
    # Create boxplot
    plt.figure(figsize=(10, 6))
    plt.boxplot(data_groups, labels=labels)
    plt.ylabel('ANI Score (%)', fontsize=16)
    plt.yticks(fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=15)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    # Save BEFORE showing, with transparent background
    plt.savefig(results_folder + "ani_boxplot_pres.png", transparent=True, dpi=300)
    plt.show()
    


plot_ani_boxplot(files, comparison_labels)

# filter for lines where ANI is below 99 in Infantis comparison
# Define the input and output file names
input_file = ref3  # This should be the file you want to filter (e.g., all vs infantis)
output_file = "/home/senekowitsch/Thesis/QC/01_ANI/results/fastANI_all_vs_infantis_reference_output_below_99.txt"

def filter_ani_below_threshold(file_path, threshold=99.0):
    with open(file_path, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Split the line by whitespace (tabs or spaces)
            columns = line.split()
            # Check if there are enough columns
            if len(columns) >= 3:
                try:
                    # Convert the 3rd column (index 2) to a float
                    ani_percentage = float(columns[2])
                    
                    # Filter for values below the threshold
                    if ani_percentage < threshold:
                        outfile.write(line)
                except ValueError:
                    # Handle the case where conversion to float fails (e.g., header or malformed line)
                    continue
    print(f"Filtering complete. Results saved to {output_file}")

filter_ani_below_threshold(input_file)

