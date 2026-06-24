# !/usr/bin/env python3

"""
CheckM Quality Report Statistical Analysis and Visualization.

This script parses the output from CheckM to provide a comprehensive 
statistical summary of genome assembly quality. It evaluates core metrics 
including Completeness, Contamination, and Contig N50 to characterize the 
overall health of the genomic dataset.

Key Features:
    1.  Dual Logging: Utilizes a custom Logger class to simultaneously 
        output results to the terminal and a permanent text log.
    2.  Descriptive Statistics: Calculates mean, median, standard deviation, 
        and range for quality scores.
    3.  Distribution Visualization: Generates a multi-panel histogram 
        showing the frequency distributions of Completeness, Contamination, 
        and N50 values.

Input:
    - /path/to/quality_report.tsv: The tab-separated output from CheckM.

Output:
    - results/checkM_quality_report_analysis_log.txt: Text file containing 
      summary statistics.
    - results/checkM_quality_report_distribution.png: Three-panel histogram 
      plot for visual inspection of dataset quality.

Dependencies:
    - pandas: For robust TSV parsing and statistical calculations.
    - matplotlib: For generating distribution histograms.
    - sys, os: For stream redirection and file path management.
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import sys
import numpy as np

input_file = "/home/senekowitsch/Thesis/QC/03_checkM/results/quality_report.tsv"
results_folder = "/home/senekowitsch/Thesis/QC/03_checkM/results/"
output_log_file = results_folder + "checkM_quality_report_analysis_log.txt"

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
# calculate mean, median and standard deviation for completeness and contamination
def analyze_checkm_quality_report(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
        completeness = df['Completeness']
        contamination = df['Contamination']

        print("CheckM Quality Report Analysis:")
        print(f"Number of samples: {len(df)}")
        print("\nCompleteness:")
        print(f"Mean: {completeness.mean():.2f}")
        print(f"Median: {completeness.median():.2f}")
        print(f"Standard Deviation: {completeness.std():.2f}")
        print(f"Minimum: {completeness.min():.2f}")
        print(f"Maximum: {completeness.max():.2f}")

        print("\nContamination:")
        print(f"Mean: {contamination.mean():.2f}")
        print(f"Median: {contamination.median():.2f}")
        print(f"Standard Deviation: {contamination.std():.2f}")
        print(f"Minimum: {contamination.min():.2f}")
        print(f"Maximum: {contamination.max():.2f}")

    except Exception as e:
        print(f"Error analyzing CheckM quality report: {e}")


# plot distribution of completeness and contamination
def plot_checkm_quality_report(file_path, output_folder):
    try:
        df = pd.read_csv(file_path, sep='\t')
        completeness = df['Completeness']
        contamination = df['Contamination']
        N50 = df['Contig_N50']

        plt.figure(figsize=(12, 5))

        plt.subplot(2, 2, 1)
        plt.hist(completeness, bins=20, color='blue', alpha=0.7)
        plt.axvline(90, color='red', linestyle='--', linewidth=2, label='Threshold (90%)')
        plt.title('Distribution of Completeness', fontsize=16)
        plt.xlabel('Completeness (%)', fontsize=14)
        plt.ylabel('Frequency', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=12)

        plt.subplot(2, 2, 2)
        plt.hist(contamination, bins=20, color='red', alpha=0.7)
        plt.axvline(5, color='red', linestyle='--', linewidth=2, label='Threshold (5%)')
        plt.title('Distribution of Contamination', fontsize=16)
        plt.xlabel('Contamination (%)', fontsize=14)
        plt.ylabel('Frequency', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=12)

        plt.subplot(2, 2, 3)
        log_bins = np.logspace(np.log10(N50.min()), np.log10(N50.max()), 50)
        plt.hist(N50, bins=log_bins, color='green', alpha=0.7)
        plt.xscale('log')
        plt.axvline(150000, color='red', linestyle='--', linewidth=2, label='Threshold (150,000)')
        plt.title('Distribution of Contig N50', fontsize=16)
        plt.xlabel('Contig N50 (log scale)', fontsize=14)
        plt.ylabel('Frequency', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=12)

        plt.tight_layout()
        output_path = f"{output_folder}checkM_quality_report_distribution_presi.png"
        plt.savefig(output_path, transparent=True, dpi=300)
        print(f"Distribution plot saved to: {output_path}")

    except Exception as e:
        print(f"Error plotting CheckM quality report: {e}")

# --- Main Execution ---
if __name__ == "__main__":
    # Redirect stdout to both console and log file
    sys.stdout = Logger(output_log_file)

    # Analyze the CheckM quality report
    analyze_checkm_quality_report(input_file)

    # Plot the distribution of completeness and contamination
    plot_checkm_quality_report(input_file, results_folder)
    print(f"\nAnalysis complete. Results saved to: {output_log_file}")
    