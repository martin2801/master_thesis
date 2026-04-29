#!/usr/bin/env python3
"""
Post-Filtering ANI Validation and Visualization.

This script subsets a global FastANI identity matrix to include only the 
genomes present in a filtered directory (representatives). It then re-performs 
hierarchical clustering and visualization to verify that the genetic 
diversity of the original dataset is adequately represented in the subset.

Workflow:
    1.  Global Matrix Parsing: Loads the original all-vs-all FastANI matrix.
    2.  Subset Identification: Detects genomes in the filtered directory 
        (via symlinks or files).
    3.  Matrix Extraction: Slices the global matrix to create a symmetric 
        identity matrix of only the selected representatives.
    4.  Clustering: Re-calculates UPGMA linkage for the subset.
    5.  Visualization: 
        - Generates a high-resolution Clustermap (Heatmap) of the subset.
        - Generates a detailed Dendrogram with leaf labels to inspect individual 
          genome placement.

Input:
    - /path/to/fastANI_all_vs_all_output.txt.matrix: The original full matrix.
    - /path/to/08_filtered_data: Directory containing the selected genomes.

Output:
    - results/salmonella_ani_heatmap_filtered.png: Identity heatmap of subset.
    - results/salmonella_ani_dendrogram_filtered.png: Labeled tree of subset.

Dependencies:
    - pandas, numpy, seaborn, matplotlib: Data processing and plotting.
    - scipy.cluster.hierarchy: Linkage and dendrogram generation.
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage

# --- Configuration ---
input_file = "/home/senekowitsch/Thesis/QC/06_pairwise_ANI/results/fastANI_all_vs_all_output.txt.matrix"
filtered_dir = "/home/senekowitsch/Thesis/QC/08_filtered_data"
output_heatmap = "results/salmonella_ani_heatmap_filtered.png"
output_diagnostic = "results/salmonella_ani_dendrogram_filtered.png"

# Optional: set a threshold line on the dendrogram for visual reference.
# Set to None to skip the line.
DISTANCE_THRESHOLD = None  # e.g. 0.4

# --- 1. Parse the full FastANI matrix ---
with open(input_file, 'r') as f:
    lines = f.readlines()

num_genomes = int(lines[0].strip())
names = []
values = []

for line in lines[1:]:
    parts = line.strip().split('\t')
    names.append(parts[0].split('/')[-1])
    if len(parts) > 1:
        values.append([float(x) for x in parts[1:]])
    else:
        values.append([])

# --- 2. Reconstruct full symmetric matrix ---
matrix_data = np.full((num_genomes, num_genomes), 100.0)
for i in range(num_genomes):
    for j, val in enumerate(values[i]):
        matrix_data[i, j] = val
        matrix_data[j, i] = val

full_df = pd.DataFrame(matrix_data, index=names, columns=names)

# --- 3. Identify which genomes are in the filtered set ---
# The filtered_dir contains symlinks named after the genome files.
filtered_genomes = [
    entry for entry in os.listdir(filtered_dir)
    if os.path.islink(os.path.join(filtered_dir, entry))
       or os.path.isfile(os.path.join(filtered_dir, entry))
]

# Check all filtered genomes exist in the full matrix
missing = [g for g in filtered_genomes if g not in full_df.index]
if missing:
    raise ValueError(
        f"These genomes from {filtered_dir} were not found in the ANI matrix:\n"
        + "\n".join(missing)
    )

print(f"Full matrix: {num_genomes} genomes")
print(f"Filtered set: {len(filtered_genomes)} genomes")

# --- 4. Subset the matrix ---
filtered_df = full_df.loc[filtered_genomes, filtered_genomes]

# --- 5. Hierarchical clustering on the filtered matrix ---
row_linkage = linkage(filtered_df, method='average', metric='euclidean')

min_dist = row_linkage[:, 2].min()
max_dist = row_linkage[:, 2].max()
print(f"Filtered distances — Min: {min_dist:.4f}, Max: {max_dist:.4f}")

# --- 6. Clustermap ---
g = sns.clustermap(
    filtered_df,
    row_linkage=row_linkage,
    col_linkage=row_linkage,
    cmap="magma_r",
    vmin=99.4,
    vmax=100,
    xticklabels=False,
    yticklabels=False,
    figsize=(12, 12),
    cbar_kws={'label': 'ANI (%)'}
)

if DISTANCE_THRESHOLD is not None:
    g.ax_col_dendrogram.axhline(
        y=DISTANCE_THRESHOLD, color='red', linestyle='--', linewidth=2
    )
    g.ax_row_dendrogram.axvline(
        x=DISTANCE_THRESHOLD, color='red', linestyle='--', linewidth=2
    )
    title = f"Salmonella ANI — Filtered Genomes (Threshold: {DISTANCE_THRESHOLD})"
else:
    title = "Salmonella ANI — Filtered Genomes"

plt.suptitle(title, y=1.02)
os.makedirs("results", exist_ok=True)
plt.savefig(output_heatmap, dpi=300, bbox_inches='tight')
print(f"Heatmap saved to: {output_heatmap}")

# --- 7. Diagnostic dendrogram ---
from scipy.cluster.hierarchy import dendrogram as _dendrogram

plt.figure(figsize=(max(10, len(filtered_genomes) * 0.15), 7))
_dendrogram(
    row_linkage,
    labels=filtered_df.index.tolist(),
    leaf_rotation=90,
    color_threshold=DISTANCE_THRESHOLD if DISTANCE_THRESHOLD else 0,
)
if DISTANCE_THRESHOLD is not None:
    plt.axhline(
        y=DISTANCE_THRESHOLD, color='r', linestyle='--',
        label=f'Threshold: {DISTANCE_THRESHOLD}'
    )
    plt.legend()
plt.title("Dendrogram — Filtered Genomes")
plt.xlabel("Genomes")
plt.ylabel("Euclidean Distance")
plt.tight_layout()
plt.savefig(output_diagnostic, dpi=300, bbox_inches='tight')
print(f"Dendrogram saved to: {output_diagnostic}")
plt.show()
