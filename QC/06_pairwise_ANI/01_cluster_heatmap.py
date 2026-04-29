#!/usr/bin/env python3

"""
Combined Hierarchical Clustering and Network-Based ANI Analysis.

This script processes an all-vs-all Average Nucleotide Identity (ANI) matrix 
to characterize genomic relationships. It uses two complementary approaches:
1.  Hierarchical Clustering: Generates a Clustermap using UPGMA (Average 
    Linkage) to visualize global population structure.
2.  Network Component Analysis: Constructs undirected graphs at increasing 
    identity thresholds (99.9% to 99.995%) to identify connected components 
    (clusters) and assess dataset fragmentation.

Methodology:
    - Clustering: Euclidean distance on identity values with Average Linkage.
    - Network: Nodes represent genomes; edges are formed if pairwise ANI 
      meets or exceeds the specified threshold.

Input:
    - results/fastANI_all_vs_all_output.txt.matrix: FastANI matrix file.

Output:
    - results/salmonella_ani_heatmap.png: High-resolution clustered heatmap.
    - Console Output: Detailed breakdown of cluster counts and the size of 
      the largest connected component per threshold.

Dependencies:
    - pandas, numpy, seaborn, matplotlib: Data processing and visualization.
    - scipy.cluster.hierarchy: Hierarchical clustering linkage.
    - networkx: Graph construction and connected component analysis.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage

# --- Configuration ---
input_file = "results/fastANI_all_vs_all_output.txt.matrix"
output_heatmap = "results/salmonella_ani_heatmap.png"

# --- 1. Custom Parser for FastANI Matrix ---
# We cannot use pd.read_csv because of the header line and triangular structure.
with open(input_file, 'r') as f:
    lines = f.readlines()

# Line 1 is the number of genomes [cite: 520]
num_genomes = int(lines[0].strip())
names = []
values = []

# Lines 2+ contain the filename followed by ANI values
for line in lines[1:]:
    parts = line.strip().split('\t')
    # Clean name: remove path and keep filename
    names.append(parts[0].split('/')[-1]) 
    if len(parts) > 1:
        values.append([float(x) for x in parts[1:]])
    else:
        values.append([])

# --- 2. Reconstruct Full Symmetric Matrix ---
# Initialize with 100% identity on the diagonal
matrix_data = np.full((num_genomes, num_genomes), 100.0)

for i in range(num_genomes):
    for j, val in enumerate(values[i]):
        matrix_data[i, j] = val
        matrix_data[j, i] = val # Ensure the matrix is symmetric

# Convert to DataFrame for easier handling in Seaborn
matrix_df = pd.DataFrame(matrix_data, index=names, columns=names)
print(matrix_df)

# --- 3. Hierarchical Clustering (Paper Methodology) ---
# Parameters: Euclidean distance and Average Linkage [cite: 141, 524]
row_linkage = linkage(matrix_df, method='average', metric='euclidean')

# --- 4. Generate the Clustermap ---
# Replicating the visual style of Figure 1 [cite: 138, 140]
plt.figure(figsize=(12, 10))
g = sns.clustermap(
    matrix_df,
    row_linkage=row_linkage,
    col_linkage=row_linkage, # Symmetric clustering
    cmap="magma_r",          # Darker colors = higher identity
    vmin=99.4,               # Adjusted for high-ANI data
    vmax=100,
    annot=False,              # Set to True to see values in the cells
    xticklabels=False, 
    yticklabels=False,
    figsize=(10, 10),
    cbar_kws={'label': 'ANI (%)'}
)

# Title and Layout
# g.ax_heatmap.set_title("Salmonella ANI Clustering\n(Euclidean Dist + Average Linkage)")
plt.savefig(output_heatmap, dpi=300, bbox_inches='tight')
print(f"Heatmap saved to: {output_heatmap}")
plt.show()


import networkx as nx

for thresh in [99.9, 99.95, 99.97, 99.99, 99.995]:
    G = nx.Graph()
    G.add_nodes_from(matrix_df.index)

    for i in range(len(matrix_df)):
        for j in range(i + 1, len(matrix_df)):
            if matrix_df.iat[i, j] >= thresh:
                G.add_edge(matrix_df.index[i], matrix_df.index[j])

    n_clusters = nx.number_connected_components(G)
    biggest = max(len(c) for c in nx.connected_components(G))

    print(f"{thresh}% -> clusters: {n_clusters}, largest cluster: {biggest}")