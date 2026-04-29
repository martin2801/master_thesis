#!/usr/bin/env python3

"""
Salmonella Genomic Representative Selection via ANI-based Clustering.

This script identifies clusters within a Salmonella dataset using FastANI 
all-vs-all identity matrices and reduces the dataset to a representative 
subset. It prioritizes genetic diversity by selecting the most distant 
members within each cluster (lowest ANI) plus random representatives.

Workflow:
    1.  Parse & Symmetrize: Converts FastANI output into a full identity matrix.
    2.  Clustering: Performs UPGMA (Average Linkage) hierarchical clustering.
    3.  User Interaction: Prompts for a distance threshold to define clusters.
    4.  Visualization: 
        - Generates a Clustermap with cluster-coded sidebars and ANI gradients.
        - Generates a Dendrogram showing the 'cut-off' line for cluster definition.
    5.  Organization: Creates symbolic links for all members grouped by cluster ID.
    6.  Strategic Filtering: 
        - Clusters <= 5 members: Keep all.
        - Clusters > 5 members: Select the "Furthest Pair" (lowest ANI) + 3 random samples.

Output:
    - results/salmonella_ani_heatmap.png: Visual cluster matrix.
    - results/salmonella_ani_dendrogram.png: Tree structure with cut-off line.
    - clusters/: Directory containing subfolders for every detected cluster.
    - filtered_data/: Directory containing symlinks to selected representative genomes.

Dependencies:
    - scipy.cluster.hierarchy: For linkage and tree-cutting logic.
    - seaborn, matplotlib: For genomic heatmaps and diagnostic plots.
    - sklearn: For data structure handling.
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.metrics import silhouette_score
import random
import itertools
from shutil import copy2

# --- Configuration ---
input_file = "/home/senekowitsch/Thesis/QC/06_pairwise_ANI/results/fastANI_all_vs_all_output.txt.matrix"
output_heatmap = "results/salmonella_ani_heatmap.png"
output_diagnostic = "results/salmonella_ani_dendrogram.png"
filtered_dir = "/home/senekowitsch/Thesis/QC/08_filtered_data"
# read from manual input from the comand line after checking the tree
DISTANCE_THRESHOLD = float(input("Enter the distance threshold (e.g., 0.4): "))

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

# --- 3. Hierarchical Clustering ---
# Parameters: Euclidean distance and Average Linkage [cite: 141, 524]
row_linkage = linkage(matrix_df, method='average', metric='euclidean')

# --- 4. Extract Clusters and Map Colors ---
# Assign each genome to a cluster ID based on the distance threshold
cluster_labels = fcluster(row_linkage, t=DISTANCE_THRESHOLD, criterion='distance')
matrix_df['cluster'] = cluster_labels

# Create a color palette for the clusters to visualize them on the heatmap
unique_clusters = np.unique(cluster_labels)
palette = sns.color_palette("hls", len(unique_clusters))
cluster_color_map = dict(zip(unique_clusters, palette))
row_colors = pd.Series(cluster_labels, index=matrix_df.index).map(cluster_color_map)

print(f"Detected {len(unique_clusters)} clusters at distance {DISTANCE_THRESHOLD}")


# Check the range of merge distances
min_dist = row_linkage[:, 2].min()
max_dist = row_linkage[:, 2].max()
median_dist = np.median(row_linkage[:, 2])

print(f"Distances in tree: Min={min_dist:.2f}, Max={max_dist:.2f}, Median={median_dist:.2f}")

# --- 5. Visualization: Clustermap with Cluster Sidebar ---
g = sns.clustermap(
    matrix_df.drop(columns=['cluster']), # Don't cluster the label column
    row_linkage=row_linkage,
    col_linkage=row_linkage,
    row_colors=row_colors,    # Adds the colored cluster bar
    cmap="magma_r",          #
    vmin=99.4,               
    vmax=100,
    xticklabels=False, 
    yticklabels=False,
    figsize=(12, 12),
    cbar_kws={'label': 'ANI (%)'}
)

# Add the threshold line directly to the Clustermap's dendrogram
# Note: In clustermap, the dendrograms are on separate axes
g.ax_col_dendrogram.axhline(y=DISTANCE_THRESHOLD, color='red', linestyle='--', linewidth=2)
g.ax_row_dendrogram.axvline(x=DISTANCE_THRESHOLD, color='red', linestyle='--', linewidth=2)

plt.suptitle(f"Salmonella ANI Clustermap (Threshold: {DISTANCE_THRESHOLD})", y=1.02)
plt.savefig(output_heatmap, dpi=300, bbox_inches='tight')


# --- 6. Diagnostic Plot: Dendrogram with Cut Line ---
plt.figure(figsize=(10, 7))
dendrogram(row_linkage, no_labels=True, color_threshold=DISTANCE_THRESHOLD)
plt.axhline(y=DISTANCE_THRESHOLD, color='r', linestyle='--', label=f'Threshold: {DISTANCE_THRESHOLD}')
plt.title("Hierarchical Clustering Dendrogram")
plt.xlabel("Genomes")
plt.ylabel("Euclidean Distance")
plt.legend()
plt.savefig(output_diagnostic, dpi=300, bbox_inches='tight')
plt.show()

# --- 7. Summary Table ---
cluster_summary = matrix_df['cluster'].value_counts().sort_index()
print("\nCluster Membership Summary:")
print(cluster_summary)


# Example of how to see which samples are in Cluster 1
#print(matrix_df[matrix_df['cluster'] == 1].index.tolist())

# --- 8. Separate Data by Clusters ---
# Create a list of the samples in each cluster
cluster_samples = {}
for cluster_id in sorted(matrix_df['cluster'].unique()):
    cluster_samples[cluster_id] = matrix_df[matrix_df['cluster'] == cluster_id].index.tolist()
print("\nSamples in each cluster:")
for cluster_id, samples in cluster_samples.items():
    print(f"Cluster {cluster_id}: {samples}")

# remove the clusters folder if it exists
if os.path.exists("clusters"):
    import shutil
    shutil.rmtree("clusters")

# check if folder clusters exists, if not create it
if not os.path.exists("clusters"):
    os.makedirs("clusters")

# make a folder for each cluster and create symbolic links to the original files in the clusters folder
for cluster_id, samples in cluster_samples.items():
    cluster_folder = f"clusters/cluster_{cluster_id}"
    if not os.path.exists(cluster_folder):
        os.makedirs(cluster_folder)
    for sample in samples:
        original_file = f"/home/senekowitsch/Thesis/QC/00_data/{sample}"
        link_name = os.path.join(cluster_folder, sample)
        if not os.path.exists(link_name):
            os.symlink(original_file, link_name)

# --- 9. Filter Data ---
# Refresh the filtered_data directory
if os.path.exists(filtered_dir):
    import shutil
    shutil.rmtree(filtered_dir)
os.makedirs(filtered_dir)

print("\n--- Filtering Clusters ---")

for cluster_id, samples in cluster_samples.items():
    selected_samples = []
    
    if len(samples) <= 5:
        # Take everyone if the group is small
        selected_samples = samples
        print(f"Cluster {cluster_id}: Keeping all {len(samples)} samples.")
    else:
        # 1. Find the two most distant samples (Lowest ANI)
        # We use the pre-calculated matrix_df (excluding the 'cluster' column)
        cluster_matrix = matrix_df.loc[samples, samples]
        
        # Find the pair with the minimum ANI value
        # stack() turns the matrix into a Series with MultiIndex (Sample A, Sample B)
        # we filter out 100.0 (the diagonal) to find the true furthest pair
        pairwise_ani = cluster_matrix.stack()
        pairwise_ani = pairwise_ani[pairwise_ani < 99.999] # ignore self-comparisons
        # print(pairwise_ani.head(30))
        if not pairwise_ani.empty:
            furthest_pair = pairwise_ani.idxmin() # returns (sample1, sample2)
            min_ani_value = pairwise_ani.min()      # The actual ANI % value of the furthest pair
            selected_samples.extend(list(furthest_pair))
            # 2. Pick 3 additional random samples from the remainder
            remaining_pool = [s for s in samples if s not in selected_samples]
            random_picks = random.sample(remaining_pool, min(3, len(remaining_pool)))
            selected_samples.extend(random_picks)
            # --- VERIFICATION CHECK ---
            print(f"--- Cluster {cluster_id} Selection Summary ---")
            print(f"  Total samples in cluster: {len(samples)}")
            print(f"  Furthest Pair: {furthest_pair[0]} and {furthest_pair[1]}")
            print(f"  Minimum ANI: {min_ani_value:.4f}%")
            print(f"  Randomly Selected: {', '.join(random_picks)}")
            print("-" * 40)

        else:
            # Fallback if all samples are 100% identical
            selected_samples.extend(samples[:2])
        
        print(f"Cluster {cluster_id}: Filtered {len(samples)} down to {len(selected_samples)} (Furthest pair + random).")

    # 3. Create symbolic links in the filtered_data folder
    for sample in selected_samples:
        original_file = f"/home/senekowitsch/Thesis/QC/00_data/{sample}"
        link_name = os.path.join(filtered_dir, sample)
        
        # Ensure we don't try to link the same file twice (if it was picked twice)
        if not os.path.exists(link_name):
            try:
                os.symlink(original_file, link_name)
            except OSError as e:
                print(f"Error linking {sample}: {e}")

print("-" * 40)
print("-" * 40)
print(f"\nFiltering complete. Representative genomes are in: {filtered_dir}")
print("-" * 40)
print("-" * 40)