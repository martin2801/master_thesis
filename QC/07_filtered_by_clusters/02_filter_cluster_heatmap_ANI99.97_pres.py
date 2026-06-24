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
import random
from itertools import combinations
from scipy.stats import mannwhitneyu

# --- Configuration ---
input_file = "/home/senekowitsch/Thesis/QC/06_pairwise_ANI/results/fastANI_all_vs_all_output.txt.matrix"
output_heatmap = "results/salmonella_ani_heatmap_pres.png"
output_diagnostic = "results/salmonella_ani_dendrogram.png"
filtered_dir = "/home/senekowitsch/Thesis/QC/08_filtered_data"
# read from manual input from the comand line after checking the tree
DISTANCE_THRESHOLD = float(input("Enter the distance threshold (e.g., 0.4): "))
ANI_CEILING = 99.98
random.seed(42)

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
g.ax_col_dendrogram.axhline(y=DISTANCE_THRESHOLD, color='red', linestyle='--', linewidth=2)
g.ax_row_dendrogram.axvline(x=DISTANCE_THRESHOLD, color='red', linestyle='--', linewidth=2)

# Enlarge the colorbar (the "legend" in the top-left)
g.ax_cbar.set_ylabel('ANI (%)', fontsize=18)
g.ax_cbar.tick_params(labelsize=14)

g.savefig(output_heatmap, dpi=300, bbox_inches='tight', transparent=True)


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

# save a text file with the cluster membership
with open("results/cluster_membership.txt", "w") as f:
    for cluster_id, samples in cluster_samples.items():
        f.write(f"Cluster {cluster_id}:\n")
        for sample in samples:
            f.write(f"  {sample}\n")
        f.write("\n")
#print("\nSamples in each cluster:")
#for cluster_id, samples in cluster_samples.items():
#    print(f"Cluster {cluster_id}: {samples}")

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

# --- 9. Visualize Cluster within vs. between ANI ---
# Create a boxplot to show the distribution of ANI values within and between clusters
within_ani = []
between_ani = []

# Get list of all genome names
genome_names = matrix_df.index.tolist()

# Iterate through every unique pair of genomes (i, j)
for s1, s2 in combinations(genome_names, 2):
    ani_val = matrix_df.at[s1, s2]
    
    # Check if they belong to the same cluster ID
    if cluster_labels[genome_names.index(s1)] == cluster_labels[genome_names.index(s2)]:
        within_ani.append(ani_val)
    else:
        between_ani.append(ani_val)

# Calculate Statistics
stats = {
    "Within-Cluster": {"Mean": np.mean(within_ani), "Median": np.median(within_ani)},
    "Between-Cluster": {"Mean": np.mean(between_ani), "Median": np.median(between_ani)}
}

print("\n--- ANI Distribution Summary ---")
for group, values in stats.items():
    print(f"{group}: Mean={values['Mean']:.4f}%, Median={values['Median']:.4f}%")

# Check for statistical significance (Mann-Whitney U test)
stat, p_value = mannwhitneyu(within_ani, between_ani, alternative='two-sided')

# Visualization: Boxplot
plt.figure(figsize=(8, 6))
plt.boxplot([within_ani, between_ani], tick_labels=['Within Cluster', 'Between Cluster'])
plt.ylabel('ANI (%)')
plt.title('Genomic Identity Distribution: Within vs. Between Clusters')
plt.grid(axis='y', linestyle='--', alpha=0.7)
# add mean and median values as text annotations
offset = 0.25
plt.text(1 - offset, stats["Within-Cluster"]["Mean"], f"Mean: {stats['Within-Cluster']['Mean']:.4f}%", ha='center', va='bottom', fontsize=9, color='blue')
plt.text(1 - offset, stats["Within-Cluster"]["Median"], f"Median: {stats['Within-Cluster']['Median']:.4f}%", ha='center', va='top', fontsize=9, color='blue')
plt.text(2 + offset, stats["Between-Cluster"]["Mean"], f"Mean: {stats['Between-Cluster']['Mean']:.4f}%", ha='center', va='bottom', fontsize=9, color='orange')
plt.text(2 + offset, stats["Between-Cluster"]["Median"], f"Median: {stats['Between-Cluster']['Median']:.4f}%", ha='center', va='top', fontsize=9, color='orange')
plt.savefig("results/ani_distribution_boxplot.png", dpi=300)
plt.show()

# --- 10. Filter Data (With Logging and 99.97% Threshold) ---
log_file_path = "results/filtering_log.txt"

if os.path.exists(filtered_dir):
    import shutil
    shutil.rmtree(filtered_dir)
os.makedirs(filtered_dir)

print(f"\n--- Filtering Clusters (Threshold: {ANI_CEILING}%) ---")

with open(log_file_path, "w") as log_file:
    log_file.write("Salmonella Filtering Log\n")
    log_file.write(f"Identity Ceiling: {ANI_CEILING}%\n")
    log_file.write("="*40 + "\n")

    for cluster_id, samples in cluster_samples.items():
        # A. Redundancy Filter
        cluster_matrix = matrix_df.loc[samples, samples]
        to_remove = set()
        
        for i, s1 in enumerate(samples):
            if s1 in to_remove: continue
            for s2 in samples[i+1:]:
                if s2 in to_remove: continue
                if cluster_matrix.at[s1, s2] >= ANI_CEILING:
                    to_remove.add(s2)
                    log_file.write(f"CLUSTER {cluster_id}: REMOVED: {s2} (Too similar to {s1}: {cluster_matrix.at[s1, s2]}%)\n")
        
        final_pool = [s for s in samples if s not in to_remove]
        selected_samples = []
        
        # B. Selection Logic
        if len(final_pool) <= 5:
            selected_samples = final_pool
            log_file.write(f"CLUSTER {cluster_id}: Kept all {len(final_pool)} unique samples.\n")
        else:
            clean_matrix = matrix_df.loc[final_pool, final_pool]
            pairwise_ani = clean_matrix.stack()
            pairwise_ani = pairwise_ani[pairwise_ani < 99.999] 
            
            if not pairwise_ani.empty:
                furthest_pair = list(pairwise_ani.idxmin())
                selected_samples.extend(furthest_pair)
                
                remaining = [s for s in final_pool if s not in selected_samples]
                random_picks = random.sample(remaining, min(3, len(remaining)))
                selected_samples.extend(random_picks)
                
                log_file.write(f"CLUSTER {cluster_id}: Selected furthest pair ({furthest_pair}) + {len(random_picks)} random.\n")
            else:
                selected_samples = final_pool[:1]
                log_file.write(f"CLUSTER {cluster_id}: All remaining samples identical, kept one.\n")

        # C. Create Symlinks
        for sample in selected_samples:
            os.symlink(f"/home/senekowitsch/Thesis/QC/00_data/{sample}", os.path.join(filtered_dir, sample))

print(f"Done! Check {log_file_path} to see which genomes were excluded.")



# Final Summary

print("-" * 40)
print(f"\nMann-Whitney U Test: U={stat}, p-value={p_value:.4e} (two-sided)")
if p_value < 0.05:
    print("The difference in ANI distributions is statistically significant.")
else:
    print("The difference in ANI distributions is NOT statistically significant.")
print("-" * 40)
print(f"\nFiltering complete. Representative genomes are in: {filtered_dir}")
print("-" * 40)

print("-" * 40)


