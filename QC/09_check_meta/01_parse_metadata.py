#!/usr/bin/env python3

import pandas as pd

# Set up file structure
metadata_file = "/home/senekowitsch/Thesis/QC/infantis_clean_genome_collection.csv"
filtered_genomes_file = "/home/senekowitsch/Thesis/QC/07_filtered_by_clusters/results/filtered_genomes.txt"
output_file = "/home/senekowitsch/Thesis/QC/09_check_meta/results/filtered_metadata.csv"
clusters_file = "/home/senekowitsch/Thesis/QC/07_filtered_by_clusters/results/cluster_membership.txt"

# Load metadata and filtered genome list
metadata = pd.read_csv(metadata_file)
with open(filtered_genomes_file, 'r') as f:
    filtered_genomes = set(line.strip() for line in f)

# Change the filtered_genomes to only contain the genome IDs (remove everything after the second underscore)
filtered_genomes = set(genome_id.split('_')[1].split('.')[0] for genome_id in filtered_genomes)
print(filtered_genomes)
# Check how many genomes are in the filtered list
# print(f"Number of genomes in the filtered list: {len(filtered_genomes)}")
# Check how many unique genome IDs are in the filtered list
unique_genome_ids = filtered_genomes
print(f"Number of unique genome IDs in the filtered list: {len(unique_genome_ids)}")


# Transform the metadata column to match the new format
metadata['ShortID'] = metadata['AssemblyAccession'].str.split('_').str[1].str.split('.').str[0]

# Filter metadata using the new ShortID column
filtered_metadata = metadata[metadata['ShortID'].isin(filtered_genomes)]
print(filtered_metadata)

# Save the filtered metadata to a new CSV file
filtered_metadata.to_csv(output_file, index=False)
print(f"Filtered metadata saved to {output_file}")

# Check how many unique GeoLocations are in the filtered metadata
unique_geolocations = filtered_metadata['GeoLocation'].nunique()
print(f"Number of unique GeoLocations in the filtered metadata: {unique_geolocations}")
# Then list them
#print("Unique GeoLocations:")
#print(filtered_metadata['GeoLocation'].unique())

# Check how many unique CollectionDate are in the filtered metadata
unique_collection_dates = filtered_metadata['CollectionDate'].nunique()
print(f"Number of unique CollectionDate in the filtered metadata: {unique_collection_dates}")
# Then list them
#print("Unique CollectionDates:")
#print(filtered_metadata['CollectionDate'].unique())

# Create the combinations of GeoLocation and CollectionDate and check how many unique combinations there are
filtered_metadata['GeoLocation_CollectionDate'] = filtered_metadata['GeoLocation'] + '_' + filtered_metadata['CollectionDate']
unique_combinations = filtered_metadata['GeoLocation_CollectionDate'].nunique()
not_unique_combinations = filtered_metadata['GeoLocation_CollectionDate'].duplicated()
print(f"Number of unique combinations of GeoLocation and CollectionDate in the filtered metadata: {unique_combinations}")
print(f"Number of non-unique combinations of GeoLocation and CollectionDate in the filtered metadata: {not_unique_combinations.sum()}")
# Then list them
print("Unique combinations of GeoLocation and CollectionDate:")
print(filtered_metadata['GeoLocation_CollectionDate'].unique())
print("Non-unique combinations of GeoLocation and CollectionDate:")
print(filtered_metadata[not_unique_combinations]['GeoLocation_CollectionDate'].unique())

# Create the combinations of GeoLocation, CollectionDate and Source and check how many unique combinations there are
filtered_metadata['GeoLocation_CollectionDate_Source'] = filtered_metadata['GeoLocation'] + '_' + filtered_metadata['CollectionDate'] + '_' + filtered_metadata['Source']
unique_combinations_source = filtered_metadata['GeoLocation_CollectionDate_Source'].nunique()
not_unique_combinations_source = filtered_metadata['GeoLocation_CollectionDate_Source'].duplicated()
print(f"Number of unique combinations of GeoLocation, CollectionDate and Source in the filtered metadata: {unique_combinations_source}")
print(f"Number of non-unique combinations of GeoLocation, CollectionDate and Source in the filtered metadata: {not_unique_combinations_source.sum()}")
# Then list them
print("Unique combinations of GeoLocation, CollectionDate and Source:")
print(filtered_metadata['GeoLocation_CollectionDate_Source'].unique())
print("Non-unique combinations of GeoLocation, CollectionDate and Source:")
print(filtered_metadata[not_unique_combinations_source]['GeoLocation_CollectionDate_Source'].unique())

# Create a set of all accessions present in the CSV
csv_accessions = set(metadata['ShortID'])

# Find the IDs that are in your filtered list but MISSING from the CSV
missing_ids = filtered_genomes - csv_accessions

print(f"\nDiscrepancy found: {len(missing_ids)} IDs are in the text file but NOT in the CSV.")
print("Missing IDs:")
print(sorted(list(missing_ids)))


# Find the IDs that are missing in the complete list of genomes
all_genomes_file = "/home/senekowitsch/Thesis/QC/00_data/genome_list.txt"
with open(all_genomes_file, 'r') as f:
    all_genomes = set(line.strip().split('_')[1].split('.')[0] for line in f)

# Compare using the transformed metadata IDs
missing_in_all_genomes = metadata[~metadata['ShortID'].isin(all_genomes)]
print(f"\nDiscrepancy found: {len(missing_in_all_genomes)} IDs are in the metadata file but NOT in the complete list of genomes.")
#print("Missing IDs in complete list of genomes:")
#print(sorted(list(missing_in_all_genomes['AssemblyAccession'])))

missing_in_metadata = all_genomes - set(metadata['ShortID'])
print(f"\nDiscrepancy found: {len(missing_in_metadata)} IDs are in the complete list of genomes but NOT in the metadata file.")
#print("Missing IDs in metadata file:")
#print(sorted(list(missing_in_metadata)))
