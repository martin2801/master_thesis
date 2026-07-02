##set working directory
setwd("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/PopCoGenomeS/Output_500_0.5")

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

# Load data
data_raw <- read.delim("salmonella_outgroups_rerun.length_bias_500.txt", sep = "\t", header = TRUE)
data_filtered <- read.delim("salmonella_outgroups_rerun.length_bias.filtered.txt", sep = "\t", header = TRUE)
metadata <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/infantis_clean_genome_collection_with_bioproject.csv", sep = ",", header = TRUE)
metadata_snpclusters <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/PopCoGenomeS/Output_500_0.5/isolates_infantis_env_filtered.tsv", sep = ",", header = TRUE)

genomes_in_analysis <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/Recombination/filtered_samples.txt", header = FALSE)
genomes_in_analysis$id <- gsub("\\.fna$", "", genomes_in_analysis$V1)
genomes_in_analysis <- genomes_in_analysis %>%
  mutate(id = str_extract(id, "(?<=GCA_)\\d+"))
genomes_in_analysis <- genomes_in_analysis %>%
  mutate(id = str_remove(id, "^0+"))


# Read in clusters_table.tsv
# Contains genomes after filtering plus to which cluster they belong (from QC 07)
clusters <- read.table("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/Clock/Trees/clusters_table.tsv", 
                       header = FALSE, 
                       sep = "\t", 
                       col.names = c("Genome", "Cluster"),
                       stringsAsFactors = FALSE)
clusters$Cluster <- as.factor(clusters$Cluster)
clusters$Genome <- gsub("\\.fna$", "", clusters$Genome)
clusters_short <- clusters %>%
  mutate(Genome = str_extract(Genome, "(?<=GCA_)\\d+"))
clusters_short <- clusters_short %>%
  mutate(id = str_remove(Genome, "^0+"))

# --- Genomes in each sweep (hard-coded, 'out_' prefix removed) ---

genomes_sweep1 <- c(
  "32036575", "32037335", "31832435"
)
# n = 3

genomes_sweep2 <- c(
  "45808965", "45808925", "45808885", "45808905"
)
# n = 4

genomes_sweep3 <- c(
  "11452785", "8046645", "11453035", "8091485", "11452675", "11452375",
  "11452255"
)
# n = 7

genomes_sweep4 <- c(
  "20671385", "20671465", "20671355", "20671165"
)
# n = 4

genomes_sweep5 <- c(
  "8927125", "8928825", "8026935"
)
# n = 3

genomes_sweep6 <- c(
  "8637495", "8484365", "26733835"
)
# n = 3

genomes_sweep7 <- c(
  "45176625", "9318515", "31904705", "29413815", "27606475",
  "11168375", "11595885", "11588445", "32316375", "25671015", "28916195",
  "28916085", "31765335", "24242935", "24000255", "32498505", "32498225",
  "31895175", "14051005", "24910175", "8713855", "8813075", "22363375",
  "11740395", "32773835", "17181415", "32858445", "16626965", "41292945",
  "7823845", "8778775", "14599635", "10805445", "7794975", "17337285",
  "10134355", "9164925", "19205865", "27171365", "9536415", "18277005",
  "18276945", "44464205", "23514945", "32498565", "14768535", "14075525",
  "33637285", "31725495", "23416215", "44468645", "31765255", "10895095",
  "33877195", "7792295", "21308675", "9207225", "27606455", "16459765",
  "8559715", "7824425", "17689515", "15165115", "26044135", "14600175",
  "20607955", "19562115", "8765535", "21761305", "10141235", "9167145",
  "20516485", "11519185", "29769425", "31725395", "20081485", "16626985",
  "7883985", "33632345", "23837735", "22905815", "32643285", "21180765",
  "7792575", "33015095", "8222005", "15376085", "8461385", "10142775",
  "10142075", "9046505", "31904685", "9474185", "26252495", "11518205",
  "9342185", "16058255", "9523365", "9026405", "31059355", "14550555",
  "20964875", "11156655", "8666235", "9579385", "16976055", "9124805",
  "33688975", "33015535", "31895395", "32858425", "20653895", "18276585",
  "11528025", "8715535", "22617875", "9401885", "11169795", "8461765",
  "20808085", "46616815", "32450575", "19357875", "32323275", "10902475",
  "25806255", "19267755", "14599955", "10833965", "9206385", "8649995",
  "8769955", "10142695", "31895675", "15165015", "12310975", "40741165",
  "8989965", "20487425", "7918635", "29323955", "19233445", "16058695",
  "9167845", "26370835", "11582185", "21007735", "14581275", "8444725",
  "9476285", "14053045", "9123725", "9365255", "24851145", "8770075",
  "24001045", "8714695", "9131505", "33015275", "31059995", "32498545",
  "31770105", "31184075", "31027665", "26044315", "31722655", "24688825",
  "17659675", "31059475", "11600345", "22330065", "11592475", "11551395",
  "33757895", "32759455", "11464235", "16349885", "22331385", "15176205",
  "8547575", "7795135", "32858385", "33687925", "25597515", "15165415",
  "9217935", "10892735", "14045005", "10133815", "14581255", "33014935",
  "33631805", "31895515", "25671005", "28547265", "20658025", "19562235",
  "15169455", "7736995", "9317375", "8461165", "20864745", "16121955",
  "14738275", "11163695", "10133415"
)
# n = 206


metadata <- metadata %>%
  mutate(
    id = str_extract(AssemblyAccession, "(?<=_)[0-9]+") %>%
      str_remove("^0+")
  )


# Filter by exact genome-ID membership using the hard-coded sweep vectors
sweep_metadata1 <- metadata %>% filter(id %in% genomes_sweep1)
sweep_metadata2 <- metadata %>% filter(id %in% genomes_sweep2)
sweep_metadata3 <- metadata %>% filter(id %in% genomes_sweep3)
sweep_metadata4 <- metadata %>% filter(id %in% genomes_sweep4)
sweep_metadata5 <- metadata %>% filter(id %in% genomes_sweep5)
sweep_metadata6 <- metadata %>% filter(id %in% genomes_sweep6)
sweep_metadata7 <- metadata %>% filter(id %in% genomes_sweep7)

sweep_clusters1 <- clusters_short %>% filter(id %in% genomes_sweep1)
sweep_clusters2 <- clusters_short %>% filter(id %in% genomes_sweep2)
sweep_clusters3 <- clusters_short %>% filter(id %in% genomes_sweep3)
sweep_clusters4 <- clusters_short %>% filter(id %in% genomes_sweep4)
sweep_clusters5 <- clusters_short %>% filter(id %in% genomes_sweep5)
sweep_clusters6 <- clusters_short %>% filter(id %in% genomes_sweep6)
sweep_clusters7 <- clusters_short %>% filter(id %in% genomes_sweep7)

  
# Define regions in sweeps
regions1 <- c(	
  "usa:ga")
regions_pattern1 <- paste(regions1, collapse = "|")

regions2 <- c(	
  "china:nan tong")
regions_pattern2 <- paste(regions2, collapse = "|")

regions3 <- c(	
  "	usa:co")
regions_pattern3 <- paste(regions3, collapse = "|")

regions4 <- c(	
  "canada:guelph")
regions_pattern4 <- paste(regions4, collapse = "|")

regions5 <- c(	
  "usa:oh", "usa:oh")
regions_pattern5 <- paste(regions5, collapse = "|")

regions6 <- c(	
  "usa:al", "usa:oh", "usa:tn")
regions_pattern6 <- paste(regions6, collapse = "|")

regions7 <- c(	
  "canada:ontario", "saudi arabia", "usa:al", "usa:ar", "usa:ca", 
  "usa:de", "usa:ga", "usa:ia", "usa:il", "usa:in", "usa:ky", "usa:la",
  "usa:md", "usa:mi", "usa:mn", "usa:mo", "usa:ms", "usa:nc", "usa:ne",
  "usa:ny", "usa:nj", "usa:oh", "usa:ok", "usa:pa", "usa:sc", "usa:sd", 
  "usa:tn", "usa:tx", "usa:va", "usa:wv")
regions_pattern7 <- paste(regions7, collapse = "|")

# Define dates in sweeps
dates1 <- c(
  "2022-10-12", "2022-09-08", "2022-09-08")
dates_pattern1 <- paste(dates1, collapse = "|")

dates2 <- c(
  "2019-09-26", "2019-08-20")
dates_pattern2 <- paste(dates2, collapse = "|")

dates3 <- c(
  "2007-01-09")
dates_pattern3 <- paste(dates3, collapse = "|")

dates4 <- c(
  "2012-07-11", "2012-07-12", "2012-08-08")
dates_pattern4 <- paste(dates4, collapse = "|")

dates5 <- c(
  "2012-07-26", "2012")
dates_pattern5 <- paste(dates5, collapse = "|")

dates6 <- c(
  "2022", "2015", "2017")
dates_pattern6 <- paste(dates6, collapse = "|")

# define clusters in sweep
clusters1 <- c(51)
clusters2 <- c(32)
clusters3 <- c(38, 39, 43, 45)
clusters4 <- c(38, 40)
clusters5 <- c(3)
clusters6 <- c(35)
clusters7 <- c(15, 16, 17, 18, 31)
clusters_pattern1 <- paste(clusters1, collapse = "|")
clusters_pattern2 <- paste(clusters2, collapse = "|")
clusters_pattern3 <- paste(clusters3, collapse = "|")
clusters_pattern4 <- paste(clusters4, collapse = "|")
clusters_pattern5 <- paste(clusters5, collapse = "|")
clusters_pattern6 <- paste(clusters6, collapse = "|")
clusters_pattern7 <- paste(clusters7, collapse = "|")

# Filter metadata by sweep parameters
metadata_by_region1 <- metadata %>%
  filter(str_detect(GeoLocation, regions_pattern1))
metadata_by_region2 <- metadata %>%
  filter(str_detect(GeoLocation, regions_pattern2))
metadata_by_region3 <- metadata %>%
  filter(str_detect(GeoLocation, regions_pattern3))
metadata_by_region4 <- metadata %>%
  filter(str_detect(GeoLocation, regions_pattern4))
metadata_by_region5 <- metadata %>%
  filter(str_detect(GeoLocation, regions_pattern5))
metadata_by_region6 <- metadata %>%
  filter(str_detect(GeoLocation, regions_pattern6))
metadata_by_region7 <- metadata %>%
  filter(str_detect(GeoLocation, regions_pattern7))



metadata_by_date1 <- metadata %>%
  filter(str_detect(CollectionDate, dates_pattern1))
metadata_by_date2 <- metadata %>%
  filter(str_detect(CollectionDate, dates_pattern2))
metadata_by_date3 <- metadata %>%
  filter(str_detect(CollectionDate, dates_pattern3))
metadata_by_date4 <- metadata %>%
  filter(str_detect(CollectionDate, dates_pattern4))
metadata_by_date5 <- metadata %>%
  filter(str_detect(CollectionDate, dates_pattern5))
metadata_by_date6 <- metadata %>%
  filter(str_detect(CollectionDate, dates_pattern6))

# Filter genomes in clusters by clusters in sweeps
genomes_in_sweep_clusters1 <- clusters_short %>%
  filter(str_detect(Cluster, clusters_pattern1))
genomes_in_sweep_clusters2 <- clusters_short %>%
  filter(str_detect(Cluster, clusters_pattern2))
genomes_in_sweep_clusters3 <- clusters_short %>%
  filter(str_detect(Cluster, clusters_pattern3))
genomes_in_sweep_clusters4 <- clusters_short %>%
  filter(str_detect(Cluster, clusters_pattern4))
genomes_in_sweep_clusters5 <- clusters_short %>%
  filter(str_detect(Cluster, clusters_pattern5))
genomes_in_sweep_clusters6 <- clusters_short %>%
  filter(str_detect(Cluster, clusters_pattern6))
genomes_in_sweep_clusters7 <- clusters_short %>%
  filter(str_detect(Cluster, clusters_pattern7))

genomes_in_sweep_clustersanddataset1 <- genomes_in_sweep_clusters1 %>%
  filter(id %in% genomes_in_analysis$id)
genomes_in_sweep_clustersanddataset2 <- genomes_in_sweep_clusters2 %>%
  filter(id %in% genomes_in_analysis$id)
genomes_in_sweep_clustersanddataset3 <- genomes_in_sweep_clusters3 %>%
  filter(id %in% genomes_in_analysis$id)
genomes_in_sweep_clustersanddataset4 <- genomes_in_sweep_clusters4 %>%
  filter(id %in% genomes_in_analysis$id)
genomes_in_sweep_clustersanddataset5 <- genomes_in_sweep_clusters5 %>%
  filter(id %in% genomes_in_analysis$id)
genomes_in_sweep_clustersanddataset6 <- genomes_in_sweep_clusters6 %>%
  filter(id %in% genomes_in_analysis$id)
genomes_in_sweep_clustersanddataset7 <- genomes_in_sweep_clusters7 %>%
  filter(id %in% genomes_in_analysis$id)

genomes_not_in_sweep_clustersbutdataset1 <- genomes_in_sweep_clustersanddataset1 %>%
  filter(!id %in% genomes_sweep1)
genomes_not_in_sweep_clustersbutdataset2 <- genomes_in_sweep_clustersanddataset2 %>%
  filter(!id %in% genomes_sweep2)
genomes_not_in_sweep_clustersbutdataset3 <- genomes_in_sweep_clustersanddataset3 %>%
  filter(!id %in% genomes_sweep3)
genomes_not_in_sweep_clustersbutdataset4 <- genomes_in_sweep_clustersanddataset4 %>%
  filter(!id %in% genomes_sweep4)
genomes_not_in_sweep_clustersbutdataset5 <- genomes_in_sweep_clustersanddataset5 %>%
  filter(!id %in% genomes_sweep5)
genomes_not_in_sweep_clustersbutdataset6 <- genomes_in_sweep_clustersanddataset6 %>%
  filter(!id %in% genomes_sweep6)
genomes_not_in_sweep_clustersbutdataset7 <- genomes_in_sweep_clustersanddataset7 %>%
  filter(!id %in% genomes_sweep7)






library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
# scales comes with ggplot2; ggalluvial is optional (Part 2, flow diagram)

# ---- 0.  Denominator choice -----------------------------------------
# TRUE  : "total cluster" = cluster genomes that are in the recombination
#         dataset (i.e. genomes that *could* have been in a sweep).
#         -> matches your *...anddataset* tables.
# FALSE : "total cluster" = the full QC07 ANI cluster (all genomes).
restrict_to_dataset <- TRUE

# ---- 1.  One tidy lookup:  genome id -> sweep number ----------------
sweep_lookup <- bind_rows(
  data.frame(id = genomes_sweep1, sweep = 1L),
  data.frame(id = genomes_sweep2, sweep = 2L),
  data.frame(id = genomes_sweep3, sweep = 3L),
  data.frame(id = genomes_sweep4, sweep = 4L),
  data.frame(id = genomes_sweep5, sweep = 5L),
  data.frame(id = genomes_sweep6, sweep = 6L),
  data.frame(id = genomes_sweep7, sweep = 7L)
)
# one genome -> one sweep, so there should be no duplicate ids
stopifnot(!any(duplicated(sweep_lookup$id)))

# ---- 2.  Master table:  genome | cluster | sweep --------------------
genome_master <- clusters_short %>%
  select(id, Cluster) %>%
  mutate(Cluster = as.character(Cluster))

if (restrict_to_dataset) {
  genome_master <- genome_master %>% filter(id %in% genomes_in_analysis$id)
}

genome_master <- genome_master %>%
  left_join(sweep_lookup, by = "id") %>%
  mutate(sweep_lab = ifelse(is.na(sweep), "none", paste0("Sweep ", sweep)))

# warn about any sweep genome that isn't in the cluster/dataset master
missing_ids <- setdiff(sweep_lookup$id, genome_master$id)
if (length(missing_ids))
  message(length(missing_ids),
          " sweep genome(s) not found in clusters/dataset master: ",
          paste(missing_ids, collapse = ", "))

# clusters that are involved in at least one sweep
clusters_per_sweep <- list(
  `1` = clusters1, `2` = clusters2, `3` = clusters3, `4` = clusters4,
  `5` = clusters5, `6` = clusters6, `7` = clusters7
)
sweep_clusters_all <- sort(unique(unlist(clusters_per_sweep)))
cl_levels <- as.character(sweep_clusters_all)

# consistent colours for sweeps (+ grey for "none")
sweep_cols <- setNames(scales::hue_pal()(7), paste0("Sweep ", 1:7))
sweep_cols <- c(sweep_cols, none = "grey85")
sweep_levels <- c(paste0("Sweep ", 1:7), "none")


# =====================================================================
#  PART 1 :  per-sweep pie  ("how much of the cluster is in the sweep")
# =====================================================================
pie_data <- bind_rows(lapply(names(clusters_per_sweep), function(s) {
  cls <- as.character(clusters_per_sweep[[s]])
  sub <- genome_master %>% filter(Cluster %in% cls)
  in_sweep <- sum(sub$sweep == as.integer(s), na.rm = TRUE)
  data.frame(
    sweep    = paste0("Sweep ", s),
    category = c("In this sweep", "Rest of cluster(s)"),
    count    = c(in_sweep, nrow(sub) - in_sweep)
  )
}))

pie_data <- pie_data %>%
  group_by(sweep) %>%
  mutate(pct   = count / sum(count),
         label = ifelse(count > 0,
                        paste0(count, "\n(", round(100 * pct), "%)"), "")) %>%
  ungroup()

# Easier-to-compare alternative: faceted 100% stacked bar (same numbers)
p_pies_bar <- ggplot(pie_data, aes(x = sweep, y = count, fill = category)) +
  geom_col(position = "fill", colour = "white") +
  geom_text(aes(label = ifelse(count > 0, count, "")),
            position = position_fill(vjust = 0.5), size = 3) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("In this sweep"      = "#2c7fb8",
                               "Rest of cluster(s)" = "#d9d9d9")) +
  labs(title = "Proportion of each sweep's cluster(s) captured by the sweep",
       x = NULL, y = "Share of cluster genomes", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_pies_bar)


# =====================================================================
#  PART 2 :  which ANI clusters go into which sweeps, in what proportion
# =====================================================================
cluster_composition <- genome_master %>%
  filter(Cluster %in% cl_levels) %>%
  count(Cluster, sweep_lab, name = "n") %>%
  group_by(Cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(Cluster   = factor(Cluster, levels = cl_levels),
         sweep_lab = factor(sweep_lab, levels = sweep_levels))

# 2a. per-cluster composition (proportions) -- the main "where does each
#     cluster go" plot.  Cluster 38 will appear split across sweeps 3 & 4.
p_cluster_prop <- ggplot(cluster_composition,
                         aes(x = Cluster, y = prop, fill = sweep_lab)) +
  geom_col(colour = "white") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = sweep_cols, drop = FALSE) +
  labs(title = "Composition of each sweep-associated ANI cluster",
       subtitle = "Fraction of each cluster's genomes falling into each sweep",
       x = "ANI cluster", y = "Share of cluster genomes", fill = NULL) +
  theme_minimal()
print(p_cluster_prop)

# 2b. same but absolute counts (shows cluster sizes too)
p_cluster_count <- ggplot(cluster_composition,
                          aes(x = Cluster, y = n, fill = sweep_lab)) +
  geom_col(colour = "white") +
  scale_fill_manual(values = sweep_cols, drop = FALSE) +
  labs(title = "ANI cluster size and sweep membership (counts)",
       x = "ANI cluster", y = "Number of genomes", fill = NULL) +
  theme_minimal()
print(p_cluster_count)

# 2c. heatmap: cluster x sweep, fill = share of the cluster.
#     If your hypothesis holds, each row lights up in ~one column.
p_heat <- ggplot(cluster_composition,
                 aes(x = sweep_lab, y = Cluster, fill = prop)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(prop > 0, paste0(round(100 * prop), "%"), "")),
            size = 2.8) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c",
                      labels = scales::percent_format()) +
  labs(title = "ANI cluster x sweep membership",
       x = NULL, y = "ANI cluster", fill = "Share of\ncluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_heat)

# 2d. (optional) alluvial flow: cluster -> sweep, ribbon width = genomes
if (requireNamespace("ggalluvial", quietly = TRUE)) {
  library(ggalluvial)
  allu <- genome_master %>%
    filter(Cluster %in% cl_levels) %>%
    count(Cluster, sweep_lab, name = "n") %>%
    mutate(Cluster   = factor(Cluster, levels = cl_levels),
           sweep_lab = factor(sweep_lab, levels = sweep_levels))
  p_allu <- ggplot(allu,
                   aes(axis1 = Cluster, axis2 = sweep_lab, y = n)) +
    geom_alluvium(aes(fill = Cluster), alpha = 0.8) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("ANI cluster", "Sweep"), expand = c(.1, .1)) +
    labs(title = "Flow of ANI clusters into sweeps",
         y = "Number of genomes") +
    theme_minimal() +
    theme(legend.position = "none")
  print(p_allu)
} else {
  message("Install ggalluvial for the flow diagram: install.packages('ggalluvial')")
}


# =====================================================================
#  PART 3 :  quantify the correlation
# =====================================================================
# Purity: for each cluster, of the genomes that ARE in some sweep, what
# fraction sit in that cluster's single dominant sweep (1 = perfectly clean).
purity <- cluster_composition %>%
  filter(sweep_lab != "none") %>%
  group_by(Cluster) %>%
  summarise(genomes_in_any_sweep = sum(n),
            dominant_sweep_share = max(n) / sum(n),
            .groups = "drop") %>%
  arrange(desc(dominant_sweep_share))
print(purity)

# One overall association number (Cramer's V) for cluster x sweep.
# Note: chi-square approximation is rough on sparse tables -> heuristic only.
sub_for_test <- genome_master %>%
  filter(sweep_lab != "none", Cluster %in% cl_levels) %>%
  mutate(Cluster   = droplevels(factor(Cluster)),
         sweep_lab = droplevels(factor(sweep_lab)))
tab <- table(sub_for_test$Cluster, sub_for_test$sweep_lab)
chi <- suppressWarnings(chisq.test(tab))
cramers_v <- sqrt(as.numeric(chi$statistic) /
                    (sum(tab) * (min(dim(tab)) - 1)))
cat("Cramer's V (ANI cluster vs sweep):", round(cramers_v, 3),
    "  (0 = independent, 1 = each cluster maps to one sweep)\n")














# ── SNP cluster analysis per sweep ──────────────────────────────────────────

# load snp clusters metadata
metadata_snpclusters <- metadata_snpclusters %>%
  mutate(
    id = str_extract(Assembly, "(?<=_)[0-9]+") %>%
      str_remove("^0+")
  )


# Filter metadata_snp_clusters by sweep gene patterns
snp_sweep1 <- metadata_snpclusters %>%
  filter(str_detect(id, gene_pattern1))
snp_sweep2 <- metadata_snpclusters %>%
  filter(str_detect(id, gene_pattern2))
snp_sweep3 <- metadata_snpclusters %>%
  filter(str_detect(id, gene_pattern3))
snp_sweep4 <- metadata_snpclusters %>%
  filter(str_detect(id, gene_pattern4))

# Add sweep label to each and combine
snp_by_sweep <- bind_rows(
  snp_sweep1 %>% mutate(Sweep = "Sweep 1"),
  snp_sweep2 %>% mutate(Sweep = "Sweep 2"),
  snp_sweep3 %>% mutate(Sweep = "Sweep 3"),
  snp_sweep4 %>% mutate(Sweep = "Sweep 4")
)

# Summary: count of genomes per SNP cluster per sweep
snp_cluster_summary <- snp_by_sweep %>%
  group_by(Sweep, SNP.cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Sweep, desc(n))

print(snp_cluster_summary)

# Plot: SNP cluster distribution per sweep
ggplot(snp_cluster_summary, aes(x = factor(SNP.cluster), y = n, fill = Sweep)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Sweep, scales = "free_x") +
  labs(
    title = "SNP Cluster Distribution per Sweep",
    x = "SNP Cluster",
    y = "Number of Genomes"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# ── Find all samples in sweep-associated SNP clusters ───────────────────────

# Extract the unique SNP clusters found in each sweep
snp_clusters_sweep1 <- unique(snp_sweep1$SNP.cluster)
snp_clusters_sweep2 <- unique(snp_sweep2$SNP.cluster)
snp_clusters_sweep3 <- unique(snp_sweep3$SNP.cluster)
snp_clusters_sweep4 <- unique(snp_sweep4$SNP.cluster)

# Find all samples in metadata_snpclusters that share those SNP clusters
# (including samples not matched by the gene pattern)
all_in_snpcluster1 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep1)
all_in_snpcluster2 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep2)
all_in_snpcluster3 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep3)
all_in_snpcluster4 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep4)

# ── Find all samples in sweep-associated SNP clusters THAT ARE IN THE ANALYSIS ──

all_in_snpcluster_and_dataset1 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep1) %>%
  filter(id %in% genomes_in_analysis$id)
all_in_snpcluster_and_dataset2 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep2) %>%
  filter(id %in% genomes_in_analysis$id)
all_in_snpcluster_and_dataset3 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep3) %>%
  filter(id %in% genomes_in_analysis$id)
all_in_snpcluster_and_dataset4 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep4) %>%
  filter(id %in% genomes_in_analysis$id)

# Update the summary function to reflect the dataset filter
summarise_snpcluster_expansion <- function(sweep_label, sweep_df, all_df) {
  n_sweep  <- nrow(sweep_df)
  n_total  <- nrow(all_df)
  n_extra  <- n_total - n_sweep
  additional <- all_df %>% filter(!id %in% sweep_df$id)
  
  cat("\n──", sweep_label, "──\n")
  cat("  SNP clusters:                        ", paste(unique(all_df$SNP.cluster), collapse = ", "), "\n")
  cat("  Genomes in sweep:                    ", n_sweep, "\n")
  cat("  Total in SNP cluster & in analysis:  ", n_total, "\n")
  cat("  Additional samples (same cluster):   ", n_extra, "\n")
  
  return(additional)
}

additional_sweep1 <- summarise_snpcluster_expansion("Sweep 1", snp_sweep1, all_in_snpcluster_and_dataset1)
additional_sweep2 <- summarise_snpcluster_expansion("Sweep 2", snp_sweep2, all_in_snpcluster_and_dataset2)
additional_sweep3 <- summarise_snpcluster_expansion("Sweep 3", snp_sweep3, all_in_snpcluster_and_dataset3)
additional_sweep4 <- summarise_snpcluster_expansion("Sweep 4", snp_sweep4, all_in_snpcluster_and_dataset4)

# Combined view of all additional samples across sweeps
additional_al <- bind_rows(
  additional_sweep1 %>% mutate(Sweep = "Sweep 1"),
  additional_sweep2 %>% mutate(Sweep = "Sweep 2"),
  additional_sweep3 %>% mutate(Sweep = "Sweep 3"),
  additional_sweep4 %>% mutate(Sweep = "Sweep 4")
)

martrix <- read.delim("/Users/martinsenekowitsch/Downloads/snp_matrix.tab", sep = "\t", header = TRUE)
