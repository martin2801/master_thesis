##set working directory
setwd("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/PopCoGenomeS/Output_500_0.5")

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

# Load data
data_raw           <- read.delim("salmonella_outgroups_rerun.length_bias_500.txt", sep = "\t", header = TRUE)
data_filtered      <- read.delim("salmonella_outgroups_rerun.length_bias.filtered.txt", sep = "\t", header = TRUE)
metadata           <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/infantis_clean_genome_collection_with_bioproject.csv", sep = ",", header = TRUE)
metadata_snpclusters <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/PopCoGenomeS/Output_500_0.5/isolates_infantis_env_filtered.tsv", sep = ",", header = TRUE)

genomes_in_analysis <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/Recombination/filtered_samples.txt", header = FALSE)
genomes_in_analysis$id <- gsub("\\.fna$", "", genomes_in_analysis$V1)
genomes_in_analysis <- genomes_in_analysis %>%
  mutate(id = str_extract(id, "(?<=GCA_)\\d+"))
genomes_in_analysis <- genomes_in_analysis %>%
  mutate(id = str_remove(id, "^0+"))

# Read in clusters_table.tsv
clusters <- read.table("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/Clock/Trees/clusters_table.tsv",
                       header = FALSE, sep = "\t",
                       col.names = c("Genome", "Cluster"),
                       stringsAsFactors = FALSE)
clusters$Cluster <- as.factor(clusters$Cluster)
clusters$Genome  <- gsub("\\.fna$", "", clusters$Genome)
clusters_short   <- clusters %>%
  mutate(Genome = str_extract(Genome, "(?<=GCA_)\\d+")) %>%
  mutate(id = str_remove(Genome, "^0+"))

# ── Define genomes in each sweep (numeric IDs) ──────────────────────────────

# Sweep 1  (n=3)
genomes_sweep1 <- c(
  "32036575", "32037335", "31832435"
)

# Sweep 2  (n=4)
genomes_sweep2 <- c(
  "45808965", "45808925", "45808885", "45808905"
)

# Sweep 3  (n=7)
genomes_sweep3 <- c(
  "11452785", "8046645", "11453035", "8091485", "11452675", "11452375",
  "11452255"
)

# Sweep 4  (n=4)
genomes_sweep4 <- c(
  "20671385", "20671465", "20671355", "20671165"
)

# Sweep 5  (n=3)
genomes_sweep5 <- c(
  "8927125", "8928825", "8026935"
)

# Sweep 6  (n=3)
genomes_sweep6 <- c(
  "8637495", "8484365", "26733835"
)

# Sweep 7  (n=206)
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

# ── Build regex patterns ─────────────────────────────────────────────────────
gene_pattern1 <- paste(genomes_sweep1, collapse = "|")
gene_pattern2 <- paste(genomes_sweep2, collapse = "|")
gene_pattern3 <- paste(genomes_sweep3, collapse = "|")
gene_pattern4 <- paste(genomes_sweep4, collapse = "|")
gene_pattern5 <- paste(genomes_sweep5, collapse = "|")
gene_pattern6 <- paste(genomes_sweep6, collapse = "|")
gene_pattern7 <- paste(genomes_sweep7, collapse = "|")

# ── Extract metadata per sweep ───────────────────────────────────────────────
metadata <- metadata %>%
  mutate(id = str_extract(AssemblyAccession, "(?<=_)[0-9]+") %>% str_remove("^0+"))

sweep_metadata1 <- metadata %>% filter(str_detect(id, gene_pattern1))
sweep_metadata2 <- metadata %>% filter(str_detect(id, gene_pattern2))
sweep_metadata3 <- metadata %>% filter(str_detect(id, gene_pattern3))
sweep_metadata4 <- metadata %>% filter(str_detect(id, gene_pattern4))
sweep_metadata5 <- metadata %>% filter(str_detect(id, gene_pattern5))
sweep_metadata6 <- metadata %>% filter(str_detect(id, gene_pattern6))
sweep_metadata7 <- metadata %>% filter(str_detect(id, gene_pattern7))

# ── Extract cluster membership per sweep ────────────────────────────────────
sweep_clusters1 <- clusters_short %>% filter(str_detect(id, gene_pattern1))
sweep_clusters2 <- clusters_short %>% filter(str_detect(id, gene_pattern2))
sweep_clusters3 <- clusters_short %>% filter(str_detect(id, gene_pattern3))
sweep_clusters4 <- clusters_short %>% filter(str_detect(id, gene_pattern4))
sweep_clusters5 <- clusters_short %>% filter(str_detect(id, gene_pattern5))
sweep_clusters6 <- clusters_short %>% filter(str_detect(id, gene_pattern6))
sweep_clusters7 <- clusters_short %>% filter(str_detect(id, gene_pattern7))

# ── SNP cluster analysis per sweep ──────────────────────────────────────────
metadata_snpclusters <- metadata_snpclusters %>%
  mutate(id = str_extract(Assembly, "(?<=_)[0-9]+") %>% str_remove("^0+"))

snp_sweep1 <- metadata_snpclusters %>% filter(str_detect(id, gene_pattern1))
snp_sweep2 <- metadata_snpclusters %>% filter(str_detect(id, gene_pattern2))
snp_sweep3 <- metadata_snpclusters %>% filter(str_detect(id, gene_pattern3))
snp_sweep4 <- metadata_snpclusters %>% filter(str_detect(id, gene_pattern4))
snp_sweep5 <- metadata_snpclusters %>% filter(str_detect(id, gene_pattern5))
snp_sweep6 <- metadata_snpclusters %>% filter(str_detect(id, gene_pattern6))
snp_sweep7 <- metadata_snpclusters %>% filter(str_detect(id, gene_pattern7))

snp_by_sweep <- bind_rows(
  snp_sweep1 %>% mutate(Sweep = "Sweep 1"),
  snp_sweep2 %>% mutate(Sweep = "Sweep 2"),
  snp_sweep3 %>% mutate(Sweep = "Sweep 3"),
  snp_sweep4 %>% mutate(Sweep = "Sweep 4"),
  snp_sweep5 %>% mutate(Sweep = "Sweep 5"),
  snp_sweep6 %>% mutate(Sweep = "Sweep 6"),
  snp_sweep7 %>% mutate(Sweep = "Sweep 7")
)

snp_cluster_summary <- snp_by_sweep %>%
  group_by(Sweep, SNP.cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Sweep, desc(n))

print(snp_cluster_summary)

ggplot(snp_cluster_summary, aes(x = factor(SNP.cluster), y = n, fill = Sweep)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Sweep, scales = "free_x") +
  labs(title = "SNP Cluster Distribution per Sweep", x = "SNP Cluster", y = "Number of Genomes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# ── Find all samples in sweep-associated SNP clusters ───────────────────────
snp_clusters_sweep1 <- unique(snp_sweep1$SNP.cluster)
snp_clusters_sweep2 <- unique(snp_sweep2$SNP.cluster)
snp_clusters_sweep3 <- unique(snp_sweep3$SNP.cluster)
snp_clusters_sweep4 <- unique(snp_sweep4$SNP.cluster)
snp_clusters_sweep5 <- unique(snp_sweep5$SNP.cluster)
snp_clusters_sweep6 <- unique(snp_sweep6$SNP.cluster)
snp_clusters_sweep7 <- unique(snp_sweep7$SNP.cluster)

all_in_snpcluster_and_dataset1 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep1) %>% filter(id %in% genomes_in_analysis$id)
all_in_snpcluster_and_dataset2 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep2) %>% filter(id %in% genomes_in_analysis$id)
all_in_snpcluster_and_dataset3 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep3) %>% filter(id %in% genomes_in_analysis$id)
all_in_snpcluster_and_dataset4 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep4) %>% filter(id %in% genomes_in_analysis$id)
all_in_snpcluster_and_dataset5 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep5) %>% filter(id %in% genomes_in_analysis$id)
all_in_snpcluster_and_dataset6 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep6) %>% filter(id %in% genomes_in_analysis$id)
all_in_snpcluster_and_dataset7 <- metadata_snpclusters %>%
  filter(SNP.cluster %in% snp_clusters_sweep7) %>% filter(id %in% genomes_in_analysis$id)

# ── Summarise SNP cluster expansion per sweep ────────────────────────────────
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
additional_sweep5 <- summarise_snpcluster_expansion("Sweep 5", snp_sweep5, all_in_snpcluster_and_dataset5)
additional_sweep6 <- summarise_snpcluster_expansion("Sweep 6", snp_sweep6, all_in_snpcluster_and_dataset6)
additional_sweep7 <- summarise_snpcluster_expansion("Sweep 7", snp_sweep7, all_in_snpcluster_and_dataset7)

additional_all <- bind_rows(
  additional_sweep1 %>% mutate(Sweep = "Sweep 1"),
  additional_sweep2 %>% mutate(Sweep = "Sweep 2"),
  additional_sweep3 %>% mutate(Sweep = "Sweep 3"),
  additional_sweep4 %>% mutate(Sweep = "Sweep 4"),
  additional_sweep5 %>% mutate(Sweep = "Sweep 5"),
  additional_sweep6 %>% mutate(Sweep = "Sweep 6"),
  additional_sweep7 %>% mutate(Sweep = "Sweep 7")
)


# Bubble plot
library(forcats)

snp_cluster_summary <- snp_by_sweep %>%
  group_by(Sweep, SNP.cluster) %>%
  summarise(n = n(), .groups = "drop")

ggplot(snp_cluster_summary,
       aes(x = factor(SNP.cluster), y = fct_rev(Sweep), size = n)) +
  geom_point(colour = "steelblue", alpha = 0.8) +
  geom_text(aes(label = n), size = 2.8, vjust = -1.3) +
  scale_size_area(max_size = 12) +
  labs(x = "SNP cluster", y = NULL, size = "Genomes",
       title = "Sweep membership across SNP clusters") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Barchart but in % like above
snp_cluster_summary <- snp_by_sweep %>%
  group_by(Sweep, SNP.cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sweep) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  arrange(Sweep, desc(pct))

print(snp_cluster_summary)

ggplot(snp_cluster_summary, aes(x = factor(SNP.cluster), y = pct, fill = Sweep)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Sweep, scales = "free_x") +
  labs(title = "SNP Cluster Distribution per Sweep",
       x = "SNP Cluster", y = "% of sweep genomes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# extra samples in SNP clusters
sweep_dfs <- list(snp_sweep1, snp_sweep2, snp_sweep3, snp_sweep4,
                  snp_sweep5, snp_sweep6, snp_sweep7)
all_dfs   <- list(all_in_snpcluster_and_dataset1, all_in_snpcluster_and_dataset2,
                  all_in_snpcluster_and_dataset3, all_in_snpcluster_and_dataset4,
                  all_in_snpcluster_and_dataset5, all_in_snpcluster_and_dataset6,
                  all_in_snpcluster_and_dataset7)

expansion_summary <- tibble(
  Sweep    = paste("Sweep", seq_along(sweep_dfs)),
  in_sweep = sapply(sweep_dfs, nrow),
  total    = sapply(all_dfs, nrow)
) %>%
  mutate(additional = total - in_sweep) %>%
  pivot_longer(c(in_sweep, additional),
               names_to = "Type", values_to = "Genomes")

print(expansion_summary)

ggplot(expansion_summary,
       aes(x = Genomes, y = fct_rev(Sweep), fill = Type)) +
  geom_col() +
  scale_fill_manual(
    values = c(in_sweep = "#2c7fb8", additional = "#a6bddb"),
    labels = c(in_sweep = "In sweep",
               additional = "Same SNP cluster, not in sweep")) +
  labs(title = "Sweep size vs. SNP-cluster expansion",
       x = "Number of genomes", y = NULL, fill = NULL) +
  theme_bw()
