.libPaths("~/R/library")

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(tidyverse)


# install.packages("tidyverse")

# set working directory
setwd("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/BLAST")

# load data
# Data related to enriched genes
genes_enr <- read.delim("enriched_genes.tsv", sep = "\t", header = TRUE)
genes_dep <- read.delim("depleted_genes.tsv", sep = "\t", header = TRUE)
go_enr <- read.delim("enriched_go.tsv", sep = "\t", header = TRUE)
go_dep <- read.delim("depleted_go.tsv", sep = "\t", header = TRUE)
cog_sig <- read.delim("sig_cog_results.tsv", sep = "\t", header = TRUE)

# Data related to pident
sum_pident <- read.delim("pident_summary.tsv", sep = "\t", header = TRUE)
int_pident <- read.delim("pident_interesting.tsv", sep = "\t", header = TRUE)
spec_pident <- read.delim("pident_sweep_specific.tsv", sep = "\t", header = TRUE)

int_pident_cog <- read.delim("pident_interesting_cog.tsv", sep = "\t", header = TRUE)

# Data related to pESI presence 
pESI_presence <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/pESI/sweep_labels_pESI.txt", sep = "\t", header = TRUE)

# Data related to other plasmids
plasmids_present_matrix <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/pESI/plasmidfinder_summary.tsv", sep = "\t", header = TRUE)
plasmids_inventory <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/pESI/plasmid_inventory.tsv", sep = "\t", header = TRUE)

# sweep identity
sweep_labels <- read.delim("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/pESI/genome_sweep_labels.txt", sep = "\t", header = TRUE)


# ------------------------------------
# COG
cog_sig_filtered <- cog_sig %>%
  filter(abs(perc_diff) >= 5)

# Genes
genes_enr_sweep1 <- genes_enr %>%
  filter(group == "sweep_1")
genes_enr_sweep2 <- genes_enr %>%
  filter(group == "sweep_2")
genes_enr_sweep3 <- genes_enr %>%
  filter(group == "sweep_3")
genes_enr_sweep4 <- genes_enr %>%
  filter(group == "sweep_4")
genes_enr_sweep5 <- genes_enr %>%
  filter(group == "sweep_5")
genes_enr_sweep6 <- genes_enr %>%
  filter(group == "sweep_6")
genes_enr_sweep7 <- genes_enr %>%
  filter(group == "sweep_7")

write.table(genes_enr_sweep7, "sweep7_enr.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------------------
# pESI
plot_data <- pESI_presence %>%
  mutate(pESI_status = ifelse(pESI == "YES", "pESI +", "pESI -")) %>%
  count(sweep, pESI_status)

ggplot(plot_data, aes(x = pESI_status, y = n, fill = pESI_status)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.3) +
  facet_wrap(~ sweep) +
  scale_fill_manual(values = c("pESI +" = "#2c7fb8", "pESI -" = "#bdbdbd")) +
  labs(x = NULL, y = "Number of genomes", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

# Plasmids
plasmids_filtered <- plasmids_inventory %>%
  # filter(replicons != "IncFIB(pN55391)_1") %>%
  left_join(sweep_labels, by = "genome")

plasmids_filtered_split <- split(plasmids_filtered, plasmids_filtered$replicons)

library(dplyr)
library(ggplot2)

plasmids_filtered %>%
  count(replicons, sort = TRUE) %>%
  ggplot(aes(x = reorder(replicons, n), y = n)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Replicon",
    y = "Count",
    title = "Number of genomes per replicon"
  )

plasmids_filtered %>%
  count(replicons, sweep) %>%
  ggplot(aes(x = replicons, y = n, fill = sweep)) +
  geom_col() +
  coord_flip()

plasmids_filtered %>%
  count(replicons, sweep) %>%
  group_by(replicons) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = replicons, y = prop, fill = sweep)) +
  geom_col() +
  coord_flip()






# Total genomes per sweep group
total_genomes <- sweep_labels %>%
  count(sweep, name = "total_genomes")

# Number of genomes carrying each replicon per sweep group
replicon_counts <- plasmids_filtered %>%
  distinct(genome, sweep, replicons) %>%   # count a replicon once per genome
  count(sweep, replicons, name = "n")

# Calculate proportion
plot_data <- replicon_counts %>%
  left_join(total_genomes, by = "sweep") %>%
  mutate(prop = n / total_genomes * 100)

ggplot(plot_data,
       aes(x = replicons, y = prop, fill = replicons)) +
  geom_col() +
  geom_text(aes(label = paste0(round(prop, 1), "%")),
            vjust = -0.3,
            size = 3) +
  facet_wrap(~ sweep) +
  labs(
    x = "Replicon",
    y = "Genomes carrying replicon (%)",
    fill = "Replicon"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
