.libPaths("~/R/library")


#install.packages("patchwork", lib="~/R/library", repos="https://cloud.r-project.org/")
#remotes::install_version("GlobalOptions", version="0.1.2", lib="~/R/library", repos="https://cloud.r-project.org/")
#BiocManager::install("ComplexHeatmap", lib="~/R/library")

library(dplyr)
library(tidyr)
library(data.table)
library(GO.db)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(remotes)
library(patchwork)

#set working directory
setwd("/home/senekowitsch/Thesis/Functional/05_analyze_BLAST")

# load the data
all_data <- fread("/home/senekowitsch/Thesis/Functional/04_clean_BLAST/all_data.tsv.gz")
head(all_data)
# -------------------------------------------------
# Filter: remove qseqid × query_genome combos that
# hit ALL 472 genomes with pident == 100 everywhere
# -------------------------------------------------

# Since the same qseqid label can appear across different query genomes,
# we group by the combination of query_genome + qseqid to treat each
# blast query independently.

all_data_filtered <- all_data %>%
  group_by(query_genome, qseqid) %>%
  filter(
    !(n_distinct(sseqid) >= 472 & all(pident == 100))
  ) %>%
  ungroup()

# Quick sanity check
cat("Rows before filtering:", nrow(all_data), "\n")
cat("Rows after filtering: ", nrow(all_data_filtered), "\n")
cat("Rows removed:         ", nrow(all_data) - nrow(all_data_filtered), "\n")

# How many unique qseqids were dropped entirely?
dropped_qseqids <- setdiff(
  paste(all_data$query_genome, all_data$qseqid),
  paste(all_data_filtered$query_genome, all_data_filtered$qseqid)
)
cat("query_genome+qseqid combos removed:", length(dropped_qseqids), "\n")

# --------------------------------------------------------------------------------------------------
# -------------------------------------------------
# Check genes in sweep groups
# -------------------------------------------------
# get unique Preferred_names per query genome, excluding "-"
gene_lists <- all_data_filtered %>%
  filter(Preferred_name != "-") %>%
  group_by(query_group, query_genome) %>%
  summarise(genes = list(unique(Preferred_name)), .groups = "drop")

# function to compare the 3 genomes within a sweep group
compare_sweep_genes <- function(group_name) {
  g <- gene_lists %>% filter(query_group == group_name)
  sets <- setNames(g$genes, g$query_genome)
  
  # genes in all 3
  core <- Reduce(intersect, sets)
  # genes in at least 1
  total_union <- Reduce(union, sets)
  # genes not in all 3
  not_shared <- setdiff(total_union, core)
  
  cat("\n====", group_name, "====\n")
  cat("Genes in all 3 genomes (core):", length(core), "\n")
  cat("Genes in union:", length(total_union), "\n")
  cat("Genes not shared by all 3:", length(not_shared), "\n")
  
  # which genome has unique genes
  for (genome in names(sets)) {
    unique_to_this <- setdiff(sets[[genome]], Reduce(union, sets[names(sets) != genome]))
    cat("  Unique to", genome, ":", length(unique_to_this), "\n")
  }
}

# run for all sweep groups
for (grp in c("sweep_1", "sweep_2", "sweep_3", "sweep_4", "sweep_5", "sweep_6", "sweep_7", "no_sweep")) {
  compare_sweep_genes(grp)
}

# check how many unique Preferred_names there are in total (excluding "-")
all_data_filtered %>%
  filter(Preferred_name != "-") %>%
  pull(Preferred_name) %>%
  unique() %>%
  length()

# get union of all Preferred_names across all query genomes, excluding "-"
all_genes <- all_data_filtered %>%
  filter(Preferred_name != "-") %>%
  pull(Preferred_name) %>%
  unique()

cat("Total unique annotated genes:", length(all_genes), "\n")

# build presence/absence matrix
# for each gene, which target genomes have at least one hit
pa_matrix <- all_data_filtered %>%
  filter(Preferred_name != "-") %>%
  group_by(Preferred_name, sseqid) %>%
  summarise(present = 1, .groups = "drop") %>%
  tidyr::pivot_wider(names_from = sseqid, 
                     values_from = present, 
                     values_fill = 0)

cat("Matrix dimensions:", nrow(pa_matrix), "genes x", ncol(pa_matrix)-1, "genomes\n")
pa_matrix[1:5, 1:5]

# how many genes are present in all 472 genomes (universal genes)
universal <- pa_matrix %>%
  filter(rowSums(dplyr::select(., -Preferred_name)) == 472) %>%
  pull(Preferred_name)
cat("Universal genes (in all 472 genomes):", length(universal), "\n")

# how many genes are present in only 1 genome (unique genes)
unique_genes <- pa_matrix %>%
  filter(rowSums(dplyr::select(., -Preferred_name)) == 1) %>%
  pull(Preferred_name)
cat("Genes in only 1 genome:", length(unique_genes), "\n")

# distribution of how many genomes each gene is found in
presence_counts <- rowSums(dplyr::select(pa_matrix, -Preferred_name))
hist(presence_counts, breaks = 50, 
     main = "Gene presence across 472 genomes",
     xlab = "Number of genomes", ylab = "Number of genes")

# load sweep labels if not already loaded
sweep_labels <- read.delim("/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt", 
                           sep = "\t", header = TRUE)
sweep_labels$genome <- as.character(sweep_labels$genome)

# count genomes per group
group_sizes <- sweep_labels %>%
  group_by(sweep) %>%
  summarise(n = n())

print(group_sizes)

# get genome IDs for each group
get_genomes <- function(group) {
  sweep_labels %>% filter(sweep == group) %>% pull(genome)
}

no_sweep_genomes <- get_genomes("no_sweep")
sweep_1_genomes  <- get_genomes("sweep_1")
sweep_2_genomes  <- get_genomes("sweep_2")
sweep_3_genomes  <- get_genomes("sweep_3")
sweep_4_genomes  <- get_genomes("sweep_4")
sweep_5_genomes  <- get_genomes("sweep_5")
sweep_6_genomes  <- get_genomes("sweep_6")
sweep_7_genomes  <- get_genomes("sweep_7")

# calculate frequency per group
gene_freq <- pa_matrix %>%
  mutate(
    no_sweep_freq = rowSums(dplyr::select(., all_of(no_sweep_genomes)))  / 242,
    sweep_1_freq  = rowSums(dplyr::select(., all_of(sweep_1_genomes)))   / 3,
    sweep_2_freq  = rowSums(dplyr::select(., all_of(sweep_2_genomes)))   / 4,
    sweep_3_freq  = rowSums(dplyr::select(., all_of(sweep_3_genomes)))   / 7,
    sweep_4_freq  = rowSums(dplyr::select(., all_of(sweep_4_genomes)))   / 4,
    sweep_5_freq  = rowSums(dplyr::select(., all_of(sweep_5_genomes)))   / 3,
    sweep_6_freq  = rowSums(dplyr::select(., all_of(sweep_6_genomes)))   / 3,
    sweep_7_freq  = rowSums(dplyr::select(., all_of(sweep_7_genomes)))   / 206
  ) %>%
  dplyr::select(Preferred_name, ends_with("_freq"))

head(gene_freq)

# remove genes present in >90% of all 472 genomes (universal genes)
gene_freq_filtered <- gene_freq %>%
  filter(!(no_sweep_freq > 0.9 & 
             sweep_1_freq  > 0.9 & 
             sweep_2_freq  > 0.9 & 
             sweep_3_freq  > 0.9 & 
             sweep_4_freq  > 0.9 & 
             sweep_5_freq  > 0.9 & 
             sweep_6_freq  > 0.9 & 
             sweep_7_freq  > 0.9))

cat("Genes before filtering:", nrow(gene_freq), "\n")
cat("Genes after filtering:", nrow(gene_freq_filtered), "\n")

# -------------------------------------------------
# Fishers exact test for each gene and each sweep group vs all others
# -------------------------------------------------
# convert pa_matrix to long format with sweep labels attached
pa_long <- pa_matrix %>%
  tidyr::pivot_longer(cols = -Preferred_name, 
                      names_to = "genome", 
                      values_to = "present") %>%
  left_join(sweep_labels, by = "genome")

# function to run Fisher's test for one gene x one group
fisher_test <- function(gene_name, group_name) {
  gene_data <- pa_long %>% filter(Preferred_name == gene_name)
  
  a <- sum(gene_data$present[gene_data$sweep == group_name])           # has gene, in group
  b <- sum(!gene_data$present[gene_data$sweep == group_name])          # no gene, in group  
  c <- sum(gene_data$present[gene_data$sweep != group_name])           # has gene, other groups
  d <- sum(!gene_data$present[gene_data$sweep != group_name])          # no gene, other groups
  
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  
  data.frame(
    Preferred_name = gene_name,
    group          = group_name,
    n_in_group     = a + b,
    present_in_group = a,
    present_other  = c,
    odds_ratio     = ft$estimate,
    p_value        = ft$p.value
  )
}

# run for all genes x all groups
groups <- c("no_sweep", "sweep_1", "sweep_2", "sweep_3",
            "sweep_4", "sweep_5", "sweep_6", "sweep_7")

all_genes <- gene_freq_filtered$Preferred_name

cat("Running Fisher's tests:", length(all_genes), "genes x", length(groups), "groups\n")

fisher_results <- bind_rows(lapply(groups, function(grp) {
  bind_rows(lapply(all_genes, function(gene) {
    fisher_test(gene, grp)
  }))
}))

# correct for multiple testing
fisher_results <- fisher_results %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))

# look at significant results
sig_results <- fisher_results %>%
  filter(p_adjusted < 0.05) %>%
  arrange(p_adjusted)

cat("Significant results:", nrow(sig_results), "\n")
head(sig_results)

# enriched in group
enriched <- sig_results %>%
  filter(odds_ratio > 1) %>%
  arrange(desc(odds_ratio))

# depleted in group  
depleted <- sig_results %>%
  filter(odds_ratio < 1) %>%
  arrange(odds_ratio)

cat("Enriched:", nrow(enriched), "\n")
cat("Depleted:", nrow(depleted), "\n")

head(enriched)
head(depleted)

enriched <- enriched %>%
  mutate(
    freq_in_group   = present_in_group / n_in_group,
    freq_out_group  = present_other / (472 - n_in_group)
  ) %>%
  dplyr::select(Preferred_name, group, odds_ratio, freq_in_group, freq_out_group, 
         n_in_group, present_in_group, present_other, p_adjusted) %>%
  arrange(desc(odds_ratio))

depleted <- depleted %>%
  mutate(
    freq_in_group  = present_in_group / n_in_group,
    freq_out_group = present_other / (472 - n_in_group)
  ) %>%
  dplyr::select(Preferred_name, group, odds_ratio, freq_in_group, freq_out_group,
         n_in_group, present_in_group, present_other, p_adjusted) %>%
  arrange(odds_ratio)

# add COG category information
# get COG category per Preferred_name (most common non-NA value)
cog_lookup <- all_data_filtered %>%
  filter(Preferred_name != "-", !is.na(COG_category)) %>%
  group_by(Preferred_name) %>%
  summarise(COG_category = names(sort(table(COG_category), decreasing = TRUE))[1],
            .groups = "drop")

# join to enriched and depleted
enriched <- enriched %>% left_join(cog_lookup, by = "Preferred_name")
depleted <- depleted %>% left_join(cog_lookup, by = "Preferred_name")

head(enriched)
head(depleted)

write.table(enriched, "enriched_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(depleted, "depleted_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# --------------------------------------------------------------------------------------------------
# -------------------------------------------------
# Check COG categories in sweep groups
# -------------------------------------------------
# what COG categories do you have?
table(all_data_filtered$COG_category, useNA = "always")
# -> many multi letter COGs

# how many have multi-letter COGs?
all_data_filtered %>%
  filter(COG_category != "-", !is.na(COG_category)) %>%
  mutate(n_cats = nchar(COG_category)) %>%
  count(n_cats)

# count unique genes per COG category per target genome
# distinct(sseqid, Preferred_name, COG_split) ensures each gene is only counted once
# per target genome, even if multiple query genomes contributed that hit
cog_counts <- all_data_filtered %>%
  filter(COG_category != "-", !is.na(COG_category)) %>%
  mutate(COG_split = strsplit(COG_category, "")) %>%
  tidyr::unnest(COG_split) %>%
  distinct(sseqid, Preferred_name, COG_split) %>%
  group_by(sseqid, COG_split) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(sseqid = as.character(sseqid)) %>%
  left_join(sweep_labels, by = c("sseqid" = "genome"))

# sanity check - median count per COG category should be in the range of tens to hundreds
cog_counts %>%
  group_by(COG_split) %>%
  summarise(median_count = median(count)) %>%
  arrange(desc(median_count))

# fill in 0s for genomes that have no genes in a given COG category
cog_counts_complete <- cog_counts %>%
  complete(sseqid, COG_split, fill = list(count = 0)) %>%
  mutate(sseqid = as.character(sseqid)) %>%
  left_join(sweep_labels, by = c("sseqid" = "genome")) %>%
  mutate(sweep = coalesce(sweep.x, sweep.y)) %>%
  dplyr::select(-sweep.x, -sweep.y)

# Wilcoxon test for each COG category x sweep group vs all others combined
# = Mann Whitney U test, instead of t-test when we cant assume normality
# valid here becasue non-negative integers, heavily right skewed

cog_groups <- unique(cog_counts_complete$COG_split)

sweep_groups <- c("sweep_1", "sweep_2", "sweep_3", "sweep_4",
                  "sweep_5", "sweep_6", "sweep_7")

wilcox_results <- bind_rows(lapply(sweep_groups, function(grp) {
  bind_rows(lapply(cog_groups, function(cog) {
    
    in_group  <- cog_counts_complete %>%
      filter(COG_split == cog, sweep == grp) %>%
      pull(count)
    out_group <- cog_counts_complete %>%
      filter(COG_split == cog, sweep != grp) %>%  # all others combined
      pull(count)
    
    wt <- wilcox.test(in_group, out_group)
    
    data.frame(
      COG_split    = cog,
      group        = grp,
      median_sweep = median(in_group),
      median_other = median(out_group),
      mean_sweep  = mean(in_group),
      mean_other  = mean(out_group),
      perc_diff = (median(in_group) - median(out_group)) / median(out_group) * 100,
      p_value      = wt$p.value
    )
  }))
}))

wilcox_results <- wilcox_results %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"),
         direction  = ifelse(median_sweep > median_other, "enriched", "depleted"))

sig_cog <- wilcox_results %>%
  filter(p_adjusted < 0.05) %>%
  arrange(p_adjusted)

cat("Significant COG results:", nrow(sig_cog), "\n")
head(sig_cog)

write.table(sig_cog, "sig_cog_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

sig_cog %>%
  dplyr::select(COG_split, group, median_sweep, median_other, direction, p_adjusted) %>%
  arrange(group, direction, p_adjusted)

# Column names:
# COG_split — the single-letter COG functional category (e.g. V = defense mechanisms, N = cell motility)
# group — the sweep group being tested
# median_sweep — the median number of unique genes in that COG category per genome within the sweep group
# median_other — the median number of unique genes in that COG category per genome outside the sweep group (all others combined)
# direction — enriched means sweep genomes have more genes in that category than others; depleted means fewer
# p_adjusted — Benjamini-Hochberg corrected p-value from the Wilcoxon test
#   Controls false discovery rate. At p<0.05 we would expect at most 5% false positives
#   alternative would be Bonferroni but more strict, here we can tolerate small number of
#   false positives. We want to find candidates to investigate further.

# --------------------------------------------------------------------------------------------------
# -------------------------------------------------
# Check GO terms in sweep groups
# -------------------------------------------------
# how does the GO column look?
head(all_data_filtered$GOs)

# function to count genes meeting the threshold for each sweep group
check_thresholds <- function(freq_df, sweep_col, other_cols, sweep_thresh, other_thresh) {
  freq_df %>%
    filter(!!sym(sweep_col) >= sweep_thresh) %>%
    filter(if_all(all_of(other_cols), ~ . <= other_thresh)) %>%
    nrow()
}

other_cols <- c("no_sweep_freq", "sweep_1_freq", "sweep_2_freq", "sweep_3_freq",
                "sweep_4_freq", "sweep_5_freq", "sweep_6_freq", "sweep_7_freq")

# try different thresholds for each sweep group
for (grp in c("sweep_1_freq", "sweep_2_freq", "sweep_3_freq", "sweep_4_freq",
              "sweep_5_freq", "sweep_6_freq", "sweep_7_freq", "no_sweep_freq")) {
  
  others <- setdiff(other_cols, grp)
  
  cat("\n====", grp, "====\n")
  for (st in c(0.8, 0.9, 0.95)) {
    for (ot in c(0.1, 0.2, 0.3)) {
      n <- check_thresholds(gene_freq, grp, others, st, ot)
      cat("  sweep >=", st, "& others <=", ot, ":", n, "genes\n")
    }
  }
}

# what are the actual genes and their GO terms for sweep_1?
gene_freq %>%
  filter(sweep_1_freq >= 0.8) %>%
  filter(if_all(all_of(setdiff(other_cols, "sweep_1_freq")), ~ . <= 0.2)) %>%
  left_join(
    all_data_filtered %>% distinct(Preferred_name, GOs),
    by = "Preferred_name"
  ) %>%
  dplyr::select(Preferred_name, GOs)

# sweep_2
gene_freq %>%
  filter(sweep_2_freq >= 0.8) %>%
  filter(if_all(all_of(setdiff(other_cols, "sweep_2_freq")), ~ . <= 0.2)) %>%
  left_join(
    all_data_filtered %>% distinct(Preferred_name, GOs),
    by = "Preferred_name"
  ) %>%
  dplyr::select(Preferred_name, GOs)

# sweep_3
gene_freq %>%
  filter(sweep_3_freq >= 0.8) %>%
  filter(if_all(all_of(setdiff(other_cols, "sweep_3_freq")), ~ . <= 0.2)) %>%
  left_join(
    all_data_filtered %>% distinct(Preferred_name, GOs),
    by = "Preferred_name"
  ) %>%
  dplyr::select(Preferred_name, GOs)

# sweep_4
gene_freq %>%
  filter(sweep_4_freq >= 0.8) %>%
  filter(if_all(all_of(setdiff(other_cols, "sweep_4_freq")), ~ . <= 0.2)) %>%
  left_join(
    all_data_filtered %>% distinct(Preferred_name, GOs),
    by = "Preferred_name"
  ) %>%
  dplyr::select(Preferred_name, GOs)

# sweep_5
gene_freq %>%
  filter(sweep_5_freq >= 0.8) %>%
  filter(if_all(all_of(setdiff(other_cols, "sweep_5_freq")), ~ . <= 0.2)) %>%
  left_join(
    all_data_filtered %>% distinct(Preferred_name, GOs),
    by = "Preferred_name"
  ) %>%
  dplyr::select(Preferred_name, GOs)

# sweep_6
gene_freq %>%
  filter(sweep_6_freq >= 0.8) %>%
  filter(if_all(all_of(setdiff(other_cols, "sweep_6_freq")), ~ . <= 0.2)) %>%
  left_join(
    all_data_filtered %>% distinct(Preferred_name, GOs),
    by = "Preferred_name"
  ) %>%
  dplyr::select(Preferred_name, GOs)

# sweep_7
gene_freq %>%
  filter(sweep_7_freq >= 0.8) %>%
  filter(if_all(all_of(setdiff(other_cols, "sweep_7_freq")), ~ . <= 0.2)) %>%
  left_join(
    all_data_filtered %>% distinct(Preferred_name, GOs),
    by = "Preferred_name"
  ) %>%
  dplyr::select(Preferred_name, GOs)

# -> not conclusive, we need different approach to look at GO terms, maybe enrichment analysis


# Look at which GO terms we have total and if there is an enrichment of certain GO terms in sweep groups
# split GO terms and build long format
go_long <- all_data_filtered %>%
  filter(GOs != "-", !is.na(GOs)) %>%
  distinct(sseqid, Preferred_name, GOs) %>%
  mutate(GO = strsplit(GOs, ",")) %>%
  tidyr::unnest(GO) %>%
  mutate(GO = trimws(GO))

# then build presence/absence matrix with GO terms as rows
go_pa_matrix <- go_long %>%
  distinct(GO, sseqid) %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = sseqid,
                     values_from = present,
                     values_fill = 0)

cat("GO matrix dimensions:", nrow(go_pa_matrix), "terms x", ncol(go_pa_matrix)-1, "genomes\n")
write.table(go_pa_matrix, "go_pa_matrix.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# convert to long format with sweep labels
go_pa_long <- go_pa_matrix %>%
  tidyr::pivot_longer(cols = -GO,
                      names_to = "genome",
                      values_to = "present") %>%
  left_join(sweep_labels, by = "genome")

# fisher test function for GO terms
fisher_test_go <- function(go_term, group_name) {
  gene_data <- go_pa_long %>% filter(GO == go_term)
  
  a <- sum(gene_data$present[gene_data$sweep == group_name])
  b <- sum(!gene_data$present[gene_data$sweep == group_name])
  c <- sum(gene_data$present[gene_data$sweep != group_name])
  d <- sum(!gene_data$present[gene_data$sweep != group_name])
  
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  
  data.frame(
    GO         = go_term,
    group      = group_name,
    n_in_group = a + b,
    present_in_group = a,
    present_other    = c,
    odds_ratio = ft$estimate,
    p_value    = ft$p.value
  )
}

# run for all GO terms x all groups
all_go_terms <- go_pa_matrix$GO
groups <- c("no_sweep", "sweep_1", "sweep_2", "sweep_3",
            "sweep_4", "sweep_5", "sweep_6", "sweep_7")

cat("Running Fisher's tests:", length(all_go_terms), "GO terms x", length(groups), "groups\n")

go_fisher_results <- bind_rows(lapply(groups, function(grp) {
  bind_rows(lapply(all_go_terms, function(go) {
    fisher_test_go(go, grp)
  }))
}))

# correct for multiple testing
go_fisher_results <- go_fisher_results %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))

# significant results
go_sig <- go_fisher_results %>%
  filter(p_adjusted < 0.05)

go_enriched <- go_sig %>% filter(odds_ratio > 1) %>% arrange(desc(odds_ratio))
go_depleted <- go_sig %>% filter(odds_ratio < 1) %>% arrange(odds_ratio)

cat("GO terms enriched:", nrow(go_enriched), "\n")
cat("GO terms depleted:", nrow(go_depleted), "\n")


# get GO term descriptions
go_terms <- AnnotationDbi::select(GO.db, 
                                   keys = unique(c(go_enriched$GO, go_depleted$GO)),
                                   columns = c("GOID", "TERM", "ONTOLOGY"),
                                   keytype = "GOID")

# join to results
go_enriched <- go_enriched %>% left_join(go_terms, by = c("GO" = "GOID"))
go_depleted <- go_depleted %>% left_join(go_terms, by = c("GO" = "GOID"))

write.table(go_enriched, "enriched_go.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(go_depleted, "depleted_go.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# --------------------------------------------------------------------------------------------------
# -------------------------------------------------
# Plotting results
# -------------------------------------------------
setwd("/home/senekowitsch/Thesis/Functional/05_analyze_BLAST")

# --------------------------------------------------------------------------------------------------
# Heatmap, use all genes from  filtered gene frequency table
# These are the genes that survived your >90% universal filter
all_accessory_genes <- gene_freq$Preferred_name

# 2. Prepare the matrix using all accessory genes
mat_all <- pa_matrix %>% 
  filter(Preferred_name %in% all_accessory_genes) %>% 
  tibble::column_to_rownames("Preferred_name") %>%
  as.matrix()

# 3. Ensure the column (genome) order matches your annotation
# We use the full sweep_labels here
anno_df_all <- sweep_labels %>% 
  filter(genome %in% colnames(mat_all)) %>% 
  arrange(sweep)

mat_all <- mat_all[, anno_df_all$genome]

# 4. Create the annotation with the full color palette
sweep_colors <- c(
  "sweep_1"  = "#E41A1C", "sweep_2"  = "#377EB8", 
  "sweep_3"  = "#4DAF4A", "sweep_4"  = "#984EA3", 
  "sweep_5"  = "#FF7F00", "sweep_6"  = "#FFFF33", 
  "sweep_7"  = "#A65628", "no_sweep" = "#999999"
)

col_anno_all <- HeatmapAnnotation(
  Sweep = anno_df_all$sweep, 
  col = list(Sweep = sweep_colors)
)

draw(Heatmap(mat_all, 
             name = "Presence", 
             col = c("0" = "white", "1" = "black"),
             top_annotation = col_anno_all,
             show_column_names = FALSE,  # Too many genomes (472) to show names
             cluster_columns = FALSE, 
             cluster_rows = TRUE,        # This will group genes with similar P/A patterns
             show_row_names = FALSE,     # Hide names if there are >1000 genes to avoid overlap
             row_title = paste(nrow(mat_all), "Accessory Genes"),
             column_title = "Complete Accessory Pangenome Signatures in S. Infantis"))

# 5. Save the large Heatmap to PDF
# Note: Increased height to 20 inches to give rows more room to breathe
pdf("Sweeps_All_Genes_Heatmap.pdf", width = 12, height = 20)

draw(Heatmap(mat_all, 
             name = "Presence", 
             col = c("0" = "white", "1" = "black"),
             top_annotation = col_anno_all,
             show_column_names = FALSE,  # Too many genomes (472) to show names
             cluster_columns = FALSE, 
             cluster_rows = TRUE,        # This will group genes with similar P/A patterns
             show_row_names = FALSE,     # Hide names if there are >1000 genes to avoid overlap
             row_title = paste(nrow(mat_all), "Accessory Genes"),
             column_title = "Complete Accessory Pangenome Signatures in S. Infantis"))

dev.off()

# --------------------------------------------------------------------------------------------------
# Heatmap of presence/absence of top enriched genes across genomes, grouped by ALL sweeps
# 1. Select the top 100 genes most enriched in ANY sweep 
# this will be basically all genes that are enriched in at least one sweep from the fishers test before
# and then look for those genes in the genomes, check if present and plot
top_genes <- enriched %>% group_by(group) %>% slice_max(odds_ratio, n = 100) %>% pull(Preferred_name) %>% unique()

# 2. Prepare the matrix (Genes as rows, Genomes as columns)
mat <- pa_matrix %>% filter(Preferred_name %in% top_genes) %>% tibble::column_to_rownames("Preferred_name")
mat <- as.matrix(mat)

# 3. Prepare column annotations (the sweep labels)
anno_df <- sweep_labels %>% filter(genome %in% colnames(mat)) %>% arrange(sweep)
mat <- mat[, anno_df$genome] # Sort matrix columns to match annotation
col_anno <- HeatmapAnnotation(Sweep = anno_df$sweep, 
                               col = list(Sweep = c("sweep_1" = "red", "sweep_7" = "blue", "no_sweep" = "gray")))

# Define a full color palette for all your groups
sweep_colors <- c(
  "sweep_1"  = "#E41A1C", 
  "sweep_2"  = "#377EB8", 
  "sweep_3"  = "#4DAF4A", 
  "sweep_4"  = "#984EA3", 
  "sweep_5"  = "#FF7F00", 
  "sweep_6"  = "#FFFF33", 
  "sweep_7"  = "#A65628", 
  "no_sweep" = "#999999"
)

# Update the annotation with the full list
col_anno <- HeatmapAnnotation(
  Sweep = anno_df$sweep, 
  col = list(Sweep = sweep_colors)
)

# Now run the Heatmap command
Heatmap(mat, 
        name = "Presence", 
        col = c("0" = "white", "1" = "black"),
        top_annotation = col_anno,
        show_column_names = FALSE, 
        cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 8),
        column_title = "Sweep-Specific Gene Signatures in S. Infantis")

pdf("Sweeps_all_Heatmap.pdf", width = 10, height = 8)
draw(Heatmap(mat, 
      name = "Presence", 
      col = c("0" = "white", "1" = "black"),
      top_annotation = col_anno,
      show_column_names = FALSE, 
      cluster_columns = FALSE, 
      row_names_gp = gpar(fontsize = 8),
      column_title = "Sweep-Specific Gene Signatures in S. Infantis")
)
dev.off()

# --------------------------------------------------------------------------------------------------
# Heatmap of presence/absence of top enriched genes across genomes, grouped by ALL sweep
# 1. Select the top 10 genes most enriched in ANY sweep 
# this will select top genes that are enriched in at least one sweep from the fishers test before
# and use this list to plot the heatmap
top_genes_sub <- enriched %>% group_by(group) %>% slice_max(odds_ratio, n = 10) %>% pull(Preferred_name) %>% unique()

# 2. Prepare the matrix (Genes as rows, Genomes as columns)
mat <- pa_matrix %>% filter(Preferred_name %in% top_genes_sub) %>% tibble::column_to_rownames("Preferred_name")
mat <- as.matrix(mat)

# 3. Prepare column annotations (the sweep labels)
anno_df <- sweep_labels %>% filter(genome %in% colnames(mat)) %>% arrange(sweep)
mat <- mat[, anno_df$genome] # Sort matrix columns to match annotation
col_anno <- HeatmapAnnotation(Sweep = anno_df$sweep, 
                               col = list(Sweep = c("sweep_1" = "red", "sweep_7" = "blue", "no_sweep" = "gray")))

# Define a full color palette for all your groups
sweep_colors <- c(
  "sweep_1"  = "#E41A1C", 
  "sweep_2"  = "#377EB8", 
  "sweep_3"  = "#4DAF4A", 
  "sweep_4"  = "#984EA3", 
  "sweep_5"  = "#FF7F00", 
  "sweep_6"  = "#FFFF33", 
  "sweep_7"  = "#A65628", 
  "no_sweep" = "#999999"
)

# Update the annotation with the full list
col_anno <- HeatmapAnnotation(
  Sweep = anno_df$sweep, 
  col = list(Sweep = sweep_colors)
)

# Now run the Heatmap command
Heatmap(mat, 
        name = "Presence", 
        col = c("0" = "white", "1" = "black"),
        top_annotation = col_anno,
        show_column_names = FALSE, 
        cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 8),
        column_title = "Sweep-Specific Gene Signatures in S. Infantis")

pdf("Sweeps_top10_Heatmap.pdf", width = 10, height = 8)
draw(Heatmap(mat, 
      name = "Presence", 
      col = c("0" = "white", "1" = "black"),
      top_annotation = col_anno,
      show_column_names = FALSE, 
      cluster_columns = FALSE, 
      row_names_gp = gpar(fontsize = 8),
      column_title = "Sweep-Specific Gene Signatures in S. Infantis")
)
dev.off()

# --------------------------------------------------------------------------------------------------
# Heatmap for only sweeps 1-6, with more genes (not just top 50), but still manageable 
# Define the target sweeps
target_sweeps <- c("sweep_1", "sweep_2", "sweep_3", "sweep_4", "sweep_5", "sweep_6" )

# Filter the annotation data for only these sweeps
anno_subset <- sweep_labels %>% 
  filter(sweep %in% target_sweeps) %>% 
  arrange(sweep)

# Subset the presence/absence matrix to these genomes
# We use the Preferred_name from your existing pa_matrix
mat_subset <- pa_matrix %>% 
  filter(Preferred_name %in% top_genes) %>% # Using the top genes identified earlier
  tibble::column_to_rownames("Preferred_name") %>%
  dplyr::select(all_of(anno_subset$genome)) %>%
  as.matrix()

# 4. Filter out any genes that might be 0 across all these specific genomes
mat_subset <- mat_subset[rowSums(mat_subset) > 0, ]

# 5. Define colors for only these 6 sweeps
sweep_colors_subset <- c(
  "sweep_1" = "#E41A1C", 
  "sweep_2" = "#377EB8", 
  "sweep_3" = "#4DAF4A", 
  "sweep_4" = "#984EA3", 
  "sweep_5" = "#FF7F00", 
  "sweep_6" = "#FFFF33"
)

# 6. Create the new annotation
col_anno_subset <- HeatmapAnnotation(
  Sweep = anno_subset$sweep, 
  col = list(Sweep = sweep_colors_subset),
  show_annotation_name = TRUE
)

# 7. Draw the refined Heatmap
Heatmap(mat_subset, 
        name = "Presence", 
        col = c("0" = "white", "1" = "black"),
        top_annotation = col_anno_subset,
        show_column_names = TRUE, # Enabled because there are only ~24 genomes now
        column_names_gp = gpar(fontsize = 7),
        cluster_columns = FALSE, # Keep them grouped by sweep
        cluster_rows = TRUE,     # Cluster genes to see functional blocks
        row_names_gp = gpar(fontsize = 8),
        column_title = "Gene Signatures Specific to Sweeps 1-6 (Salmonella Infantis)")

# Save as PDF
pdf("Sweeps_1-6_Heatmap.pdf", width = 10, height = 12)
draw(Heatmap(mat_subset, 
             name = "Presence", 
             col = c("0" = "white", "1" = "black"),
             top_annotation = col_anno_subset,
             show_column_names = TRUE, 
             cluster_columns = FALSE, 
             cluster_rows = TRUE,
             row_names_gp = gpar(fontsize = 8),
             column_title = "Gene Signatures Specific to Sweeps 1-6"))
dev.off()
# --------------------------------------------------------------------------------------------------






# 1. Define the function to create a labeled volcano plot for a single group
make_volcano <- function(grp_name, data) {
  
  # Filter and prepare data for this specific group
  plot_df <- data %>% 
    filter(group == grp_name) %>%
    mutate(
      log10p = -log10(p_adjusted),
      log2OR = log2(as.numeric(odds_ratio) + 0.01),
      significance = ifelse(p_adjusted < 0.05 & abs(log2OR) > 1, "Significant", "Not Significant")
    )
  
  # Identify top 5 genes for labeling (to keep the composite plot clean)
  top_labels <- plot_df %>%
    filter(significance == "Significant") %>%
    slice_max(log2OR, n = 5)
  
  # Build the plot
  p <- ggplot(plot_df, aes(x = log2OR, y = log10p, color = significance)) +
    geom_point(alpha = 0.3, size = 0.8) +
    theme_minimal(base_size = 8) + # Smaller text for a composite plot
    scale_color_manual(values = c("gray80", "firebrick"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "blue") +
    geom_text_repel(
      data = top_labels,
      aes(label = Preferred_name),
      size = 2,
      box.padding = 0.3,
      segment.size = 0.2
    ) +
    labs(title = grp_name, x = NULL, y = NULL)
  
  return(p)
}
x
# 2. Get the list of all groups
all_groups <- c("sweep_1", "sweep_2", "sweep_3", "sweep_4", 
                "sweep_5", "sweep_6", "sweep_7", "no_sweep")

# 3. Generate all plots using lapply
plot_list <- lapply(all_groups, make_volcano, data = fisher_results)

# 4. Use patchwork to concatenate them
# wrap_plots combines the list, and we specify the layout
composite_volcano <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(
    title = "Comparative Functional Enrichment across Salmonella Infantis Sweeps",
    subtitle = "Red points indicate significant genes (p_adj < 0.05, Log2OR > 1)",
    caption = "Data source: Prokka + EggNOG-mapper + BLAST against 472 genomes"
  )

# 5. Display and Save
print(composite_volcano)

# Save as a large PDF for your thesis
ggsave("Combined_Volcano_Plots.pdf", composite_volcano, width = 10, height = 14)

