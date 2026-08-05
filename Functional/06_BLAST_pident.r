.libPaths("~/R/library")

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

setwd("/home/senekowitsch/Thesis/Functional/05_analyze_BLAST")

# load filtered data
all_data_filtered <- fread("all_data_filtered.tsv.gz")
n_distinct(all_data_filtered$sseqid) # 472 distinct target genomes with hits after filtering

# load sweep labels to get sweep group sizes
sweep_labels <- read.delim("/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt",
                           sep = "\t", header = TRUE)
sweep_labels$genome <- as.character(sweep_labels$genome)

# total number of genomes
n_total <- nrow(sweep_labels)  # 472
n_total

# size of each sweep group (used as denominators for freq calculations)
sweep_sizes <- sweep_labels %>%
  group_by(sweep) %>%
  summarise(n = n(), .groups = "drop")

setwd("/home/senekowitsch/Thesis/Functional/06_BLAST_pident")
# -------------------------------------------------
# Classify each gene hit as inside or outside the own sweep
# -------------------------------------------------
all_data_filtered <- all_data_filtered %>%
  mutate(hit_location = ifelse(target_sweep == query_group, "inside", "outside"))

# -------------------------------------------------
# Exclude unannotated genes
# -------------------------------------------------
data_annotated <- all_data_filtered %>%
  filter(
    !is.na(COG_category),
    COG_category != "",
    COG_category != "-"
  )
head(data_annotated)

# -------------------------------------------------
# Overall inside summary:
# median pident + distinct genomes hit inside sweep
# -------------------------------------------------
# keeps only "inside" hits, groups by gene (Preferred_name) and the query's sweep (query_group)
# for each gene/sweep combination computes the median/mean/sd of percent identity
# n_genomes_hit_inside: the number of distinct target genomes that gene matched within its own sweep
# then joins sweep_sizes to get n, the total genome count in that sweep, and divides to get 
# fraction of genomes in the sweep that carry a hit for this gene.
# So each row answers: "within sweep X, how identical and how prevalent is gene G?"
# -------------------------------------------------
inside_summary <- data_annotated %>%
  # keep only hits inside the sweep for this summary
  filter(hit_location == "inside") %>%
  # group by gene annotation + query sweep group
  # so keep 1 chunk for gene1 in sweepA, 1 chunk for gene1 in sweepB, etc.
  group_by(Preferred_name, query_group) %>%
  # collapse each chunk to one row containing:
  summarise(
    median_pident_inside = median(pident),
    mean_pident_inside   = mean(pident),
    sd_pident_inside     = sd(pident),
    # how many unique target genomes have a hit, because we have 3 query genomes per sweep
    n_genomes_hit_inside = n_distinct(sseqid),  # distinct genomes, not raw hits
    .groups = "drop"
  ) %>%
  # join sweep size (from outside summary) to compute freq_inside
  left_join(sweep_sizes, by = c("query_group" = "sweep")) %>%
  # number of distinct genomes hit by genomes in this sweep group
  mutate(freq_inside = n_genomes_hit_inside / n) %>%
  dplyr::select(-n)

# -------------------------------------------------
# Overall outside summary:
# median pident + distinct genomes hit outside sweep
# -------------------------------------------------
# for hits landing outside the query's own sweep, lumped together regardless of which other sweep they hit.
# Same median/mean/sd/pident logic, same distinct-genome counting. 
# since "outside" means "everywhere else." This gives one overall outside-identity number per gene/sweep
# ignoring which specific sweep the outside hits came from.
# -------------------------------------------------
outside_summary <- data_annotated %>%
  filter(hit_location == "outside") %>%
  group_by(Preferred_name, query_group) %>%
  summarise(
    median_pident_outside = median(pident),
    mean_pident_outside   = mean(pident),
    sd_pident_outside     = sd(pident),
    n_genomes_hit_outside = n_distinct(sseqid),
    .groups = "drop"
  ) %>%
  # denominator = total genomes minus the query sweep size
  left_join(sweep_sizes, by = c("query_group" = "sweep")) %>%
  mutate(freq_outside = n_genomes_hit_outside / (n_total - n)) %>%
  # remove n column from join (comes from sweep_sizes)
  dplyr::select(-n)

# -------------------------------------------------
# Per-target-sweep outside breakdown:
# median pident + freq hit for each outside sweep group
# ------------------------------------------------
#nthis is the more granular version of the outside summary
# instead of collapsing all outside sweeps into one number, it groups by gene, query sweep
# and target_sweep, so a gene from sweep_1 hitting sweep_2 is tracked separately from that same gene hitting sweep_3.
# -------------------------------------------------
per_sweep_outside <- data_annotated %>%
  filter(hit_location == "outside") %>%
  # one row per Preferred_name + query_group + target_sweep
  # so gene1 from sweepA hitting sweepB is separate from gene1 from sweepA hitting sweepC, etc.
  group_by(Preferred_name, query_group, target_sweep) %>%
  summarise(
    median_pident    = median(pident),
    n_genomes_hit    = n_distinct(sseqid),
    .groups = "drop"
  ) %>%
  # join target sweep size for freq denominator
  # here with target_sweep because we want to know the size of the sweep we're hitting, not the query sweep
  left_join(sweep_sizes, by = c("target_sweep" = "sweep")) %>%
  mutate(freq = n_genomes_hit / n) %>%
  dplyr::select(-n) %>%
  # rotates target_sweep values into columns
  # one row per gene and query_group, with separate columns for each target_sweep's median_pident, freq, and n_genomes_hit
  pivot_wider(
    names_from  = target_sweep,
    values_from = c(median_pident, freq, n_genomes_hit),
    # controls how the new column names are formed, e.g. median_pident_sweepA, freq_sweepA, etc.
    names_glue  = "{.value}_{target_sweep}"
  )

# -------------------------------------------------
# Join everything together
# -------------------------------------------------
pident_full <- inside_summary %>%
  left_join(outside_summary,   by = c("Preferred_name", "query_group")) %>%
  left_join(per_sweep_outside, by = c("Preferred_name", "query_group"))

# -------------------------------------------------
# Save full table
# -------------------------------------------------
write.table(pident_full, "pident_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Full summary written: pident_summary.tsv\n")
cat("Dimensions:", nrow(pident_full), "rows x", ncol(pident_full), "cols\n")

# -------------------------------------------------
# Filter for interesting cases:
# inside pident == 100, freq_inside >= 0.9,
# outside pident < 100
# -------------------------------------------------
interesting <- pident_full %>%
  filter(
    median_pident_inside >= 99.999,
    freq_inside >= 0.9,
    !is.na(median_pident_outside),
    median_pident_outside < 100
  ) %>%
  arrange(median_pident_outside) #sorts by median_pident_outside ascending

cat("Interesting cases (inside=100, freq_inside>=0.9, outside<100):", nrow(interesting), "\n")
head(interesting)

write.table(interesting, "pident_interesting.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Interesting cases written: pident_interesting.tsv\n")

# -------------------------------------------------
# Genes with NO hits outside the sweep
# (truly sweep-specific at the BLAST level)
# also require freq_inside >= 0.9
# -------------------------------------------------
sweep_specific <- pident_full %>%
  filter(
    median_pident_inside >= 99.999,
    freq_inside >= 0.9,
    is.na(median_pident_outside)
  )

cat("Genes with hits only inside sweep (no outside hits):", nrow(sweep_specific), "\n")

write.table(sweep_specific, "pident_sweep_specific.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Sweep-specific genes written: pident_sweep_specific.tsv\n")

# -------------------------------------------------
# Quick sanity checks
# -------------------------------------------------

# distribution of inside pident - should be overwhelmingly 100
cat("\nInside pident distribution:\n")
print(table(cut(pident_full$median_pident_inside,
                breaks = c(0, 90, 95, 99, 99.9, 100),
                include.lowest = TRUE)))

# distribution of freq_inside - most genes should be high
cat("\nInside frequency distribution:\n")
print(table(cut(pident_full$freq_inside,
                breaks = c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0),
                include.lowest = TRUE)))

# distribution of outside pident for interesting cases
cat("\nOutside pident distribution (interesting cases only):\n")
print(table(cut(interesting$median_pident_outside,
                breaks = c(0, 70, 80, 90, 95, 99, 99.9, 100),
                include.lowest = TRUE)))

# distribution of freq_outside for interesting cases
cat("\nOutside frequency distribution (interesting cases only):\n")
print(table(cut(interesting$freq_outside,
                breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0),
                include.lowest = TRUE)))

# -------------------------------------------------
# COG category breakdown of interesting cases
# -------------------------------------------------
 
# extract distinct Preferred_name x COG_split mapping
# split multi-letter COG annotations (e.g. "LV") into one row per letter
# so a gene annotated as "LV" appears as both "L" and "V"
cog_map <- data_annotated %>%
  filter(!is.na(COG_category), COG_category != "", COG_category != "-") %>%
  distinct(Preferred_name, COG_category) %>%
  mutate(COG_split = strsplit(COG_category, "")) %>%
  tidyr::unnest(COG_split) %>%
  distinct(Preferred_name, COG_split)
 
# join COG onto interesting cases
# genes with multiple COG letters will get one row per letter
interesting_cog <- interesting %>%
  left_join(cog_map, by = "Preferred_name") %>%
  arrange(COG_split, median_pident_outside)
 
write.table(interesting_cog, "pident_interesting_cog.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Interesting cases with COG written: pident_interesting_cog.tsv\n")
 
# summary: how many interesting genes per COG_split x query_group?
# useful to see if certain functional categories are overrepresented
cog_counts_interesting <- interesting_cog %>%
  filter(!is.na(COG_split)) %>%
  group_by(COG_split, query_group) %>%
  summarise(
    n_genes                       = n_distinct(Preferred_name),
    median_pident_outside_overall = median(median_pident_outside, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(COG_split, query_group)
 
write.table(cog_counts_interesting, "pident_cog_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("COG summary written: pident_cog_summary.tsv\n")
 
# same for sweep_specific genes
sweep_specific_cog <- sweep_specific %>%
  left_join(cog_map, by = "Preferred_name") %>%
  arrange(COG_split)
 
write.table(sweep_specific_cog, "pident_sweep_specific_cog.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Sweep-specific genes with COG written: pident_sweep_specific_cog.tsv\n")
 
# quick overview of COG distribution in interesting vs sweep_specific
cat("\nCOG category counts in interesting cases:\n")
print(sort(table(interesting_cog$COG_split), decreasing = TRUE))
 
cat("\nCOG category counts in sweep-specific genes:\n")
print(sort(table(sweep_specific_cog$COG_split), decreasing = TRUE))

# -------------------------------------------------
# Visualization: pident inside vs outside per COG
# -------------------------------------------------
# for sweep_1
# filter to sweep_1 and pivot to long format for plotting
# so each gene has two rows: one for inside, one for outside
plot_data <- interesting_cog %>%
  filter(query_group == "sweep_1", !is.na(COG_split)) %>%
  dplyr::select(Preferred_name, COG_split, median_pident_inside, median_pident_outside) %>%
  pivot_longer(
    cols      = c(median_pident_inside, median_pident_outside),
    names_to  = "location",
    values_to = "median_pident"
  ) %>%
  mutate(location = recode(location,
    "median_pident_inside"  = "inside",
    "median_pident_outside" = "outside"
  ))

# plot
p <- ggplot(plot_data, aes(x = location, y = median_pident, colour = Preferred_name)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ COG_split, scales = "free_y") +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(
    title    = "Median pident inside vs outside sweep_1 per COG category",
    x        = NULL,
    y        = "Median pident (%)",
    colour   = "Gene"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  )
 
ggsave("pident_sweep1_cog.pdf", plot = p, width = 14, height = 10)
cat("Plot saved: pident_sweep1_cog.pdf\n")

# -------------------------------------------------
# for sweep_2
# filter to sweep_2 and pivot to long format for plotting
# so each gene has two rows: one for inside, one for outside
plot_data <- interesting_cog %>%
  filter(query_group == "sweep_2", !is.na(COG_split)) %>%
  dplyr::select(Preferred_name, COG_split, median_pident_inside, median_pident_outside) %>%
  pivot_longer(
    cols      = c(median_pident_inside, median_pident_outside),
    names_to  = "location",
    values_to = "median_pident"
  ) %>%
  mutate(location = recode(location,
    "median_pident_inside"  = "inside",
    "median_pident_outside" = "outside"
  ))

# plot
p <- ggplot(plot_data, aes(x = location, y = median_pident, colour = Preferred_name)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ COG_split, scales = "free_y") +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(
    title    = "Median pident inside vs outside sweep_2 per COG category",
    x        = NULL,
    y        = "Median pident (%)",
    colour   = "Gene"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  )

ggsave("pident_sweep2_cog.pdf", plot = p, width = 14, height = 10)
cat("Plot saved: pident_sweep2_cog.pdf\n")

# -------------------------------------------------
# for sweep_3
# filter to sweep_3 and pivot to long format for plotting
# so each gene has two rows: one for inside, one for outside
plot_data <- interesting_cog %>%
  filter(query_group == "sweep_3", !is.na(COG_split)) %>%
  dplyr::select(Preferred_name, COG_split, median_pident_inside, median_pident_outside) %>%
  pivot_longer(
    cols      = c(median_pident_inside, median_pident_outside),
    names_to  = "location",
    values_to = "median_pident"
  ) %>%
  mutate(location = recode(location,
    "median_pident_inside"  = "inside",
    "median_pident_outside" = "outside"
  ))

# plot
p <- ggplot(plot_data, aes(x = location, y = median_pident, colour = Preferred_name)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ COG_split, scales = "free_y") +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(
    title    = "Median pident inside vs outside sweep_3 per COG category",
    x        = NULL,
    y        = "Median pident (%)",
    colour   = "Gene"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  )
 
ggsave("pident_sweep3_cog.pdf", plot = p, width = 14, height = 10)
cat("Plot saved: pident_sweep3_cog.pdf\n")

# -------------------------------------------------
# for sweep_4
# filter to sweep_4 and pivot to long format for plotting
# so each gene has two rows: one for inside, one for outside
plot_data <- interesting_cog %>%
  filter(query_group == "sweep_4", !is.na(COG_split)) %>%
  dplyr::select(Preferred_name, COG_split, median_pident_inside, median_pident_outside) %>%
  pivot_longer(
    cols      = c(median_pident_inside, median_pident_outside),
    names_to  = "location",
    values_to = "median_pident"
  ) %>%
  mutate(location = recode(location,
    "median_pident_inside"  = "inside",
    "median_pident_outside" = "outside"
  ))

# plot
p <- ggplot(plot_data, aes(x = location, y = median_pident, colour = Preferred_name)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ COG_split, scales = "free_y") +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(
    title    = "Median pident inside vs outside sweep_4 per COG category",
    x        = NULL,
    y        = "Median pident (%)",
    colour   = "Gene"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  )

ggsave("pident_sweep4_cog.pdf", plot = p, width = 14, height = 10)
cat("Plot saved: pident_sweep4_cog.pdf\n")

# -------------------------------------------------
# for sweep_5
# filter to sweep_5 and pivot to long format for plotting
# so each gene has two rows: one for inside, one for outside
plot_data <- interesting_cog %>%
  filter(query_group == "sweep_5", !is.na(COG_split)) %>%
  dplyr::select(Preferred_name, COG_split, median_pident_inside, median_pident_outside) %>%
  pivot_longer(
    cols      = c(median_pident_inside, median_pident_outside),
    names_to  = "location",
    values_to = "median_pident"
  ) %>%
  mutate(location = recode(location,
    "median_pident_inside"  = "inside",
    "median_pident_outside" = "outside"
  ))

# plot
p <- ggplot(plot_data, aes(x = location, y = median_pident, colour = Preferred_name)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ COG_split, scales = "free_y") +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(
    title    = "Median pident inside vs outside sweep_5 per COG category",
    x        = NULL,
    y        = "Median pident (%)",
    colour   = "Gene"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  )

ggsave("pident_sweep5_cog.pdf", plot = p, width = 14, height = 10)
cat("Plot saved: pident_sweep5_cog.pdf\n")

# -------------------------------------------------
# for sweep_6
# filter to sweep_6 and pivot to long format for plotting
# so each gene has two rows: one for inside, one for outside
plot_data <- interesting_cog %>%
  filter(query_group == "sweep_6", !is.na(COG_split)) %>%
  dplyr::select(Preferred_name, COG_split, median_pident_inside, median_pident_outside) %>%
  pivot_longer(
    cols      = c(median_pident_inside, median_pident_outside),
    names_to  = "location",
    values_to = "median_pident"
  ) %>%
  mutate(location = recode(location,
    "median_pident_inside"  = "inside",
    "median_pident_outside" = "outside"
  ))

# plot
p <- ggplot(plot_data, aes(x = location, y = median_pident, colour = Preferred_name)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ COG_split, scales = "free_y") +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(
    title    = "Median pident inside vs outside sweep_6 per COG category",
    x        = NULL,
    y        = "Median pident (%)",
    colour   = "Gene"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  )

ggsave("pident_sweep6_cog.pdf", plot = p, width = 14, height = 10)
cat("Plot saved: pident_sweep6_cog.pdf\n")

# -------------------------------------------------
# for sweep_7
# filter to sweep_7 and pivot to long format for plotting
# so each gene has two rows: one for inside, one for outside
plot_data <- interesting_cog %>%
  filter(query_group == "sweep_7", !is.na(COG_split)) %>%
  dplyr::select(Preferred_name, COG_split, median_pident_inside, median_pident_outside) %>%
  pivot_longer(
    cols      = c(median_pident_inside, median_pident_outside),
    names_to  = "location",
    values_to = "median_pident"
  ) %>%
  mutate(location = recode(location,
    "median_pident_inside"  = "inside",
    "median_pident_outside" = "outside"
  ))

# plot
p <- ggplot(plot_data, aes(x = location, y = median_pident, colour = Preferred_name)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ COG_split, scales = "free_y") +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(
    title    = "Median pident inside vs outside sweep_7 per COG category",
    x        = NULL,
    y        = "Median pident (%)",
    colour   = "Gene"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  )

ggsave("pident_sweep7_cog.pdf", plot = p, width = 14, height = 10)
cat("Plot saved: pident_sweep7_cog.pdf\n")


# -------------------------------------------------
# for no_sweep
# filter to no_sweep and pivot to long format for plotting
# so each gene has two rows: one for inside, one for outside
plot_data <- interesting_cog %>%
  filter(query_group == "no_sweep", !is.na(COG_split)) %>%
  dplyr::select(Preferred_name, COG_split, median_pident_inside, median_pident_outside) %>%
  pivot_longer(
    cols      = c(median_pident_inside, median_pident_outside),
    names_to  = "location",
    values_to = "median_pident"
  ) %>%
  mutate(location = recode(location,
    "median_pident_inside"  = "inside",
    "median_pident_outside" = "outside"
  ))

# plot
p <- ggplot(plot_data, aes(x = location, y = median_pident, colour = Preferred_name)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ COG_split, scales = "free_y") +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(
    title    = "Median pident inside vs outside no_sweep per COG category",
    x        = NULL,
    y        = "Median pident (%)",
    colour   = "Gene"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  )

ggsave("pident_no_sweep_cog.pdf", plot = p, width = 14, height = 10)
cat("Plot saved: pident_no_sweep_cog.pdf\n")



sweep1_genes <- interesting_cog %>% filter(query_group == "sweep_1") %>% pull(Preferred_name) %>% unique()
sweep2_genes <- interesting_cog %>% filter(query_group == "sweep_2") %>% pull(Preferred_name) %>% unique()
sweep3_genes <- interesting_cog %>% filter(query_group == "sweep_3") %>% pull(Preferred_name) %>% unique()
sweep4_genes <- interesting_cog %>% filter(query_group == "sweep_4") %>% pull(Preferred_name) %>% unique()
sweep5_genes <- interesting_cog %>% filter(query_group == "sweep_5") %>% pull(Preferred_name) %>% unique()
sweep6_genes <- interesting_cog %>% filter(query_group == "sweep_6") %>% pull(Preferred_name) %>% unique()
sweep7_genes <- interesting_cog %>% filter(query_group == "sweep_7") %>% pull(Preferred_name) %>% unique()

all_sweep_genes <- list(
  sweep_1 = sweep1_genes,
  sweep_2 = sweep2_genes,
  sweep_3 = sweep3_genes,
  sweep_4 = sweep4_genes,
  sweep_5 = sweep5_genes,
  sweep_6 = sweep6_genes,
  sweep_7 = sweep7_genes
)

# all unique genes across all sweeps
all_genes_union <- unique(unlist(all_sweep_genes))

# presence/absence matrix
pa <- data.frame(
  Preferred_name = all_genes_union
)

for (grp in names(all_sweep_genes)) {
  pa[[grp]] <- as.integer(all_genes_union %in% all_sweep_genes[[grp]])
}

# how many sweeps is each gene interesting in?
pa$n_sweeps <- rowSums(pa[, names(all_sweep_genes)])

# sort by how broadly shared they are
pa <- pa %>% arrange(desc(n_sweeps))

head(pa)

# genes interesting in all 7 sweeps
pa %>% filter(n_sweeps == 7) %>% pull(Preferred_name)

# genes unique to sweep_7 only
pa %>% filter(n_sweeps == 1, sweep_1 == 1) %>% pull(Preferred_name)


# =================================================
# Visualization: pident inside vs outside per COG
# Replaces the per-sweep copy-paste blocks.
# One function, looped over every query group.
# Shared, fixed y-axis so all plots are comparable.
# =================================================

# Shared y-axis range, computed ONCE across all groups and both
# locations, so every output PDF uses identical limits.
y_min <- interesting_cog %>%
  filter(!is.na(COG_split)) %>%
  summarise(m = min(c(median_pident_inside, median_pident_outside),
                    na.rm = TRUE)) %>%
  pull(m)

y_limits <- c(floor(y_min), 100)   # same on every plot

plot_pident_by_cog <- function(group,
                               data = interesting_cog,
                               ylim = y_limits) {

  pd <- data %>%
    filter(query_group == group, !is.na(COG_split)) %>%
    dplyr::select(Preferred_name, COG_split,
                  median_pident_inside, median_pident_outside) %>%
    pivot_longer(
      cols      = c(median_pident_inside, median_pident_outside),
      names_to  = "location",
      values_to = "median_pident"
    ) %>%
    mutate(location = recode(location,
      "median_pident_inside"  = "inside",
      "median_pident_outside" = "outside"
    ))

  # some groups may have no interesting genes -> skip instead of erroring
  if (nrow(pd) == 0) {
    message("No interesting genes for ", group, " - skipping plot")
    return(invisible(NULL))
  }

  p <- ggplot(pd, aes(x = location, y = median_pident,
                      colour = Preferred_name)) +
    geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
    facet_wrap(~ COG_split) +                 # fixed scales -> facets comparable
    scale_y_continuous(limits = ylim) +       # identical y-axis on every plot
    labs(
      title  = paste0("Median pident inside vs outside ", group,
                      " per COG category"),
      x      = NULL,
      y      = "Median pident (%)",
      colour = "Gene"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      strip.text      = element_text(face = "bold")
    )

  outfile <- paste0("pident_", group, "_cog.pdf")
  ggsave(outfile, plot = p, width = 14, height = 10)
  cat("Plot saved:", outfile, "\n")
  invisible(p)
}

# Run for every query group
groups <- c("sweep_1", "sweep_2", "sweep_3", "sweep_4",
            "sweep_5", "sweep_6", "sweep_7", "no_sweep")

for (g in groups) plot_pident_by_cog(g)