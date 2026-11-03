#!/bin/bash
source /home/senekowitsch/miniconda3/etc/profile.d/conda.sh
THREADS=30

# Variables
base_dir='/home/senekowitsch/Thesis/Functional_2/03_clade_timeing'
output_base="${base_dir}/output"

# Inputs
assemblies_infantis="/home/senekowitsch/Thesis/Functional_2/00_data"  # *.fna
metadata_NCBI="/home/senekowitsch/Thesis/Functional_2/03_clade_timeing/infantis_clean_genome_collection.csv"
sampling_years_all="/home/senekowitsch/Thesis/Functional_2/03_clade_timeing/tips_years_all.tsv"

# Group/lineage lists (reuse your existing ones)
LIST_SWEEP7="/home/senekowitsch/Thesis/Functional_2/sweep7_genomes.txt"
# for tree visualization
infantis_tree="/home/senekowitsch/Thesis/Functional_2/01_plasmid_id/output/tree/full_tree_renamed.treefile"
infantis_tree_pesi="/home/senekowitsch/Salmonella/AGES/Functional/03_pESI_tree/output/infantis/tree/infantis_pESI_tree.treefile"

# ------------------------------------------
# Create directories
# ------------------------------------------
mkdir -p "${base_dir}"
mkdir -p "${output_base}"

cd "${base_dir}"


# Output: tab-separated "tip<TAB>year" (year = NA if missing)
# infantis_clean_genome_collection.csv (quoted CSV; AssemblyAccession = col 1, CollectionDate = col 4)
tail -n +2 $metadata_NCBI | tr -d '\r' |
awk '{
  n=0; f=""; q=0                                  # quote-aware CSV split
  for (i=1; i<=length($0); i++) {
    c=substr($0,i,1)
    if (c=="\"") q=!q
    else if (c=="," && !q) { a[++n]=f; f="" }
    else f=f c
  }
  a[++n]=f
  acc=a[1]; sub(/^GC[AF]_0*/,"",acc); sub(/\.[0-9]+$/,"",acc)   # GCA_009397805.1 -> 9397805
  y=(a[4] ~ /^[0-9][0-9][0-9][0-9]/) ? substr(a[4],1,4) : "NA"  # YYYY / YYYY-MM / YYYY-MM-DD
  print acc "\t" y
}' | sort -u > tips_years_collection.tsv

cat tips_years_collection.tsv > tips_years_all.tsv


# ------------------------------------------
# Sampling year: Sweep7 vs. other NCBI tips in the tree
# ------------------------------------------

# Tip names from the tree
grep -oE '[(,][^(),:;]+:' "$infantis_tree" | tr -d '(,:' | sort -u > tree_tips.txt

# Normalise Sweep7 IDs (drops path, .fna/.fasta, GCA_/GCF_ + zeros, version)
sed -E 's#.*/##; s/\r//; s/\.(fna|fa|fasta|gff)$//; s/^GC[AF]_0*//; s/\.[0-9]+$//' \
  "$LIST_SWEEP7" | sort -u > sweep7_ids.txt

# NCBI tips that are in the tree, with year and group
awk -F'\t' -v OFS='\t' '
  FILENAME==ARGV[1] {tree[$1]; next}
  FILENAME==ARGV[2] {s7[$1];   next}
  ($1 in tree) && $2!="NA" { print $1, $2, (($1 in s7) ? "Sweep7" : "Other") }
' tree_tips.txt sweep7_ids.txt tips_years_collection.tsv > ncbi_tree_years_groups.tsv

echo "Group sizes:"; cut -f3 ncbi_tree_years_groups.tsv | sort | uniq -c
echo "Sweep7 IDs not found among NCBI tree tips:"
cut -f1 ncbi_tree_years_groups.tsv | sort | comm -23 sweep7_ids.txt - | head

# ------------------------------------------
# Statistics (printed + saved)
# ------------------------------------------
export output_base
Rscript - <<'EOF' 2>&1 | tee "${output_base}/sampling_year_sweep7_vs_other.txt"
d <- read.delim("ncbi_tree_years_groups.tsv", header = FALSE,
                col.names = c("tip", "year", "group"),
                colClasses = c("character", "integer", "character"))
print(tapply(d$year, d$group, summary))
w <- wilcox.test(year ~ group, data = d, conf.int = TRUE)
k <- ks.test(d$year[d$group == "Sweep7"], d$year[d$group == "Other"])
print(w); print(k)
 
# Compact table of the key numbers
res <- data.frame(
  n_Sweep7      = sum(d$group == "Sweep7"),
  n_Other       = sum(d$group == "Other"),
  median_Sweep7 = median(d$year[d$group == "Sweep7"]),
  median_Other  = median(d$year[d$group == "Other"]),
  shift_years   = round(w$estimate, 2),    # Other minus Sweep7
  CI_low        = round(w$conf.int[1], 2),
  CI_high       = round(w$conf.int[2], 2),
  wilcox_W      = w$statistic,
  wilcox_p      = signif(w$p.value, 3),
  ks_D          = round(k$statistic, 3),
  ks_p          = signif(k$p.value, 3))
write.table(res, file.path(Sys.getenv("output_base"), "sampling_year_sweep7_vs_other_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
EOF


# ------------------------------------------
# Boxplot: sampling year Sweep7 vs. Other
# ------------------------------------------
Rscript - <<'EOF'
suppressPackageStartupMessages(library(ggplot2))
d <- read.delim("ncbi_tree_years_groups.tsv", header = FALSE,
                col.names = c("tip", "year", "group"),
                colClasses = c("character", "integer", "character"))
d$group <- factor(d$group, levels = c("Other", "Sweep7"))
w <- wilcox.test(year ~ group, data = d, conf.int = TRUE)
n <- table(d$group)
pval <- w$p.value
lab  <- if (pval < 2.2e-16) "p < 2.2e-16" else paste("p =", signif(pval, 2))
# stars <- cut(pval, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "ns"))  # alternative label
yb <- max(d$year) + 1      # height of the significance bracket
p <- ggplot(d, aes(group, year, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6) +
  geom_jitter(width = 0.15, height = 0.25, size = 1, alpha = 0.4) +
  scale_fill_manual(values = c(Other = "grey70", Sweep7 = "#D55E00"), guide = "none") +
  scale_x_discrete(labels = sprintf("%s\n(n = %d)", names(n), n)) +
  scale_y_continuous(breaks = seq(2000, 2026, 2)) +
  annotate("segment", x = c(1, 1, 2), xend = c(2, 1, 2),
           y = c(yb, yb, yb), yend = c(yb, yb - 0.4, yb - 0.4)) +
  annotate("text", x = 1.5, y = yb + 0.8, label = lab, size = 3.8) +
  labs(x = NULL, y = "Sampling year", title = "NCBI genomes in the phylogeny",
       caption = "Wilcoxon rank-sum test") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(size = 12, hjust = 0.5))
out <- file.path(Sys.getenv("output_base"), "sampling_year_sweep7_vs_other_boxplot")
ggsave(paste0(out, ".pdf"), p, width = 4.5, height = 6)
ggsave(paste0(out, ".png"), p, width = 4.5, height = 6, dpi = 300)
EOF





