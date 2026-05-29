.libPaths("~/R/library")
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
blast_dir  <- "/home/senekowitsch/Thesis/Functional/03_blast/output"
annot_dir  <- "/home/senekowitsch/Thesis/Functional/04_clean_BLAST/annotations"
labels_file <- "/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt"
output_file <- "/home/senekowitsch/Thesis/Functional/04_clean_BLAST/all_data.rds"

# ── Column definitions ─────────────────────────────────────────────────────────
blast_cols <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen",
                "qcovs", "evalue", "bitscore")

annot_cols <- c("qseqid", "seed_ortholog", "evalue_annot", "score",
                "eggNOG_OGs", "max_annot_lvl", "COG_category",
                "Description", "Preferred_name", "GOs", "EC",
                "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
                "KEGG_Reaction", "KEGG_rclass", "BRITE",
                "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")

keep_cols <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs",
               "evalue", "bitscore", "COG_category", "Description",
               "Preferred_name", "GOs", "KEGG_ko", "KEGG_Pathway", "CAZy",
               "PFAMs", "target_sweep", "query_genome", "query_group")

# ── Load sweep labels ──────────────────────────────────────────────────────────
sweep_labels <- read.delim(labels_file, sep = "\t", header = TRUE)
sweep_labels$genome <- as.character(sweep_labels$genome)

# ── Load and process all files ─────────────────────────────────────────────────
blast_files <- list.files(blast_dir, pattern = "_blast\\.tsv$", full.names = TRUE)

all_data <- lapply(blast_files, function(f) {
  
  fname        <- basename(f)
  genome_id    <- sub(".*_(\\d+)_blast\\.tsv$", "\\1", fname)
  query_group  <- sub("_\\d+_blast\\.tsv$", "", fname)
  
  cat("Processing:", fname, "\n")
  
  # load BLAST
  blast <- read.delim(f, sep = "\t", header = FALSE, col.names = blast_cols)
  
  # load annotation
  annot_file <- file.path(annot_dir, paste0(genome_id, ".emapper.annotations"))
  annot <- read.delim(annot_file, sep = "\t", header = FALSE, skip = 4,
                      col.names = annot_cols)
  annot <- head(annot, -3)  # remove last 3 summary lines
  
  # merge BLAST + annotations
  merged <- merge(blast, annot, by = "qseqid", all.x = TRUE)
  
  # clean sseqid to just the number
  merged$sseqid <- sub("\\|.*", "", merged$sseqid)
  
  # join sweep labels for target genome
  merged <- merge(merged, sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
  merged <- rename(merged, target_sweep = sweep)
  
  # add query genome identity
  merged$query_genome <- genome_id
  merged$query_group  <- query_group
  
  # subset to keep_cols only
  merged <- merged[, keep_cols]
  
  merged
})

# ── Combine and save ───────────────────────────────────────────────────────────
cat("Combining all data...\n")
all_data <- do.call(rbind, all_data)

cat("Saving to", output_file, "\n")
saveRDS(all_data, output_file)

cat("Done! Final dimensions:", nrow(all_data), "rows x", ncol(all_data), "cols\n")