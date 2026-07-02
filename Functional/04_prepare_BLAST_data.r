.libPaths("~/R/library")

library(dplyr)
library(data.table)

##set working directory
setwd("/home/senekowitsch/Thesis/Functional/03_blast/output")

# -------------------------------------------------
# Context
# -------------------------------------------------
# Identified genomes that were in sweeps and outside of sweeps
# Took 3 random genomes from each sweep and from outside of sweep
# Used Prokka to annotate those 24 genomes
# Merged all 473 genomes into one large file and used this to create blastdb
# Blasted all Prokka identified genes from the 24 genomes against this db
# Extracted these blast results, the prokka annotations (they have the same identifier)
# Extracted list with sweeps

# -------------------------------------------------
# Load data
# -------------------------------------------------
blast_cols <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "evalue", "bitscore")

# Load BLAST
no_sweep_11080585 <- read.delim("no_sweep_11080585_blast.tsv", sep = "\t", header = FALSE, col.names = blast_cols)
no_sweep_14738655 <- read.delim("no_sweep_14738655_blast.tsv", sep = "\t", header = FALSE, col.names = blast_cols)
no_sweep_8866805  <- read.delim("no_sweep_8866805_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_1_31832435  <- read.delim("sweep_1_31832435_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_1_32036575  <- read.delim("sweep_1_32036575_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_1_32037335  <- read.delim("sweep_1_32037335_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_2_45808885  <- read.delim("sweep_2_45808885_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_2_45808905  <- read.delim("sweep_2_45808905_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_2_45808925  <- read.delim("sweep_2_45808925_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_3_11453035  <- read.delim("sweep_3_11453035_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_3_8046645   <- read.delim("sweep_3_8046645_blast.tsv",   sep = "\t", header = FALSE, col.names = blast_cols)
sweep_3_8091485   <- read.delim("sweep_3_8091485_blast.tsv",   sep = "\t", header = FALSE, col.names = blast_cols)
sweep_4_20671165  <- read.delim("sweep_4_20671165_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_4_20671355  <- read.delim("sweep_4_20671355_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_4_20671385  <- read.delim("sweep_4_20671385_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_5_8026935   <- read.delim("sweep_5_8026935_blast.tsv",   sep = "\t", header = FALSE, col.names = blast_cols)
sweep_5_8927125   <- read.delim("sweep_5_8927125_blast.tsv",   sep = "\t", header = FALSE, col.names = blast_cols)
sweep_5_8928825   <- read.delim("sweep_5_8928825_blast.tsv",   sep = "\t", header = FALSE, col.names = blast_cols)
sweep_6_26733835  <- read.delim("sweep_6_26733835_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_6_8484365   <- read.delim("sweep_6_8484365_blast.tsv",   sep = "\t", header = FALSE, col.names = blast_cols)
sweep_6_8637495   <- read.delim("sweep_6_8637495_blast.tsv",   sep = "\t", header = FALSE, col.names = blast_cols)
sweep_7_16459765  <- read.delim("sweep_7_16459765_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_7_33015275  <- read.delim("sweep_7_33015275_blast.tsv",  sep = "\t", header = FALSE, col.names = blast_cols)
sweep_7_8813075   <- read.delim("sweep_7_8813075_blast.tsv",   sep = "\t", header = FALSE, col.names = blast_cols)

# Load annotation files
setwd("/home/senekowitsch/Thesis/Functional/04_clean_BLAST/annotations")
no_sweep_11080585_annot <- read.delim("11080585.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
no_sweep_14738655_annot <- read.delim("14738655.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
no_sweep_8866805_annot  <- read.delim("8866805.emapper.annotations",   sep = "\t", skip = 4, header = TRUE)
sweep_1_31832435_annot  <- read.delim("31832435.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_1_32036575_annot  <- read.delim("32036575.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_1_32037335_annot  <- read.delim("32037335.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_2_45808885_annot  <- read.delim("45808885.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_2_45808905_annot  <- read.delim("45808905.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_2_45808925_annot  <- read.delim("45808925.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_3_11453035_annot  <- read.delim("11453035.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_3_8046645_annot   <- read.delim("8046645.emapper.annotations",   sep = "\t", skip = 4, header = TRUE)
sweep_3_8091485_annot   <- read.delim("8091485.emapper.annotations",   sep = "\t", skip = 4, header = TRUE)
sweep_4_20671165_annot  <- read.delim("20671165.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_4_20671355_annot  <- read.delim("20671355.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_4_20671385_annot  <- read.delim("20671385.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_5_8026935_annot   <- read.delim("8026935.emapper.annotations",   sep = "\t", skip = 4, header = TRUE)
sweep_5_8927125_annot   <- read.delim("8927125.emapper.annotations",   sep = "\t", skip = 4, header = TRUE)
sweep_5_8928825_annot   <- read.delim("8928825.emapper.annotations",   sep = "\t", skip = 4, header = TRUE)
sweep_6_26733835_annot  <- read.delim("26733835.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_6_8484365_annot   <- read.delim("8484365.emapper.annotations",   sep = "\t", skip = 4, header = TRUE)
sweep_6_8637495_annot   <- read.delim("8637495.emapper.annotations",   sep = "\t", skip = 4, header = TRUE)
sweep_7_16459765_annot  <- read.delim("16459765.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_7_33015275_annot  <- read.delim("33015275.emapper.annotations",  sep = "\t", skip = 4, header = TRUE)
sweep_7_8813075_annot   <- read.delim("8813075.emapper.annotations",   sep = "\t", skip = 4, header = TRUE)

# Remove last 3 lines
no_sweep_11080585_annot <- head(no_sweep_11080585_annot, -3)
no_sweep_14738655_annot <- head(no_sweep_14738655_annot, -3)
no_sweep_8866805_annot  <- head(no_sweep_8866805_annot,  -3)
sweep_1_31832435_annot  <- head(sweep_1_31832435_annot,  -3)
sweep_1_32036575_annot  <- head(sweep_1_32036575_annot,  -3)
sweep_1_32037335_annot  <- head(sweep_1_32037335_annot,  -3)
sweep_2_45808885_annot  <- head(sweep_2_45808885_annot,  -3)
sweep_2_45808905_annot  <- head(sweep_2_45808905_annot,  -3)
sweep_2_45808925_annot  <- head(sweep_2_45808925_annot,  -3)
sweep_3_11453035_annot  <- head(sweep_3_11453035_annot,  -3)
sweep_3_8046645_annot   <- head(sweep_3_8046645_annot,   -3)
sweep_3_8091485_annot   <- head(sweep_3_8091485_annot,   -3)
sweep_4_20671165_annot  <- head(sweep_4_20671165_annot,  -3)
sweep_4_20671355_annot  <- head(sweep_4_20671355_annot,  -3)
sweep_4_20671385_annot  <- head(sweep_4_20671385_annot,  -3)
sweep_5_8026935_annot   <- head(sweep_5_8026935_annot,   -3)
sweep_5_8927125_annot   <- head(sweep_5_8927125_annot,   -3)
sweep_5_8928825_annot   <- head(sweep_5_8928825_annot,   -3)
sweep_6_26733835_annot  <- head(sweep_6_26733835_annot,  -3)
sweep_6_8484365_annot   <- head(sweep_6_8484365_annot,   -3)
sweep_6_8637495_annot   <- head(sweep_6_8637495_annot,   -3)
sweep_7_16459765_annot  <- head(sweep_7_16459765_annot,  -3)
sweep_7_33015275_annot  <- head(sweep_7_33015275_annot,  -3)
sweep_7_8813075_annot   <- head(sweep_7_8813075_annot,   -3)

setwd("/home/senekowitsch/Thesis/Functional/03_blast/output")

# Merge
no_sweep_11080585_merged <- merge(no_sweep_11080585, no_sweep_11080585_annot, by.x = "qseqid", by.y = "X.query", all.x = TRUE)
no_sweep_14738655_merged <- merge(no_sweep_14738655, no_sweep_14738655_annot, by.x = "qseqid", by.y = "X.query", all.x = TRUE)
no_sweep_8866805_merged  <- merge(no_sweep_8866805,  no_sweep_8866805_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_1_31832435_merged  <- merge(sweep_1_31832435,  sweep_1_31832435_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_1_32036575_merged  <- merge(sweep_1_32036575,  sweep_1_32036575_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_1_32037335_merged  <- merge(sweep_1_32037335,  sweep_1_32037335_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_2_45808885_merged  <- merge(sweep_2_45808885,  sweep_2_45808885_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_2_45808905_merged  <- merge(sweep_2_45808905,  sweep_2_45808905_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_2_45808925_merged  <- merge(sweep_2_45808925,  sweep_2_45808925_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_3_11453035_merged  <- merge(sweep_3_11453035,  sweep_3_11453035_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_3_8046645_merged   <- merge(sweep_3_8046645,   sweep_3_8046645_annot,   by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_3_8091485_merged   <- merge(sweep_3_8091485,   sweep_3_8091485_annot,   by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_4_20671165_merged  <- merge(sweep_4_20671165,  sweep_4_20671165_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_4_20671355_merged  <- merge(sweep_4_20671355,  sweep_4_20671355_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_4_20671385_merged  <- merge(sweep_4_20671385,  sweep_4_20671385_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_5_8026935_merged   <- merge(sweep_5_8026935,   sweep_5_8026935_annot,   by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_5_8927125_merged   <- merge(sweep_5_8927125,   sweep_5_8927125_annot,   by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_5_8928825_merged   <- merge(sweep_5_8928825,   sweep_5_8928825_annot,   by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_6_26733835_merged  <- merge(sweep_6_26733835,  sweep_6_26733835_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_6_8484365_merged   <- merge(sweep_6_8484365,   sweep_6_8484365_annot,   by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_6_8637495_merged   <- merge(sweep_6_8637495,   sweep_6_8637495_annot,   by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_7_16459765_merged  <- merge(sweep_7_16459765,  sweep_7_16459765_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_7_33015275_merged  <- merge(sweep_7_33015275,  sweep_7_33015275_annot,  by.x = "qseqid", by.y = "X.query", all.x = TRUE)
sweep_7_8813075_merged   <- merge(sweep_7_8813075,   sweep_7_8813075_annot,   by.x = "qseqid", by.y = "X.query", all.x = TRUE)

# Split sseqid to keep only the sample ID
no_sweep_11080585_merged$sseqid <- sub("\\|.*", "", no_sweep_11080585_merged$sseqid)
no_sweep_14738655_merged$sseqid <- sub("\\|.*", "", no_sweep_14738655_merged$sseqid)
no_sweep_8866805_merged$sseqid  <- sub("\\|.*", "", no_sweep_8866805_merged$sseqid)
sweep_1_31832435_merged$sseqid  <- sub("\\|.*", "", sweep_1_31832435_merged$sseqid)
sweep_1_32036575_merged$sseqid  <- sub("\\|.*", "", sweep_1_32036575_merged$sseqid)
sweep_1_32037335_merged$sseqid  <- sub("\\|.*", "", sweep_1_32037335_merged$sseqid)
sweep_2_45808885_merged$sseqid  <- sub("\\|.*", "", sweep_2_45808885_merged$sseqid)
sweep_2_45808905_merged$sseqid  <- sub("\\|.*", "", sweep_2_45808905_merged$sseqid)
sweep_2_45808925_merged$sseqid  <- sub("\\|.*", "", sweep_2_45808925_merged$sseqid)
sweep_3_11453035_merged$sseqid  <- sub("\\|.*", "", sweep_3_11453035_merged$sseqid)
sweep_3_8046645_merged$sseqid   <- sub("\\|.*", "", sweep_3_8046645_merged$sseqid)
sweep_3_8091485_merged$sseqid   <- sub("\\|.*", "", sweep_3_8091485_merged$sseqid)
sweep_4_20671165_merged$sseqid  <- sub("\\|.*", "", sweep_4_20671165_merged$sseqid)
sweep_4_20671355_merged$sseqid  <- sub("\\|.*", "", sweep_4_20671355_merged$sseqid)
sweep_4_20671385_merged$sseqid  <- sub("\\|.*", "", sweep_4_20671385_merged$sseqid)
sweep_5_8026935_merged$sseqid   <- sub("\\|.*", "", sweep_5_8026935_merged$sseqid)
sweep_5_8927125_merged$sseqid   <- sub("\\|.*", "", sweep_5_8927125_merged$sseqid)
sweep_5_8928825_merged$sseqid   <- sub("\\|.*", "", sweep_5_8928825_merged$sseqid)
sweep_6_26733835_merged$sseqid  <- sub("\\|.*", "", sweep_6_26733835_merged$sseqid)
sweep_6_8484365_merged$sseqid   <- sub("\\|.*", "", sweep_6_8484365_merged$sseqid)
sweep_6_8637495_merged$sseqid   <- sub("\\|.*", "", sweep_6_8637495_merged$sseqid)
sweep_7_16459765_merged$sseqid  <- sub("\\|.*", "", sweep_7_16459765_merged$sseqid)
sweep_7_33015275_merged$sseqid  <- sub("\\|.*", "", sweep_7_33015275_merged$sseqid)
sweep_7_8813075_merged$sseqid   <- sub("\\|.*", "", sweep_7_8813075_merged$sseqid)

# Load sweep annotation
# load the labels file
sweep_labels <- read.delim("/home/senekowitsch/Thesis/Functional/01_prokka/genome_sweep_labels.txt", sep = "\t", header = TRUE)
sweep_labels$genome <- as.character(sweep_labels$genome)  # make sure it's character to match sseqid

no_sweep_11080585_merged <- merge(no_sweep_11080585_merged, sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
no_sweep_14738655_merged <- merge(no_sweep_14738655_merged, sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
no_sweep_8866805_merged  <- merge(no_sweep_8866805_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_1_31832435_merged  <- merge(sweep_1_31832435_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_1_32036575_merged  <- merge(sweep_1_32036575_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_1_32037335_merged  <- merge(sweep_1_32037335_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_2_45808885_merged  <- merge(sweep_2_45808885_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_2_45808905_merged  <- merge(sweep_2_45808905_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_2_45808925_merged  <- merge(sweep_2_45808925_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_3_11453035_merged  <- merge(sweep_3_11453035_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_3_8046645_merged   <- merge(sweep_3_8046645_merged,   sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_3_8091485_merged   <- merge(sweep_3_8091485_merged,   sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_4_20671165_merged  <- merge(sweep_4_20671165_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_4_20671355_merged  <- merge(sweep_4_20671355_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_4_20671385_merged  <- merge(sweep_4_20671385_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_5_8026935_merged   <- merge(sweep_5_8026935_merged,   sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_5_8927125_merged   <- merge(sweep_5_8927125_merged,   sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_5_8928825_merged   <- merge(sweep_5_8928825_merged,   sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_6_26733835_merged  <- merge(sweep_6_26733835_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_6_8484365_merged   <- merge(sweep_6_8484365_merged,   sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_6_8637495_merged   <- merge(sweep_6_8637495_merged,   sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_7_16459765_merged  <- merge(sweep_7_16459765_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_7_33015275_merged  <- merge(sweep_7_33015275_merged,  sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)
sweep_7_8813075_merged   <- merge(sweep_7_8813075_merged,   sweep_labels, by.x = "sseqid", by.y = "genome", all.x = TRUE)

# add query genome identity to each
no_sweep_11080585_merged$query_genome <- "11080585"
no_sweep_11080585_merged$query_group  <- "no_sweep"
no_sweep_14738655_merged$query_genome <- "14738655"
no_sweep_14738655_merged$query_group  <- "no_sweep"
no_sweep_8866805_merged$query_genome  <- "8866805"
no_sweep_8866805_merged$query_group   <- "no_sweep"
sweep_1_31832435_merged$query_genome  <- "31832435"
sweep_1_31832435_merged$query_group   <- "sweep_1"
sweep_1_32036575_merged$query_genome  <- "32036575"
sweep_1_32036575_merged$query_group   <- "sweep_1"
sweep_1_32037335_merged$query_genome  <- "32037335"
sweep_1_32037335_merged$query_group   <- "sweep_1"
sweep_2_45808885_merged$query_genome  <- "45808885"
sweep_2_45808885_merged$query_group   <- "sweep_2"
sweep_2_45808905_merged$query_genome  <- "45808905"
sweep_2_45808905_merged$query_group   <- "sweep_2"
sweep_2_45808925_merged$query_genome  <- "45808925"
sweep_2_45808925_merged$query_group   <- "sweep_2"
sweep_3_11453035_merged$query_genome  <- "11453035"
sweep_3_11453035_merged$query_group   <- "sweep_3"
sweep_3_8046645_merged$query_genome   <- "8046645"
sweep_3_8046645_merged$query_group    <- "sweep_3"
sweep_3_8091485_merged$query_genome   <- "8091485"
sweep_3_8091485_merged$query_group    <- "sweep_3"
sweep_4_20671165_merged$query_genome  <- "20671165"
sweep_4_20671165_merged$query_group   <- "sweep_4"
sweep_4_20671355_merged$query_genome  <- "20671355"
sweep_4_20671355_merged$query_group   <- "sweep_4"
sweep_4_20671385_merged$query_genome  <- "20671385"
sweep_4_20671385_merged$query_group   <- "sweep_4"
sweep_5_8026935_merged$query_genome   <- "8026935"
sweep_5_8026935_merged$query_group    <- "sweep_5"
sweep_5_8927125_merged$query_genome   <- "8927125"
sweep_5_8927125_merged$query_group    <- "sweep_5"
sweep_5_8928825_merged$query_genome   <- "8928825"
sweep_5_8928825_merged$query_group    <- "sweep_5"
sweep_6_26733835_merged$query_genome  <- "26733835"
sweep_6_26733835_merged$query_group   <- "sweep_6"
sweep_6_8484365_merged$query_genome   <- "8484365"
sweep_6_8484365_merged$query_group    <- "sweep_6"
sweep_6_8637495_merged$query_genome   <- "8637495"
sweep_6_8637495_merged$query_group    <- "sweep_6"
sweep_7_16459765_merged$query_genome  <- "16459765"
sweep_7_16459765_merged$query_group   <- "sweep_7"
sweep_7_33015275_merged$query_genome  <- "33015275"
sweep_7_33015275_merged$query_group   <- "sweep_7"
sweep_7_8813075_merged$query_genome   <- "8813075"
sweep_7_8813075_merged$query_group    <- "sweep_7"

colnames(no_sweep_11080585_merged)

# Slim down dataset by keeping only important columns
keep_cols <- c("sseqid", "qseqid", "pident", "length", "qlen", "slen", "qcovs", 
               "evalue.x", "bitscore", "COG_category", 
               "Preferred_name", "GOs", "KEGG_ko", "KEGG_Pathway", "CAZy",
               "PFAMs", "sweep", "query_genome", "query_group", "Description")

no_sweep_11080585_merged <- no_sweep_11080585_merged[, keep_cols]
no_sweep_14738655_merged <- no_sweep_14738655_merged[, keep_cols]
no_sweep_8866805_merged  <- no_sweep_8866805_merged[,  keep_cols]
sweep_1_31832435_merged  <- sweep_1_31832435_merged[,  keep_cols]
sweep_1_32036575_merged  <- sweep_1_32036575_merged[,  keep_cols]
sweep_1_32037335_merged  <- sweep_1_32037335_merged[,  keep_cols]
sweep_2_45808885_merged  <- sweep_2_45808885_merged[,  keep_cols]
sweep_2_45808905_merged  <- sweep_2_45808905_merged[,  keep_cols]
sweep_2_45808925_merged  <- sweep_2_45808925_merged[,  keep_cols]
sweep_3_11453035_merged  <- sweep_3_11453035_merged[,  keep_cols]
sweep_3_8046645_merged   <- sweep_3_8046645_merged[,   keep_cols]
sweep_3_8091485_merged   <- sweep_3_8091485_merged[,   keep_cols]
sweep_4_20671165_merged  <- sweep_4_20671165_merged[,  keep_cols]
sweep_4_20671355_merged  <- sweep_4_20671355_merged[,  keep_cols]
sweep_4_20671385_merged  <- sweep_4_20671385_merged[,  keep_cols]
sweep_5_8026935_merged   <- sweep_5_8026935_merged[,   keep_cols]
sweep_5_8927125_merged   <- sweep_5_8927125_merged[,   keep_cols]
sweep_5_8928825_merged   <- sweep_5_8928825_merged[,   keep_cols]
sweep_6_26733835_merged  <- sweep_6_26733835_merged[,  keep_cols]
sweep_6_8484365_merged   <- sweep_6_8484365_merged[,   keep_cols]
sweep_6_8637495_merged   <- sweep_6_8637495_merged[,   keep_cols]
sweep_7_16459765_merged  <- sweep_7_16459765_merged[,  keep_cols]
sweep_7_33015275_merged  <- sweep_7_33015275_merged[,  keep_cols]
sweep_7_8813075_merged   <- sweep_7_8813075_merged[,   keep_cols]

# Check if the number of rows after merge is greater than the BLAST hits
nrow(sweep_1_32036575_merged) > nrow(sweep_1_32036575)

# Combine to one df
all_data <- rbind(
  no_sweep_11080585_merged, no_sweep_14738655_merged, no_sweep_8866805_merged,
  sweep_1_31832435_merged,  sweep_1_32036575_merged,  sweep_1_32037335_merged,
  sweep_2_45808885_merged,  sweep_2_45808905_merged,  sweep_2_45808925_merged,
  sweep_3_11453035_merged,  sweep_3_8046645_merged,   sweep_3_8091485_merged,
  sweep_4_20671165_merged,  sweep_4_20671355_merged,  sweep_4_20671385_merged,
  sweep_5_8026935_merged,   sweep_5_8927125_merged,   sweep_5_8928825_merged,
  sweep_6_26733835_merged,  sweep_6_8484365_merged,   sweep_6_8637495_merged,
  sweep_7_16459765_merged,  sweep_7_33015275_merged,  sweep_7_8813075_merged
)

# rename sweep to target_sweep"
all_data <- rename(all_data, target_sweep = sweep)

head(all_data)


# -------------------------------------------------
# Column names
# -------------------------------------------------
# There is a separate file for each blast that was done! We have to keep this in mind somehow when we do the analysis
# qseqid: Id of the gene that was used to blast against the db. This can be found in the annotation
# sseqid: Id of the genome that the query was found in
# target_sweep: which sweep was the genome of sseqid in
# query_genome: which genome the query gene came from
# query_group: which sweep group the query genome belongs to (sweep_1, sweep_2, ..., no_sweep)

head(no_sweep_11080585_merged[, c("sseqid", "target_sweep")])
table(no_sweep_11080585_merged$target_sweep, useNA = "always")  # how many hits per sweep group, and NAs (sseqids not in the label file)
table(no_sweep_11080585_merged$sseqid, useNA = "always") # how many hits per sseqid
table(no_sweep_11080585_merged$qseqid, useNA = "always") # how many hits per qseqid
table(no_sweep_11080585_merged$query_genome, useNA = "always")
table(no_sweep_11080585_merged$query_group, useNA = "always")

# -------------------------------------------------
# Basic checks of the combined data frame
# -------------------------------------------------
dim(all_data)          # how many rows and columns
str(all_data)          # column types
summary(all_data)      # basic stats for each column

# check for NAs - important to know which genes lack annotations
colSums(is.na(all_data))

# how many unique query genomes
n_distinct(all_data$query_genome)  # should be 24

# how many unique target genomes
n_distinct(all_data$sseqid)        # should be ~472

# hit distribution across sweep groups
table(all_data$target_sweep)
table(all_data$query_group)

# -------------------------------------------------
# Safe data frame for easyer loading later
# -------------------------------------------------
fwrite(all_data, "all_data_2.tsv.gz", sep = "\t")

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

# Drop columns that are not needed for the analysis to save memory and simplify downstream work
keep_cols <- c("sseqid", "qseqid", "pident",
               "COG_category", "Preferred_name", "GOs", 
               "KEGG_ko", "KEGG_Pathway", "CAZy", "PFAMs", 
               "target_sweep", "query_genome", "query_group",
               "Description")

all_data_filtered <- all_data_filtered[, keep_cols]

# -------------------------------------------------
# Safe data frame for easyer loading later
# -------------------------------------------------
fwrite(all_data_filtered, "all_data_filtered_2.tsv.gz", sep = "\t")
object.size(all_data_filtered) |> format(units = "GB")
