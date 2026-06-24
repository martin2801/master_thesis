#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# =============================================================================
# recombination_pipeline.sh
# =============================================================================
#
# DESCRIPTION:
#   Recombination detection and characterisation pipeline for a collection of
#   Salmonella Infantis genome assemblies.  The script benchmarks three
#   independent approaches — parsnp/ClonalFrameML, parsnp/Gubbins, and
#   snippy/Gubbins+ClonalFrameML — and quantifies their agreement via
#   BEDTools interval arithmetic.
#
# WORKFLOW OVERVIEW:
#   1.  Align genomes with parsnp (reference-free, random reference)
#   2.  Run ClonalFrameML on the parsnp alignment (with & without reference)
#   3.  Run Gubbins on the parsnp alignment (without reference only)
#   4.  Repeat steps 1–3 using a fixed, named reference genome
#   5.  Align genomes against a named reference with snippy
#   6.  Run Gubbins and ClonalFrameML on the snippy core alignment
#   7.  Compare Gubbins vs ClonalFrameML predictions (snippy branch) using
#       BEDTools: compute overlap, unique regions, and Jaccard index
#
# INPUT:
#   - 117 Salmonella Infantis genome assemblies in FASTA (.fna) format
#   - One reference genome: GCF_000506925.1 (S. Infantis SI119944)
#
# OUTPUT (under $output_base):
#   parsnp/                   parsnp alignment (random ref), tree, XMFA, FASTA
#   parsnp_with_ref/          same, but with the named reference forced in
#   clonalframeML/            ClonalFrameML results (parsnp, no named ref)
#   clonalframeML_with_ref/   ClonalFrameML results (parsnp, named ref)
#   clonalframeML_snippy/     ClonalFrameML results (snippy core alignment)
#   gubbins/                  Gubbins results (parsnp alignment, no named ref)
#   gubbins_with_ref/         Gubbins results (parsnp alignment, named ref)
#   gubbins_snippy/           Gubbins results (snippy core alignment)
#   snippy/                   per-sample snippy outputs + core alignment
#   *.bed                     BEDTools intermediate and comparison files
#   overlapping_regions.txt   Regions predicted by BOTH tools
#   gubbins_only.bed          Regions unique to Gubbins
#   cfml_only.bed             Regions unique to ClonalFrameML
#
# DEPENDENCIES (each activated via conda):
#   parsnp, harvesttools  -> conda env: parsnp
#   ClonalFrameML         -> conda env: clonalframeML
#   Gubbins               -> conda env: gubbins
#   snippy                -> conda env: snippy_env
#   IQ-TREE               -> conda env: phylo_pipeline
#   BioPython             -> conda env: fastANI  (used for tree pruning)
#   BEDTools              -> conda env: bedtools
#   Python/matplotlib     -> conda env: ani_heatmap
#
# AUTHOR:    [Martin Senekowitsch]
# DATE:      [16.4.2025]
# VERSION:   1.0
# =============================================================================

# =============================================================================
# CONFIGURATION
# Variables are defined once here so the script is easy to adapt to new
# datasets without touching the tool-invocation sections.
# =============================================================================

# Variables
# Root output directory — all tool-specific sub-directories are created here
base_dir='/home/senekowitsch/Thesis/Sweeps/01_recombination_checks'
output_base="$base_dir/output"

# Tool-specific output sub-directories (created with mkdir -p before use)
output_parsnp="${output_base}/parsnp"
output_snippy="${output_base}/snippy"
output_clonalframeML="${output_base}/clonalframeML"
output_gubbins="${output_base}/gubbins"

# Directory holding the .fna genome assemblies used as the query set
raw_data='/home/senekowitsch/Thesis/PopCoGenomeS/00_data_alt' # 117 genome assemblies in .fna format

# Named reference genome (GCF_000506925.1, S. Infantis SI119944).
# Used explicitly in the "with_ref" parsnp run and as the snippy reference. 
reference_genome='/home/senekowitsch/Thesis/QC/01_ANI/Genomes/Reference/infantis/GCF_000506925.1/ncbi_dataset/data/GCF_000506925.1/GCF_000506925.1_SI119944_genomic.fna'

# Number of CPU threads passed to all multi-threaded tools
threads=30

# ==========================================================================================================================================================
# Module 1: PARSNP + ClonalFrameML + Gubbins (random reference)
# ==========================================================================================================================================================


# =============================================================================
# SECTION 1: PARSNP — WHOLE-GENOME ALIGNMENT (RANDOM REFERENCE)
# =============================================================================
# parsnp is run with '-r !' which tells it to pick a reference genome at random
# from the input set.  The resulting alignment, tree, and VCF are used as
# input for both ClonalFrameML and Gubbins in the sections that follow.
#
# Key flags:
#   -r !              randomly select reference from the input directory
#   -d                directory of input FASTA assemblies
#   -o                output directory
#   -p                number of threads
#   -e                enable extra sensitivity (MUMi-based filtering relaxed)
#   --extend-ani-cutoff 0.80
#                     include genomes with at least 80% ANI to the reference
#   --vcf             output SNP calls in VCF format
#   --fo              force overwrite of existing output
#   -c                use all sequences (don't stop at MUMi threshold)
# =============================================================================
conda activate parsnp
mkdir -p "$output_parsnp"

# Run parsnp
parsnp -r ! \
       -d "$raw_data" \
       -o "$output_parsnp" \
       -p "$threads" \
       -e \
       --extend-ani-cutoff 0.80 \
       --vcf \
       --fo \
       -c

# =============================================================================
# SECTION 2: CLONALFRAMEML — PARSNP ALIGNMENT, RANDOM REFERENCE (FULL SET)
# =============================================================================
# harvesttools converts the parsnp XMFA alignment into a multi-FASTA file
# that ClonalFrameML can consume.
#
# ClonalFrameML is first run on the full alignment (which includes the parsnp-
# chosen reference sequence) to serve as a baseline before the reference is
# removed in Section 3.
# =============================================================================

mkdir -p "$output_clonalframeML"

# Run harvesttools to convert parsnp output to clonalframeML input (fasta format)
harvesttools -x "$output_parsnp"/parsnp.xmfa -M "$output_parsnp"/parsnp.fasta

# Run clonalframeML with the reference included
conda activate clonalframeML
ClonalFrameML "$output_parsnp/parsnp.tree" \
               "$output_parsnp/parsnp.fasta"\
               "$output_clonalframeML/infantis_cfml"


# =============================================================================
# SECTION 3: REMOVE REFERENCE SEQUENCE & PRUNE TREE
# =============================================================================
# parsnp appends ".ref" to the header of whichever sequence it chose as the
# reference.  Because the reference is typically over-represented in the
# alignment (it defines the coordinate space), many pipelines exclude it from
# downstream analyses to avoid biasing recombination inference.
#
# Step A: identify the reference sequence name in the FASTA header.
# Step B: write a new FASTA without the reference sequence (but rename the
#         ".ref" header variant so sample names are consistent).
# Step C: prune the same tip from the parsnp NJ tree using BioPython so that
#         the tree and alignment stay in sync.
# =============================================================================

# Step A — Identify the parsnp-chosen reference sequence name.
# The reference header looks like: >19357715.fna.ref
# REF_FULL_NAME captures the full name including ".ref" suffix.
# REF_CLEAN_NAME strips the ".ref" suffix (used for renaming inside the FASTA).
REF_FULL_NAME=$(grep ">" "$output_parsnp/parsnp.fasta" | grep "\.ref" | head -n 1 | sed 's/>//')
REF_CLEAN_NAME=$(echo "$REF_FULL_NAME" | sed 's/\.ref//')

echo "Reference found: $REF_FULL_NAME"

# Step B — Remove the ".ref" duplicate entry; rename ".ref" header to clean name.
# The awk script works in record-separator mode (RS=">") so each FASTA entry
# is one record.  Records whose first field matches the clean name (the
# ordinary copy already present) are skipped; the ".ref" copy is kept but
# renamed; all other records are passed through unchanged.
awk -v ref_full="$REF_FULL_NAME" -v ref_clean="$REF_CLEAN_NAME" '
BEGIN {RS=">"; ORS=""} 
$1 == ref_clean {next} 
$1 == ref_full {sub(ref_full, ref_clean); print ">"$0; next} 
NF {print ">"$0}' "$output_parsnp/parsnp.fasta" > "$output_parsnp/alignment_unique.fasta"


# Step C — Prune the reference tip from the parsnp NJ tree.
# The resulting tree (parsnp_clean.tree) has the same leaf set as
# alignment_unique.fasta and is used by ClonalFrameML and Gubbins below.
conda activate fastANI
python3 << EOF
from Bio import Phylo
tree = Phylo.read("$output_parsnp/parsnp.tree", "newick")
tree.prune("$REF_FULL_NAME")
Phylo.write(tree, "$output_parsnp/parsnp_clean.tree", "newick")
print("Done - tips remaining:", tree.count_terminals())
EOF


# =============================================================================
# SECTION 4: CLONALFRAMEML — PARSNP ALIGNMENT, REFERENCE EXCLUDED
# =============================================================================
# Re-run ClonalFrameML on the reference-free alignment and pruned tree.
# This is the primary ClonalFrameML result for the parsnp-based branch.
# =============================================================================

conda activate clonalframeML
ClonalFrameML "$output_parsnp/parsnp_clean.tree" \
              "$output_parsnp/alignment_unique.fasta"\
              "$output_clonalframeML/infantis_cfml_unique"


# =============================================================================
# SECTION 5: GUBBINS — PARSNP ALIGNMENT, REFERENCE EXCLUDED
# =============================================================================
# Gubbins identifies recombinant regions by iteratively building a phylogeny
# and detecting spatial clustering of SNPs that is inconsistent with clonal
# descent.  It takes a whole-genome alignment as input (no separate tree
# required — it generates and refines the tree internally).
#
# Key flags:
#   --prefix     output file prefix
#   --threads    number of CPU threads
#   --min-snps 3 minimum SNPs required to call a recombinant block
# =============================================================================

mkdir -p "$output_gubbins"
cd "$output_gubbins"
conda activate gubbins
run_gubbins.py \
  --prefix infantis_gubbins \
  --threads "$threads" \
  --min-snps 3 \
  "$output_parsnp/alignment_unique.fasta"


# ==========================================================================================================================================================
# Module 2: PARSNP + ClonalFrameML + Gubbins (external named reference)
# ==========================================================================================================================================================

# =============================================================================
# SECTION 6: PARSNP — WHOLE-GENOME ALIGNMENT (NAMED REFERENCE)
# =============================================================================
# Repeat the parsnp alignment using the biologically meaningful reference
# genome (GCF_000506925.1) instead of a randomly chosen one.  This allows
# us to assess whether the choice of reference affects the downstream
# recombination predictions.
#
# The "_with_ref" suffix is appended to all output directories in this section
# to distinguish them from the random-reference results above.
# =============================================================================

conda activate parsnp
mkdir -p "${output_parsnp}_with_ref"

# Run parsnp with the reference included
parsnp -r "$reference_genome" \
       -d "$raw_data" \
       -o "${output_parsnp}_with_ref" \
       -p "$threads" \
       -e \
       --extend-ani-cutoff 0.80 \
       --vcf \
       --fo \
       -c

# Run harvesttools to convert parsnp output to clonalframeML input (fasta format)
harvesttools -x "${output_parsnp}_with_ref/parsnp.xmfa" -M "${output_parsnp}_with_ref/parsnp.fasta"

# Run clonalframeML with the reference included
mkdir -p "${output_clonalframeML}_with_ref"
conda activate clonalframeML
ClonalFrameML "${output_parsnp}_with_ref/parsnp.tree" \
               "${output_parsnp}_with_ref/parsnp.fasta"\
               "${output_clonalframeML}_with_ref/infantis_cfml_with_ref"

# Run Gubbins on the alignment with the reference included
mkdir -p "${output_gubbins}_with_ref"
cd "${output_gubbins}_with_ref"
conda activate gubbins
run_gubbins.py \
  --prefix infantis_gubbins_with_ref \
  --threads "$threads" \
  --min-snps 3 \
  "${output_parsnp}_with_ref/parsnp.fasta"


# ==========================================================================================================================================================
# Module 3: Snippy + ClonalFrameML + Gubbins (external named reference)
# ==========================================================================================================================================================

# =============================================================================
# SECTION 7: SNIPPY — READ-FREE VARIANT CALLING USING ASSEMBLIES
# =============================================================================
# snippy is used here to call SNPs from assembled contigs against the named
# reference genome.  Although snippy was designed for short-read data, the
# --ctgs flag accepts contigs (assemblies), making it applicable here.
#
# A separate snippy output directory is created for each of the 117 genomes.
# snippy-core is then used to merge all per-sample VCFs into a single core
# SNP alignment (core.full.aln) that is suitable for whole-genome recombination
# detection.
#
# The core alignment provides an alternative starting point to the parsnp-
# based alignments above, allowing us to assess the impact of the alignment
# strategy on recombination inference.
# =============================================================================

conda activate snippy_env
mkdir -p "$output_snippy"

# --- Per-sample snippy runs ---
# Each genome is processed independently against the named reference.
# Output is written to a sample-specific subdirectory (out_<sample_name>).
for genome in "$raw_data"/*.fna; do
    # Extract the filename without the path and extension for the folder name
    sample_name=$(basename "$genome" .fna)
    
    echo "Processing $sample_name..."

    # Key flags:
    #   --cpus   number of CPU threads
    #   --outdir per-sample output directory
    #   --ref    reference genome (FASTA or GenBank)
    #   --ctgs   input assembled contigs (instead of reads)
    snippy --cpus "$threads" \
           --outdir "$output_snippy/out_$sample_name" \
           --ref "$reference_genome" \
           --ctgs "$genome"
done

# --- Core alignment ---
# snippy-core merges all per-sample VCFs and produces a multi-sample alignment.
# core.full.aln includes ALL reference positions (not just polymorphic sites),
# which is required by Gubbins (it needs invariant sites to model substitution
# rates correctly).
snippy-core --prefix core --ref "$reference_genome" $output_snippy/*
mv core.* "$output_snippy/"

# Remove the snippy folders
rm -rf "$output_snippy"/out_*

# Quick sanity-check stats
echo "Core SNP alignment done"
echo "Number of SNPs: $(grep -v '#' $output_snippy/core.vcf | wc -l)"
echo "Alignment length: $(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{print length(seq)}' $output_snippy/core.full.aln | sort -nu)"


# =============================================================================
# SECTION 8: GUBBINS — SNIPPY CORE ALIGNMENT
# =============================================================================
# Gubbins is run on the snippy core alignment.  The full alignment
# (core.full.aln) is preferred over core.aln because Gubbins requires sites
# that are invariant across samples to correctly model the substitution process.
# =============================================================================

mkdir -p "${output_gubbins}_snippy"
cd "${output_gubbins}_snippy"
conda activate gubbins
run_gubbins.py \
  --prefix infantis_gubbins_snippy \
  --threads $threads \
  --min-snps 3 \
  "$output_snippy/core.full.aln"

echo "Gubbins complete"
echo "Recombinant regions found: $(grep -v '^#' ${output_gubbins}_snippy/infantis_gubbins_snippy.recombination_predictions.gff | wc -l)"


# =============================================================================
# SECTION 9: IQTREE + CLONALFRAMEML — SNIPPY CORE ALIGNMENT
# =============================================================================
# ClonalFrameML requires a pre-computed input tree.  Unlike Gubbins, it does
# not build its own tree.  IQ-TREE is used here to infer a maximum-likelihood
# tree from the snippy core alignment under the GTR+Gamma model, which is
# then passed to ClonalFrameML.
#
# IQ-TREE flags:
#   -s     input alignment
#   -m     substitution model (GTR+G = General Time Reversible + Gamma rates)
#   -T     threads
#   -B 1000 UFBoot ultrafast bootstrap replicates
#   --redo  overwrite any existing checkpoint files
#
# ClonalFrameML models recombination as a Poisson process on the phylogeny,
# estimating: import length (δ), import rate (ρ/θ), and divergence of
# imported DNA (1/ν) per branch.
# =============================================================================

conda activate phylo_pipeline

# create a tree from the snippy core alignment using iqtree
mkdir -p "${output_clonalframeML}_snippy/iqtree"
cd "${output_clonalframeML}_snippy/iqtree"
iqtree \
    -s "$output_snippy/core.full.aln" \
    --prefix "${output_clonalframeML}_snippy/iqtree/raw_tree" \
    -m GTR+G \
    -T "$threads" \
    -B 1000 \
    --redo

# Run ClonalFrameML using the IQ-TREE ML tree
conda activate clonalframeML
ClonalFrameML "${output_clonalframeML}_snippy/iqtree/raw_tree.treefile" \
               "${output_snippy}/core.full.aln"\
               "${output_clonalframeML}_snippy/infantis_cfml_snippy"

echo "------------------------------------------------------------"
echo "ClonalFrameML complete"
echo "Recombinant regions found: $(wc -l < ${output_clonalframeML}_snippy/infantis_cfml_snippy.importation_status.txt) (including header)"
echo "------------------------------------------------------------"


# =============================================================================
# SECTION 10: BEDTOOLS — COMPARE GUBBINS VS CLONALFRAMEML (SNIPPY BRANCH)
# =============================================================================
# Both tools output genomic intervals predicted to be recombinant.  BEDTools
# is used to convert and compare these interval sets.
#
# Gubbins output:  GFF3 file (<prefix>.recombination_predictions.gff)
#   Columns used:  seq_name, start, end, attributes
#
# ClonalFrameML output:  TSV file (<prefix>.importation_status.txt)
#   Columns used:  Node, Beg, End  (1-based coordinates)
#
# Both sets are converted to 3-column BED format (chrom, start, end) and
# harmonised to the same chromosome label ("core") before comparison.
# =============================================================================

echo "Comparing Gubbins and ClonalFrameML outputs for the snippy core alignment"

cd "${output_base}"

# --- Convert Gubbins GFF to BED ---
# GFF fields: seqname(1), source(2), feature(3), start(4), end(5), ..., attrs(9)
# BED output: seqname, start, end, attributes  (sorted by start position)
grep -v "^#" ${output_gubbins}_snippy/infantis_gubbins_snippy.recombination_predictions.gff \
  | awk '{print $1"\t"$4"\t"$5"\t"$9}' \
  | sort -k2,2n > gubbins_regions.bed

# --- Convert ClonalFrameML importation_status to BED ---
# Header line is skipped (tail -n +2).
# A dummy chromosome name "core" is prepended to match the Gubbins output.
tail -n +2 ${output_clonalframeML}_snippy/infantis_cfml_snippy.importation_status.txt \
  | awk '{print "core\t"$2"\t"$3"\t"$1}' \
  | sort -k2,2n > cfml_regions.bed

# --- Diagnostic checks ---
# Verify that chromosome labels and region counts look sensible before
# running the actual set operations below.
cut -f1 gubbins_regions.bed | sort -u   # Should show the Gubbins seq name
cut -f1 cfml_regions.bed | sort -u      # Should show "core" as the CFML seq name
wc -l gubbins_regions.bed               # Total gubbins regions
wc -l cfml_regions.bed                  # Total CFML regions
cat gubbins_regions.bed | head -5       # Preview first 5 Gubbins entries
cat cfml_regions.bed | head -5          # Preview first 5 CFML entries

# --- Normalise chromosome names ---
# Gubbins uses the sequence name from the alignment (e.g., "SEQUENCE" or the
# sample name); ClonalFrameML already says "core".  Rename the Gubbins
# chromosome to "core" so bedtools can perform set operations across both files.
awk '{$1="core"; print}' OFS="\t" gubbins_regions.bed > gubbins_fixed.bed

# CFML is already labelled "core" — just copy to keep naming consistent
cp cfml_regions.bed cfml_fixed.bed

# --- BEDTools set operations ---
conda activate bedtools
# Regions predicted by BOTH tools (intersection with overlap size appended)
bedtools intersect -a gubbins_fixed.bed -b cfml_fixed.bed -wo > overlapping_regions.txt
# Regions unique to Gubbins (not overlapping any CFML region)
bedtools subtract -a gubbins_fixed.bed -b cfml_fixed.bed > gubbins_only.bed
# Regions unique to ClonalFrameML (not overlapping any Gubbins region)
bedtools subtract -a cfml_fixed.bed -b gubbins_fixed.bed > cfml_only.bed

# --- Quick summary ---
echo "Overlapping regions: $(wc -l < overlapping_regions.txt)"
echo "Gubbins only: $(wc -l < gubbins_only.bed)"
echo "CFML only: $(wc -l < cfml_only.bed)"

# Size distribution of Gubbins-unique fragments (useful to judge whether they
# are likely artefacts or genuine biology)
awk '{print $3-$2}' gubbins_only.bed | sort -n | uniq -c | head -20

# Total base-pair coverage for Gubbins-only vs all Gubbins predictions
awk '{sum+=$3-$2} END{print "Gubbins-only bp:", sum}' gubbins_only.bed
awk '{sum+=$3-$2} END{print "Total Gubbins bp:", sum}' gubbins_fixed.bed

# Check whether Gubbins-unique fragments cluster at particular genomic loci
# (e.g. mobile elements, rRNA operons)
awk '{print $1"\t"$2"\t"$3}' gubbins_only.bed | sort -k2,2n | uniq -c | sort -rn | head -20

# Spot-check two genomic windows: do CFML predictions cover these coordinates
# at all, even if the boundaries differ from Gubbins?
echo -e "core\t1340000\t1400000" | bedtools intersect -a cfml_fixed.bed -b - 
echo -e "core\t4750000\t4800000" | bedtools intersect -a cfml_fixed.bed -b -


# =============================================================================
# SECTION 11: PERFORMANCE COMPARISON SUMMARY
# =============================================================================
# Compute a quantitative summary of the agreement between Gubbins and
# ClonalFrameML for the snippy-based branch.
#
# Metrics reported:
#   - Number of predicted recombinant regions per tool
#   - Total base pairs (bp) in predicted regions (after merging within each tool)
#   - Shared bp (intersection), tool-specific bp, and union bp
#   - Jaccard index (bp-based): J = |intersection| / |union|
#     (0 = no overlap; 1 = identical predictions)
#   - Percentage of each tool's predicted bp that overlaps the other tool
# =============================================================================

echo ""
echo "============================================================"
echo "  GUBBINS vs CLONALFRAMEML — PERFORMANCE COMPARISON SUMMARY"
echo "============================================================"
 
# --- Region counts ---
gubbins_n=$(wc -l < gubbins_fixed.bed)
cfml_n=$(wc -l < cfml_fixed.bed)
echo ""
echo "Number of predicted recombinant regions:"
echo "  Gubbins:      $gubbins_n"
echo "  ClonalFrameML: $cfml_n"
 
# --- Total bp per tool ---
# Merge overlapping regions within each tool first to avoid double-counting
gubbins_bp=$(bedtools merge -i gubbins_fixed.bed \
  | awk '{sum += $3 - $2} END {printf "%d", sum+0}')
cfml_bp=$(bedtools merge -i cfml_fixed.bed \
  | awk '{sum += $3 - $2} END {printf "%d", sum+0}')
echo ""
echo "Total size of predicted recombinant regions (bp, after merging overlaps within each tool):"
echo "  Gubbins:       $gubbins_bp bp"
echo "  ClonalFrameML: $cfml_bp bp"
 
# --- Overlap ---
# Pre-sort merged BED files to disk (avoids stdin sorting errors in bedtools)
bedtools merge -i gubbins_fixed.bed | sort -k1,1 -k2,2n > gubbins_merged.bed
bedtools merge -i cfml_fixed.bed    | sort -k1,1 -k2,2n > cfml_merged.bed
 
# Intersection: bp covered by BOTH tools
overlap_bp=$(bedtools intersect -a gubbins_merged.bed -b cfml_merged.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  | awk '{sum += $3 - $2} END {printf "%d", sum+0}')
 
# Union: bp covered by EITHER tool (cat + sort + merge collapses all)
union_bp=$(cat gubbins_fixed.bed cfml_fixed.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  | awk '{sum += $3 - $2} END {printf "%d", sum+0}')
 
# bp unique to each tool = their total minus the shared intersection
gubbins_only_bp=$(( gubbins_bp - overlap_bp ))
cfml_only_bp=$(( cfml_bp - overlap_bp ))
 
# Jaccard index: intersection / union (should be 0–1)
jaccard=$(awk -v o="$overlap_bp" -v u="$union_bp" \
  'BEGIN {if (u > 0) printf "%.4f", o/u; else print "N/A"}')
 
echo ""
echo "Overlap between tools:"
echo "  Shared bp (intersection):  $overlap_bp bp"
echo "  Gubbins-only bp:           $gubbins_only_bp bp"
echo "  ClonalFrameML-only bp:     $cfml_only_bp bp"
echo "  Union bp:                  $union_bp bp"
echo "  Jaccard index (bp-based):  $jaccard  (0 = no overlap, 1 = identical)"
 
# --- Fraction of each tool's regions that are shared ---
gubbins_pct=$(awk -v o="$overlap_bp" -v t="$gubbins_bp" \
  'BEGIN {if (t > 0) printf "%.1f", 100*o/t; else print "N/A"}')
cfml_pct=$(awk -v o="$overlap_bp" -v t="$cfml_bp" \
  'BEGIN {if (t > 0) printf "%.1f", 100*o/t; else print "N/A"}')
echo ""
echo "Fraction of each tool's bp that overlap with the other:"
echo "  Gubbins bp covered by CFML: ${gubbins_pct}%"
echo "  CFML bp covered by Gubbins: ${cfml_pct}%"
 
echo ""
echo "============================================================"
echo ""


# =============================================================================
# SECTION 12: VISUALISATION
# =============================================================================
# recombination_plot.py produces summary figures from the BEDTools output and
# tool-specific result files.  The script is expected to reside in the working
# directory.  Output format (PNG/SVG/PDF) is defined within the Python script.
# =============================================================================

conda activate ani_heatmap
python3 ../recombination_plot.py \
    -g gubbins_fixed.bed \
    -f cfml_fixed.bed \
    -o recombination_comparison.png \
    -l 4900000

conda activate ani_heatmap
python3 ../recombination_plot_pres.py \
    -g gubbins_fixed.bed \
    -f cfml_fixed.bed \
    -o recombination_comparison_pres.png \
    -l 5000000