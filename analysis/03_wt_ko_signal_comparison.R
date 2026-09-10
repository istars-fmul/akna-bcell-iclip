#!/usr/bin/env Rscript

# Comparison of WT and Akna KO crosslink signal at reproducible binding sites
# -------------------------------------------------------------------------
# Quantify the same WT-defined, annotated intervals in all five libraries and
# append per-group means and WT/KO ratios to the annotation workbook.
#
# Inputs from analysis 02: annotated CSV (1-based inclusive coordinates) and
# its matching workbook. BigWigs must contain strand-specific signal normalized
# per library, with explicit WT1-3 and KO1-2 paths in the order below.
# Output: a workbook retaining the annotated sites and adding KO1, KO2, WT_mean,
# KO_mean, fold_change and pass_threshold. Both passing and failing sites remain
# in the workbook; analysis 04 selects the WT-enriched genes from this flag.
#
# Usage:
# Rscript analysis/03_wt_ko_signal_comparison.R ANNOTATED.csv INPUT.xlsx \
#   OUTPUT.xlsx WT1.plus.bw WT1.minus.bw WT2.plus.bw WT2.minus.bw \
#   WT3.plus.bw WT3.minus.bw KO1.plus.bw KO1.minus.bw KO2.plus.bw KO2.minus.bw

## Load libraries
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(openxlsx)
})

## Load data: explicit sample paths prevent filesystem-order-dependent labels
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13) {
  stop(
    "Expected 13 positional arguments. See the usage comment at the top of ",
    "analysis/03_wt_ko_signal_comparison.R."
  )
}

annotated_csv <- args[1]
input_workbook <- args[2]
output_workbook <- args[3]
track_paths <- list(
  WT1 = c(plus = args[4], minus = args[5]),
  WT2 = c(plus = args[6], minus = args[7]),
  WT3 = c(plus = args[8], minus = args[9]),
  KO1 = c(plus = args[10], minus = args[11]),
  KO2 = c(plus = args[12], minus = args[13])
)

# Run from the repository root; AKNA_UTILS_R can override the shared helper path.
utils_r <- Sys.getenv("AKNA_UTILS_R", "workflow/utils/utils.r")
source(utils_r)

# Import each strand as run-length encoded coverage for interval-wise summation.
tracks <- lapply(track_paths, function(paths) {
  list(
    plus = import(paths[["plus"]], "BigWig", as = "Rle"),
    minus = import(paths[["minus"]], "BigWig", as = "Rle")
  )
})

# Load reproducible sites with gene and region annotations from the CSV.
# Its coordinates are already 1-based; unlike BED import, no start shift is needed.
bs_df <- read.csv(annotated_csv, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("seqnames", "start", "end", "strand", "gene_name")
missing_columns <- setdiff(required_columns, names(bs_df))
if (length(missing_columns) > 0) {
  stop("Annotated CSV is missing: ", paste(missing_columns, collapse = ", "))
}

bs <- makeGRangesFromDataFrame(
  bs_df,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

## Sum signal over each binding site using only its matching strand
# The sum covers the whole 9-nt interval, not just its maximum-signal position.
bs.p <- bs[strand(bs) == "+"]
bs.m <- bs[strand(bs) == "-"]
for (sample in names(tracks)) {
  mcols(bs.p)[[sample]] <- sum_track_over_ranges(tracks[[sample]]$plus, bs.p)
  mcols(bs.m)[[sample]] <- sum_track_over_ranges(tracks[[sample]]$minus, bs.m)
}
# Recombine the strand-specific sites before calculating group-level summaries.
bs_combined <- c(bs.p, bs.m)

## Calculate WT/KO signal ratios
# Average site sums separately over three WT and two KO libraries. Add 0.01
# to both means to define ratios at zero signal; use the strict >1.5 threshold.
# This is a descriptive enrichment filter, not a differential-binding test.
pseudocount <- 0.01
signal_data <- as.data.frame(bs_combined) %>%
  rowwise() %>%
  mutate(
    WT_mean = mean(c(WT1, WT2, WT3), na.rm = TRUE),
    KO_mean = mean(c(KO1, KO2), na.rm = TRUE),
    fold_change = (WT_mean + pseudocount) / (KO_mean + pseudocount),
    pass_threshold = fold_change > 1.5
  ) %>%
  ungroup()

# Select the new signal columns to append to the existing annotation workbook.
new_columns <- signal_data %>%
  select(
    seqnames, start, end, strand, KO1, KO2,
    WT_mean, KO_mean, fold_change, pass_threshold
  )

## Update and export the workbook
# Read all input sheets and join by chromosome, coordinates and strand.
# Analysis 02 supplies the primary all_reproducible_binding_sites sheet.
sheet_names <- getSheetNames(input_workbook)
sheet_data <- setNames(
  lapply(sheet_names, function(sheet) read.xlsx(input_workbook, sheet = sheet)),
  sheet_names
)
updated_data <- lapply(sheet_data, function(sheet) {
  left_join(sheet, new_columns, by = c("seqnames", "start", "end", "strand"))
})
if (!"all_reproducible_binding_sites" %in% names(updated_data)) {
  names(updated_data)[1] <- "all_reproducible_binding_sites"
}

# Save a separate signal-comparison workbook, preserving the input annotations.
dir.create(dirname(output_workbook), recursive = TRUE, showWarnings = FALSE)
write.xlsx(updated_data, file = output_workbook, overwrite = TRUE)

message("Strict WT/KO ratio >1.5 sites: ", sum(signal_data$pass_threshold))
message(
  "Genes represented by strict sites: ",
  length(unique(signal_data$gene_name[signal_data$pass_threshold]))
)
