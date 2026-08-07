#!/usr/bin/env Rscript

# Usage:
# Rscript analysis/03_wt_ko_signal_comparison.R ANNOTATED.csv INPUT.xlsx \
#   OUTPUT.xlsx WT1.plus.bw WT1.minus.bw WT2.plus.bw WT2.minus.bw \
#   WT3.plus.bw WT3.minus.bw KO1.plus.bw KO1.minus.bw KO2.plus.bw KO2.minus.bw

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(openxlsx)
})

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

utils_r <- Sys.getenv("AKNA_UTILS_R", "workflow/utils/utils.r")
source(utils_r)

tracks <- lapply(track_paths, function(paths) {
  list(
    plus = import(paths[["plus"]], "BigWig", as = "Rle"),
    minus = import(paths[["minus"]], "BigWig", as = "Rle")
  )
})

# Analysis 02 deliberately provides this metadata-rich CSV separately from its
# standards-compliant BED export. The CSV contains gene and region annotations.
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

bs.p <- bs[strand(bs) == "+"]
bs.m <- bs[strand(bs) == "-"]
for (sample in names(tracks)) {
  mcols(bs.p)[[sample]] <- sum_track_over_ranges(tracks[[sample]]$plus, bs.p)
  mcols(bs.m)[[sample]] <- sum_track_over_ranges(tracks[[sample]]$minus, bs.m)
}
bs_combined <- c(bs.p, bs.m)

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

new_columns <- signal_data %>%
  select(
    seqnames, start, end, strand, KO1, KO2,
    WT_mean, KO_mean, fold_change, pass_threshold
  )

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

dir.create(dirname(output_workbook), recursive = TRUE, showWarnings = FALSE)
write.xlsx(updated_data, file = output_workbook, overwrite = TRUE)

message("Strict WT/KO ratio >1.5 sites: ", sum(signal_data$pass_threshold))
message(
  "Genes represented by strict sites: ",
  length(unique(signal_data$gene_name[signal_data$pass_threshold]))
)
