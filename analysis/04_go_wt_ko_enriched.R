#!/usr/bin/env Rscript

# Gene Ontology enrichment of WT-enriched AKNA target genes
# -------------------------------------------------------
# Use the signal-comparison workbook from analysis 03 to select target genes,
# then test GO enrichment against protein-coding genes from filtered GENCODE M25.
# The unit of this analysis is a unique gene, not a binding site. Multiple sites
# in one gene therefore contribute only one foreground entry.
#
# Each ontology (BP, CC, MF) is tested with topGO weight01 and Fisher's exact
# test. BH correction uses the test scores across all tested terms
# within that ontology. A curated subset of top significant terms is used for the
# compact manuscript figure; the complete tested-term tables are also exported.
#
# Outputs in OUT_DIR: all tested GO terms, compact figure source table, mapping
# snapshot, an Excel results workbook, and the compact plot as PDF/PNG/SVG.
# An optional previous plot-term CSV produces a comparison table only; it does
# not affect gene selection, GO testing, or the corrected plot.
#
# Usage:
# Rscript analysis/04_go_wt_ko_enriched.R SIGNAL.xlsx GENCODE.gtf OUT_DIR \
#   [PREVIOUS_PLOT_TERMS.csv]

## Load libraries and input/output arguments
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(biomaRt)
  library(dplyr)
  library(forcats)
  library(GenomicFeatures)
  library(ggplot2)
  library(GO.db)
  library(readr)
  library(readxl)
  library(rtracklayer)
  library(stringr)
  library(svglite)
  library(tidyr)
  library(topGO)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3, 4)) {
  stop(
    "Usage: 04_go_wt_ko_enriched.R SIGNAL.xlsx GENCODE.gtf OUT_DIR ",
    "[PREVIOUS_PLOT_TERMS.csv]"
  )
}
primary_xlsx <- args[1]
annotation_file <- args[2]
out_dir <- args[3]
old_plot_terms_file <- if (length(args) == 4) args[4] else NA_character_

stopifnot(file.exists(primary_xlsx), file.exists(annotation_file))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Select the WT-enriched foreground from the reproducible-site workbook
# The count checks below verify that the expected manuscript dataset was supplied.
bs <- read_excel(primary_xlsx, sheet = "all_reproducible_binding_sites") %>%
  mutate(pass_threshold = as.logical(pass_threshold))

if (sum(bs$pass_threshold %in% TRUE) != 2412) {
  stop(
    "Expected 2,412 reproducible sites with strict WT/KO ratio >1.5, but found ",
    sum(bs$pass_threshold %in% TRUE), "."
  )
}

# Collapse multiple binding sites to unique, nonempty gene symbols.
foreground_genes_unfiltered <- sort(unique(bs$gene_name[bs$pass_threshold %in% TRUE]))
foreground_genes_unfiltered <- foreground_genes_unfiltered[
  !is.na(foreground_genes_unfiltered) & foreground_genes_unfiltered != ""
]

## Construct the protein-coding background using the analysis-02 annotation filters
message("Loading and filtering GENCODE vM25 annotation")
anno <- rtracklayer::import(annotation_file, format = "GTF")
anno <- anno[anno$level != 3]
# Internal 0 retains records with a missing TSL attribute; it is not a GENCODE
# support category. Literal "NA" is mapped to 10 and excluded; levels 1-3 remain.
anno$transcript_support_level[is.na(anno$transcript_support_level)] <- 0
anno$transcript_support_level[anno$transcript_support_level == "NA"] <- 10
anno <- anno[anno$transcript_support_level %in% c(0, 1, 2, 3)]

# Build the transcript database, then recover gene symbols/types from the GTF.
anno_db <- GenomicFeatures::makeTxDbFromGRanges(anno)
gns <- GenomicFeatures::genes(anno_db)
idx <- match(gns$gene_id, anno$gene_id)
elementMetadata(gns) <- cbind(elementMetadata(gns), elementMetadata(anno)[idx, ])

background_genes <- sort(unique(gns$gene_name[gns$gene_type == "protein_coding"]))
background_genes <- background_genes[
  !is.na(background_genes) & background_genes != ""
]
foreground_genes <- sort(intersect(foreground_genes_unfiltered, background_genes))

if (length(foreground_genes) != 1243) {
  stop(
    "Expected 1,243 WT/KO-enriched target genes, but found ",
    length(foreground_genes), "."
  )
}

## Retrieve gene-to-GO annotation through Ensembl BioMart
# Try mirrors if the primary host is unavailable and record the successful host.
message("Connecting to Ensembl BioMart")
mart_hosts <- c(
  "https://www.ensembl.org",
  "https://useast.ensembl.org",
  "https://uswest.ensembl.org",
  "https://asia.ensembl.org"
)
mart <- NULL
mart_host_used <- NA_character_
for (mart_host in mart_hosts) {
  mart <- tryCatch(
    biomaRt::useMart(
      "ENSEMBL_MART_ENSEMBL",
      dataset = "mmusculus_gene_ensembl",
      host = mart_host
    ),
    error = function(e) NULL
  )
  if (!is.null(mart)) {
    mart_host_used <- mart_host
    break
  }
}
if (is.null(mart)) stop("Could not connect to Ensembl BioMart or its mirrors.")

message("Retrieving GO mappings for the fixed GENCODE protein-coding background")
go_mapping <- biomaRt::getBM(
  attributes = c("go_id", "external_gene_name", "namespace_1003"),
  filters = "external_gene_name",
  values = background_genes,
  mart = mart
) %>%
  filter(go_id != "", external_gene_name != "") %>%
  distinct(external_gene_name, go_id, namespace_1003) %>%
  arrange(external_gene_name, go_id)

# topGO expects a named list of GO identifiers for each gene. Genes without
# mapped GO terms cannot contribute to the ontology-specific enrichment tests.
gene_2_GO <- unstack(go_mapping[, c("go_id", "external_gene_name")])
foreground_genes_with_go <- sort(intersect(foreground_genes, names(gene_2_GO)))

## Run one ontology: foreground membership, GO graph, testing, then BH correction
run_topgo <- function(ontology) {
  gene_list <- factor(as.integer(background_genes %in% foreground_genes_with_go))
  names(gene_list) <- background_genes

  go_data <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = gene_list,
    annot = annFUN.gene2GO,
    gene2GO = gene_2_GO,
    nodeSize = 20
  )
  weight_fisher <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
  result <- GenTable(
    go_data,
    weightFisher = weight_fisher,
    orderBy = "weightFisher",
    topNodes = length(usedGO(go_data))
  )
  # GenTable formats P values for display. Retrieve numeric scores from the test
  # object instead, so rounding or strings such as "< 1e-30" do not affect BH.
  exact_p <- score(weight_fisher)

  result %>%
    transmute(
      ontology = recode(
        ontology,
        BP = "Biological process",
        CC = "Cellular component",
        MF = "Molecular function"
      ),
      GO.ID,
      Term,
      Annotated_background_genes = as.numeric(Annotated),
      Observed_AKNA_target_genes = as.numeric(Significant),
      Expected_AKNA_target_genes = as.numeric(Expected),
      # Fold enrichment is observed / expected target-gene count, not signal ratio.
      Fold_enrichment = Observed_AKNA_target_genes / Expected_AKNA_target_genes,
      weightFisher_P = unname(exact_p[GO.ID])
    ) %>%
    mutate(
      BH_adjusted_P = p.adjust(weightFisher_P, method = "BH"),
      minus_log10_adjusted_P = -log10(
        pmax(BH_adjusted_P, .Machine$double.xmin)
      )
    ) %>%
    arrange(BH_adjusted_P, weightFisher_P)
}

# Each call adjusts its own ontology before the three result tables are combined.
message("Running topGO for BP, CC, and MF")
go_all <- bind_rows(lapply(c("BP", "CC", "MF"), run_topgo))
# Recover full GO term names rather than the abbreviated labels in GenTable.
go_term_lookup <- AnnotationDbi::mapIds(
  GO.db::GO.db,
  keys = unique(go_all$GO.ID),
  column = "TERM",
  keytype = "GOID",
  multiVals = "first"
)
go_all <- go_all %>%
  mutate(Term = coalesce(unname(go_term_lookup[GO.ID]), Term))
go_significant <- go_all %>% filter(BH_adjusted_P <= 0.05)

## Select representative functional themes for the manuscript plot
# This display list is applied after testing; it does not restrict the tests or
# the multiple-testing correction. The complete tables retain other GO terms.
go_theme_terms <- tribble(
  ~functional_theme, ~GO.ID,
  "Translation / ribosome", "GO:0022627",
  "Translation / ribosome", "GO:0042274",
  "Translation / ribosome", "GO:0045727",
  "Translation / ribosome", "GO:0002183",
  "Translation / ribosome", "GO:0006446",
  "Translation / ribosome", "GO:0003743",
  "Translation / ribosome", "GO:0043022",
  "Translation / ribosome", "GO:0032040",
  "Nuclear / chromatin", "GO:0005654",
  "Nuclear / chromatin", "GO:0005730",
  "Nuclear / chromatin", "GO:0016607",
  "Nuclear / chromatin", "GO:0000785",
  "Nuclear / chromatin", "GO:0000122",
  "Nuclear / chromatin", "GO:0045944",
  "Nuclear / chromatin", "GO:0006606",
  "Nuclear / chromatin", "GO:0003682",
  "B-cell / immune function", "GO:0050853",
  "B-cell / immune function", "GO:0019886",
  "B-cell / immune function", "GO:0001782",
  "B-cell / immune function", "GO:0030183",
  "B-cell / immune function", "GO:0042611",
  "B-cell / immune function", "GO:0050869",
  "B-cell / immune function", "GO:0030888",
  "Cytoskeleton / microtubule", "GO:0072686",
  "Cytoskeleton / microtubule", "GO:0008608",
  "Cytoskeleton / microtubule", "GO:0005813",
  "Cytoskeleton / microtubule", "GO:0015629",
  "Cytoskeleton / microtubule", "GO:0034314",
  "Cytoskeleton / microtubule", "GO:0030027",
  "Cytoskeleton / microtubule", "GO:0051015",
  "Protein homeostasis / complexes", "GO:0050821",
  "Protein homeostasis / complexes", "GO:0031625",
  "Protein homeostasis / complexes", "GO:0044877",
  "Protein homeostasis / complexes", "GO:0032991",
  "Protein homeostasis / complexes", "GO:0032436",
  "Protein homeostasis / complexes", "GO:0051082",
  "Protein homeostasis / complexes", "GO:0140534"
)

# Generic binding terms are omitted from the compact display only.
generic_go_ids <- c(
  "GO:0003729", "GO:0019843", "GO:0003723",
  "GO:0003730", "GO:0005515", "GO:0042802"
)

# Within each curated theme, choose up to three significant terms by raw P value.
go_functional <- go_significant %>%
  filter(!GO.ID %in% generic_go_ids) %>%
  inner_join(go_theme_terms, by = "GO.ID") %>%
  filter(Expected_AKNA_target_genes > 0) %>%
  group_by(functional_theme) %>%
  slice_min(weightFisher_P, n = 3, with_ties = FALSE) %>%
  ungroup()

# Of those representative terms, display only enrichment >=3; wrap labels for
# the figure and preserve single-line labels in the exported source table.
go_compact <- go_functional %>%
  filter(Fold_enrichment >= 3) %>%
  arrange(Fold_enrichment) %>%
  mutate(
    Term_for_plot = str_wrap(Term, width = 52),
    term_label_export = str_replace_all(Term_for_plot, fixed("\n"), " "),
    Term_for_plot = factor(Term_for_plot, levels = Term_for_plot)
  )

go_compact_export <- go_compact %>%
  arrange(desc(Fold_enrichment), weightFisher_P) %>%
  transmute(
    functional_theme, ontology, GO.ID, Term,
    Term_for_plot = term_label_export,
    Annotated_background_genes,
    Observed_AKNA_target_genes,
    Expected_AKNA_target_genes,
    Fold_enrichment,
    weightFisher_P,
    BH_adjusted_P,
    minus_log10_adjusted_P
  )

# Also export every significant >=3-fold term, independent of the display list.
go_ge3 <- go_significant %>%
  filter(Fold_enrichment >= 3) %>%
  arrange(desc(Fold_enrichment), BH_adjusted_P)

## Record selection, background, test settings and software in the results workbook
analysis_summary <- tibble(
  parameter = c(
    "Site selection", "WT/KO threshold",
    "Reproducible binding sites passing threshold",
    "Foreground target genes before annotation intersection",
    "Foreground target genes in GENCODE background",
    "Foreground target genes with GO annotation",
    "Protein-coding GENCODE background genes",
    "topGO algorithm", "Statistical test", "Minimum GO node size",
    "Multiple-testing correction", "Display fold-enrichment threshold",
    "BioMart host used", "topGO version", "biomaRt version", "R version"
  ),
  value = c(
    "Reproducible WT-defined sites passing WT/KO signal filtering",
    "> 1.5", as.character(sum(bs$pass_threshold %in% TRUE)),
    as.character(length(foreground_genes_unfiltered)),
    as.character(length(foreground_genes)),
    as.character(length(foreground_genes_with_go)),
    as.character(length(background_genes)),
    "weight01", "Fisher exact test", "20",
    "Benjamini-Hochberg within each ontology across all tested terms",
    ">= 3", mart_host_used, as.character(packageVersion("topGO")),
    as.character(packageVersion("biomaRt")), R.version.string
  )
)

## Optional comparison with a previous plot's term table
# Without that input, the comparison table simply identifies the new plot terms.
old_plot_terms <- if (
  !is.na(old_plot_terms_file) && file.exists(old_plot_terms_file)
) {
  read_csv(old_plot_terms_file, show_col_types = FALSE) %>%
    transmute(
      GO.ID, Term_old = Term,
      Observed_old = Observed_AKNA_target_genes,
      Expected_old = Expected_AKNA_target_genes,
      Fold_enrichment_old = Fold_enrichment,
      BH_adjusted_P_old = BH_adjusted_P_for_plot,
      old_plot = TRUE
    )
} else {
  tibble(
    GO.ID = character(), Term_old = character(), Observed_old = numeric(),
    Expected_old = numeric(), Fold_enrichment_old = numeric(),
    BH_adjusted_P_old = numeric(), old_plot = logical()
  )
}

new_plot_terms <- go_compact_export %>%
  transmute(
    GO.ID, Term_new = Term, Observed_new = Observed_AKNA_target_genes,
    Expected_new = Expected_AKNA_target_genes,
    Fold_enrichment_new = Fold_enrichment,
    BH_adjusted_P_new = BH_adjusted_P, new_plot = TRUE
  )

plot_comparison <- full_join(old_plot_terms, new_plot_terms, by = "GO.ID") %>%
  mutate(
    old_plot = replace_na(old_plot, FALSE),
    new_plot = replace_na(new_plot, FALSE),
    plot_membership = case_when(
      old_plot & new_plot ~ "retained",
      old_plot ~ "old plot only",
      new_plot ~ "corrected plot only"
    ),
    change_in_observed_genes = Observed_new - Observed_old,
    percent_change_in_fold_enrichment = 100 *
      (Fold_enrichment_new - Fold_enrichment_old) / Fold_enrichment_old
  ) %>%
  arrange(plot_membership, GO.ID)

target_genes <- tibble(gene_name = foreground_genes) %>%
  mutate(has_GO_annotation = gene_name %in% foreground_genes_with_go)

## Export numerical results and the gene-to-GO mapping used for this run
write_csv(go_all, file.path(out_dir, "go_enrichment_all_WT_KO_enriched_corrected.csv"))
write_csv(
  go_compact_export,
  file.path(out_dir, "go_compact_enrichment_ge3_WT_KO_enriched_plot_terms.csv")
)
write_csv(
  plot_comparison,
  file.path(out_dir, "go_compact_previous_vs_WT_KO_enriched_comparison.csv")
)
write_csv(
  go_mapping,
  file.path(out_dir, "go_gene_to_GO_mapping_WT_KO_enriched_run.csv")
)

write_xlsx(
  list(
    analysis_summary = analysis_summary,
    compact_plot_terms = go_compact_export,
    comparison_with_previous_plot = plot_comparison,
    significant_GO_terms_ge3 = go_ge3,
    significant_GO_terms = go_significant,
    all_tested_GO_terms = go_all,
    foreground_target_genes = target_genes
  ),
  file.path(out_dir, "AKNA_noKO_GO_enrichment_WT_KO_enriched_corrected.xlsx")
)

## Plot the compact enrichment summary
# X: observed/expected target genes; point area: observed target-gene count;
# colour: -log10(BH-adjusted P).
p_go_compact <- go_compact %>%
  ggplot(aes(x = Fold_enrichment, y = Term_for_plot)) +
  geom_point(
    aes(size = Observed_AKNA_target_genes, color = minus_log10_adjusted_P),
    alpha = 0.9
  ) +
  scale_color_gradient(
    low = "#DDEAF6", high = "#08306B",
    name = expression(-log[10](adjusted ~ italic(P)))
  ) +
  scale_size_area(
    max_size = 5.8, breaks = c(10, 20, 30),
    name = "Observed AKNA\ntarget genes"
  ) +
  labs(
    x = "Observed / expected AKNA target genes", y = NULL,
    title = "GO terms enriched >=3-fold among AKNA target mRNAs"
  ) +
  theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0),
    axis.line.y = element_blank(), axis.ticks.y = element_blank()
  )

plot_base <- file.path(
  out_dir,
  "figure_go_compact_enrichment_ge3_WT_KO_enriched_corrected"
)
ggsave(paste0(plot_base, ".pdf"), p_go_compact, width = 5.2, height = 3.8,
       useDingbats = FALSE)
ggsave(paste0(plot_base, ".png"), p_go_compact, width = 5.2, height = 3.8,
       dpi = 300)
ggsave(paste0(plot_base, ".svg"), p_go_compact, width = 5.2, height = 3.8,
       device = svglite::svglite)

message("Corrected GO analysis complete")
message("Foreground genes: ", length(foreground_genes))
message("Significant terms: ", nrow(go_significant))
message("Significant terms with fold enrichment >=3: ", nrow(go_ge3))
message("Terms in compact plot: ", nrow(go_compact_export))
