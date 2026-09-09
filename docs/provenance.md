# Provenance

## Source revisions

The publication code was assembled from committed content at these immutable
source revisions:

| Source repository | Branch | Commit |
|---|---|---|
| [istars-fmul/proj-rbp-akna](https://github.com/istars-fmul/proj-rbp-akna) | `figures` | `5688ea25bde5b2844057043996251160e1dce40d` (`5688ea2`) |
| [istars-fmul/iclipseq_pipeline](https://github.com/istars-fmul/iclipseq_pipeline) | `post_peak_calling_read_density_validation` | `8ba100eb5b29e8fd07514581ec34c4570671641d` (`8ba100e`) |

Analysis sources from `proj-rbp-akna` were:

| Publication file | Source path |
|---|---|
| `analysis/01_pureclip_postprocessing_noKO.ipynb` | `notebooks/bcells/peakCalling_withoutKOinput_signalComparison/PureCLIP_Postprocessing_noKO.Rmd` |
| `analysis/02_reproducibility_annotation.ipynb` | `notebooks/bcells/peakCalling_withoutKOinput_signalComparison/akna-pureclip-reproducibility-genomic-analysis_allPeaks_noKO.Rmd` |
| `analysis/03_wt_ko_signal_comparison.R` | `notebooks/bcells/peakCalling_withoutKOinput_signalComparison/signal_comparison.R` |
| `analysis/04_go_wt_ko_enriched.R` | `notebooks/bcells/paper_figures/rerun_go_wt_ko_enriched.R` |
| `analysis/05_manuscript_figures.ipynb` | `notebooks/bcells/paper_figures/akna_bcell_noKO_paper_figures.Rmd` |

The workflow is based on `pipeline/Snakefile` and its two retained Python helper
scripts at the `iclipseq_pipeline` revision above.

## Deliberate adjustments

- Replaced machine-specific paths with repository-relative configuration,
  notebook configuration cells, or documented R-script arguments.
- Converted analyses 01, 02, and 05 from R Markdown to Jupyter notebooks
  using IRkernel. The conversion preserves the publication R Markdown analysis
  chunks; only parameter initialization and notebook-specific path setup change.
  Source paths above refer to the original R Markdown analyses, which remain
  available in the source repository.
- Replaced Snakemake's deprecated `singularity:` rule directive with the
  runtime-agnostic `container:` directive without changing the images.
- Trimmed the workflow default to the coherent preprocessing, diagnostics, and
  no-input-control PureCLIP path used by the manuscript. Broken read-density
  comparison and mismatched old downstream workflow rules were not copied.
- Fixed the singular/plural sequencing-run configuration mismatch by using
  `sequencing_runs` consistently.
- Changed Flexbar adapter handling from the single-sequence `--adapter-seq`
  argument to `--adapters`, because the configured input is an adapter FASTA.
- Pinned UMI-tools 1.1.2 by immutable image digest instead of a mutable
  `latest` tag.
- Pinned the unchanged Kent utilities image by immutable digest instead of its
  mutable `latest` tag.
- Replaced shared mutable diagnostics output with one CSV per sample followed
  by deterministic collation; the reported metrics are unchanged.
- Replaced filesystem-order-dependent BigWig discovery with explicit WT1,
  WT2, WT3, KO1, and KO2 plus/minus mappings.
- Corrected WT3 reproducibility status to use `aa$cutoff[3]`, rather than the
  WT2 value at `aa$cutoff[2]`.
- Kept the annotated BED and metadata-rich CSV as distinct outputs. Signal
  comparison reads the CSV rather than incorrectly treating the BED as CSV.
- Named the primary workbook sheet `all_reproducible_binding_sites` so the
  corrected GO and manuscript-figure analyses consume the preceding stage
  without a manual sheet rename.
- Kept the WT/KO formula `(WT_mean + 0.01) / (KO_mean + 0.01)` and the strict
  `> 1.5` threshold.
- Made manuscript figure and exploratory STRING sections configurable rather
  than deleting them. STRING network generation remains external.
- Made diagnostics sample-local before final collation to remove concurrent
  writes to one temporary CSV.

No source data, generated result, figure, barcode sequence, adapter sequence,
or reference file was copied.

## Rerun comparison

After a complete rerun, compare outputs to the manuscript checks of 2,853
reproducible annotated 9-nt sites / 1,351 genes and 2,412 strict WT/KO-ratio
sites / 1,243 genes. These counts document the source analysis and have not
been regenerated merely by assembling this repository.
