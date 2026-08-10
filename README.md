# AKNA B-cell iCLIP2 analysis

This repository contains the computational analysis of AKNA iCLIP2 in mouse
in vitro-induced germinal center B (iGB) cells. The experiment comprises three
wild-type (WT1, WT2, and WT3) and two *Akna* knockout (KO1 and KO2) replicates.


## Analysis stages

1. `workflow/Snakefile` performs raw-read QC, barcode-region quality
   filtering, demultiplexing and adapter trimming, STAR mapping, UMI-aware
   deduplication, WT/KO group merging, crosslink-track generation, library
   diagnostics, iCLIPro diagnostics, and no-input-control PureCLIP calling.
   The default PureCLIP targets are WT1, WT2, WT3, and pooled WT
   (`replicates`). KO samples are processed into tracks but are not PureCLIP
   input controls or default PureCLIP targets.
2. `analysis/01` through `analysis/05` perform PureCLIP postprocessing,
   replicate reproducibility and annotation, WT/KO signal comparison, GO analysis, and manuscript figure exploration.

## Required inputs

Place inputs at the paths in `config/config.example.yaml`, or copy that file
and change the paths. The example expects:

| File | Purpose |
|---|---|
| `data/raw/akna.fastq.gz` | Raw AKNA iCLIP2 sequencing run |
| `data/raw/barcodes.fasta` | Exact experiment-specific Flexbar barcode definitions, including UMI positions |
| `data/raw/adapter.fasta` | Exact experiment-specific adapter sequence in FASTA format |
| `data/reference/GRCm38.primary_assembly.genome.fa` | Mouse GRCm38 primary-assembly genome FASTA |
| `data/reference/gencode.vM25.annotation.gtf` | GENCODE mouse release M25 annotation |

The barcode and adapter FASTA files are user-supplied and ignored by Git. The
repository intentionally contains neither their exact sequences nor example
substitutes. GRCm38 and GENCODE M25 must be downloaded from their respective
providers. The workflow creates the FASTA index and STAR index if absent.

Analyses `04` and `05` contact Ensembl BioMart. Analysis `05` also accepts the
optional external metagene and STRING-summary files shown in its YAML
parameters; these are not distributed.

## Software

Workflow tools run from rule-specific containers through Apptainer. The
versions used for this analysis are:

| Software | Version |
|---|---:|
| Snakemake | 8.10.7 |
| FastQC | 0.11.5 |
| FASTX-Toolkit | 0.0.14 |
| Seqtk | 1.3 |
| Flexbar | 3.4.0 |
| STAR | 2.5.4b |
| SAMtools | 1.9 |
| BEDTools | 2.27.1 |
| UMI-tools | 1.1.2 |
| iCLIPro | 0.1.1_patch |
| PureCLIP | 1.3.1 |
| R | 4.4.0 |
| topGO | 2.56.0 |
| biomaRt | 2.60.1 |

See [`docs/software.md`](docs/software.md) for container and R-package notes.

## Run preprocessing

From the repository root, inspect and edit the example configuration, then
run:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/config.example.yaml \
  --software-deployment-method apptainer \
  --cores 8
```

Snakemake's older `--use-singularity` flag is an alias for enabling the same
Singularity/Apptainer deployment backend. Additional bind arguments may be
needed when configured inputs or outputs are outside the repository.

The default outputs are under `results/preprocessing/akna/` and include FastQC
reports, filtered/demultiplexed reads, mapped and deduplicated BAMs, pooled WT
and KO BAMs, strand-specific RPM BigWigs, diagnostics, and PureCLIP BED files.
All are ignored by Git.

## Run downstream analyses

Render R Markdown files from the repository root. Parameters shown below are
the defaults and can be overridden with `params`.

```bash
Rscript -e 'rmarkdown::render("analysis/01_pureclip_postprocessing_noKO.Rmd")'

Rscript -e 'rmarkdown::render("analysis/02_reproducibility_annotation.Rmd")'
```

To render with another input or output, pass named parameters, for example:

```bash
Rscript -e 'rmarkdown::render(
  "analysis/01_pureclip_postprocessing_noKO.Rmd",
  params = list(output_bed = "results/analysis/01/custom_curated.bed")
)'
```

Run the WT/KO comparison with the annotated CSV from analysis `02`, its
workbook, an output workbook, and explicit plus/minus tracks in this fixed
order:

```bash
Rscript analysis/03_wt_ko_signal_comparison.R \
  results/analysis/02_reproducibility/bs_clean_noKO.csv \
  results/analysis/02_reproducibility/reproducible_binding_sites_groups_noKO.xlsx \
  results/analysis/03_signal_comparison/reproducible_binding_sites_groups_noKO_signal_comparison.xlsx \
  results/preprocessing/akna/crosslinked_nucleotides/WT1.plus.bw \
  results/preprocessing/akna/crosslinked_nucleotides/WT1.minus.bw \
  results/preprocessing/akna/crosslinked_nucleotides/WT2.plus.bw \
  results/preprocessing/akna/crosslinked_nucleotides/WT2.minus.bw \
  results/preprocessing/akna/crosslinked_nucleotides/WT3.plus.bw \
  results/preprocessing/akna/crosslinked_nucleotides/WT3.minus.bw \
  results/preprocessing/akna/crosslinked_nucleotides/KO1.plus.bw \
  results/preprocessing/akna/crosslinked_nucleotides/KO1.minus.bw \
  results/preprocessing/akna/crosslinked_nucleotides/KO2.plus.bw \
  results/preprocessing/akna/crosslinked_nucleotides/KO2.minus.bw
```

Run the corrected GO analysis:

```bash
Rscript analysis/04_go_wt_ko_enriched.R \
  results/analysis/03_signal_comparison/reproducible_binding_sites_groups_noKO_signal_comparison.xlsx \
  data/reference/gencode.vM25.annotation.gtf \
  results/analysis/04_go
```

Render the manuscript exploration after checking the YAML parameters at the
top of the notebook:

```bash
Rscript -e 'rmarkdown::render("analysis/05_manuscript_figures.Rmd")'
```

Analyses `01`, `02`, and `05` declare their configurable paths in the YAML
headers. Analyses `03` and `04` document positional arguments in their file
headers and fail on an incorrect argument count.

## Expected manuscript checks

The historical manuscript analysis reported:

| Check | Sites | Genes |
|---|---:|---:|
| Reproducible annotated 9-nt sites | 2,853 | 1,351 |
| Strict WT/KO ratio `> 1.5` | 2,412 | 1,243 |

These values are validation targets from the source analysis, not results newly
regenerated in this repository. Compare them after a complete rerun. Analysis
`04` stops if its input does not contain the expected 2,412 strict sites and
1,243 GENCODE-background genes.

## Repository layout

```text
config/                 Example experiment configuration
workflow/Snakefile      Preprocessing and no-input PureCLIP workflow
workflow/scripts/       Retained workflow helper scripts
workflow/utils/         Shared R helpers
analysis/01-05          Ordered downstream analyses
docs/provenance.md      Source revisions and deliberate adjustments
docs/software.md        Software and container details
```

Generated data and reports belong under ignored `data/` and `results/`
directories and must not be committed.

## Citation

Please cite the associated manuscript when its publication details become
available. The repository URL, release tag, and sequencing-data accession will
be added to the manuscript data-availability statement before publication.

## Contact

For questions about the workflow or analysis code, open an issue in this
repository.
