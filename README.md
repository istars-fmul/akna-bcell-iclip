# AKNA B-cell iCLIP2 analysis

This repository contains the computational analysis of AKNA iCLIP2 in mouse
in vitro-induced germinal center B (iGB) cells. The experiment comprises three wild-type (WT1, WT2, and WT3) and two *Akna* knockout (KO1 and KO2) replicates.


## Analysis stages

1. `workflow/Snakefile` performs raw-read QC, barcode-region quality
   filtering, demultiplexing and adapter trimming, STAR mapping, UMI-aware
   deduplication, WT/KO group merging, crosslink-track generation, library
   diagnostics, iCLIPro diagnostics, and no-input-control PureCLIP calling.
   The default PureCLIP targets are WT1, WT2, WT3, and pooled WT (`replicates`). KO samples are processed into tracks but are not PureCLIP
   input controls or default PureCLIP targets.
2. `analysis/01` through `analysis/05` perform PureCLIP postprocessing,
   replicate reproducibility and annotation, WT/KO signal comparison, GO analysis, binding-site region plotting, and IGV preparation.

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

Analysis `04` contacts Ensembl BioMart. Analysis `02` uses local GENCODE M25 annotation and does not perform GO enrichment. Analysis `05` uses the
signal-comparison workbook from analysis `03` and the replicate tracks for
IGV preparation; it does not require GO mappings or reference annotation.

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

See [`docs/software.md`](docs/software.md) for container details and
[R-package installation instructions](docs/software.md#r-packages).

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

Snakemake's older `--use-singularity` flag is an alias for enabling the same Singularity/Apptainer deployment backend. Additional bind arguments may be needed when configured inputs or outputs are outside the repository.

The default outputs will be under `results/preprocessing/akna/` and include FastQC reports, filtered/demultiplexed reads, mapped and deduplicated BAMs, pooled WT and KO BAMs, strand-specific RPM BigWigs, diagnostics, and PureCLIP BED files.

## Run downstream analyses

Analyses `01`, `02`, and `05` are Jupyter notebooks using the **R (IRkernel)** kernel. See [R-kernel setup](docs/software.md#jupyter-and-the-r-kernel) before opening them. Start JupyterLab from the repository root:

```bash
jupyter lab
```

Open and run `analysis/01_pureclip_postprocessing_noKO.ipynb`, followed by
`analysis/02_reproducibility_annotation.ipynb`. Select the **R** kernel, edit the first code cell (`params`) to specify input/output paths, then restart the kernel and run all cells in order. Relative paths are resolved from the repository root whether the kernel starts there or in `analysis/`.

The default analysis `01` inputs are pooled-WT tracks and PureCLIP positions.
Its output is the default input to analysis `02`. If you change an output path, update the corresponding input in the following analysis.

Analysis `02` exports an annotated BED, a metadata-rich CSV, and an Excel
workbook with one sheet named `all_reproducible_binding_sites`. These contain the reproducible annotated sites before WT/KO enrichment filtering.

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

Run the corrected GO analysis, which also generates the manuscript GO figure:

```bash
Rscript analysis/04_go_wt_ko_enriched.R \
  results/analysis/03_signal_comparison/reproducible_binding_sites_groups_noKO_signal_comparison.xlsx \
  data/reference/gencode.vM25.annotation.gtf \
  results/analysis/04_go
```

After analyses `03` and `04` finish, open
`analysis/05_manuscript_figures.ipynb`, check its input/output configuration cell, and run all cells with the R kernel. Analyses `03` and `04` remain standalone R scripts with positional arguments documented in their headers.

## Inspect loci in IGV

Analysis `05` writes the following files under its configured output directory (default: `results/analysis/05_manuscript_figures/`):

- `igv_loci.csv`: nine manually curated loci with 1-based, inclusive navigation coordinates. 
- `final_WT_KO_enriched_binding_sites.bed`: the final WT-enriched site set in BED6 format (0-based start, exclusive end).
- `WT.average.plus.bw` and `WT.average.minus.bw`: per-position arithmetic means of the three individually normalized WT replicates.
- `KO.average.plus.bw` and `KO.average.minus.bw`: corresponding means of the two individually normalized KO replicates.

Set `wt_tracks_dir` and `ko_tracks_dir` in the notebook's first cell. Both default to the workflow crosslink-track directory; they can differ for existing datasets.
The averaging step requires all five replicates, with both strands. Within each group and strand, it exports chromosomes shared by the input replicates and checks that their lengths agree. Scaffolds absent from any input are reported and omitted, matching the existing average-track convention. Uncovered positions on shared chromosomes contribute zero to the average.

In IGV:

1. Select the mouse **mm10** genome.
2. Use **File → Load from File** to open the four average BigWigs and the BED.
3. Optionally load the GENCODE M25 GTF for the annotation used in the analysis.
4. Paste an `igv_locus` value from `igv_loci.csv` into the search box.
5. Compare WT and KO on the same strand using a shared vertical scale; adjust zoom and scaling interactively as needed.


## Repository layout

```text
config/                 Example experiment configuration
workflow/Snakefile      Preprocessing and no-input PureCLIP workflow
workflow/scripts/       Retained workflow helper scripts
workflow/utils/         Shared R helpers
analysis/01-05          Ordered downstream analyses
docs/software.md        Software and container details
```

Generated data and reports belong under ignored `data/` and `results/`
directories.

## Citation

If this repository contributes to your research, please acknowledge our work by citing the associated publication.

## Contact

For questions about the workflow or analysis code, open an issue in this repository.
