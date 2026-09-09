# Software environment

## Verified versions

| Component | Version | Deployment in the workflow |
|---|---:|---|
| Snakemake | 8.10.7 | Host workflow manager |
| FastQC | 0.11.5 | `biocontainers/fastqc:v0.11.5_cv4` |
| FASTX-Toolkit | 0.0.14 | `biocontainers/fastx-toolkit:v0.0.14-6-deb_cv1` |
| Seqtk | 1.3 | `biocontainers/seqtk:v1.3-1-deb_cv1` |
| Flexbar | 3.4.0 | `biocontainers/flexbar:v1-3.4.0-2-deb_cv1` |
| STAR | 2.5.4b | `istarsfmul/star:2.5.4b` |
| SAMtools | 1.9 | `biocontainers/samtools:v1.9-4-deb_cv1` |
| BEDTools | 2.27.1 | `biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1` |
| UMI-tools | 1.1.2 | Digest-pinned image shown below |
| iCLIPro | 0.1.1_patch | `istarsfmul/iclipro:0.1.1_patch` |
| PureCLIP | 1.3.1 | `istarsfmul/pureclip:1.3.1` |
| R | 4.4.0 | Downstream analysis environment |
| topGO | 2.56.0 | R/Bioconductor package |
| biomaRt | 2.60.1 | R/Bioconductor package |

The immutable UMI-tools reference is:

```text
docker://flomicsbiotech/umitools_last_version@sha256:d8e4554ece88a0cededd1cbe024cdf767df15c47a7c6a90ed778ab2ada03bd58
```

The workflow also uses Python helper containers and the following immutable
Kent utilities image for `bedGraphToBigWig`:

```text
docker://abralab/kentutils@sha256:2cad15e98b9c4dd7da4ee4c9504ba9ec26a5cd0322fcd4454e4321d3259512d1
```

These supporting containers do not alter the named manuscript tool versions.

## R packages

The downstream analyses load packages explicitly in each file. In addition to
`topGO` and `biomaRt`, they require packages including `IRkernel`,
`rtracklayer`, `GenomicRanges`, `GenomicFeatures`, `AnnotationDbi`, `GO.db`,
`Rgraphviz`, `dplyr`, `tidyr`, `readr`, `readxl`, `openxlsx`, `writexl`,
`ggplot2`, `UpSetR`, `reshape2`, `forcats`, `stringr`, `scales`, and `svglite`.

Use R 4.4.0 with a matching Bioconductor release and verify package versions
with `sessionInfo()`. Analysis `05` prints session information in its final
notebook. BioMart-backed analyses require network access and may vary if the
live Ensembl annotation service changes; the GENCODE M25 protein-coding
background and manuscript count checks constrain the analysis input.

## Jupyter and the R kernel

Install JupyterLab in your chosen Python environment:

```bash
python -m pip install jupyterlab
```

With that environment active and `jupyter` on `PATH`, use the R installation
containing the analysis packages to install and register IRkernel:

```r
install.packages("IRkernel", repos = "https://cloud.r-project.org")
IRkernel::installspec(user = TRUE, name = "ir", displayname = "R")
```

Check registration with `jupyter kernelspec list`, then launch `jupyter lab`
from the repository root. Open an analysis notebook and select **R**. Check
`R.version.string` and `sessionInfo()` in the kernel to confirm that it uses
the intended R installation and packages. IRkernel provides the notebook
interface; the R packages listed above are still required for the analyses.

The notebooks do not require R Markdown rendering or Pandoc. They are stored
with empty outputs and execution counts; run cells in order after setting the
paths in the first code cell.

## Apptainer

Snakemake 8.10.7 enables containers with:

```bash
--software-deployment-method apptainer
```

The older `--use-singularity` spelling is retained by Snakemake as an alias.
When input/output paths are moved outside the repository, configure Apptainer
bind mounts so every rule can see the configured paths.
