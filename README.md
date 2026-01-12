# Mammary microbiota under heat stress 🐄🔥🧫

**Multi-omics and in vitro validation of heat stress–associated mammary inflammation in dairy cows**

To elucidate the mechanisms underlying heat stress–induced mammary inflammation in dairy cows, we integrated **proteomic analysis of blood from the mammary vein**, **time-series metagenomic analysis of milk samples**, and **in vitro co-culture of mammary epithelial cells and microorganisms**.

---

## Overview

This repository provides analysis scripts and workflows used to investigate host–microbiota interactions in the bovine mammary gland under heat stress conditions. The study combines:

- Proteomics of mammary-vein blood  
- Time-series milk metagenomics  
- Mammary epithelial cell RNA-seq  
- In vitro co-culture experiments with live or heat-killed microorganisms  

---

## Highlights

- **Time-series milk metagenomics workflow**  
  QC → host depletion → taxonomic profiling

- **Transcriptomics workflow (RNA-seq)**  
  QC / trimming → quantification (kallisto) and/or alignment-based analysis

- **Figure-ready R scripts**  
  For reproducible visualization and manuscript-quality figures

- **Shell-first, HPC-friendly design**  
  Easy to run with `nohup`, `sbatch`, or GNU parallel

---

## Repository contents

Key scripts included in this repository:

- `trimmomatic.sh` — Read trimming and quality control (short reads)
- `hostfree_bowtie.sh` — Host read removal using Bowtie2
- `metaphlan-rawdata.sh` — MetaPhlAn taxonomic profiling
- `run_RNAseq.sh` — RNA-seq processing pipeline (wrapper)
- `run_kallisto.sh` — RNA-seq quantification using kallisto
- `codes for Figure plot.R` — R scripts for plotting and visualization

> **Tip**  
> If you later add directories such as `data/`, `results/`, or `figures/`, this README will remain valid; only file paths in the examples need to be updated.

---

## Software requirements

The workflows were designed for **Linux / HPC environments**.

### Required tools

- **Trimmomatic**  
  Adapter removal and quality trimming of sequencing reads

- **Bowtie2**  
  Host read removal by mapping reads to the host reference genome

- **MetaPhlAn**  
  Taxonomic profiling of host-depleted metagenomic reads

- **kallisto**  
  Alignment-free quantification of RNA-seq data

- **R**  
  Downstream analysis and figure generation  
  Required packages include (but are not limited to):
  - `tidyverse`
  - `ggplot2`

### Reproducibility note

> **Note**  
> For full reproducibility, we strongly recommend recording exact software versions.  
> This can be achieved by:
>
> ```bash
> conda env export > env.yml
> ```
>
> or by using a containerized environment (Docker or Singularity).

---

## Workflow A — Milk metagenomics (time-series)

### A1. Quality control and trimming

Trim adapters and low-quality bases using `trimmomatic.sh`.

```bash
bash trimmomatic.sh \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o results/01_trimmed/
```
Expected output
- Cleaned paired-end FASTQ files for downstream profiling

### A2. Host read removal (host depletion)
Remove host-derived reads by mapping to the host reference genome.

```bash
bash hostfree_bowtie.sh \
  -1 results/01_trimmed/sample.clean.R1.fastq.gz \
  -2 results/01_trimmed/sample.clean.R2.fastq.gz \
  -x /path/to/host_bowtie2_index \
  -o results/02_hostfree/
```
Expected output
- Host-depleted paired FASTQ files (recommended input for MetaPhlAn)
- Mapping statistics and logs for quality assessment

### A3. Taxonomic profiling (MetaPhlAn)
Generate taxonomic profiles using MetaPhlAn.

```bash
bash metaphlan-rawdata.sh \
  -i results/02_hostfree/ \
  -o results/03_metaphlan/
```
Expected output
- Per-sample MetaPhlAn taxonomic tables
- Optional merged abundance tables (depending on script settings)

---

## Workflow B — Mammary epithelial cell RNA-seq (in vitro)
### B1. RNA-seq processing
Run the RNA-seq workflow wrapper.

```bash
bash run_RNAseq.sh \
  -i rnaseq_fastq/ \
  -o results/rnaseq/
```
### B2. Quantification with kallisto (optional / recommended)
Perform pseudoalignment-based RNA-seq quantification.

```bash 
bash run_kallisto.sh \
  -i rnaseq_fastq/ \
  -x /path/to/kallisto_index \
  -o results/04_kallisto/
```
Expected output
- Transcript- or gene-level abundance estimates
- Log files for quality control and troubleshooting

---

## Workflow C — Figures and visualization
Reproduce manuscript-style figures using the R plotting scripts.

```bash
Rscript "codes for Figure plot.R"
```
If the plotting script expects specific input filenames, consider documenting the input table format or adding a small config.R file.

Suggested project structure (optional)
```text
.
├── data/
│   ├── metagenome_fastq/
│   ├── rnaseq_fastq/
│   └── metadata/
├── results/
│   ├── 01_trimmed/
│   ├── 02_hostfree/
│   ├── 03_metaphlan/
│   ├── rnaseq/
│   └── figures/
└── scripts/
    ├── trimmomatic.sh
    ├── hostfree_bowtie.sh
    ├── metaphlan-rawdata.sh
    ├── run_RNAseq.sh
    ├── run_kallisto.sh
    └── codes for Figure plot.R
```

---

## Reproducibility checklist
- Record software versions (conda env export, sessionInfo() in R)

- Maintain a complete sample metadata table (metadata.tsv)

- Save host depletion statistics and sequencing depth summaries

- Generate merged abundance matrices for downstream analyses (alpha/beta diversity, longitudinal modeling)

---

## Citation
If you use this workflow or code in your research, please cite the associated manuscript (DOI to be added upon publication).
Until then, please cite this GitHub repository.

---

## License
This project is released under the GPL-3.0 License.

---

## Contact
Shengtao Gao
📧 shengtaogao@163.com




