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
