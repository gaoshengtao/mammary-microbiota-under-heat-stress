# mammary-microbiota-under-heat-stress
To elucidate the mechanisms underlying heat stress-induced mammary inflammation in dairy cows, we integrated proteomic analysis of blood from the mammary vein, time-series metagenomic analysis of milk samples, and in vitro co-culture of mammary epithelial cells and microorganisms. 
# Mammary microbiota under heat stress 🐄🔥🧫  
**Multi-omics + in vitro validation workflow for heat stress–associated mammary inflammation in dairy cows**

> To elucidate the mechanisms underlying heat stress–induced mammary inflammation in dairy cows, this project integrates  
> **(i) proteomics of mammary-vein blood**, **(ii) time-series milk metagenomics**, and **(iii) in vitro co-culture of mammary epithelial cells and microorganisms**.  

---

## Highlights
- **Time-series milk metagenomics workflow**: QC → host depletion → taxonomic profiling  
- **Transcriptomics workflow (RNA-seq)**: QC/trim → quantification (kallisto) and/or alignment-based pipeline  
- **Figure-ready R scripts** for reproducible plotting and downstream visualization  
- Shell-first, HPC-friendly (easy to wrap with `nohup`, `sbatch`, GNU parallel, etc.)

---

## Repository contents
Key scripts included in this repository:
- `trimmomatic.sh` — read trimming / quality control (short reads)
- `hostfree_bowtie.sh` — host read removal (Bowtie2-based host depletion)
- `metaphlan-rawdata.sh` — MetaPhlAn profiling from (host-depleted) reads
- `run_RNAseq.sh` — RNA-seq processing pipeline (alignment-based or wrapper)
- `run_kallisto.sh` — RNA-seq quantification with kallisto
- `codes for Figure plot.R` — plotting / figure generation in R

> Tip: If you later add folders (e.g., `data/`, `results/`, `figures/`), this README will still work—just update paths in the examples.

---

## Quick start (recommended order)
### 0) Prepare software
You’ll typically need:
- **Trimmomatic**
- **Bowtie2** (+ host genome index)
- **MetaPhlAn** (and its database configured)
- **kallisto** (optional, for RNA-seq quantification)
- **R** (tidyverse/ggplot2 etc., depending on your plotting script)

> ✅ Best practice: record exact versions with `conda env export > env.yml` (or a Docker/Singularity recipe) for full reproducibility.

---

## Workflow A — Milk metagenomics (time-series)
### A1) QC / trimming
Use `trimmomatic.sh` to trim adapters and low-quality bases.

**Example (pseudo-usage):**
```bash
bash trimmomatic.sh \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -o results/01_trimmed/
Expected output

cleaned paired FASTQ files for downstream profiling

A2) Host read removal (host depletion)
Use hostfree_bowtie.sh to remove host-derived reads by mapping to the host reference.

Example (pseudo-usage):

bash
复制代码
bash hostfree_bowtie.sh \
  -1 results/01_trimmed/sample.clean.R1.fastq.gz \
  -2 results/01_trimmed/sample.clean.R2.fastq.gz \
  -x /path/to/host_bowtie2_index \
  -o results/02_hostfree/
Expected output

host-depleted paired FASTQ files (recommended input for MetaPhlAn)

mapping logs/statistics (useful for QC reporting)

A3) Taxonomic profiling (MetaPhlAn)
Use metaphlan-rawdata.sh to generate taxonomic profiles.

Example (pseudo-usage):

bash
复制代码
bash metaphlan-rawdata.sh \
  -i results/02_hostfree/ \
  -o results/03_metaphlan/
Expected output

MetaPhlAn taxonomic tables (per-sample)

optional merged abundance tables (depending on your script settings)

Workflow B — Mammary epithelial cell RNA-seq (in vitro)
B1) RNA-seq processing (wrapper)
run_RNAseq.sh typically handles trimming + alignment/counting OR provides an orchestrated workflow.

Example (pseudo-usage):

bash
复制代码
bash run_RNAseq.sh \
  -i rnaseq_fastq/ \
  -o results/rnaseq/
B2) Quantification with kallisto (optional / recommended for speed)
Use run_kallisto.sh for pseudoalignment-based quantification.

Example (pseudo-usage):

bash
复制代码
bash run_kallisto.sh \
  -i rnaseq_fastq/ \
  -x /path/to/kallisto_index \
  -o results/04_kallisto/
Expected output

transcript/gene abundance estimates

logs for QC and troubleshooting

Workflow C — Figures & visualization
Use codes for Figure plot.R to reproduce manuscript-style figures.

Example

bash
复制代码
Rscript "codes for Figure plot.R"
If your plotting script expects specific input filenames, consider adding a small config.R or documenting the expected input table format (see “Reproducibility” below).

Suggested project structure (optional but clean)
If you want a tidy layout:

text
复制代码
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
Reproducibility checklist (highly recommended)
 Record software versions (conda env export, sessionInfo() in R)

 Keep sample metadata (metadata.tsv) with time points, THI/heat-stress labels, etc.

 Save host depletion stats + sequencing depth summaries

 Make merged abundance matrices for downstream stats (alpha/beta diversity, longitudinal models)

Citation
If you use this workflow/code in your research, please cite the associated manuscript (add DOI here once available).
For now, you can cite this repository.

License
This repository is released under the GPL-3.0 license.

Contact
shengtaogao@163.com
Issues and pull requests are welcome.
