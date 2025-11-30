# Distinct Gastrointestinal and Reproductive Microbial Patterns in Female Holobiont of Infertility

[![DOI](https://img.shields.io/badge/DOI-10.3390%2Fmicroorganisms12050989-blue)](https://doi.org/10.3390/microorganisms12050989)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![QIIME2](https://img.shields.io/badge/QIIME2-2022.8-green)](https://qiime2.org/)
[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://www.r-project.org/)

Analysis code for the study of gastrointestinal (GIT) and female reproductive tract (FRT) microbiome patterns in women with infertility conditions, including endometriosis.

## Overview

This repository contains the complete bioinformatics workflow for 16S rRNA gene sequencing analysis from Ion Torrent platform, including:

- **QIIME2 pipeline** for sequence processing and taxonomic classification
- **R analysis** for diversity metrics, statistical comparisons, and visualization

The study investigates the "female holobiont" concept—the complex microbial balance across multiple body sites (oral, fecal, vaginal, endometrial) and its relationship with infertility conditions.

## Study Design

| Parameter | Description |
|-----------|-------------|
| Cohort | 21 women with infertility (8 endometriosis, 13 other conditions) |
| Body sites | Oral, feces, vagina, endometrial fluid, endometrium |
| Total samples | 105 (5 sites × 21 patients) |
| Sequencing | Ion Torrent (Ion 16S Metagenomics Kit) |
| Target region | 16S rRNA gene (V2-V9 regions) |

The Ion 16S Metagenomics Kit was selected for its broad coverage of multiple hypervariable regions (V2, V3, V4, V6-V7, V8, V9), which improves species-level identification across diverse body sites and is validated for human microbiome studies in reproductive and gastrointestinal samples.

Taxonomic classification was performed using Greengenes 13_8 (99% OTUs), a widely validated reference database with consistent nomenclature that enables direct comparison with published microbiome studies.

## Key Findings

- **Distinct GIT vs FRT microbiomes**: FRT dominated by *Lactobacillus* (67% average), while GIT showed higher diversity with *Streptococcus* and *Ruminococcus*
- **Significant diversity differences**: GIT richness (119.8) significantly higher than FRT (44.7), p = 6.8 × 10⁻⁸
- **Endometriosis-specific patterns**: Lower bacterial diversity (p = 0.027) and elevated *Lactobacillus* dominance (>82% in vaginal and endometrial sites)
- **GIT biomarker**: *Haemophilus* identified as a taxon associated with GIT in endometriosis cases

## Repository Structure

```
16S-female-microbiome-infertility/
├── README.md
├── LICENSE
├── scripts/
│   ├── 01_qiime2_iontorrent_pipeline.sh    # QIIME2 processing pipeline
│   └── 02_statistical_analysis.Rmd          # R analysis and visualization
└── data/
    ├── README_data.md                       # Data access instructions
    └── metadata.xlsx                        # Clinical metadata (anonymized)
```

## Methods

### Sequence Processing (QIIME2)

| Step | Tool/Method | Parameters |
|------|-------------|------------|
| Import | SingleEndFastqManifestPhred33V2 | Demultiplexed reads |
| Denoising | DADA2 denoise-pyro | trim-left=15, trunc-q=20 |
| Taxonomy | VSEARCH + Greengenes 13_8 | 99% OTUs |
| Phylogeny | MAFFT + FastTree | Midpoint rooted |
| Rarefaction | core-metrics-phylogenetic | depth=3400 |

### Statistical Analysis (R)

- **Alpha diversity**: Richness, Shannon index, Faith's PD
- **Beta diversity**: Bray-Curtis, Weighted/Unweighted UniFrac
- **Statistics**: Wilcoxon test (BH correction), PERMANOVA (Bonferroni)
- **Nestedness analysis**: Jaccard-based turnover and nestedness components

## Requirements

### QIIME2 Environment

```bash
conda activate qiime2-2022.8
```

### R Packages

```r
library(phyloseq)
library(vegan)
library(qiimer)
library(betapart)
library(tidyverse)
library(ggplot2)
library(ggpubr)
```

## Quick Start

```bash
# Run QIIME2 pipeline
bash scripts/01_qiime2_iontorrent_pipeline.sh
```

```r
# Run R analysis
rmarkdown::render("scripts/02_statistical_analysis.Rmd")
```

## Data Access

- **Raw sequences**: NCBI SRA [PRJNA1031834](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1031834)
- **Original code**: Zenodo [10.5281/zenodo.10362432](https://doi.org/10.5281/zenodo.10362432)

## Citation

```bibtex
@article{Marcos2024,
  title     = {Distinct Gastrointestinal and Reproductive Microbial Patterns
               in Female Holobiont of Infertility},
  author    = {Marcos, Ana T. and Rus, Maria J. and Areal-Quecuty, Victoria
               and Simon-Soro, Aurea and Navarro-Pando, Jos{\'e} Manuel},
  journal   = {Microorganisms},
  volume    = {12},
  number    = {5},
  pages     = {989},
  year      = {2024},
  publisher = {MDPI},
  doi       = {10.3390/microorganisms12050989}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Maria J. Rus**
[![ORCID](https://img.shields.io/badge/ORCID-0000--0003--3659--2821-green)](https://orcid.org/0000-0003-3659-2821)
Email: marjimrus@gmail.com

## Acknowledgments

- INEBIR (Instituto para el Estudio de la Biología de la Reproducción Humana)
- Universidad de Sevilla - Program for Emerging Research (Ref. 2022/00000333)
