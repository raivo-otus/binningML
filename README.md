# Quantile Binning for Microbiome ML
 
EuroBioC2026 poster project — Rasmus Hindström, Juho Pelto, Tuomas Borman & Leo Lahti  
Department of Computing, University of Turku
 
## Overview
 
Benchmarks quantile binning against standard microbiome transformations (rCLR, TSS, presence/absence) for ML classification across curated case/control datasets. Implemented as a reproducible `{targets}` pipeline.
 
## Requirements
 
- R 4.4+ / Bioconductor 3.20+
- Nix + direnv (environment managed via `default.nix` / `.envrc`)
 
```bash
direnv allow   # activates the Nix environment
```
 
## Usage
 
```r
source("_targets.R")        # Source the pipeline
targets::tar_make()         # run the pipeline
targets::tar_visnetwork()   # inspect the dependency graph
```
 
The Quarto report in `analysis/report.qmd` reads results from the targets store via `tar_read()`.
 
## Structure
 
```
R/              # Pipeline functions (data, transformations, modeling)
analysis/       # Quarto report and saved plots
_targets.R      # Pipeline definition
```
