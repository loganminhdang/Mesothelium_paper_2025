# Inference of regulatory logic in the embryonic mesothelium

> Repository description

This repository contains the analysis script for the manuscript  
"Organ-specific and conserved regulatory logic orchestrates gene expression in the embryonic mesothelium" (Dang et al, 2025)

---

## 📂 Repository Map
```text
PROJECT-NAME/


├── R_scripts/
│   ├── ATAC-seq TF-IDF normalisation and visualisation.Rmd               # Perform normalisation and visualisation of ATAC-seq tracks
│   └── Bulk ATAC-seq TMM normalisation.Rmd                               # Perform normalisation of epicardium and epicardium-derived cell ATAC-seq data
|   └── Differential chromatin accessibility analysis.Rmd                 # Differential chromatin accessibility analysis of epicardium and epicardium-derived cell ATAC-seq data
|   └── hdWGCNA analysis of mesothelia.Rmd                                # Detection of organ mesothelial gene co-expression modules
├── notebooks/
│   ├── Analysis of E13.5 heart scRNA-seq.ipynb                            # E13.5 mouse heart scRNA-seq analysis
│   └── Analysis of lung and pancreas scRNA-seq datasets in early development.ipynb # Lung and pancreas scRNA-seq analysis
│   └── scRNA-seq pathway enrichment.ipynb                                 # scRNA-seq pathway enrichment analysis
│   └── Metacell differential expression analysis.ipynb                    # scRNA-seq differential analysis
│   └── Epicardium EMT reconstruction via scRNA-seq.ipynb                  # Single-cell reconstruction of mouse epicardial EMT
│   └── Analysis of human foetal heart scRNA-seq datasets.ipynb            # Human heart scRNA-seq analysis
│   └── Benchmarking of human fetal heart scRNA-seq integration.ipynb      # Single-cell integration of the developing human heart
│   └── Integration and benchmarking - mouse embryonic heart scRNA-seq     # Single-cell integration of the developing mouse heart
│   └── Integration of mouse and human epicardium.ipynb                    # Single-cell integration of the epicardium across species
│   └── scATACseq analysis of embryonic mesothelia.ipynb                   # Embryonic mesothelial scATAC-seq analysis
│   └── scATAC-seq peak calling and topic modelling.ipynb                  # Embryonic mesothelial scATAC-seq peak-calling and topic modelling (for GRN inference)
├── unix_scripts/
│   ├── Code for ATAC-seq and Cut&Run-seq signal visualisation.txt       
│   ├── Code for MACS2 peak filtering.txt             # ATAC-seq peak processing
│   └── Code for digital footprinting.txt             # ATAC-seq footprinting
│   └── Code for E-P prediction.txt                   # Prediction of Enhancer-promote connections
└── README.md              # You are here.

```
---
## 📊 Outputs

- **Figures**: `figures/` (700 dpi PNG).  
- **Single-cell integration benchmarking results**: `integrations/` (PNG).  
