# Organ-specific and conserved regulatory logic orchestrates gene expression in the embryonic mesothelium

[![GitHub release](https://img.shields.io/github/v/release//loganminhdang/Mesothelium_paper_2025)](https://github.com//loganminhdang/Mesothelium_paper_2025/releases)

> One-sentence description of what this repository does.

This repository contains the analysis script for the manuscript  
"Organ-specific and conserved regulatory logic orchestrates gene expression in the embryonic mesothelium" (Dang et al, 2025)

---

## 📂 Repository Map
```text
PROJECT-NAME/


├── R_scripts/
│   ├── ATAC-seq TF-IDF normalisation and visualisation.Rmd               # Perform normalisation and visualisation of ATAC-seq tracks
│   └── Bulk ATAC-seq TMM normalisation.Rmd           # Perform normalisation of epicardium and epicardium-derived cell ATAC-seq data
|   └── Differential chromatin accessibility analysis.Rmd           # Differential chromatin accessibility analysis of epicardium and epicardium-derived cell ATAC-seq data
|   └── hdWGCNA analysis of mesothelia.Rmd           # Detection of organ mesothelial gene co-expression modules
├── notebooks/
│   ├── 01-exploratory.ipynb
│   └── 02-model-validation.ipynb
├── scripts/
│   ├── 00-setup.sh        # One-command install (macOS + Linux).
│   ├── 01-preprocess.sh
│   └── 02-run_models.R
├── Rmd/
│   ├── paper-figures.Rmd
│   └── supplementary.Rmd
├── renv.lock              # R package snapshot (use with renv::restore()).
├── requirements.txt       # Python packages for the notebooks.
├── Makefile               # Optional: make all reproduces the whole study.
└── README.md              # You are here.

```
---

## 📓 Interactive Notebooks
| Notebook                              | Purpose                  | Run in Browser                                                                                                                                         |
| ------------------------------------- | ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `notebooks/01-exploratory.ipynb`      | EDA & sanity checks      | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/YOUR-USERNAME/PROJECT-NAME/HEAD?filepath=notebooks%2F01-exploratory.ipynb) |
| `notebooks/02-model-validation.ipynb` | Cross-validation metrics | same link as above                                                                                                                                     |

## 📊 Outputs

- **Figures**: saved to `outputs/figures/` (vector PDF + 300 dpi PNG).  
- **Supplementary Tables**: auto-written to `outputs/tables/` as both `.csv` and LaTeX `.tex`.  
- **Rendered R Markdown**:  
  - `Rmd/paper-figures.Rmd` → `docs/paper-figures.html` (manuscript plots).  
  - `Rmd/supplementary.Rmd` → `docs/supplementary.html`.
