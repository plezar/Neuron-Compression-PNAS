# Mechanical Compression Induces Neuronal Apoptosis, Reduces Synaptic Activity, and Promotes Glial Neuroinflammation in Mice and Humans

This repository contains all code and accompanying files required to reproduce the analyses and figures from:

**Zarodniuk, M., Wenninger, A., et al.**  
*Mechanical compression induces neuronal apoptosis, reduces synaptic activity, and promotes glial neuroinflammation in mice and humans.* **PNAS (2025).**

---

## Overview

This repository provides reproducible R- and MATLAB-based pipelines used to generate all figures presented in the manuscript. Each figure is accompanied by its own `.Rmd` file (and additional scripts where applicable), enabling full regeneration of results. All analyses use publicly available or soon-to-be–public datasets. Pre-rendered HTML files are also included for convenient browsing without re-running code.

---

## Repository Structure

Only top-level figure directories are described below. Each contains its own `README` file, analysis scripts, data subdirectories, and rendered outputs.

### **`Figure_1/`**
Contains scripts and data supporting quantification of neuronal viability, nuclear morphology, and glial reactivity following mechanical compression. Includes raw CSV measurements and rendered output figures.

### **`Figure_2/`**
Contains MATLAB and Python scripts for quantitative calcium imaging analysis, including fluorescence extraction, spike detection, and spectral characterization of neuronal activity.

### **`Figure_3/`**
Includes all workflows for bulk RNA-seq differential expression, GSEA/ORA, HOMER motif analysis, sample preprocessing, and supporting datasets for both mouse and human compressed neuron models.

### **`Figure_4/`**
Contains analyses of immunocytochemistry (ICC), immunohistochemistry (IHC), morphology quantification, qPCR datasets, and associated R Markdown files supporting molecular and phenotypic responses to compression.

### **`Figure_5/`**
Includes preprocessing, QC, Seurat workflows, multi-dataset integration, BayesPrism outputs, and all scripts used to analyze and visualize single-cell and spatially informed human GBM datasets relevant to mechanical stress signatures.

### **Other files**
- **`Patzke-Datta.Rproj`** — RStudio project file.
- **`structure.txt`** — Overview of directory layout and dependencies.

---

## Data Availability

**GEO accession numbers for the two primary datasets used in this study will be added once available.**  
A placeholder for the final accession numbers will be updated upon public release.

---

## License

This project is released under the **MIT License**.  
See the `LICENSE` file for details.