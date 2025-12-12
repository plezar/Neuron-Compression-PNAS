## Figure 5 scripts

### `01_qc`
- Performs QC of the GBM-CARE scRNA-seq data
- Follows the procedures described in Nomura et al, 2025 (*Nature Genetics*)
- Most of the code is borrowed from `dravishays/GBM-CARE-WT`.

### `02_create_seurat.r`
- Creates sample-level Seurat objects
- Most of the code is borrowed from `dravishays/GBM-CARE-WT`.

### `03_gbm_care_figure_5.r`
- Figure 5L
- Figure 5M
- Figure S6E

### `04_figure_s5`

#### `01_tiff_nd2_analysis.py`
- Processes multi-channel TIFF and ND2 imaging data using fixed DAPI/GFAP thresholds to segment nuclei, classify astrocytes vs neurons, quantify FosB intensity across z-stacks for control and experimental replicates, and save per-slice metrics and summary statistics to disk.

#### `02_analyze.py`
- Loads slice-level FosB quantification data, extracts the middle 20 z-slices from each ND2 series, computes per-series averaged FosB ratios, performs group-wise statistics, and generates summary plots and CSV outputs.

### `05_figure_5.Rmd`
- Figure 5C
- Figure 5D
- Figure 5E
- Figure 5G
- Figure 5H
- Figure 5J
- Figure 5K
- Figure S6C
- Figure S6D