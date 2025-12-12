## Figure 3 scripts

### `01_deseq2.R`
- Performs differential expression analysis (DESeq2) for human induced neurons (iN) and mouse glia (mG) under control vs compression conditions, including per–cell-line and combined analyses.
- Implements sequencing depth correction via downsampling and exports DE gene tables and processed DESeq2 objects for downstream GSEA and visualization.

### `02_gsea.R`
- Performs GO and Hallmark GSEA for each iN and mGlia cell line using DESeq2 differential expression results, generating both readable tables and merged summary files.
- Converts gene identifiers, filters for overlapping genes where needed, runs `clusterProfiler` GSEA, and saves merged outputs for downstream visualization and figure generation.

### `03_figure3.Rmd`
- Figure S1C
- Figure S1D
- Figure 3C
- Figure S3D
- Figure 3D
- Figure 3E
- Figure 3A

### `04_supplement.R`
- Figure S3B
- Figure S3C

### `util.R`
- Miscellaneous utility functions

### `preprocessing`
- A folder containing scripts for RNA-seq data preprocessing
- Splits human and mouse reads using `bbmap`
- Aligns and counts the data