# Transcriptomic Analysis Scripts — Resveratrol–SIRT1–SNCA Axis in Parkinson's Disease

R scripts for the integrative transcriptomic analysis of Parkinson's disease across 33 GEO datasets (postmortem, iPSC, and cell-based models), supporting the manuscript submitted to *Molecular Omics*.

## Repository contents

| Script | Purpose |
|---|---|
| `RNA-seq_script_21_07_2026.R` | DESeq2 pipeline: per-dataset differential expression for RNA-seq count data |
| `microarray_script_21_07_2026.R` | GEOquery + limma pipeline: per-dataset differential expression for microarray data |
| `Transcriptomic_figures_scripts.R` | Combines DE results from all 33 datasets into the manuscript's summary figures |

## Script → figure/table mapping

| Output | Script | Section |
|---|---|---|
| Per-dataset PCA, heatmaps, volcano plots, marker gene panels (RNA-seq datasets) | `RNA-seq_script_21_07_2026.R` | `#PCA`, `#heatmap`, `#volcanoplot`, `#marker tidyplot` |
| Per-dataset PCA, heatmaps, volcano plots, marker gene panels (microarray datasets) | `microarray_script_21_07_2026.R` | corresponding sections |
| Per-dataset GO/KEGG enrichment tables and dot plots | both scripts | `#GO KEGG` |
| Per-dataset STRING PPI networks and centrality tables | both scripts | `#string` |
| Per-dataset SIRT1–SNCA Spearman correlation | both scripts | `#correlation` |
| **Figure 1** — resveratrol target gene × dataset dot plot | `Transcriptomic_figures_scripts.R` | `##FIGURE 1 SCRIPT` |
| **Figure 2** — SIRT1–SNCA pattern classification across all 33 datasets | `Transcriptomic_figures_scripts.R` | `##FIGURE 2 SCRIPT` |
| **Figure 3** — drug–gene–dataset chord diagram | `Transcriptomic_figures_scripts.R` | `#FIGURE 3 SCRIPT` |

## Datasets

The analysis integrates 33 GEO series spanning postmortem, iPSC, and cell-based Parkinson's disease models. All datasets are publicly available from the [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/); accession numbers are listed in the manuscript and its supplementary tables.

## Run order

The three scripts are not standalone — they must be run in this order:

1. **Per-dataset DE analysis.** For each of the 33 datasets, run either `RNA-seq_script_21_07_2026.R` (for RNA-seq count data) or `microarray_script_21_07_2026.R` (for microarray intensity data), depending on the platform. Update the input paths and the `keep_samples` / group definitions for each dataset before running. Each run produces:
   - a `deseq_<GSE_ID>` object/table (log2FC, padj, gene symbol)
   - a `df_vsd_<GSE_ID>` object/table (normalized expression matrix)
   - per-dataset figures (PCA, heatmap, volcano, marker panel, GO/KEGG, STRING network) and their supporting `.csv`/`.xlsx` files
2. **Drug-gene interaction filtering** (upstream step, not included in this repository) produces the `All_Approved_Drug_Interactions_Combined_*.xlsx` files and the `all_resveratrol` table consumed by Figure 1.
3. **Combined figures.** Once all 33 `deseq_*` / `df_vsd_*` objects and the drug-interaction tables exist in the environment, run `Transcriptomic_figures_scripts.R` in order (Figure 1 → Figure 2 → Figure 3) to produce the manuscript's summary figures.

## Requirements

R ≥ 4.2 (developed under R 4.x on Windows) with the following packages:

**Bioconductor:** `DESeq2`, `limma`, `GEOquery`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `DOSE`, `EnhancedVolcano`, `apeglm`, `glmpca`, `genefilter`, `BSgenome.Hsapiens.UCSC.hg38`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, `GenomicFeatures`, `tximport`, `msigdbr`, `PoiClaClu`

**CRAN:** `dplyr`, `tidyr`, `stringr`, `ggplot2`, `tibble`, `readxl`, `purrr`, `pheatmap`, `RColorBrewer`, `circlize`, `patchwork`, `scales`, `httr`, `jsonlite`, `igraph`, `ggraph`, `tidygraph`, `openxlsx`, `ggbeeswarm`, `gridExtra`, `reshape2`

STRING and DGIdb queries require an internet connection (public REST APIs; STRING calls use `required_score` thresholds of 400–700 depending on script).

For exact reproducibility, run `sessionInfo()` after installing the above and include its output alongside the analysis, since specific package versions are not pinned in these scripts.

## Input data (not included)

- Raw count matrices / GEO series matrices for each of the 33 datasets (downloadable from GEO using the accessions listed above; `GEOquery::getGEO()` is used directly for microarray series)
- Sample metadata files (`coldata_*.csv`) with group assignments per dataset
- Consolidated drug-gene interaction tables (`All_Approved_Drug_Interactions_Combined_*.xlsx`), produced by an upstream DGIdb-filtering step not included here

## Usage notes

- File paths (`setwd()`, input/output paths) are hardcoded to the original working environment and must be updated before running.
- Each per-dataset script is organized into named, collapsible code blocks (`{ }` sections such as `#PCA`, `#heatmap`, `#volcanoplot`, `#GO KEGG`, `#string`, `#correlation`) that can be run independently once the upstream objects exist.
- Outputs are written to a local `figures/` folder as `.tiff`/`.png` (600 dpi for publication figures) alongside `.csv`/`.xlsx` supplementary tables.

## License

No license file is currently included. Add a `LICENSE` file (e.g. MIT) if you intend this code to be reused or adapted by others.

## Citation

These scripts support the manuscript on the resveratrol–SIRT1–SNCA axis in Parkinson's disease, submitted to *Molecular Omics*. Please cite the associated publication if you use or adapt this code.
