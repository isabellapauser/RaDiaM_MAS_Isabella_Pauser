# Melanoma DIABLO Multi-Omics Analysis

This repository contains an R Markdown workflow for integrating transcriptomics, proteomics, and surface-enhanced Raman spectroscopy (SERS) data from the melanoma cell lines MCM-1G and MCM-DLN. The analysis uses DIABLO from the `mixOmics` package to identify a correlated multi-omics signature that distinguishes the non-metastatic MCM-1G and metastatic MCM-DLN cell models.

The workflow retains six matched biological samples:

- `1G_A`, `1G_B`, and `1G_C`
- `DLN_A`, `DLN_B`, and `DLN_C`

## Main script

`03_DIABLO_SGL_GA_WS.Rmd`

The Rmd can produce HTML or PDF output and supports three alternative methods for reducing the SERS feature space:

| Method | Description |
| --- | --- |
| `SGL` | Sparse group lasso applied to contiguous fixed-width Raman windows. It can select both complete groups and individual wavenumbers within groups. |
| `GA` | Genetic algorithm that selects between 2 and 10 fixed-width Raman windows using replicate-aware cross-validation and a penalty for larger solutions. |
| `WS` | Unsupervised variable-width window construction. Adjacent wavenumbers are joined when their Pearson correlation reaches the selected threshold. |

All methods summarize selected or constructed spectral regions by their mean intensity before the spectra are aggregated to the six biological samples and supplied to DIABLO.

## Project structure

The script expects the following directory structure:

```text
DIABLO_Melanoma_Pauser/
├── 03_DIABLO_SGL_GA_WS.Rmd
├── Transcriptome/
│   ├── 1G-1_STAR_WGF_sorted.rmdup.featurecount
│   ├── 1G-2_STAR_WGF_sorted.rmdup.featurecount
│   ├── 1G-3_STAR_WGF_sorted.rmdup.featurecount
│   ├── DLN-2D-1_STAR_WGF_sorted.rmdup.featurecount
│   ├── DLN-2D-2_STAR_WGF_sorted.rmdup.featurecount
│   └── DLN-2D-3_STAR_WGF_sorted.rmdup.featurecount
├── Proteome/
│   └── proteinGroups_results_run1_mqpar_2025__FH_250515_Mayr_extracts_mcm_new_merged.txt.txt
└── SERS/
    └── melamonma_g_dln_sers_rerun_annamasterthesis_preprocessed.csv
```

The SERS file must contain the metadata columns `label` and `file` together with numeric spectral variables. Raman shifts are extracted from column names containing formats such as `SERS_X802.5`, `SERS_802.5`, or `X802.5`.

## Requirements

The workflow uses R and the following packages:

```r
# CRAN packages
install.packages(c(
  "dplyr",
  "GA",
  "igraph",
  "knitr",
  "msgl",
  "doParallel"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "AnnotationDbi",
  "BiocParallel",
  "DESeq2",
  "edgeR",
  "mixOmics",
  "org.Hs.eg.db"
))
```

The Rmd uses `set.seed(123)` for reproducibility.

## Selecting the SERS method

Choose the method in the YAML header of the Rmd:

```yaml
params:
  sers_selection_method: "SGL"
  ws_correlation_threshold: 0.99
```

Valid values are `SGL`, `GA`, and `WS`.

### SGL

```yaml
sers_selection_method: "SGL"
```

The spectrum is divided into fixed 10 cm⁻¹ windows. Sparse group lasso is evaluated using three replicate-aware folds that hold out replicate A, B, or C across both cell lines. Among models with the minimum cross-validation error, the workflow selects the model containing the fewest SERS variables.

### GA

```yaml
sers_selection_method: "GA"
```

The genetic algorithm searches among the same fixed 10 cm⁻¹ windows. Its fitness combines mean balanced accuracy across the replicate-aware folds with a penalty for the number of selected regions. The current settings are:

```text
Selected regions: 2–10
Population size: 100
Maximum iterations: 300
Early stopping: 50 generations without improvement
Crossover probability: 0.8
Mutation probability: 0.1
Region penalty: 0.002
```

### WS

```yaml
sers_selection_method: "WS"
ws_correlation_threshold: 0.99
```

WS calculates Pearson correlations between consecutive wavenumbers across all spectra. Adjacent variables are assigned to the same window when their correlation is greater than or equal to the threshold. WS does not use the MCM-1G/MCM-DLN labels during window construction.

A lower threshold produces fewer and broader windows, whereas a higher threshold produces more and narrower windows. If the script reports that fewer than two windows were produced, increase `ws_correlation_threshold`, for example from `0.99` to `0.999`.

## Analysis workflow

1. Read and harmonize the six RNA-seq, proteomics, and SERS samples.
2. Remove proteins with missing values; no proteomics imputation is performed.
3. Construct or select SERS regions using SGL, GA, or WS.
4. Average SERS measurements within each selected region and biological sample.
5. Filter RNA features with `edgeR::filterByExpr` and apply variance-stabilizing transformation with DESeq2.
6. Map Ensembl identifiers to gene symbols where possible.
7. Retain up to 2,000 RNA, 1,000 proteomics, and 1,000 SERS features by variance.
8. Scale each omics block and verify matching sample order.
9. Fit a fully connected DIABLO model with two starting components.
10. Evaluate the number of components and tune `keepX` using repeated three-fold cross-validation.
11. Fit the final sparse DIABLO model and export its selected features and visualizations.

## Running the analysis

Open `03_DIABLO_SGL_GA_WS.Rmd` in RStudio, select the desired parameters, and click **Knit**. Alternatively, render it from R:

```r
rmarkdown::render(
  "03_DIABLO_SGL_GA_WS.Rmd",
  params = list(
    sers_selection_method = "GA",
    ws_correlation_threshold = 0.99
  )
)
```

To compare the three approaches, render the Rmd separately with `SGL`, `GA`, and `WS`. Each run uses a separate output directory and therefore does not overwrite the other methods.

## Output

Results are named according to the selected SERS method:

```text
DIABLO_<METHOD>_results/
├── <METHOD>_selected_SERS_wavenumbers.tsv
├── <METHOD>_selected_SERS_regions.tsv
├── <METHOD>_sample_info_from_names.tsv
├── <METHOD>_DIABLO_RNA_gene_annotation.tsv
├── <METHOD>_DIABLO_selected_features.tsv
├── <METHOD>_DIABLO_final_model.rds
└── <METHOD>_DIABLO_figures/
```

Method-specific files include:

- SGL: cross-validation error path, tied best models, and fitted SGL model.
- GA: selection summary and fitted GA model.
- WS: adjacent-wavenumber correlations and window summary.

The figure directory contains DIABLO sample plots, arrow plots, variable plots, loading plots, a circos plot, a network plot, and a clustered image map.

## Interpretation and limitations

The analysis contains only six independent biological samples: three MCM-1G and three MCM-DLN replicates. The repeated spectra do not represent independent biological replicates. Replicate-aware validation is therefore used during SGL and GA selection to reduce information leakage.

Nevertheless, model tuning and performance estimates remain unstable at this sample size. Perfect or near-perfect classification should be treated as exploratory and should not be interpreted as evidence of clinical validity. Selected proteins, transcripts, and SERS regions are candidate features requiring validation in independent biological material.

