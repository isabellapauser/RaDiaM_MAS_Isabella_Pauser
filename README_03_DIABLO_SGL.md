# Melanoma DIABLO Multi-Omics Analysis

This repository contains the R Markdown workflow `03_DIABLO_SGL.Rmd`, which integrates transcriptomics, proteomics and surface-enhanced Raman spectroscopy (SERS) data from the melanoma cell lines MCM-1G and MCM-DLN.

The workflow uses Sparse Group LASSO (SGL) to select informative SERS regions and DIABLO, implemented in `mixOmics`, to identify a correlated multi-omics signature that distinguishes the non-metastatic MCM-1G and metastatic MCM-DLN cell lines.

## Study design

The analysis includes six matched biological samples:

- MCM-1G: `1G_A`, `1G_B`, `1G_C`
- MCM-DLN: `DLN_A`, `DLN_B`, `DLN_C`

Only 2D transcriptomics samples are included. Proteomic measurement replicates are supplied in an already merged MaxQuant `proteinGroups` table. Individual SERS spectra are aggregated to the same six biological samples after SGL-based region selection.

## Required directory structure

Place the R Markdown file in the `DIABLO_Melanoma_Pauser` project directory with the following structure:

```text
DIABLO_Melanoma_Pauser/
├── 03_DIABLO_SGL.Rmd
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

The two unusual filename spellings (`.txt.txt` and `melamonma`) must match the actual files or be corrected in the `paths` chunk.

## Software requirements

The workflow requires R and the following packages:

```r
mixOmics
DESeq2
edgeR
dplyr
igraph
org.Hs.eg.db
AnnotationDbi
BiocParallel
msgl
sglOptim
parallel
doParallel
foreach
knitr
```

Install Bioconductor packages through `BiocManager` and CRAN packages through `install.packages()`. Record the package versions used for the final analysis with `sessionInfo()`.

## Analysis workflow

1. Read and harmonize sample names across the three data layers.
2. Merge the six featureCounts tables and filter lowly expressed genes.
3. Apply the DESeq2 variance-stabilizing transformation to RNA-seq counts.
4. Import MaxQuant LFQ intensities, remove proteins containing missing values and log2-transform the retained intensities without imputation.
5. Divide the SERS spectrum into contiguous 10 cm⁻¹ windows.
6. Perform replicate-aware, three-fold SGL cross-validation, holding out the corresponding MCM-1G and MCM-DLN biological replicate together.
7. Aggregate selected SERS wavenumbers into spectral regions and then into six biological samples.
8. Map Ensembl gene identifiers to gene symbols where mappings are available.
9. Retain up to 2,000 transcriptomic features, 1,000 proteins and all eligible SGL-selected SERS regions based on pre-scaling variance.
10. Scale the selected features and perform pairwise sPLS comparisons.
11. Fit a fully connected DIABLO model with off-diagonal design values of 1 and diagonal values of 0.
12. Evaluate the number of components, tune `keepX`, fit the final sparse model and export its selected signature.
13. Generate loading, sample, correlation-network, circos and clustered-image-map visualizations.

The random seed is set to `123` for reproducibility. SGL and DIABLO validation use biological-replicate-aware or stratified three-fold schemes, respectively.

## Current model configuration

The current run produced the following input dimensions:

| Data block | Samples | Candidate features |
|---|---:|---:|
| Transcriptome | 6 | 2,000 |
| Proteome | 6 | 1,000 |
| SERS | 6 | 9 regions |

Performance evaluation indicated that one component was sufficient. Feature-number tuning selected five transcriptomic features, five proteins and two SERS regions for the final component. A separate two-component model is fitted only for exploratory two-dimensional visualization.

## Main outputs

The workflow writes the following files:

- `sample_info_from_names.tsv`: sample groups and biological replicates
- `DIABLO_RNA_gene_annotation.tsv`: original RNA identifiers and mapped labels
- `DIABLO_selected_features.tsv`: final DIABLO signature and loading weights
- `SERS/SGL_selected_SERS_wavenumbers.tsv`: individual SGL-selected wavenumbers
- `SERS/SGL_selected_SERS_regions.tsv`: selected 10 cm⁻¹ spectral regions
- `SERS/SGL_tied_best_models.tsv`: equally performing SGL models
- `SERS/SGL_cross_validation_error_path.tsv`: SGL error across the fitted path
- `SERS/SGL_cross_validation.rds`: SGL cross-validation object
- `SERS/SGL_final_model.rds`: final fitted SGL path
- `DIABLO_figures/`: exported DIABLO figures in PNG format

## Running the analysis

Open `03_DIABLO_SGL.Rmd` in RStudio and render the document to HTML or PDF. Confirm that the working directory resolves to `DIABLO_Melanoma_Pauser` and that all files in the `paths` and `rna-load` chunks exist.

Before running the full document, resolve these issues in the current R Markdown file:

- Remove the duplicated `selected-features` chunk.
- Remove the duplicated `visualization-model`, `plot-diablo`, `plot-individuals`, `plot-arrows` and `plot-variables` chunks.
- Change `final.diablo.model$design` to `final_diablo_model$design`.
- Ensure every chunk label is unique, because duplicated labels prevent knitting.
- Run the network export with `sane_par=FALSE` if the global plotting settings cause `figure margins too large`.

After these checks, run the chunks sequentially because later steps depend on objects created earlier in the document.

## Interpretation and limitations

This analysis is exploratory. Only six biological samples are available, and supervised SERS selection is performed before downstream DIABLO validation. Consequently, classification performance may be optimistic and does not constitute external validation. The selected features should be treated as candidate multi-omics markers requiring confirmation in independent biological samples.

Unmapped Ensembl identifiers are retained rather than discarded. Positive and negative loading signs indicate opposite directions along a latent component but should not automatically be interpreted as MCM-1G- or MCM-DLN-associated without checking the orientation of the corresponding sample scores.
