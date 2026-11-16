# Single-cell atlas of HPV-negative oral cavity lesions

Re-analysis of a public single-cell RNA-seq dataset of head and neck cancer
([GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919), Choi et al. 2023), restricted to
**HPV-negative oral cavity** samples: normal mucosa (NL), leukoplakia (LP) and oral squamous cell carcinoma (CA).

The long-term goal is to find epithelial genes that separate tumour from normal epithelium **in most
patients**, as candidate markers for checking surgical margins in oral cancer. This repository covers the first
part: quality control, integration, clustering and cell-type annotation, plus checks of the data that
affect the later marker search.

![UMAP after Harmony integration](figures/umap_after_harmony.png)

## Data

| | |
|---|---|
| Source | GEO [GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) |
| Full dataset | 54,239 cells × 20,000 genes, head and neck subsites, HPV+ and HPV− |
| Subset used | Oral cavity (`subsite == "OC"`), HPV-negative: 24,156 cells, 17 samples |
| Samples | 10 carcinoma, 4 leukoplakia, 3 normal mucosa |

## Pipeline

All steps are in [`oscc_scrnaseq_analysis.ipynb`](oscc_scrnaseq_analysis.ipynb) (Python, scanpy).

| Step | Method | Key settings |
|---|---|---|
| Load | Dense text matrix read in chunks into a sparse matrix | — |
| QC | MAD-based outlier detection **per sample** | 5 MADs on log total counts, log genes, top-20-gene %; mitochondrial: 3 MADs (logit scale, upper) **and** > 8% |
| Doublets | Scrublet per sample | flagged, not removed (see *Limitations*) |
| Normalization | `normalize_total` + `log1p`; raw counts kept in `layers["counts"]` | — |
| HVGs | Seurat flavour, selected within each sample | 3,000 genes; mitochondrial and ribosomal genes excluded |
| PCA | on HVGs | 50 PCs, 40 used downstream |
| Integration | Harmony on `sample.id` | used for clustering only |
| Clustering | Leiden on the Harmony neighbour graph | resolutions 0.2 / 0.5 / 1.0; 0.5 chosen |
| Annotation | Manual, from marker dotplots and Wilcoxon marker genes, then compared with the authors' labels | — |

QC kept **21,935 of 24,156 cells** (91%).

## Main observations

**1. Harmony over-integrates tumour cells.** Each tumour has its own copy-number changes, so malignant cells
should partly separate by patient. After Harmony, the share of malignant cells among the 15 nearest
neighbours of normal epithelial cells rose from **0.16 to 0.37** (complete mixing = 0.78). So Harmony (by sample)
is used for cell-type clustering only. Malignancy calls and differential expression will not use the integrated
embedding.

| Before Harmony | After Harmony |
|---|---|
| ![](figures/umap_before_harmony.png) | ![](figures/umap_after_harmony.png) |

**2. Ambient RNA differs between tumour and normal samples.** Malignant cells carry a median of 12 IGKC
UMIs (an antibody gene they do not make), while normal epithelial cells carry 0. The amount follows the
plasma-cell share of each tumour sample. A naive tumour-vs-normal comparison would therefore report
plasma-cell genes as "tumour markers".

<img src="figures/ambient_igkc.png" width="380" alt="Ambient IGKC vs plasma-cell share">

**3. Keratins.** KRT5 (basal keratin) is detected in 98% of malignant cells, while KRT13 (a differentiation
keratin of normal oral mucosa) is mostly absent from them. Loss of a normal marker may be as useful for a
margin test as gain of a tumour marker.

**4. Few normal epithelial cells.** Normal mucosa dissociates poorly, so normal samples give mostly immune and
stromal cells. Epithelial cells per sample: NL 32–85 (3 samples), LP 23–69 (4 samples), CA 10 samples, of which
3 have fewer than 20. A pseudobulk comparison is possible but small, so results will be candidate markers that
need validation.

## Cell-type annotation

17 Leiden clusters (resolution 0.5), annotated from markers before comparing with the authors' labels.

![Marker dotplot](figures/marker_dotplot.png)

| Cluster | Label | Cells |
|---|---|---|
| 0 | Epithelial | 1,410 |
| 1 | Lymphatic endothelial | 371 |
| 2 | Fibroblasts (matrix) | 2,273 |
| 3 | Pericytes | 615 |
| 4 | Fibroblasts (resting) | 4,647 |
| 5 | Myeloid | 2,862 |
| 6 | CD4 T cells | 1,003 |
| 7 | T cells (ribo-high, check QC) | 2,752 |
| 8 | CD8 T / NK | 2,232 |
| 9 | B cells | 575 |
| 10 | Endothelial (blood) | 1,817 |
| 11 | Plasma cells (IgG) | 498 |
| 12 | Mature DC (LAMP3+) | 182 |
| 13 | Mast cells | 125 |
| 14 | Plasmacytoid DC | 368 |
| 15 | Low quality? (check QC) | 111 |
| 16 | Skeletal muscle | 94 |

The epithelial cluster is a broad label. Malignant vs normal epithelium will be decided by CNV inference, not by
keratins or by position in the integrated embedding.

## Limitations

- Scrublet predicted 126 doublets (0.6%). They are flagged in `obs["predicted_doublet"]` but **not removed** in
  this version; they will be removed before differential expression.
- Cluster 15 is unresolved, and cluster 7 needs a final label.
- Only 3 normal mucosa samples, so statistical power for tumour-vs-normal comparisons is limited.

## Next steps

1. CNV inference on epithelial cells to call malignant vs normal epithelium.
2. Pseudobulk differential expression (per sample, raw counts) between tumour and normal epithelium,
   accounting for ambient RNA.
3. Validation of candidate markers in independent data (TCGA, spatial transcriptomics), later by IHC/RNAscope.

## How to run

1. Download `GSE181919_Barcode_metadata.txt.gz` and `GSE181919_UMI_counts.txt.gz` from the
   [GEO page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) into `data/`
   (gzipped or extracted, both work).
2. Create the environment:
   ```bash
   conda create -n oscc-sc python=3.13
   conda activate oscc-sc
   pip install -r requirements.txt
   ```
3. Run the notebook from top to bottom. Reading the count matrix takes about 6 minutes. Checkpoints are written to `results/` and figures to `figures/`.

```
.
├── oscc_scrnaseq_analysis.ipynb   # full analysis
├── figures/                       # figures used in this README
├── requirements.txt
├── data/                          # GEO files (not tracked)
└── results/                       # .h5ad checkpoints (not tracked)
```

## References

- Choi J-H, Lee B-S, Jang JY, et al. Single-cell transcriptome profiling of the stepwise progression of head
  and neck cancer. *Nat Commun* (2023). https://doi.org/10.1038/s41467-023-36691-x
- Heumos L, Schaar AC, Lance C, et al. Best practices for single-cell analysis across modalities.
  *Nat Rev Genet* 24, 550–572 (2023). https://doi.org/10.1038/s41576-023-00586-w
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis.
  *Genome Biol* 19, 15 (2018).
- Korsunsky I, et al. Fast, sensitive and accurate integration of single-cell data with Harmony.
  *Nat Methods* 16, 1289–1296 (2019).
- Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of cell doublets in single-cell
  transcriptomic data. *Cell Syst* 8, 281–291 (2019).
