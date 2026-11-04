# Supplementary Material Package
## Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping
*Book chapter, Academic Press (Elsevier). Package built on 4 November 2026.*

This folder holds the reproducible supplementary-material package for the chapter. It contains the cleaned
public datasets with data dictionaries, an independent check of the headline results, Tables S1–S96 and
Figures S1–S30 in publication formats, the final Supplementary Material (DOCX and PDF) with automatic captions
and cross-references, and the supporting reports and code.

> **Scope.** The chapter's statistical analysis was completed and independently audited on 31 October 2026
> (folders `Updated_Research_Analysis`, `Independent_Audit`). This package does **not** re-run or alter that
> analysis. It cleans and documents the data and recomputes the headline results from the cleaned data
> (75 of 75 checks reproduced). It also corrects the text of the supplement and packages everything for release.

---

## 1. Chapter context and objectives

| # | Objective | Research question |
|---|---|---|
| 1 | Urban microbial diversity | How do microbial composition and diversity vary across urban environments, cities and seasons? |
| 2 | AMR patterns and microbiome associations | Which ARGs, classes and mechanisms occur, how abundant are they, and how do they relate to community characteristics? |
| 3 | Environmental and public-health implications | Which environmental factors are associated with AMR patterns? |
| 4 | Predictive and spatial assessment | Do the data support predictive modelling and spatial pattern identification? |
| 5 | Food systems and policy | What does the evidence say about environmental AMR and food systems? |

The evidence is analysed in **three streams that are never pooled**, because they use different pipelines,
databases and abundance units:

* **Stream 1:** Indian open-drain wastewater (DS02). 447 samples from 19 sites in 4 cities. CARD/RGI annotation; unit RPKM.
* **Stream 2:** the curated global ResFinder resource (Martiny et al.). Urban matrices plus food-related contextual units. KMA/ResFinder annotation; unit FPKM.
* **Stream 3:** 51 public runs reprocessed uniformly (16 from Dhaka, plus South Asian comparators and calibration runs). Each run is subsampled to 4 M read pairs and processed with fastp, Kraken 2/SILVA 138 and KMA/ResFinder 4.7.2; unit FPKM-16S.

## 2. Folder structure

```text
00_Project_Inventory/        Source_File_Inventory.xlsx (3,036 files, SHA-256, integrity, duplicates, roles),
                             Data_Quality_Report.docx, Quality_Control_Checklist.xlsx,
                             Public_Release_Screen.csv, Processing_Log.txt
01_Raw_Data_Archive/         read-only copies of the original source data + SHA-256 manifest (not published)
02_Cleaned_Datasets/         Cleaned_CSV (12 datasets), Cleaned_Excel, Data_Dictionaries, Cleaning_Operations_Log.csv
03_Data_Analysis/            Descriptive_Statistics, Statistical_Analysis, Comparative_Analysis,
                             Analytical_Results (independent verification of supplementary results)
04_Tables/                   CSV_Tables (Table_S01..S96), Excel_Tables (formatted workbook), Table_Index_S1-S96.csv
05_Figures/                  PNG (300 dpi), TIFF (600 dpi LZW), PDF, Editable_Figures (plotting code), Figure_Index
06_Research_Documents/       Supplementary_Material_Final.docx / .pdf
                             Supplementary_Appendices/ (text-corrections log, build report)
07_Code_and_Scripts/         Python (workflow steps 00-09, run_all.py), Jupyter_Notebooks (20 executed analysis
                             notebooks), Original_Analysis_Code (archival), Requirements.txt
08_Methodology_and_Reports/  Data_Cleaning_Methodology.docx, Statistical_Analysis_Methodology.docx, Technical_Report.docx
```

## 3. Datasets

| Dataset | Content | Rows |
|---|---|---|
| DS_A_sample_master | Harmonized metadata for all 9,734 screened samples (8,053 included urban shotgun metagenomes) | 9,734 |
| DS_B_run_master | Sequencing runs (platform, reads, bases, URLs) | 11,866 |
| DS_C_amr_genes_long (.csv.gz) | Producer-processed AMR gene detections (CARD/RGI DS02; ResFinder/KMA resource) | 316,983 |
| DS_D_taxonomy_long | DS02 Kraken 2/Bracken phylum profiles | 75,309 |
| DS_E_diversity_long | Diversity indices from the ResFinder resource | 3,011 |
| DS_F_env_covariates_long | Source-reported design/environmental covariates | 19,441 |
| DS_G … DS_K | Stream 3 reprocessing: run QC, ResFinder hits, Kraken 2 clades, ARG summaries, comparison variables | 45–25,446 |
| DS_L_dhaka_geographic_verification | Dhaka coordinate provenance and point-in-polygon checks | 17 |

**Cleaning (non-destructive):**
* No record was removed; row counts are asserted equal before and after cleaning.
* Variable names are in snake_case.
* Numeric columns are typed, and missing-value codes are kept in `<variable>_missing_code` columns. Missing AMR values are never set to zero.
* Dates are classified by precision.
* Range, coordinate, duplicate, outlier (3 × IQR on log10) and referential-integrity checks are recorded as `cleaning_flags`.

The full methods are in `08_Methodology_and_Reports/Data_Cleaning_Methodology.docx`, and every variable is
described in `02_Cleaned_Datasets/Data_Dictionaries`.

## 4. Analytical methods (summary)

The methods are:
* descriptive statistics;
* Kruskal–Wallis tests with epsilon-squared, including site-level tests on 19 site medians with whole-site permutation;
* linear mixed models with site or study random intercepts;
* PERMANOVA/PERMDISP, Mantel and Procrustes tests;
* CLR-PCA, PCoA and NMDS ordination;
* Spearman correlations with Benjamini–Hochberg adjustment;
* Wilson confidence intervals;
* DerSimonian–Laird random-effects pooling;
* Moran's I and distance decay;
* exact permutation tests with Cliff's delta;
* exploratory random-forest, logistic and ridge models with study-grouped, temporal and leave-one-city-out validation.

`Statistical_Analysis_Methodology.docx` gives the equations, variables, a register of 383 tests and the limitations.

## 5. Supplementary tables and figures

* **Tables S1–S96** are in eight sections: evidence base; Stream 1 taxonomy and diversity; Stream 1 resistome; Stream 2; Stream 3 (Dhaka and South Asia); predictive modelling; spatial analysis; robustness and depth sensitivity. Values are copied programmatically from the verified analysis CSVs and checked cell by cell (`04_Tables/Table_Index_S1-S96.csv`).
* **Figures S1–S30** are the verified analysis figures: PNG at 300 dpi, TIFF at 600 dpi, and PDF (`05_Figures/Figure_Index_S1-S30.csv`). Fig. S28 is the title-corrected version embedded in the supplement.

## 6. Supplementary_Material_Final.docx — what changed

`Supplementary Material.docx` (author's version of 4 Nov 2026) was used as the base and was not modified.
The final copy:

1. **Corrects 254 text errors** that had entered the captions, notes and interpretations during paraphrasing. Each one was checked against the verified, author-approved interpretations. They comprise 57 factual errors (for example, Table S8 gave the denominator as 81,697 instead of 8,053; Fig. S16 attributed the any-hit VIM count to NDM; Fig. S22 said "seven" confirmed samples instead of 12), 114 meaning errors, 71 garbled passages and 12 sentence-style captions. Every change is listed with its reason in `Supplementary_Appendices/Supplementary_Text_Corrections_Log.xlsx`.
2. **Uses automatic Word captions** (`Fig. S{SEQ Figure}`, `Table S{SEQ Table}`) in Word's built-in Caption style.
3. **Turns every citation of a supplementary item into a cross-reference.** All 246 citations are hyperlinked REF fields, and every table and figure is cited at least once.
4. **Adds front matter and layout:** a Table of Contents, List of Figures, List of Tables, numbered section headings, a running header and "Page X of Y" footers. Microsoft Word updated all 741 fields (138 pages, 0 field errors).

**PDF status:** Word's PDF output stalled under automation on the build computer, even for a one-line test document, so `Supplementary_Material_Final.pdf` has not been generated yet. To create it, open the DOCX in Word and choose File > Save As > PDF, or run `07_update_fields_and_export_pdf.ps1` on a machine where Word's PDF export works.

## 7. Reproducing the workflow

```bash
pip install -r 07_Code_and_Scripts/Requirements.txt
cd 07_Code_and_Scripts/Python
python run_all.py            # full workflow (needs the original research folder as AMR_SOURCE_DIR)
python run_all.py --verify   # recompute the headline results and run QC from the package alone
```

`Technical_Report.docx` describes the execution order, environment and validation. The original analysis
notebooks are archived and contain the absolute paths of the original environment.

## 8. Data availability and attribution

All sequence data are public in INSDC (ENA/SRA/DDBJ) or MG-RAST and are listed by accession in
`DS_A_sample_master` and `DS_B_run_master`. Processed tables come from:
* Singh et al. 2026, Zenodo 10.5281/zenodo.17590646 (CC BY 4.0);
* Martiny et al. global resistome resource, Zenodo 10.5281/zenodo.6919377 (CC BY 4.0);
* Munk et al. 2022, *Nature Communications*, doi:10.1038/s41467-022-34312-7.

Dhaka administrative boundaries are from the BBS/OCHA COD-AB dataset (Humanitarian Data Exchange; see its
licence terms). Please cite these sources
when reusing the corresponding data. No patient-level or personal data are included.

## 9. Citation and contact

Please cite the chapter: *Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping*,
Academic Press (Elsevier), together with this repository. Contact: the chapter authors (see the repository README).
The analysis workflow was developed by Mohammad Saddam Hosen, MD Shahidul Islam Fakir and Kurshed Khan.
Parts of the analysis and packaging were prepared with AI assistance (Claude) for author review. The authors
remain responsible for verifying and interpreting the results.
