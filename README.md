# Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping

**Supplementary data, code and materials for the book chapter (Academic Press, Elsevier)**

This repository provides the public data, analysis code, supplementary tables and figures, and the
Supplementary Material document for the chapter. The chapter maps the urban microbiome and its
antimicrobial-resistance (AMR) landscape in public shotgun metagenomes:

* 8,053 included urban environmental metagenomes from 112 datasets, 106 countries and 274 source-labelled cities;
* three analytical streams kept separate by pipeline (Indian open-drain wastewater with CARD/RGI; the global ResFinder resource; and uniform reprocessing of 51 runs centred on Dhaka and South Asian comparators).

Developed by Mohammad Saddam Hosen, MD Shahidul Islam Fakir and Kurshed Khan.

---

## Contents at a glance

| What | Where |
|---|---|
| Supplementary Material (Figs. S1–S30, Tables S1–S96, with interpretations) | [`supplementary_materials/`](supplementary_materials) |
| Tables S1–S96, unabridged (CSV + one formatted Excel workbook) | [`tables/`](tables) |
| Figures S1–S30 (PNG 300 dpi, TIFF 600 dpi, PDF) | [`figures/`](figures) |
| Cleaned datasets (12) and data dictionaries | [`data/`](data) |
| Independent recomputation of headline results (75/75 reproduced) | [`analysis/`](analysis) |
| Workflow code, archived analysis code, executed notebooks | [`scripts/`](scripts), [`notebooks/`](notebooks) |
| Cleaning, statistical-analysis, reproducibility and QC reports | [`documentation/`](documentation) |
| Original data-collection package (31 Oct 2026) | repository root (see section 8) |

## 1. Project objectives

1. Describe how urban microbial composition and diversity vary across environments, cities and seasons.
2. Characterise ARGs, resistance classes and mechanisms, their abundance and diversity, and their relation to community structure.
3. Examine environmental factors associated with AMR patterns and what the evidence implies for public health.
4. Assess whether the data support predictive modelling and spatial pattern identification.
5. Summarise contextual evidence on environmental AMR and food systems.

## 2. Datasets (`data/`)

| File | Content | Rows |
|---|---|---|
| `DS_A_sample_master.csv` | Harmonized metadata, one row per repository sample (9,734 screened; 8,053 included) | 9,734 |
| `DS_B_run_master.csv` | Sequencing runs: platform, reads, bases, FASTQ URLs | 11,866 |
| `DS_C_amr_genes_long.csv.gz` | Producer-processed AMR gene detections (CARD/RGI for DS02; ResFinder/KMA resource) | 316,983 |
| `DS_D_taxonomy_long.csv` | DS02 Kraken 2/Bracken phylum profiles | 75,309 |
| `DS_E_diversity_long.csv` | Diversity indices of the ResFinder resource | 3,011 |
| `DS_F_env_covariates_long.csv` | Source-reported design/environmental covariates | 19,441 |
| `DS_G` … `DS_K` | Uniform reprocessing (Stream 3): run QC, ResFinder hits, Kraken 2 clades, ARG summaries, comparison variables | 45–25,446 |
| `DS_L_dhaka_geographic_verification.csv` | Dhaka coordinate provenance and point-in-polygon checks | 17 |

The data dictionary (`data/data_dictionary/Data_Dictionary_All_Datasets.xlsx`) lists each variable's
original name, description, type, unit, possible values, missingness and any transformation.

## 3. Data preparation

Cleaning was **non-destructive**:
* No record was removed.
* Variable names were standardised to snake_case.
* Numeric variables were typed, and the documented missing-value codes were kept in `<variable>_missing_code` columns. Missing AMR observations are never zero.
* Date precision was classified.
* Coordinate, range, duplicate, outlier (Tukey 3 × IQR on log10) and referential-integrity checks were recorded in a `cleaning_flags` column instead of being "fixed".

Every operation is logged in `data/Cleaning_Operations_Log.csv`; the methods are in
`documentation/data_cleaning_report/`.

## 4. Analytical methods

* Descriptive statistics.
* Kruskal–Wallis tests with epsilon-squared, including site-level tests on the 19 open-drain sites with whole-site permutation.
* Linear mixed models with site or study random intercepts.
* PERMANOVA/PERMDISP, Mantel and Procrustes tests.
* CLR-PCA, PCoA and NMDS.
* Spearman correlation with Benjamini–Hochberg adjustment.
* Wilson intervals.
* DerSimonian–Laird pooling.
* Moran's I and distance decay.
* Exact permutation tests with Cliff's delta.
* Exploratory random-forest, logistic and ridge models with study-grouped, temporal and leave-one-city-out validation.

Equations, variables, the register of 383 tests and the limitations are in
`documentation/statistical_analysis_report/`. Streams 1–3 use different abundance units (RPKM, FPKM,
FPKM-16S) and must not be compared in magnitude. Gene detection shows DNA presence, not phenotypic resistance.

## 5. Supplementary tables and figures

* **Tables S1–S96** (`tables/csv/Table_Sxx_<ID>.csv`, `tables/Supplementary_Data_Tables_S1-S96.xlsx`) are in eight sections: evidence base; Indian open drains (taxonomy, diversity, resistome); the global ResFinder resource; Dhaka and South Asian comparisons; predictive modelling; spatial analysis; robustness and depth sensitivity. `tables/Table_Index_S1-S96.csv` gives each table's title, size and source.
* **Figures S1–S30** (`figures/`) are listed with captions, source tables, resolution and the plotting notebook in `figures/Figure_Index_S1-S30.csv`.
* **Supplementary Material** (`supplementary_materials/Supplementary_Material_Final.docx`; the PDF version will be added) uses automatic Word captions, hyperlinked cross-references, a table of contents and lists of figures and tables. `appendices/Supplementary_Text_Corrections_Log.xlsx` documents the corrections made to the paraphrased text against the verified interpretations.

## 6. Software requirements and installation

Python ≥ 3.11 is required. Install the dependencies with:

```bash
pip install -r requirements.txt
```

Microsoft Word is needed only to regenerate the document fields and PDF (step 07, Windows). The archived
analysis notebooks also need matplotlib, seaborn, scikit-learn and scikit-bio. The raw-read
reprocessing used fastp 1.4.0, Kraken 2.17.1 with SILVA 138, and KMA 1.7.1 with ResFinder DB 4.7.2 on Linux.

## 7. Reproducing the results

Recompute the headline supplementary results from the cleaned data in this repository:

```bash
python scripts/pipeline/run_all.py --verify
```

This writes `analysis/results/Independent_Verification_of_Supplementary_Results.xlsx`, where every recomputed
value is compared with the published table value. The full packaging workflow (steps 00–10; see
`documentation/reproducibility/Technical_Report.docx`) also needs the original research folder.
`scripts/original_analysis/` and `notebooks/` archive the verified analysis (twenty executed notebooks, run on
31 October 2026). These archived scripts contain the absolute paths of the original environment, which must be
adapted before re-execution.

## 8. Original data-collection package (repository root)

| File | Content |
|---|---|
| `00_Master_Data_Index.xlsx` | Entry point: headline counts, dataset list, coverage by sample type, platform, AMR and taxonomy availability, missing metadata, analysis readiness |
| `01_Public_Dataset_Inventory.xlsx` | Source registry, 236-query search log, 1,050 deduplicated candidate records, accession register, eligibility, download manifest, exclusion log, source-verification checks |
| `02_Sample_Metadata_Master.xlsx` | 9,734 harmonized samples, 11,866 runs, geography, environment, linkage, missingness, duplicates, QC |
| `03_Analysis_Ready_Metagenomic_Data.xlsx` | 8,053 included samples; long-format taxonomy, AMR genes, diversity and covariates |
| `04_AMR_Reference_Database_Register.xlsx` | FARME DB and CARD 4.0.2 reference details, other AMR resources, versions, annotation methods |
| `Public_Data_Source_Log.csv` | Every executed query and download, with exact URL |

These files were collected on 31 October 2026 from MGnify, MG-RAST, ENA, NCBI BioProject/BioSample,
Europe PMC and Zenodo. Harmonization rules:
* The unit is the repository sample accession.
* No coordinates, dates, cities or AMR values were imputed.
* Missing-value codes are Not reported, Not available, Not yet verified, Not applicable and Restricted access.
* CARD/RGI and ResFinder/KMA results are kept apart.

Some sample types and coordinates were corrected after an independent audit; the cleaned
`data/cleaned/DS_A_sample_master.csv` includes these corrections (`sample_type_correction`,
`coordinate_provenance_note`).

## 9. Data availability statement

All sequence data analysed are publicly available in INSDC (ENA/SRA/DDBJ) or MG-RAST under the accessions
listed in `data/cleaned/DS_A_sample_master.csv` and `DS_B_run_master.csv`. Processed tables were obtained from:
* Zenodo 10.5281/zenodo.17590646 (Singh et al. 2026; CC BY 4.0);
* Zenodo 10.5281/zenodo.6919377 (Martiny et al. global resistome resource; CC BY 4.0);
* the Supplementary Data of Munk et al. 2022 (*Nature Communications*, doi:10.1038/s41467-022-34312-7).

Dhaka administrative boundaries come from the BBS/OCHA COD-AB dataset (Humanitarian Data Exchange). This
repository contains no patient-level or personal data. Third-party publications and raw downloads are not
redistributed and can be retrieved from the sources above.

## 10. Citation, licence and contact

If you use these materials, please cite the chapter *Global Insights of the Urban Microbiome and AMR Landscape
from Metagenomic Mapping* (Academic Press, Elsevier) and this repository, together with the original data
sources listed above.

**Licence:** this repository does not currently include a licence file, so the default copyright applies. Data
derived from CC BY 4.0 sources remain subject to those licences.

**Contact:** Mohammad Saddam Hosen (repository owner, via GitHub).

Parts of the data preparation and packaging were carried out with AI assistance (Claude) for author review;
the authors are responsible for the verification and interpretation of the results.
