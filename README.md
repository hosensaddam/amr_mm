# Public Metagenomic Data Collection — Data Preparation README

**Project:** Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping
**Collection date:** 31 October 2026 (all API queries and downloads)
**Prepared by:** Agentic custom AI assistant for author review. System developed by MSH, MSIF, and KK. Eligibility decisions and classifications require author confirmation.

## 1\. What this package contains

|File|Content|
|-|-|
|`00\\\\\\\_Master\\\\\\\_Data\\\\\\\_Index.xlsx`|Entry point: headline counts, dataset list, coverage by sample type, platform, AMR and taxonomy availability, missing metadata, analysis readiness|
|`01\\\\\\\_Public\\\\\\\_Dataset\\\\\\\_Inventory.xlsx`|Source registry, 236-query search log, 1,050 deduplicated candidate records, accession register, eligibility, download manifest, exclusion log, 22 source-verification checks|
|`02\\\\\\\_Sample\\\\\\\_Metadata\\\\\\\_Master.xlsx`|9,734 harmonized samples (one row per repository sample), 11,866 runs, geography, environment, linkage, missingness (live formulas), duplicates, QC|
|`03\\\\\\\_Analysis\\\\\\\_Ready\\\\\\\_Metagenomic\\\\\\\_Data.xlsx`|8,053 included samples; long-format taxonomy, AMR genes, diversity and covariates; within-pipeline summaries; readiness assessment|
|`04\\\\\\\_AMR\\\\\\\_Reference\\\\\\\_Database\\\\\\\_Register.xlsx`|FARME DB and CARD 4.0.2 reference details, other AMR resources, versions, annotation methods|
|`Public\\\\\\\_Data\\\\\\\_Source\\\\\\\_Log.csv`|Every executed query and every download, with exact URL|
|`Data\\\\\\\_Collection\\\\\\\_and\\\\\\\_QC\\\\\\\_Report.docx`|Narrative report of what was searched, verified, downloaded, excluded and what remains|

The working folder (`AMR\\\\\\\_Chapter\\\\\\\_Public\\\\\\\_Data\\\\\\\_Collection/`) also holds `scripts/`, `raw\\\\\\\_downloads/` (original files, unmodified), `harmonized/` (CSV versions of every table) and `logs/`.

## 2\. Pipeline (re-run in order)

|Script|Step|
|-|-|
|`01\\\\\\\_query\\\\\\\_search.py`|35 queries × 6 databases (MGnify v2, MG-RAST, ENA, NCBI BioProject, Europe PMC, Zenodo); saves every raw response|
|`01b\\\\\\\_mgnify\\\\\\\_adapted\\\\\\\_queries.py`|26 adapted MGnify queries (v2 `search` is a phrase match; all multi-word queries returned 0)|
|`01c\\\\\\\_recompute\\\\\\\_relevant\\\\\\\_ids.py`|Recomputes screened-in record IDs from saved responses (corrects a list-filter bug; counts unchanged)|
|`02\\\\\\\_study\\\\\\\_level\\\\\\\_verification.py`|Deduplicates candidates; retrieves ENA run reports for 593 INSDC studies|
|`02b\\\\\\\_review\\\\\\\_table.py` + `curation\\\\\\\_decisions.py`|Review table and recorded eligibility decisions|
|`03\\\\\\\_fetch\\\\\\\_sample\\\\\\\_metadata.py`|Extended run/sample metadata for 132 selected studies|
|`03b\\\\\\\_fetch\\\\\\\_biosample\\\\\\\_attributes.py`|NCBI BioSample attributes where ENA lacked location or date (2,139 samples)|
|`03c\\\\\\\_subset\\\\\\\_global\\\\\\\_resistome.py`|Extracts Zenodo 6919377 ARG and diversity rows for the 1,004 covered runs|
|`04\\\\\\\_harmonize.py`|Harmonized sample/run masters and long tables|
|`05\\\\\\\_build\\\\\\\_workbooks.py`|Workbooks, source log and summaries|

Python 3 with pandas, openpyxl, xlsxwriter and pycountry. Run with `python3 -I scripts/<name>.py` from the package root.

## 3\. Harmonization rules

* **Unit of the master table:** repository sample accession (BioSample/ENA sample or MG-RAST sample). Runs are aggregated (`Run\\\\\\\_ID` lists them). For DS01, `Biological\\\\\\\_Sample\\\\\\\_ID` groups accessions that derive from one physical sample (951 accessions listed in Munk et al. 2022 Supplementary Data 1 → 798 biological samples; the other 97 DS01 accessions are PRJEB13831 samples not used in that paper).
* **Geography and date priority:** DS01 Munk Suppl. Data 1 by run, else ENA; DS02 Zenodo `meta.tsv` (city) + SRA; DS03 NCBI BioSample; DS04 MetaSUB GitHub metadata joined on `sample\\\\\\\_alias = uuid` + SRA; all others ENA `geo\\\\\\\_loc\\\\\\\_name`, then NCBI BioSample.
* **City vs locality:** `City` is filled only where the source labels a city. Other `geo\\\\\\\_loc\\\\\\\_name` text is kept in `Locality\\\\\\\_As\\\\\\\_Reported` (625 distinct strings) and is not converted to cities.
* **No imputation:** no coordinates, dates, cities or AMR values were inferred. One documented derivation: MetaSUB country propagated from same-city SRA records when the SRA field was missing and unambiguous (395 samples, flagged).
* **Missing-value codes:** Not reported / Not available / Not yet verified / Not applicable / Restricted access. Missing AMR observations are never zero.
* **Data types:** only samples with a shotgun metagenomic run are `Included`; where such a sample also has non-shotgun runs (81 included samples), those runs are listed in `Other\\\\\\\_Run\\\\\\\_IDs` and left out of `Run\\\\\\\_ID`, URLs and size totals. Amplicon, metatranscriptomic, isolate, single-cell and synthetic samples are `Excluded`; food-chain studies are `Separate contextual dataset`.
* **Sample types:** study-level category from the review, refined per sample by documented keyword rules (`Sample\\\\\\\_Type\\\\\\\_Rule`). Source matrix text is kept in `Sample\\\\\\\_Type\\\\\\\_Source\\\\\\\_Value`.
* **AMR pipelines are kept apart:** CARD/RGI (DS02) and ResFinder/KMA (Zenodo 6919377) differ in database, naming and count definition. Counts are raw, not normalized.

## 4\. Key verified findings about the sources

1. All 951 runs in Munk et al. (2022) Supplementary Data 1 map to **seven** ENA projects, not the six named in the paper: the 2016 pilot round (227 runs) is in PRJEB27054, which also underlies MGnify MGYS00005846.
2. DS01 holds 242 distinct city labels (243 city–country pairs, because Santiago occurs in Chile and Spain) in 101 countries, consistent with the paper's 243 cities in 101 countries. Spelling variants such as Kochi/Cochin are not merged.
3. MetaSUB's SRA project holds 1,898 samples, not the 4,728 published, all listed as NovaSeq 6000 (publication: HiSeq X Ten); 86 records list Kyiv under Poland (10 more Kyiv records lack a country; all 96 set to Not yet verified); 60 carry a 1970-01-01 placeholder date.
4. Indian samples come from large open wastewater drains (publication Methods); processed archive MD5 matches Zenodo; the CARD release used is not stated anywhere accessed.
5. Proposal target cities are present: Dhaka (16 included samples across 3 datasets), Delhi (202, incl. 100 wastewater), Karachi (8 sewage samples, DS01 only). MetaSUB contains no South Asian city.

## 5\. What was downloaded and what was not

* **Downloaded (15 files, \~870 MB):** Zenodo 17590646 processed archive; Zenodo 6919377 metadata, ResFinder annotation, diversity and ARG tables; Munk et al. supplementary package and full texts; MetaSUB metadata; FARME and CARD reference files. Plus metadata-only API responses (ENA, NCBI BioSample, MG-RAST).
* **Not downloaded:** raw FASTQ for any dataset (≈24.1 TB compressed across included samples; e.g., DS01 6.2 TB, DS02 1.2 TB, DS04 1.4 TB); Zenodo 6919377 `rRNA.tsv` (6.4 GB); MGnify analyses; Munk et al. Source Data and Singh et al. supplementary files (publisher returned HTTP 403 or not in Europe PMC package).

## 6\. Open items needing author decisions

1. Confirm or revise study-level eligibility decisions (`curation\\\\\\\_decisions.py`; 25 studies pending).
2. Approve any raw-read subset for uniform reprocessing; reprocessing is required for cross-study comparison.
3. Confirm FARME DB is the proposal's "FARME data management" resource.
4. Retrieve manually: Munk et al. Source Data; Singh et al. Supplementary Data S7 (site coordinates) and source data; CARD release used by Singh et al.
5. Decide whether to standardize the 625 locality strings to cities (manual lookup table).
6. Confirm the Dhaka PRJNA1399457 publication link (none found in repository records).

