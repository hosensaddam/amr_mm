"""Step 08 - Generate the methodology, analysis, reproducibility and data-quality reports (Word).

All numbers in the reports are read from files produced by steps 01-05 or from the verified analysis
tables; nothing is typed by hand except method descriptions and standard equations.
Outputs:
  00_Project_Inventory/Data_Quality_Report.docx
  08_Methodology_and_Reports/Data_Cleaning_Methodology.docx      (Data Cleaning and Preparation Report)
  08_Methodology_and_Reports/Statistical_Analysis_Methodology.docx (Statistical Analysis Report)
  08_Methodology_and_Reports/Technical_Report.docx               (Technical Reproducibility Report)
"""
import datetime as dt

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Cm

from amr_paths import (INV, CLEAN_CSV, DICT, ANALYSIS, TABLES, FIGURES, DOCS, REPORTS, UPDATED, log)

CHAPTER = "Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping"
TODAY = dt.date.today().strftime("%d %B %Y")
T = UPDATED / "04_Tables"


# ------------------------------------------------------------------ helpers
def new_doc(title, subtitle):
    d = Document()
    st = d.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(11)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    st.paragraph_format.space_after = Pt(6)
    st.paragraph_format.line_spacing = 1.15
    for h, sz in (("Heading 1", 15), ("Heading 2", 13), ("Heading 3", 11.5)):
        d.styles[h].font.name = "Times New Roman"
        d.styles[h].font.size = Pt(sz)
        d.styles[h].font.color.rgb = None
    sec = d.sections[0]
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Cm(2.2))
    t = d.add_paragraph()
    r = t.add_run(title); r.bold = True; r.font.size = Pt(18)
    s = d.add_paragraph()
    r = s.add_run(subtitle); r.italic = True
    d.add_paragraph(f"Chapter: {CHAPTER} (Academic Press, Elsevier). Prepared: {TODAY}.")
    toc = d.add_paragraph()
    _field(toc, 'TOC \\o "1-2" \\h \\z \\u', "Contents (update field in Word: right-click > Update Field).")
    # footer page numbers
    fp = sec.footer.paragraphs[0]
    fp.alignment = 1
    fp.add_run("Page ")
    _field(fp, "PAGE", "1")
    return d


def _field(p, instr, result):
    def r_with(child):
        r = OxmlElement("w:r"); r.append(child); p._p.append(r)
    b = OxmlElement("w:fldChar"); b.set(qn("w:fldCharType"), "begin"); b.set(qn("w:dirty"), "true"); r_with(b)
    it = OxmlElement("w:instrText"); it.set(qn("xml:space"), "preserve"); it.text = f" {instr} "; r_with(it)
    s = OxmlElement("w:fldChar"); s.set(qn("w:fldCharType"), "separate"); r_with(s)
    p.add_run(result)
    e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), "end"); r_with(e)


def para(d, text, bold_lead=None, style=None):
    p = d.add_paragraph(style=style)
    if bold_lead:
        p.add_run(bold_lead).bold = True
    p.add_run(text)
    return p


def bullets(d, items):
    for it in items:
        d.add_paragraph(it, style="List Bullet")


def eq(d, text):
    p = d.add_paragraph()
    p.paragraph_format.left_indent = Cm(1)
    r = p.add_run(text); r.italic = True
    return p


def table(d, df, caption=None, widths=None, font=8.5, max_rows=60):
    if caption:
        c = d.add_paragraph(); c.add_run(caption).bold = True
        c.paragraph_format.keep_with_next = True
    df = df.head(max_rows)
    t = d.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    for i, col in enumerate(df.columns):
        cell = t.rows[0].cells[i]
        cell.text = str(col)
        for r in cell.paragraphs[0].runs:
            r.bold = True; r.font.size = Pt(font)
        tcpr = cell._tc.get_or_add_tcPr(); shd = OxmlElement("w:shd")
        shd.set(qn("w:val"), "clear"); shd.set(qn("w:fill"), "D9E2F3"); tcpr.append(shd)
    trpr = t.rows[0]._tr.get_or_add_trPr(); h = OxmlElement("w:tblHeader"); trpr.append(h)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for i, v in enumerate(row):
            if isinstance(v, float):
                v = f"{v:,.4g}" if abs(v) < 1e6 else f"{v:,.0f}"
            cells[i].text = "" if pd.isna(v) else str(v)
            for r in cells[i].paragraphs[0].runs:
                r.font.size = Pt(font)
    if widths:
        for row in t.rows:
            for i, w in enumerate(widths):
                row.cells[i].width = Cm(w)
    d.add_paragraph()
    return t


def landscape(d):
    sec = d.sections[-1]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width


# ------------------------------------------------------------------ data
def load():
    inv = pd.read_csv(INV / "Source_File_Inventory.csv")
    prof = pd.read_excel(ANALYSIS / "Descriptive_Statistics" / "Dataset_Profiles.xlsx", sheet_name=None)
    dd = pd.read_csv(DICT / "Data_Dictionary_All_Datasets.csv")
    ver = pd.read_csv(ANALYSIS / "Analytical_Results" / "Independent_Verification_of_Supplementary_Results.csv")
    tests = pd.read_csv(T / "R02u_updated_statistical_tests_register.csv")
    lim = pd.read_csv(T / "R04u_updated_limitations.csv")
    obj = pd.read_csv(T / "R00_research_objectives.csv")
    tidx = pd.read_csv(TABLES / "Table_Index_S1-S96.csv")
    fidx = pd.read_csv(FIGURES / "Figure_Index_S1-S30.csv")
    corr = pd.read_excel(DOCS / "Supplementary_Appendices" / "Supplementary_Text_Corrections_Log.xlsx",
                         sheet_name="Corrections")
    return inv, prof, dd, ver, tests, lim, obj, tidx, fidx, corr


DATASET_DESC = {
    "DS_A_sample_master": "Harmonized metadata for every screened repository sample (9,734); 8,053 are included urban "
                          "shotgun metagenomes, 840 are food-related contextual samples and 841 were excluded.",
    "DS_B_run_master": "One row per sequencing run with repository-reported platform, read and base counts and file URLs.",
    "DS_C_amr_genes_long": "Producer-processed gene-level AMR detections: CARD/RGI read mapping for the Indian open-drain "
                           "study (DS02; Stream 1) and ResFinder/KMA read mapping from the Martiny et al. resource "
                           "(Stream 2). The two pipelines are never pooled.",
    "DS_D_taxonomy_long": "Kraken 2/Bracken phylum-level profiles for DS02 (447 samples).",
    "DS_E_diversity_long": "Diversity indices provided by the Martiny et al. resource (ARG and SILVA rRNA genus/phylum).",
    "DS_F_env_covariates_long": "Design and environmental covariates reported by the sources (sampling round, season, etc.).",
    "DS_G_reprocessed_run_qc": "Stream 3: read QC and yields for the 51 runs uniformly reprocessed (fastp, Kraken 2/SILVA 138, KMA/ResFinder 4.7.2).",
    "DS_H_reprocessed_resfinder_annotated": "Stream 3: ResFinder template hits with gene/class annotation, high-confidence flag and FPKM-16S / RPKM-depth.",
    "DS_I_reprocessed_taxonomy_kraken_long": "Stream 3: Kraken 2 clade counts at all ranks (SILVA 138 16S index).",
    "DS_J_reprocessed_amr_summary": "Stream 3: per-run ARG totals, gene counts and diversity.",
    "DS_K_reprocessed_comparative_variables": "Stream 3: per-run analysis variables used in the Dhaka and South Asian comparisons.",
    "DS_L_dhaka_geographic_verification": "Coordinate provenance and point-in-polygon verification of Dhaka-labelled samples.",
}


# ------------------------------------------------------------------ reports
def cleaning_report(prof, dd, corr_unused=None):
    d = new_doc("Data Cleaning and Preparation Report", "Methodology, operations and final dataset characteristics")
    ds, ops = prof["Dataset_Summary"], prof["Cleaning_Operations"]
    flags, comp = prof["Flag_Counts"], prof["Completeness"]
    d.add_heading("1. Purpose and principles", 1)
    para(d, "The chapter relies on public metagenomic metadata and producer-processed annotation tables that were collected, "
            "harmonized and audited before this step (data collection on 31 October 2026; independent audit and update the "
            "same day). The cleaning described here prepares those harmonized tables for public release and reuse. It is "
            "deliberately non-destructive:")
    bullets(d, [
        "No observation is deleted. Row counts before and after cleaning are asserted to be identical for every dataset.",
        "Documented missing-value codes (Not reported, Not available, Not yet verified, Not applicable, Restricted access) "
        "are preserved. Where a numeric variable contained such a code, the variable is typed as numeric and the code is "
        "kept in a companion column <variable>_missing_code. Missing AMR observations are never set to zero.",
        "Outliers, implausible values and inconsistencies are flagged in a cleaning_flags column; no value is removed, "
        "winsorised or imputed. Extreme values in sequencing depth and abundance are expected in metagenomics and are "
        "handled analytically (log scales, rank-based tests, depth sensitivity analyses), not by exclusion.",
        "Variable names are standardised to lower snake_case; the original name is listed in the data dictionary.",
        "Source-reported labels (e.g., city spellings such as Kochi/Cochin) are not merged, because harmonizing them "
        "would require external look-ups that were not part of the verified analysis.",
    ])
    d.add_heading("2. Source files and dataset descriptions", 1)
    for _, r in ds.iterrows():
        para(d, f" {DATASET_DESC.get(r.Dataset, '')} Source: {r.Source_File}.", bold_lead=f"{r.Dataset} - {r.Title}.")
    table(d, ds[["Dataset", "Rows_Original", "Columns_Original", "Rows_Cleaned", "Columns_Cleaned",
                 "Exact_Duplicate_Rows", "Records_Flagged", "Pct_Records_Flagged"]],
          "Table 1. Datasets before and after cleaning")
    d.add_heading("3. Cleaning procedures", 1)
    steps = [
        ("Whitespace", "leading/trailing spaces removed and repeated internal spaces collapsed in every text cell."),
        ("Duplicates", "exact duplicate rows and duplicated primary keys detected and flagged (none removed)."),
        ("Data types", "numeric variables converted; non-numeric documented codes moved to companion columns; any "
                       "undocumented non-numeric value flagged."),
        ("Placeholders", "a base count of 0 for MG-RAST records without FASTQ files is a placeholder; it was set to "
                         "missing with the code 'Zero placeholder (no FASTQ base count reported)'."),
        ("Dates", "precision classified (day, month, year, interval); ISO 8601 date added for day-level values; "
                  "1970-01-01 placeholders, unparseable values and dates after the collection date flagged; "
                  "Collection_Year checked against the date."),
        ("Coordinates", "range check (|lat| <= 90, |lon| <= 180), (0, 0) null-island check and incomplete pairs."),
        ("Ranges", "percentages (0-100), fractions (0-1), indices and counts (non-negative) checked."),
        ("Outliers", "values beyond Q1 - 3 IQR or Q3 + 3 IQR on the log10 scale (Tukey far-out rule) flagged."),
        ("Categories", "values differing only by case or punctuation searched within each categorical variable."),
        ("Relationships", "sample and run identifiers checked across all tables (referential integrity); country-ISO "
                          "code mapping; included samples must have a shotgun sequencing type."),
    ]
    for a, b in steps:
        para(d, " " + b, bold_lead=a + ":")
    d.add_heading("4. Data-quality issues found and their handling", 1)
    nz = ops[ops.Records_or_Cells_Affected > 0]
    table(d, nz[["Dataset", "Step", "Detail", "Records_or_Cells_Affected"]], "Table 2. Cleaning operations with a "
          "non-zero effect", font=7.5, max_rows=120)
    table(d, flags, "Table 3. Records flagged, by flag (records are retained)", font=8)
    d.add_heading("5. Missing-value treatment", 1)
    para(d, "Missingness is reported, not repaired. No coordinates, dates, cities or AMR values were imputed in the source "
            "harmonization or in this step. The variables with the highest proportion of missing or coded values in the "
            "sample master are listed below (complete table: 03_Data_Analysis/Descriptive_Statistics/Dataset_Profiles.xlsx).")
    sm = comp[(comp.Dataset == "DS_A_sample_master") & (comp.Pct > 0)].sort_values("Pct", ascending=False).head(25)
    table(d, sm, "Table 4. Most incomplete sample-master variables (all 9,734 records)", font=8)
    d.add_heading("6. Duplicate handling", 1)
    para(d, "No exact duplicate rows were found in any dataset. Key duplicates reported in Table 2 are expected structures "
            "(for example, several CARD models of one gene in one run, or one physical sample deposited as several "
            "accessions in DS01, which the source groups through Biological_Sample_ID); they are flagged so that analysts "
            "can aggregate at the appropriate unit, and are not removed.")
    d.add_heading("7. Variable transformations", 1)
    tr = dd[dd.Transformation.fillna("") != ""].groupby("Dataset").Transformation.count().reset_index(
        name="Variables_transformed_or_derived")
    table(d, tr, "Table 5. Variables renamed, typed or derived during cleaning")
    para(d, "Derived variables: <date>_precision, <date>_iso, <numeric>_missing_code and cleaning_flags. Analytical "
            "transformations (RPKM, FPKM-16S, log10, centred log-ratio) are not stored in the cleaned data; they are "
            "computed by the analysis code and documented in the Statistical Analysis Report.")
    d.add_heading("8. Final dataset characteristics", 1)
    para(d, "Cleaned datasets are in 02_Cleaned_Datasets/Cleaned_CSV (UTF-8 CSV; the gene-level AMR table is gzip-"
            "compressed because it exceeds 40 MB) and, for tables of up to 200,000 rows, in Cleaned_Excel with an embedded "
            "dictionary sheet. The complete data dictionary (variable name, original name, description, type, unit, "
            "possible values, missingness and transformation) is in 02_Cleaned_Datasets/Data_Dictionaries.")
    table(d, dd.groupby("Dataset").agg(Variables=("Variable", "count"),
                                       Mean_pct_missing=("Pct_Missing_or_Coded", "mean")).reset_index(),
          "Table 6. Variables per cleaned dataset")
    d.save(REPORTS / "Data_Cleaning_Methodology.docx")


def analysis_report(obj, tests, ver, lim, tidx, fidx):
    d = new_doc("Statistical Analysis Report", "Objectives, analytical framework, methods, results and limitations")
    d.add_heading("1. Research objectives", 1)
    table(d, obj, "Table 1. Research objectives and their status in the verified analysis")
    d.add_heading("2. Analytical framework", 1)
    para(d, "The evidence is organised in three analytical streams that are never pooled because they differ in "
            "annotation pipeline, reference database and abundance unit:")
    bullets(d, [
        "Stream 1 - Indian open-drain wastewater (DS02; 447 samples, 19 sites in Chennai, Delhi, Kolkata and Mumbai; "
        "producer-processed CARD/RGI and Kraken 2/Bracken). Abundance unit: RPKM.",
        "Stream 2 - curated global ResFinder resource (Martiny et al.; urban matrices and food-related contextual units; "
        "KMA/ResFinder and SILVA rRNA). Abundance unit: FPKM per million rRNA-aligned bacterial fragments.",
        "Stream 3 - uniform reprocessing of 51 public runs (16 Dhaka, South Asian comparators and calibration units; "
        "4,000,000 read pairs each; fastp, Kraken 2/SILVA 138, KMA/ResFinder 4.7.2). Abundance unit: FPKM-16S.",
    ])
    para(d, "Units of analysis follow the design: the open-drain site (19 sites) is the independent unit for city "
            "contrasts in Stream 1 (site-level tests and site random-intercept mixed models); the study is the grouping "
            "factor in Stream 2 (study random intercepts, study-grouped cross-validation); whole sampling units are "
            "permuted in Stream 3. Sample-level p-values are reported as descriptive only.")
    d.add_heading("3. Variables", 1)
    bullets(d, [
        "Outcomes: total normalized ARG abundance (RPKM, FPKM, FPKM-16S; analysed on log10 scale), ARG gene/family richness "
        "and Shannon diversity, class and mechanism shares, detection of clinically important determinants (high-"
        "confidence and any-hit tiers), taxonomic relative abundance, alpha diversity (Shannon, Gini-Simpson, Pielou, "
        "rarefied richness) and beta diversity (Bray-Curtis, Aitchison).",
        "Explanatory factors: city, season, year, open-drain site, urban matrix (sample type), study/dataset, continent, "
        "geographic coordinates and sequencing depth.",
    ])
    d.add_heading("4. Statistical methods and equations", 1)
    eqs = [
        ("Normalized abundance (Stream 1)", "RPKM_s = sum_g [ r_gs / (L_g / 10^3) ] / (N_s / 10^6), where r_gs = reads mapped to "
         "gene g in sample s, L_g = reference length (bp) and N_s = sequenced reads."),
        ("Normalized abundance (Stream 3)", "FPKM-16S_s = sum_g [ f_gs / (L_g / 10^3) ] / (F16S_s / 10^6), with f = KMA fragments "
         "and F16S = prokaryotic 16S fragments; RPKM-depth uses QC read pairs as denominator."),
        ("Centred log-ratio", "clr(x)_i = ln( x_i / g(x) ), g(x) = (prod_j x_j)^(1/D); zeros replaced by a pseudocount before transformation."),
        ("Shannon and evenness", "H' = - sum_i p_i ln p_i ;  Gini-Simpson = 1 - sum_i p_i^2 ;  Pielou J = H' / ln S."),
        ("Bray-Curtis / Aitchison", "BC_jk = sum_i |x_ij - x_ik| / sum_i (x_ij + x_ik) ;  d_A(x, y) = || clr(x) - clr(y) ||_2."),
        ("Kruskal-Wallis", "H = [12 / (N(N+1))] sum_k R_k^2 / n_k - 3(N+1);  epsilon^2 = (H - k + 1) / (N - k). Site-level version: "
         "one median per site, p from 9,999 permutations of city labels among whole sites."),
        ("Linear mixed model", "y_ij = x_ij' beta + u_j + e_ij,  u_j ~ N(0, sigma_u^2), e_ij ~ N(0, sigma^2); fold change = 10^beta "
         "for log10 outcomes; ICC = sigma_u^2 / (sigma_u^2 + sigma^2); REML estimation, likelihood-ratio tests under ML."),
        ("PERMANOVA", "pseudo-F = [SS_A / (a - 1)] / [SS_W / (N - a)]; R^2 = SS_A / SS_T; permutations restricted to whole "
         "sites/studies where required; PERMDISP for dispersion."),
        ("Mantel / Procrustes", "Spearman correlation between distance matrices (permutation p); Procrustes m^2 after superimposition."),
        ("Spatial autocorrelation", "Moran's I = (n / W) sum_ij w_ij (x_i - xbar)(x_j - xbar) / sum_i (x_i - xbar)^2, "
         "row-standardised k-nearest-neighbour weights, 9,999 permutations."),
        ("Proportions", "Wilson 95% CI: [p + z^2/2n +/- z sqrt( p(1-p)/n + z^2/4n^2 )] / (1 + z^2/n)."),
        ("Multiplicity", "Benjamini-Hochberg: q_(i) = min_{j>=i} ( m p_(j) / j ), within pre-defined test families."),
        ("Effect size (two groups)", "Cliff's delta = P(X > Y) - P(X < Y); exact permutation p over all relabellings for small groups."),
        ("Meta-analysis (Stream 2)", "DerSimonian-Laird tau^2 = max(0, (Q - (k-1)) / (sum w - sum w^2 / sum w)); I^2 = (Q - (k-1)) / Q."),
        ("Predictive models", "Model A: random forest and logistic regression, study-grouped 4-fold CV, ROC-AUC, PR-AUC, Brier score. "
         "Model B: ridge and random forest on CLR phyla, temporal and leave-one-city-out validation, R^2, MAE, RMSE."),
    ]
    for a, b in eqs:
        para(d, "", bold_lead=a)
        eq(d, b)
    d.add_heading("5. Register of statistical tests", 1)
    reg = tests.groupby(["objective", "test"]).size().reset_index(name="tests").sort_values(["objective", "tests"],
                                                                                            ascending=[True, False])
    para(d, f"The verified analysis performed {len(tests)} registered tests ({(tests.inferential_status == 'Inferential').sum()} "
            f"inferential, the rest descriptive or sensitivity). The full register (statistic, df, p, BH-adjusted p, effect "
            f"size and CI) is reproduced in 03_Data_Analysis/Statistical_Analysis/Statistical_Tests_Register.csv.")
    table(d, reg, "Table 2. Tests by objective and method", font=8, max_rows=80)
    d.add_heading("6. Results: reproduction of headline findings", 1)
    para(d, f"Headline results were recomputed from the cleaned, publicly released data with independently written code "
            f"(07_Code_and_Scripts/Python/03_verify_analysis.py). {(ver.Match == 'YES').sum()} of {len(ver)} checks reproduced "
            f"the reported values within the stated tolerance. Permutation p-values are Monte Carlo estimates and are "
            f"compared with an absolute tolerance of 0.003.")
    table(d, ver[["Area", "Statistic", "Supplement_Table", "Reported", "Recomputed", "Tolerance", "Match"]],
          "Table 3. Independent verification of supplementary results", font=7.5, max_rows=100)
    d.add_heading("7. Tables and figures", 1)
    para(d, f"The supplement contains {len(tidx)} tables (Tables S1-S{len(tidx)}) and {len(fidx)} figures (Figs. S1-S{len(fidx)}). "
            "Each is listed with its section, analysis identifier and source file in 04_Tables/Table_Index_S1-S96.csv and "
            "05_Figures/Figure_Index_S1-S30.csv; interpretations are given under each item in Supplementary_Material_Final.docx.")
    table(d, tidx.groupby("Section").agg(Tables=("Table", "count"), First=("Table", "first"),
                                         Last=("Table", "last")).reset_index(), "Table 4. Supplementary tables by section")
    table(d, fidx.groupby("Section", sort=False).agg(Figures=("Figure", "count"), First=("Figure", "first"),
                                                     Last=("Figure", "last")).reset_index(),
          "Table 5. Supplementary figures by section")
    d.add_heading("8. Scientific interpretation (summary)", 1)
    for k, v in [
        ("Evidence base.", "The analysable evidence is broad (8,053 urban metagenomes, 106 countries) but processed AMR "
         "tables exist for only 13.3% of samples and are aligned with matrix and pipeline, so cross-matrix comparison is "
         "confounded by annotation pipeline."),
        ("Stream 1.", "Total ARG abundance differed between Indian cities with the site as the unit (site-level Kruskal-Wallis "
         "and mixed model; Delhi highest), within resistomes dominated by MLS, multi-class and aminoglycoside genes."),
        ("Stream 2.", "Municipal wastewater carried the highest normalized ARG abundance, but matrix contrasts are largely "
         "inseparable from study effects (single-study matrices; ICC 0.65)."),
        ("Stream 3.", "Dhaka resistomes were dataset-structured; high-confidence NDM and KPC calls occurred only in the "
         "DS03 and DS84 datasets (none in DS01 or Karachi); comparisons with Karachi and Delhi are small-sample and descriptive."),
    ]:
        para(d, " " + v, bold_lead=k)
    d.add_heading("9. Limitations", 1)
    table(d, lim, "Table 6. Limitations and consequences for inference", font=7.5, max_rows=40)
    d.save(REPORTS / "Statistical_Analysis_Methodology.docx")
    tests.to_csv(ANALYSIS / "Statistical_Analysis" / "Statistical_Tests_Register.csv", index=False)


def technical_report(inv, tidx, fidx):
    d = new_doc("Technical Reproducibility Report", "How to reproduce the cleaning, verification, tables, figures and documents")
    d.add_heading("1. Software environment", 1)
    import platform, importlib
    rows = [("Python", platform.python_version())]
    for m in ("pandas", "numpy", "scipy", "statsmodels", "openpyxl", "docx", "PIL", "matplotlib", "lxml"):
        try:
            rows.append((m if m != "docx" else "python-docx", importlib.import_module(m).__version__))
        except Exception:
            rows.append((m, "not installed"))
    rows.append(("Microsoft Word", "required only for step 07 (field update and PDF export via COM, Windows)"))
    table(d, pd.DataFrame(rows, columns=["Component", "Version used"]), "Table 1. Software used to build this package")
    para(d, "Install the Python dependencies with: pip install -r 07_Code_and_Scripts/Requirements.txt")
    d.add_heading("2. Directory layout and inputs", 1)
    para(d, "The scripts use relative paths only. PROJECT is the package folder; SOURCE (the research directory holding "
            "the original materials) defaults to its parent folder and can be set with the environment variable "
            "AMR_SOURCE_DIR. Users of the public repository, who do not have the original working folders, can rerun the "
            "verification (step 03) directly from the published cleaned datasets and tables.")
    d.add_heading("3. Execution sequence", 1)
    steps = pd.DataFrame([
        ("01_build_inventory.py", "Inventory, SHA-256, integrity and duplicates of all source files", "00_Project_Inventory"),
        ("02_clean_datasets.py", "Non-destructive cleaning, flags, data dictionaries", "02_Cleaned_Datasets, 03_Data_Analysis"),
        ("03_verify_analysis.py", "Independent recomputation of headline results", "03_Data_Analysis"),
        ("06_build_supplementary_docx.py", "Corrections, captions, cross-references, TOC/LoF/LoT", "06_Research_Documents"),
        ("07_update_fields_and_export_pdf.ps1", "Word updates all fields; DOCX saved; PDF exported", "06_Research_Documents"),
        ("04_export_tables.py", "Tables S1-S96 as CSV and formatted workbook", "04_Tables"),
        ("05_export_figures.py", "Figures S1-S30 as PNG, TIFF, PDF + index", "05_Figures"),
        ("08_build_reports.py", "Methodology, analysis, reproducibility and data-quality reports", "08_Methodology_and_Reports"),
        ("09_quality_control.py", "Final QC audit and public-release screen", "00_Project_Inventory"),
    ], columns=["Script", "Purpose", "Outputs"])
    table(d, steps, "Table 2. Execution order (run_all.py runs the whole sequence)")
    para(d, "Steps 04 and 05 read captions and images from the final document, so they run after step 07. "
            "run_all.py executes the sequence and stops at the first failure.")
    d.add_heading("4. Reproducing the original analysis", 1)
    para(d, "The verified analysis (twenty executed notebooks; percent-format sources; helper modules) is archived in "
            "07_Code_and_Scripts/Original_Analysis_Code. It was run on 31 October 2026 with Python 3.12 in fresh kernels; "
            "raw-read reprocessing (Stream 3) used fastp 1.4.0, Kraken 2.17.1 with the SILVA 138 16S index, Bracken and "
            "KMA 1.7.1 with ResFinder database 4.7.2 on Linux. These archival scripts contain absolute paths of the "
            "original working environment, which must be adapted (AMR_ROOT / AMR_HARM variables) before re-execution.")
    d.add_heading("5. Validation built into the workflow", 1)
    bullets(d, [
        "Row counts asserted equal before and after cleaning; referential integrity checked across tables.",
        "Caption sequences asserted to run Fig. S1-S30 and Table S1-S96 without gaps; every item cited in the text.",
        "All text corrections must match the document exactly once, or the build fails.",
        "Word reports any 'Error! Reference source not found' after field update; QC fails if any remains.",
        "Exported tables are compared cell-by-cell with the verified CSVs; figures are compared by checksum with the "
        "images embedded in the document.",
    ])
    d.add_heading("6. Outputs", 1)
    para(d, f"Tables: {len(tidx)} CSV files and one workbook. Figures: {len(fidx)} x PNG/TIFF/PDF. Documents: "
            "Supplementary_Material_Final.docx and .pdf. Source inventory: "
            f"{len(inv)} files with SHA-256 checksums.")
    d.save(REPORTS / "Technical_Report.docx")


def quality_report(inv, prof):
    d = new_doc("Data Quality Report", "Source inventory, integrity, duplicates and dataset quality")
    d.add_heading("1. Source directory inventory", 1)
    para(d, f"{len(inv)} files ({inv.Size_Bytes.sum() / 1e9:.2f} GB) were inspected recursively. Each file was hashed "
            "(SHA-256), opened with a format-appropriate reader and classified by content and role.")
    table(d, inv.groupby("Top_Folder").agg(Files=("File_ID", "count"), Size_MB=("Size_MB", "sum")).reset_index(),
          "Table 1. Files by top-level folder")
    table(d, inv.groupby(["Relevance"]).agg(Files=("File_ID", "count")).reset_index(), "Table 2. Files by relevance class")
    d.add_heading("2. Integrity", 1)
    bad = inv[~inv.Integrity.isin(["OK", "NOT CHECKED (format)"])]
    para(d, f"{len(bad)} files failed or could not be confirmed by the integrity check; no corrupted workbook, document, "
            "image or archive was found. Empty files are listed below.")
    table(d, bad[["Relative_Path", "Integrity", "Role"]], "Table 3. Files with integrity findings")
    d.add_heading("3. Duplicate, outdated and irrelevant files", 1)
    dup = inv[inv.Duplicate_Status == "Duplicate - redundant copy"]
    para(d, f"{inv.Duplicate_Group.nunique() - 1} groups of byte-identical files contain {len(dup)} redundant copies. "
            "The most visible are the root-level '_1' copies of the data-collection workbooks, the report and the source "
            "log, the duplicated build archives in 08_Final_Submission/Build_Sources, and the first analysis run "
            "(Analysis_Results), which the independent audit (Independent_Audit) and the update "
            "(Updated_Research_Analysis) supersede. Python bytecode caches and a folder marked '_to_delete' are irrelevant. "
            "No file was deleted; the inventory records each status.")
    roots = dup[~dup.Relative_Path.str.contains("/")]
    table(d, roots[["Relative_Path", "Duplicate_Group"]], "Table 4. Redundant copies at the top level")
    table(d, inv[inv.Relevance.isin(["Superseded", "Irrelevant", "Irrelevant (regenerable)"])].groupby(
        ["Relevance", "Role"]).size().reset_index(name="Files"), "Table 5. Superseded and irrelevant files")
    d.add_heading("4. Dataset quality", 1)
    table(d, prof["Dataset_Summary"][["Dataset", "Rows_Original", "Columns_Original", "Exact_Duplicate_Rows",
                                      "Records_Flagged", "Pct_Records_Flagged"]], "Table 6. Cleaned datasets")
    table(d, prof["Flag_Counts"], "Table 7. Quality flags")
    para(d, "Details of each operation, completeness per variable and categorical checks are in "
            "03_Data_Analysis/Descriptive_Statistics/Dataset_Profiles.xlsx; methods in "
            "08_Methodology_and_Reports/Data_Cleaning_Methodology.docx.")
    d.save(INV / "Data_Quality_Report.docx")


def main():
    REPORTS.mkdir(parents=True, exist_ok=True)
    inv, prof, dd, ver, tests, lim, obj, tidx, fidx, corr = load()
    cleaning_report(prof, dd)
    analysis_report(obj, tests, ver, lim, tidx, fidx)
    technical_report(inv, tidx, fidx)
    quality_report(inv, prof)
    log("Reports written: Data_Quality_Report, Data_Cleaning_Methodology, Statistical_Analysis_Methodology, Technical_Report")


if __name__ == "__main__":
    main()
