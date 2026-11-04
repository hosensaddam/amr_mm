"""Build AMR_Analysis_Results_Verified.xlsx from the analysis outputs (tables/, figures/, catalog/).
Every table is written in full (no truncation) as an Excel table with filters; figures are embedded
at their native 300-dpi resolution (displayed scaled)."""
import re, json, sys
from pathlib import Path
import pandas as pd, numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo
from openpyxl.formatting.rule import CellIsRule, ColorScaleRule
from openpyxl.drawing.image import Image as XLImage
from openpyxl.chart import BarChart, Reference
from PIL import Image as PILImage
sys.path.insert(0, str(Path(__file__).resolve().parent))
from amr_common import OUT, TAB, FIG, VER, load_catalog, ANALYSIS_DATE

WB_NAME = "AMR_Analysis_Results_Verified.xlsx"
SHEETS = ["Executive_Summary", "Research_Objectives", "Method_Selection", "Dataset_Inventory", "Data_Dictionary", "Data_Quality_Audit",
          "Sample_Distribution", "Geographic_Coverage", "Metadata_Missingness", "Descriptive_Statistics", "Taxonomic_Abundance", "Alpha_Diversity",
          "Beta_Diversity", "AMR_Gene_Summary", "AMR_Resistance_Classes", "Comparative_Analysis", "PCA_Results", "Correlation_Analysis",
          "Environmental_Analysis", "Spatial_Analysis", "Predictive_Models", "Statistical_Tests", "Sensitivity_Analysis", "Figure_Gallery",
          "Results_Interpretation", "Limitations", "Verification_Report"]
DESC = {
 "Executive_Summary": "Headline counts and verified findings. Every number is traceable to a table in this workbook and to the CSV/notebook named in Verification_Report.",
 "Research_Objectives": "The five research objectives and the extent to which the available data addressed them.",
 "Method_Selection": "Objective-to-method framework: question, data, variables, methods, software, outputs, actual findings, limitations, status.",
 "Dataset_Inventory": "Files inspected (with MD5), harmonized tables analysed, and the sheets of the project workbooks.",
 "Data_Dictionary": "Column-level dictionary of the harmonized tables with observed and missing counts (source missing-value codes treated as missing, never zero).",
 "Data_Quality_Audit": "Inclusion, duplicate, linkage, geographic, sample-type, depth and annotation checks, plus the exclusion log.",
 "Sample_Distribution": "Included samples by sample type and dataset.",
 "Geographic_Coverage": "Countries, continents and cities represented; proposal target cities.",
 "Metadata_Missingness": "Missingness of core metadata fields among included samples.",
 "Descriptive_Statistics": "Evidence-base overview, platforms, read layout, collection years and annotation coverage.",
 "Taxonomic_Abundance": "DS02 domain and prokaryotic phylum relative abundance (Bracken).",
 "Alpha_Diversity": "Alpha diversity (DS02 phylum level; Martiny et al. genus level by matrix).",
 "Beta_Diversity": "PERMANOVA and PERMDISP results for taxonomic and resistome distances.",
 "AMR_Gene_Summary": "ARG gene families, loads and clinically important determinants (CARD/RGI and ResFinder/KMA analysed separately).",
 "AMR_Resistance_Classes": "Resistance-class groups and mechanisms (shares of total normalised abundance).",
 "Comparative_Analysis": "Adjusted models (OLS with HC3 SEs; linear mixed models with study random effect) and between-study heterogeneity.",
 "PCA_Results": "CLR-PCA explained variance, loadings and contributions; PCoA/NMDS; DS02 taxonomic clusters.",
 "Correlation_Analysis": "Taxonomy-resistome concordance, Spearman correlations (BH), co-occurrence network, diversity associations.",
 "Environmental_Analysis": "Seasonal and monthly patterns (DS02 longitudinal series). No external environmental covariates were linked to AMR samples.",
 "Spatial_Analysis": "Mapping coverage, regional comparison, spatial autocorrelation and the feasibility of spatial methods.",
 "Predictive_Models": "Exploratory model designs, grouped validation performance, feature importance and feasibility of predictive targets.",
 "Statistical_Tests": "Register of all tests and model metrics with BH adjustment within families, plus key test tables.",
 "Sensitivity_Analysis": "Robustness of key estimates to thresholds, denominators, aggregation, duplicate handling and seeds.",
 "Figure_Gallery": "All figures (300 dpi PNG) with number, title and caption.",
 "Results_Interpretation": "Findings register: result, denominator, statistical evidence, table/figure, interpretation and limitations.",
 "Limitations": "Limitations register with consequences for inference and remedies.",
 "Verification_Report": "Independent verification results, cross-file consistency checks, PERMANOVA cross-check and table catalogue.",
}
FIG_SHEET = {1: "Metadata_Missingness", 2: "Sample_Distribution", 3: "Descriptive_Statistics", 4: "Descriptive_Statistics", 5: "Taxonomic_Abundance", 6: "Alpha_Diversity",
             7: "Taxonomic_Abundance", 8: "Alpha_Diversity", 9: "AMR_Gene_Summary", 10: "AMR_Gene_Summary", 11: "AMR_Gene_Summary", 12: "Environmental_Analysis",
             13: "AMR_Resistance_Classes", 14: "AMR_Gene_Summary", 15: "Correlation_Analysis", 16: "Correlation_Analysis", 17: "Comparative_Analysis",
             18: "Correlation_Analysis", 19: "PCA_Results", 20: "PCA_Results", 21: "Beta_Diversity", 22: "Geographic_Coverage", 23: "Geographic_Coverage",
             24: "Spatial_Analysis", 25: "Spatial_Analysis", 26: "Spatial_Analysis", 27: "Predictive_Models", 28: "Predictive_Models", 29: "Sensitivity_Analysis",
             30: "Executive_Summary"}
HDR_FILL = PatternFill("solid", fgColor="1F4E79"); HDR_FONT = Font(bold=True, color="FFFFFF"); TITLE_FONT = Font(bold=True, size=14, color="1F4E79")
SIG_FILL = PatternFill("solid", fgColor="E2EFDA"); thin = Side(style="thin", color="BFBFBF")

def clean(v):
    if v is None: return None
    if isinstance(v, (float, np.floating)):
        return None if (np.isnan(v) or np.isinf(v)) else float(v)
    if isinstance(v, (np.integer,)): return int(v)
    if isinstance(v, (np.bool_, bool)): return bool(v)
    if isinstance(v, int): return v  # AUDIT FIX AUD-01: Python ints (from DataFrame.to_numpy(dtype=object)) were written as text
    s = str(v)
    if s.startswith("="): s = "'" + s
    return s[:32000]

def write_table(ws, df, r0, c0, tname, used):
    df = df.copy(); cols = []
    for c in df.columns:
        c = str(c) or "col"; base = c; k = 2
        while c in cols: c = f"{base}_{k}"; k += 1
        cols.append(c)
    df.columns = cols
    for j, c in enumerate(cols):
        cell = ws.cell(r0, c0 + j, c); cell.fill = HDR_FILL; cell.font = HDR_FONT; cell.alignment = Alignment(wrap_text=True, vertical="center")
    vals = df.to_numpy(dtype=object)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            ws.cell(r0 + 1 + i, c0 + j, clean(vals[i, j]))
    r1 = r0 + max(len(df), 1); cend = c0 + len(cols) - 1
    if len(df):
        name = re.sub(r"[^A-Za-z0-9_]", "_", "T_" + tname)[:200]; n = name; k = 2
        while n in used: n = f"{name}_{k}"; k += 1
        used.add(n)
        t = Table(displayName=n, ref=f"{get_column_letter(c0)}{r0}:{get_column_letter(cend)}{r1}")
        t.tableStyleInfo = TableStyleInfo(name="TableStyleLight9", showRowStripes=True)
        ws.add_table(t)
    for j, c in enumerate(cols):
        col = df.iloc[:, j]; L = get_column_letter(c0 + j)
        num = pd.api.types.is_numeric_dtype(col)
        if num and not pd.api.types.is_bool_dtype(col):
            isint = pd.api.types.is_integer_dtype(col) or (col.dropna() % 1 == 0).all()
            fmt = "#,##0" if isint else ("0.000E+00" if re.search(r"(^p$|p_value|p_BH|p_adj|_p$|^Q_p|LRT p)", c) else "0.000")
            for i in range(len(df)): ws.cell(r0 + 1 + i, c0 + j).number_format = fmt
            if re.search(r"(^p$|p_value|p_BH|p_adj)", c) and len(df):
                ws.conditional_formatting.add(f"{L}{r0+1}:{L}{r1}", CellIsRule(operator="lessThan", formula=["0.05"], fill=SIG_FILL))
            if re.search(r"pct_missing|^%|_pct$", c) and len(df) and c not in ("n",):
                ws.conditional_formatting.add(f"{L}{r0+1}:{L}{r1}", ColorScaleRule(start_type="num", start_value=0, start_color="FFFFFF", end_type="num", end_value=100, end_color="F8696B"))
        width = min(max([len(c)] + [len(str(x)) for x in col.head(300).tolist()]) + 2, 70)
        ws.column_dimensions[L].width = max(ws.column_dimensions[L].width or 0, width, 9)
        if width >= 60:
            for i in range(len(df)): ws.cell(r0 + 1 + i, c0 + j).alignment = Alignment(wrap_text=True, vertical="top")
    return r1, cend

def add_image(ws, path, anchor, display_w=760):
    with PILImage.open(path) as im: w, h = im.size
    img = XLImage(str(path)); img.width = display_w; img.height = int(h * display_w / w)
    ws.add_image(img, anchor); return img.height

def build(out_path=None):
    out_path = Path(out_path or OUT / WB_NAME)
    tabs = load_catalog("tables"); figs = load_catalog("figures").sort_values("number")
    wb = Workbook(); wb.remove(wb.active); used = set(); pos = {}; endrow = {}
    for sh in SHEETS:
        ws = wb.create_sheet(sh)
        ws["A1"] = sh.replace("_", " "); ws["A1"].font = TITLE_FONT
        ws["A2"] = DESC.get(sh, ""); ws["A2"].alignment = Alignment(wrap_text=False); ws["A2"].font = Font(italic=True, color="404040")
        ws["A3"] = f"Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping | analysis date {ANALYSIS_DATE} | source CSVs in Analysis_Results/tables"
        ws["A3"].font = Font(size=8, color="7F7F7F")
        row = 5; maxc = 1; first_hdr = None; ntab = 0
        sub = tabs[tabs.sheet == sh].sort_values("file")
        if sh == "Executive_Summary":
            sub = pd.concat([tabs[tabs.file == "R08_executive_summary.csv"], tabs[tabs.file == "D01_evidence_base_overview.csv"]])
        if sh == "Verification_Report":
            vfiles = sorted(VER.glob("V*.csv"))
            for vf in vfiles:
                df = pd.read_csv(vf)
                ws.cell(row, 1, f"Table {vf.stem}: independent verification ({vf.name})").font = Font(bold=True)
                r1, ce = write_table(ws, df, row + 1, 1, vf.stem, used); pos[vf.name] = (sh, row + 1, 1, r1, ce)
                first_hdr = first_hdr or row + 1; row = r1 + 3; maxc = max(maxc, ce); ntab += 1
        for _, t in sub.iterrows():
            df = pd.read_csv(TAB / t.file, low_memory=False)
            ws.cell(row, 1, f"Table {Path(t.file).stem}: {t.title}  [n rows = {len(df):,}; source: tables/{t.file}; notebook {t.notebook}]").font = Font(bold=True)
            r1, ce = write_table(ws, df, row + 1, 1, Path(t.file).stem, used)
            pos[t.file] = (sh, row + 1, 1, r1, ce)
            first_hdr = first_hdr or row + 1; row = r1 + 3; maxc = max(maxc, ce); ntab += 1
        endrow[sh] = row
        if ntab == 1 and first_hdr: ws.freeze_panes = ws.cell(first_hdr + 1, 1)
        else: ws.freeze_panes = "A5"
        # figures for this sheet
        frow = 5; fcol = get_column_letter(min(maxc + 2, 40))
        if sh == "Executive_Summary": fcol = "D"; frow = row + 1
        for _, f in figs[figs.number.map(FIG_SHEET) == sh].iterrows():
            ws[f"{fcol}{frow}"] = f"Figure {int(f.number)}. {f.title}"; ws[f"{fcol}{frow}"].font = Font(bold=True)
            h = add_image(ws, FIG / f.file, f"{fcol}{frow+1}", 760 if sh != "Executive_Summary" else 980)
            frow += int(h / 20) + 4
    # Figure gallery
    ws = wb["Figure_Gallery"]; row = endrow.get("Figure_Gallery", 5) + 2
    for _, f in figs.iterrows():
        ws.cell(row, 1, f"Figure {int(f.number)}. {f.title}").font = Font(bold=True, size=11)
        ws.cell(row + 1, 1, f"{f.caption}  [Objective: {f.objective}; notebook {f.notebook}; data: {f.source_data}; file: figures/{f.file}]").alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells(start_row=row + 1, start_column=1, end_row=row + 1, end_column=14); ws.row_dimensions[row + 1].height = 62
        h = add_image(ws, FIG / f.file, f"A{row+2}", 900); row += int(h / 20) + 6
    ws.column_dimensions["A"].width = 14
    # native Excel charts
    def bar(sheet, file, cat_col, val_col, title, ytitle, anchor, n=None, horizontal=True):
        if file not in pos: return
        sh_, r0, c0, r1, _ = pos[file]; ws_ = wb[sheet]
        hdr = [ws_.cell(r0, c).value for c in range(c0, c0 + 200)]
        ci, vi = hdr.index(cat_col) + c0, hdr.index(val_col) + c0
        r_end = r1 if n is None else min(r1, r0 + n)
        ch = BarChart(); ch.type = "bar" if horizontal else "col"; ch.title = title; ch.y_axis.title = ytitle; ch.legend = None
        ch.add_data(Reference(ws_, min_col=vi, min_row=r0, max_row=r_end), titles_from_data=True)
        ch.set_categories(Reference(ws_, min_col=ci, min_row=r0 + 1, max_row=r_end)); ch.height = 8; ch.width = 16
        ws_.add_chart(ch, anchor)
    bar("Sample_Distribution", "D02_by_sample_type.csv", "Sample_Type", "samples", "Included samples by sample type (native chart)", "Samples", "B" + str(pos.get("D03_by_dataset.csv", ("", 60, 0, 60))[3] + 4) if "D03_by_dataset.csv" in pos else "B200")
    bar("Geographic_Coverage", "D04_by_country.csv", "Country", "samples", "Top 15 countries by included samples (native chart)", "Samples", "H5", n=15)
    bar("AMR_Gene_Summary", "A04_ds02_arg_load_by_city.csv", "City", "median_total_RPKM", "DS02 median total ARG RPKM by city (native chart)", "RPKM", "H5", horizontal=False)
    bar("AMR_Gene_Summary", "A12_rf_arg_load_by_sample_type.csv", "Sample_Type", "median_total_FPKM", "ResFinder median total ARG FPKM by matrix (native chart)", "FPKM", "H22")
    for ws in wb.worksheets:
        ws.sheet_view.zoomScale = 90
    wb.properties.title = "AMR Analysis Results (Verified)"; wb.properties.creator = "Analysis_Results pipeline (Python/openpyxl)"
    wb.save(out_path)
    json.dump({k: list(v) for k, v in pos.items()}, open(OUT / "logs" / "workbook_table_positions.json", "w"), indent=0)
    print("saved", out_path, "| sheets:", len(wb.sheetnames), "| tables:", len(pos), "| figures:", len(figs))
    return out_path

if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else None)
