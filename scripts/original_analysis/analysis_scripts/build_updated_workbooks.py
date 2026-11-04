"""Build the updated Excel deliverables from the updated analysis outputs (tables/, figures/, catalog/):
  AMR_Updated_Final_Analysis.xlsx   - corrected global analysis + audit corrections + integration + Dhaka summary
  Dhaka_City_Advanced_Analysis.xlsx - dedicated Dhaka workbook (17 analysis sheets + README)
  Updated_Audit_Issue_Register.xlsx - all 20 audit issues with the update status and evidence
Every table is written in full as an Excel table; integers are written as numbers (AUD-01 rule);
figures are embedded from the 300-dpi PNGs; native Excel charts plot numeric cells.
Usage: python build_updated_workbooks.py <output_dir> [--with-verification]"""
import sys, re, json
from pathlib import Path
import pandas as pd, numpy as np
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.chart import BarChart, Reference
sys.path.insert(0, str(Path(__file__).resolve().parent))
from amr_common import OUT, TAB, FIG, VER, load_catalog, ANALYSIS_DATE, UPDATE_LABEL
from build_workbook import write_table, add_image, SHEETS as GSHEETS, DESC as GDESC, FIG_SHEET as GFIG, TITLE_FONT

OUTDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT
UPDATED_REPLACES = {"R02_statistical_tests_register.csv": "R02u_updated_statistical_tests_register.csv", "R03_findings_register.csv": "R03u_updated_findings_register.csv",
                    "R04_limitations.csv": "R04u_updated_limitations.csv", "R06_figure_catalog.csv": "R06u_updated_figure_catalog.csv",
                    "R07_table_catalog.csv": "R07u_updated_table_catalog.csv", "R08_executive_summary.csv": "R08u_updated_executive_summary.csv"}
DHAKA_SHEETS = ["Dhaka_Dataset_Inventory", "Dhaka_Sample_Metadata", "Dhaka_Sample_Type", "Dhaka_Data_Quality", "Dhaka_Sequencing_QC", "Dhaka_Taxonomic_Abundance",
                "Dhaka_Alpha_Diversity", "Dhaka_Beta_Diversity", "Dhaka_AMR_Genes", "Dhaka_Resistance_Classes", "Dhaka_Comparative_Analysis", "Dhaka_Statistical_Tests",
                "Dhaka_PCA_Results", "Dhaka_Spatial_Analysis", "Dhaka_Figure_Gallery", "Dhaka_Results_Interpretation", "Dhaka_Limitations"]
DDESC = {
 "Dhaka_Dataset_Inventory": "Sample, run and study accessions, platform, layout, read length, ENA total bases, download URLs, file sizes and ENA MD5 for the 16 Dhaka samples; dataset summary.",
 "Dhaka_Sample_Metadata": "Sample-level metadata after the AUD-04 correction (sample type verified against BioSample isolation_source and BioProject descriptions).",
 "Dhaka_Sample_Type": "Verified sample types by dataset (Municipal = urban sewage / sewage-treatment-plant inlet; Surface water = pond).",
 "Dhaka_Data_Quality": "Per-sample read QC and depth assessment; acquisition and QC verification checks.",
 "Dhaka_Sequencing_QC": "fastp / Kraken 2 / KMA metrics for every uniformly reprocessed run (16 Dhaka + comparators + calibration), 4,000,000 read pairs each.",
 "Dhaka_Taxonomic_Abundance": "Relative abundance (% of prokaryotic 16S fragments assigned at that rank; the share assigned at each rank is in DH10) at phylum, class, order, family and genus; dominant taxa per sample; taxonomic resolution. Species not resolvable with the 16S index.",
 "Dhaka_Alpha_Diversity": "Rarefied genus richness, Shannon, Gini-Simpson, Pielou; ARG gene richness and Shannon; depth checks.",
 "Dhaka_Beta_Diversity": "PERMANOVA/PERMDISP (descriptive), within/between-dataset dissimilarity, microbiome-resistome Mantel test.",
 "Dhaka_AMR_Genes": "Per-sample normalized ARG abundance (FPKM-16S, RPKM-depth), all detected genes with high-confidence and any-hit detection, template-level results, clinically important determinants.",
 "Dhaka_Resistance_Classes": "Antibiotic-class groups and resistance mechanisms (share of FPKM-16S), by dataset and sample type.",
 "Dhaka_Comparative_Analysis": "Within-Dhaka co-located sites and temporal series; matched-protocol South Asian sewage comparison; class composition and determinants by city group.",
 "Dhaka_Statistical_Tests": "All Dhaka comparative tests with the independent unit, group sizes, p-values (exact/Monte-Carlo permutation where groups are small), BH adjustment, effect sizes and inferential status.",
 "Dhaka_PCA_Results": "PCoA and CLR-PCA scores, explained variance and loadings for Dhaka samples.",
 "Dhaka_Spatial_Analysis": "AUD-18 geographic verification against official DNCC/DSCC boundaries; geographic status; co-located sampling points.",
 "Dhaka_Figure_Gallery": "Dhaka, comparative, integration and sensitivity figures (Figures 31-45) with captions.",
 "Dhaka_Results_Interpretation": "Dhaka findings DF01-DF10 with denominators, statistical evidence, interpretation and limitations.",
 "Dhaka_Limitations": "Limitations specific to the Dhaka reprocessing and comparisons.",
}
NEW_GLOBAL = ["Audit_Corrections", "Global_Integration", "Dhaka_City_Summary"]
NDESC = {"Audit_Corrections": "Resolution of every open audit issue (AUD-04, AUD-11 to AUD-20) with before/after values and verification evidence.",
         "Global_Integration": "Calibration of the uniform reprocessing pipeline against the ResFinder resource; CARD/RGI vs ResFinder/KMA cross-check for Delhi; updated target-city coverage.",
         "Dhaka_City_Summary": "Key Dhaka results (full detail in Dhaka_City_Advanced_Analysis.xlsx): inventory, QC, AMR summary, determinants, comparative tests and Figures 31-43."}
DH_SUMMARY_TABLES = ["DH02_dhaka_dataset_summary.csv", "DH04_dhaka_sample_type_by_dataset.csv", "DH11_dhaka_amr_sample_summary.csv", "DH14_dhaka_resistance_classes.csv",
                     "DH16_dhaka_clinical_determinants.csv", "DH28_south_asia_sewage_city_summary.csv", "DH29_dhaka_comparative_statistical_tests.csv", "DH33_dhaka_geographic_status_summary.csv"]

def header(ws, title, desc, src):
    ws["A1"] = title; ws["A1"].font = TITLE_FONT
    ws["A2"] = desc; ws["A2"].font = Font(italic=True, color="404040")
    ws["A3"] = f"Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping | {UPDATE_LABEL} | {src}"; ws["A3"].font = Font(size=8, color="7F7F7F")

def put_tables(wb, ws, files, tabs, used, pos, row=5):
    maxc = 1
    for f in files:
        t = tabs[tabs.file == f]
        title = t.title.iloc[0] if len(t) else f; nb = t.notebook.iloc[0] if len(t) else ""
        df = pd.read_csv(TAB / f, low_memory=False) if (TAB / f).exists() else pd.read_csv(VER / f, low_memory=False)
        ws.cell(row, 1, f"Table {Path(f).stem}: {title}  [n rows = {len(df):,}; source: tables/{f}; notebook {nb}]").font = Font(bold=True)
        r1, ce = write_table(ws, df, row + 1, 1, Path(f).stem, used); pos[f] = (ws.title, row + 1, 1, r1, ce)
        row = r1 + 3; maxc = max(maxc, ce)
    return row, maxc

def put_figures(ws, figs, col, row, width=760):
    for _, f in figs.iterrows():
        ws[f"{col}{row}"] = f"Figure {int(f.number)}. {f.title}"; ws[f"{col}{row}"].font = Font(bold=True)
        h = add_image(ws, FIG / f.file, f"{col}{row+1}", width); row += int(h / 20) + 4
    return row

def gallery(ws, figs, row):
    for _, f in figs.iterrows():
        ws.cell(row, 1, f"Figure {int(f.number)}. {f.title}").font = Font(bold=True, size=11)
        ws.cell(row + 1, 1, f"{f.caption}  [notebook {f.notebook}; data: {f.source_data}; file: figures/{f.file}; TIFF: figures/{f.get('tiff_file', '')}; print size {f.get('width_mm', '')} x {f.get('height_mm', '')} mm]").alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells(start_row=row + 1, start_column=1, end_row=row + 1, end_column=14); ws.row_dimensions[row + 1].height = 70
        h = add_image(ws, FIG / f.file, f"A{row+2}", 900); row += int(h / 20) + 6
    ws.column_dimensions["A"].width = 14

def bar(wb, pos, file, cat_col, val_col, title, ytitle, anchor, n=None, horizontal=True):
    if file not in pos: return False
    sh, r0, c0, r1, _ = pos[file]; ws = wb[sh]
    hdr = [ws.cell(r0, c).value for c in range(c0, c0 + 300)]
    if cat_col not in hdr or val_col not in hdr: return False
    ci, vi = hdr.index(cat_col) + c0, hdr.index(val_col) + c0
    r_end = r1 if n is None else min(r1, r0 + n)
    ch = BarChart(); ch.type = "bar" if horizontal else "col"; ch.title = title; ch.y_axis.title = ytitle; ch.legend = None
    ch.add_data(Reference(ws, min_col=vi, min_row=r0, max_row=r_end), titles_from_data=True)
    ch.set_categories(Reference(ws, min_col=ci, min_row=r0 + 1, max_row=r_end)); ch.height = 8; ch.width = 16
    ws.add_chart(ch, anchor); return True

def verification_tables():
    return sorted(p.name for p in VER.glob("V*.csv"))

def build_master(path, extra_ver=()):
    tabs = load_catalog("tables"); figs = load_catalog("figures").sort_values("number")
    tabs = tabs[~tabs.file.isin(list(UPDATED_REPLACES))]
    wb = Workbook(); wb.remove(wb.active); used, pos, endrow = set(), {}, {}
    sheets = GSHEETS[:-3] + NEW_GLOBAL + GSHEETS[-3:]
    figmap = dict(GFIG); figmap.update({n: "Dhaka_City_Summary" for n in range(31, 44)}); figmap.update({44: "Global_Integration", 45: "Sensitivity_Analysis"})
    for sh in sheets:
        ws = wb.create_sheet(sh); header(ws, sh.replace("_", " "), {**GDESC, **NDESC}.get(sh, ""), "source CSVs in Updated_Research_Analysis/04_Tables")
        if sh == "Dhaka_City_Summary": files = DH_SUMMARY_TABLES
        elif sh == "Executive_Summary": files = ["R08u_updated_executive_summary.csv", "D01_evidence_base_overview.csv"]
        else: files = sorted(tabs[(tabs.sheet == sh) & ~tabs.file.str.startswith("DH")].file)
        if sh == "Verification_Report": files = list(extra_ver) + [f for f in files]
        row, maxc = put_tables(wb, ws, files, tabs, used, pos)
        endrow[sh] = row
        ws.freeze_panes = "A5"
        fsel = figs[figs.number.map(figmap) == sh]
        if sh != "Figure_Gallery":
            put_figures(ws, fsel, "D" if sh == "Executive_Summary" else get_column_letter(min(maxc + 2, 40)), row + 1 if sh == "Executive_Summary" else 5, 980 if sh == "Executive_Summary" else 760)
    gallery(wb["Figure_Gallery"], figs, endrow["Figure_Gallery"] + 2)
    charts = [bar(wb, pos, "D02_by_sample_type.csv", "Sample_Type", "samples", "Included samples by sample type (native chart)", "Samples", "L5"),
              bar(wb, pos, "D04_by_country.csv", "Country", "samples", "Top 15 countries by included samples (native chart)", "Samples", "H5", n=15),
              bar(wb, pos, "A04_ds02_arg_load_by_city.csv", "City", "median_total_RPKM", "DS02 median total normalized ARG abundance (RPKM) by city", "RPKM", "H5", horizontal=False),
              bar(wb, pos, "A12_rf_arg_load_by_sample_type.csv", "Sample_Type", "median_total_FPKM", "ResFinder median total normalized ARG abundance (FPKM) by matrix", "FPKM", "H22"),
              bar(wb, pos, "DH11_dhaka_amr_sample_summary.csv", "Sample_ID", "total_FPKM_16S", "Dhaka: total normalized ARG abundance per sample (FPKM-16S)", "FPKM-16S", "P5")]
    for ws in wb.worksheets: ws.sheet_view.zoomScale = 90
    wb.properties.title = "AMR Updated Final Analysis"; wb.properties.creator = "Updated_Research_Analysis pipeline (Python/openpyxl)"
    wb.save(path)
    return dict(sheets=len(wb.sheetnames), tables=len(pos), figures_gallery=len(figs), charts=int(sum(charts)), positions=pos)

def build_dhaka(path):
    tabs = load_catalog("tables"); figs = load_catalog("figures").sort_values("number"); F = pd.read_csv(TAB / "R03u_updated_findings_register.csv"); L = pd.read_csv(TAB / "R04u_updated_limitations.csv")
    F[F.id.str.startswith("DF")].to_csv(TAB / "DH35_dhaka_findings.csv", index=False)
    L[L.domain.str.contains("Dhaka|reprocess|Cross-city|Pipeline integration|Taxonomic resolution", case=False)].to_csv(TAB / "DH36_dhaka_limitations.csv", index=False)
    wb = Workbook(); ws = wb.active; ws.title = "README"
    header(ws, "Dhaka City Advanced Analysis", "Dedicated workbook for the Dhaka-centred analysis (uniform reprocessing of 16 Dhaka metagenomes and South Asian comparators).", "notebooks 12-18 of Updated_Research_Analysis")
    lines = ["Contents: one sheet per analysis area (see sheet list). All values come from the CSV tables named above each table.",
             "Pipeline: first 4,000,000 read pairs per ENA run; fastp 1.4.0; Kraken 2.17.1 + SILVA 138 16S index + Bracken (genus); KMA 1.7.1 + ResFinder DB (resfinder-4.7.2).",
             "Normalized ARG abundance: FPKM-16S = aligned fragments / template kb / million prokaryotic 16S fragments; RPKM-depth = per million QC read pairs. Relative measures, not absolute ARG load.",
             "Detection: high-confidence = identity >= 90% and template coverage >= 60% (primary); any hit = unconfirmed upper bound.",
             "Inference: small groups use exact permutation tests; Dhaka vs Delhi uses sampling units (sites) as the independent unit; between-dataset Dhaka contrasts are descriptive (confounded).",
             "Gene detection is not phenotypic resistance, expression, viability, transmission or clinical infection; no risk hotspots are inferred."]
    for i, l in enumerate(lines): ws.cell(5 + i, 1, l)
    used, pos = set(), {}
    figmap = {31: "Dhaka_Sequencing_QC", 32: "Dhaka_Taxonomic_Abundance", 33: "Dhaka_Taxonomic_Abundance", 34: "Dhaka_Resistance_Classes", 35: "Dhaka_AMR_Genes", 36: "Dhaka_AMR_Genes",
              37: "Dhaka_Beta_Diversity", 38: "Dhaka_Alpha_Diversity", 39: "Dhaka_Comparative_Analysis", 40: "Dhaka_Comparative_Analysis", 41: "Dhaka_Comparative_Analysis",
              42: "Dhaka_Comparative_Analysis", 43: "Dhaka_Spatial_Analysis"}
    for sh in DHAKA_SHEETS:
        ws = wb.create_sheet(sh); header(ws, sh.replace("_", " "), DDESC.get(sh, ""), "source CSVs in Updated_Research_Analysis/04_Tables")
        if sh == "Dhaka_Results_Interpretation": files = ["DH35_dhaka_findings.csv"]
        elif sh == "Dhaka_Limitations": files = ["DH36_dhaka_limitations.csv"]
        elif sh == "Dhaka_Figure_Gallery": files = []
        else: files = sorted(tabs[tabs.sheet == sh].file)
        row, maxc = put_tables(wb, ws, files, tabs, used, pos); ws.freeze_panes = "A5"
        if sh != "Dhaka_Figure_Gallery":
            put_figures(ws, figs[figs.number.map(figmap) == sh], get_column_letter(min(maxc + 2, 40)), 5)
    gallery(wb["Dhaka_Figure_Gallery"], figs[figs.number >= 31], 5)
    charts = [bar(wb, pos, "DH11_dhaka_amr_sample_summary.csv", "Sample_ID", "total_FPKM_16S", "Total normalized ARG abundance per Dhaka sample (FPKM-16S)", "FPKM-16S", "B" + str(pos["DH11_dhaka_amr_sample_summary.csv"][3] + 4)),
              bar(wb, pos, "DH11_dhaka_amr_sample_summary.csv", "Sample_ID", "genes_high_confidence", "ARG genes detected at high confidence per Dhaka sample", "Genes", "L" + str(pos["DH11_dhaka_amr_sample_summary.csv"][3] + 4)),
              bar(wb, pos, "DH17_dhaka_alpha_diversity.csv", "Sample_ID", "rarefied_genus_richness", "Rarefied genus richness per Dhaka sample", "Genera", "B" + str(pos["DH17_dhaka_alpha_diversity.csv"][3] + 4))]
    for ws in wb.worksheets: ws.sheet_view.zoomScale = 90
    wb.properties.title = "Dhaka City Advanced Analysis"; wb.save(path)
    return dict(sheets=len(wb.sheetnames), tables=len(pos), figures=int((figs.number >= 31).sum()), charts=int(sum(charts)), positions=pos)

def build_issue_register(path, orig_register):
    src = pd.read_excel(orig_register, sheet_name="Issue_Register", header=2)
    ac = pd.read_csv(TAB / "AC00_audit_issue_resolution_register.csv").set_index("issue")
    out = src.copy()
    out["Status_before_update"] = out["Status"]
    out["Update_status_2026-10-31"] = [ac.loc[i, "status"] if i in ac.index else ("CORRECTED (independent audit; carried forward and re-verified in the updated run)" if "CORRECTED" in str(s) and "PARTIALLY" not in str(s) else s) for i, s in zip(out.ID, out.Status)]
    out["Update_correction_implemented"] = [ac.loc[i, "correction_implemented"] if i in ac.index else "Audit correction retained in the updated code (patch carried forward)" for i in out.ID]
    out["Update_before"] = [ac.loc[i, "before"] if i in ac.index else "" for i in out.ID]
    out["Update_after"] = [ac.loc[i, "after"] if i in ac.index else "" for i in out.ID]
    out["Update_verification_evidence"] = [ac.loc[i, "verification_evidence"] if i in ac.index else "Updated outputs regenerated from fresh kernels" for i in out.ID]
    sev = out.groupby("Severity").agg(issues=("ID", "size"), resolved_or_corrected=("Update_status_2026-10-31", lambda s: int((s.str.contains("RESOLVED|CORRECTED") & ~s.str.contains("PARTIAL")).sum())),
                                      partially_resolved=("Update_status_2026-10-31", lambda s: int(s.str.contains("PARTIAL").sum())),
                                      open=("Update_status_2026-10-31", lambda s: int((~s.str.contains("RESOLVED|CORRECTED")).sum()))).reindex(["Critical", "Major", "Moderate", "Minor"]).fillna(0).reset_index()
    wb = Workbook(); ws = wb.active; ws.title = "Summary"; used, pos = set(), {}
    header(ws, "Updated audit issue register - summary", "Status of all 20 issues of the Independent Research Audit after the 31 October 2026 update.", "Independent_Audit_Issue_Register.xlsx + notebook 11 (AC00)")
    write_table(ws, sev, 5, 1, "summary", used)
    ws2 = wb.create_sheet("Issue_Register"); header(ws2, "Issue register (updated)", "Original audit fields plus update status, correction, before/after values and evidence.", "notebook 11 (AC00-AC06)")
    write_table(ws2, out, 5, 1, "register", used)
    for f in ["AC01_aud04_sample_type_corrections.csv", "AC02_aud11_pseudoreplication_before_after.csv", "AC03_aud12_detection_tiers.csv", "AC04_aud17_model_leakage_before_after.csv", "AC05_aud15_figure_audit.csv", "AC06_aud19_terminology_scan.csv"]:
        w = wb.create_sheet(f.split("_")[0] + "_" + f.split("_")[1][:20]); header(w, f.replace(".csv", ""), "", f"tables/{f}")
        write_table(w, pd.read_csv(TAB / f), 5, 1, f.replace(".csv", ""), used)
    wb.save(path); return dict(sheets=len(wb.sheetnames), issues=len(out), open=int(sev.open.sum()))

def build_stat_verification(path):
    wb = Workbook(); ws = wb.active; ws.title = "Summary"; used = set()
    header(ws, "Updated statistical verification", "Independent recalculation of critical results (verify_updated.py, written separately from the analysis helpers), workbook checks and internal checks.", "verification/V04-V06")
    parts = [("Independent_Recalculation", VER / "V05_updated_independent_checks.csv"), ("Workbook_Checks", VER / "V06_workbook_checks.csv"),
             ("Raw_Read_Reproducibility", VER / "V07_raw_read_reproducibility_check.csv"), ("Separate_Agent_Recalc", VER / "V08_independent_agent_checks.csv"),
             ("Internal_Checks", TAB / "V04_updated_analysis_checks.csv"), ("AUD11_Before_After", TAB / "AC02_aud11_pseudoreplication_before_after.csv"),
             ("AUD17_Before_After", TAB / "AC04_aud17_model_leakage_before_after.csv")]
    summ = []
    for name, f in parts:
        if not f.exists(): continue
        df = pd.read_csv(f); w = wb.create_sheet(name); header(w, name.replace("_", " "), "", f"{f.parent.name}/{f.name}"); write_table(w, df, 5, 1, name, used)
        st = next((df[c] for c in ["status", "result", "verdict"] if c in df.columns), None)
        if name == "Raw_Read_Reproducibility": st = df["result"].str.startswith("IDENTICAL").map({True: "PASS", False: "FAIL"})
        summ.append(dict(sheet=name, rows=len(df), verified_or_pass=int(st.isin(["VERIFIED", "PASS", "AGREE"]).sum()) if st is not None else None,
                         mismatch_or_fail=int(st.isin(["MISMATCH", "FAIL", "DISAGREE"]).sum()) if st is not None else None,
                         note="Separate-agent verdicts are from its second pass; see status_after_final_correction for the 2 DISAGREE rows (both fixed)" if name == "Separate_Agent_Recalc" else ""))
    write_table(ws, pd.DataFrame(summ), 5, 1, "summary", used); wb.save(path); return summ

if __name__ == "__main__":
    OUTDIR.mkdir(parents=True, exist_ok=True)
    m = build_master(OUTDIR / "AMR_Updated_Final_Analysis.xlsx", extra_ver=verification_tables())
    d = build_dhaka(OUTDIR / "Dhaka_City_Advanced_Analysis.xlsx")
    orig = Path(r"C:\Users\USER\Documents\Research\AMR\Independent_Audit\Independent_Audit_Issue_Register.xlsx")
    r = build_issue_register(OUTDIR / "Updated_Audit_Issue_Register.xlsx", orig)
    sv = build_stat_verification(OUTDIR / "Updated_Statistical_Verification.xlsx")
    json.dump({"master": {k: v for k, v in m.items() if k != "positions"}, "dhaka": {k: v for k, v in d.items() if k != "positions"}, "register": r,
               "master_positions": {k: list(v) for k, v in m["positions"].items()}, "dhaka_positions": {k: list(v) for k, v in d["positions"].items()}},
              open(OUTDIR / "workbook_build_manifest.json", "w"), indent=0)
    print(json.dumps({"master": {k: v for k, v in m.items() if k != "positions"}, "dhaka": {k: v for k, v in d.items() if k != "positions"}, "register": r}))
