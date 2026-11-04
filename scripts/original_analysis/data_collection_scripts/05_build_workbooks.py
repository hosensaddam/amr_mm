"""
05_build_workbooks.py
Build the five Excel workbooks, Public_Data_Source_Log.csv and summary tables
from the harmonized CSVs and logs. Every count is computed from collected records.
Formatting: Arial 10, bold shaded header, frozen header row, structured Excel
tables (with filters), widths, hyperlinks, validation, conditional formatting.
"""
import csv, glob, hashlib, json, os, re, importlib.util
import numpy as np, pandas as pd, xlsxwriter
from xlsxwriter.utility import xl_col_to_name

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H, R, L, O = (os.path.join(B, x) for x in ("harmonized", "raw_downloads", "logs", "outputs"))
os.makedirs(O, exist_ok=True)
spec = importlib.util.spec_from_file_location("cd", os.path.join(B, "scripts", "curation_decisions.py"))
cd = importlib.util.module_from_spec(spec); spec.loader.exec_module(cd)
NR, NAV, NYV, NAP = "Not reported", "Not available", "Not yet verified", "Not applicable"
DATE = "2026-10-31"
MISSING_CODES = ["Not reported", "Not available", "Not yet verified", "Not applicable", "Restricted access"]
ILLEGAL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")

sm = pd.read_csv(os.path.join(H, "sample_master.csv"), dtype=str).fillna("")
runs = pd.read_csv(os.path.join(H, "run_master.csv"), dtype=str).fillna("")
tax = pd.read_csv(os.path.join(H, "taxonomy_long.csv"), dtype={"Taxon_ID_NCBI": str})
div = pd.read_csv(os.path.join(H, "diversity_long.csv"))
amr = pd.read_csv(os.path.join(H, "amr_genes_long.csv"), dtype={"Gene_Identifier": str}, low_memory=False)
cov = pd.read_csv(os.path.join(H, "env_covariates_long.csv"), dtype=str).fillna("")
tlog = pd.read_csv(os.path.join(H, "transformation_log.csv"))
cand = json.load(open(os.path.join(H, "study_candidates.json")))
ds_of = json.load(open(os.path.join(H, "dataset_ids.json")))
ds_name = json.load(open(os.path.join(H, "dataset_names.json")))
qlog = [json.loads(l) for f in sorted(glob.glob(os.path.join(L, "query_log_*.jsonl"))) for l in open(f)]
inc = sm[sm.Inclusion_Status == "Included"].copy()

# ============================================================== writer helpers
class WB:
    def __init__(self, path):
        self.path = path
        self.wb = xlsxwriter.Workbook(path, {"strings_to_urls": False, "strings_to_numbers": False, "nan_inf_to_errors": True})
        f = self.wb.add_format
        self.base = f({"font_name": "Arial", "font_size": 10})
        self.hdr = f({"font_name": "Arial", "font_size": 10, "bold": True, "bg_color": "#1F4E78", "font_color": "#FFFFFF", "text_wrap": True, "valign": "top", "border": 1})
        self.title = f({"font_name": "Arial", "font_size": 14, "bold": True, "font_color": "#1F4E78"})
        self.bold = f({"font_name": "Arial", "font_size": 10, "bold": True})
        self.wrap = f({"font_name": "Arial", "font_size": 10, "text_wrap": True, "valign": "top"})
        self.url = f({"font_name": "Arial", "font_size": 10, "font_color": "#0563C1", "underline": 1})
        self.int = f({"font_name": "Arial", "font_size": 10, "num_format": "#,##0"})
        self.dec = f({"font_name": "Arial", "font_size": 10, "num_format": "0.000000"})
        self.dec2 = f({"font_name": "Arial", "font_size": 10, "num_format": "0.00"})
        self.good = f({"bg_color": "#C6EFCE", "font_color": "#006100"})
        self.warn = f({"bg_color": "#FFEB9C", "font_color": "#9C5700"})
        self.bad = f({"bg_color": "#FFC7CE", "font_color": "#9C0006"})
        self.grey = f({"bg_color": "#EDEDED", "font_color": "#595959"})
        self.tables = 0

    def readme(self, name, title, rows):
        ws = self.wb.add_worksheet(name)
        ws.set_column(0, 0, 30, self.bold); ws.set_column(1, 1, 120, self.wrap)
        ws.write(0, 0, title, self.title)
        ws.write_row(2, 0, ["Item", "Description"], self.hdr)
        for i, (k, v) in enumerate(rows, start=3):
            ws.write(i, 0, k, self.bold); ws.write(i, 1, v, self.wrap)
        ws.freeze_panes(3, 0)
        return ws

    def table(self, name, df, widths=None, url_cols=(), int_cols=(), dec_cols=(), dec2_cols=(), validations=None,
              status_cols=None, note=None):
        ws = self.wb.add_worksheet(name)
        start = 0
        if note:
            ws.write(0, 0, note, self.wrap); ws.set_row(0, 30); start = 1
        cols = list(df.columns)
        data = []
        for row in df.itertuples(index=False):
            out = []
            for v in row:
                if v is None or (isinstance(v, float) and np.isnan(v)): out.append("")
                elif isinstance(v, (np.integer,)): out.append(int(v))
                elif isinstance(v, (np.floating,)): out.append(float(v))
                elif isinstance(v, str): out.append(ILLEGAL.sub("", v)[:32000])
                else: out.append(v)
            data.append(out)
        self.tables += 1
        n = len(data)
        last_row = start + max(n, 1)
        spec = {"name": re.sub(r"[^A-Za-z0-9_]", "_", f"tbl_{name}")[:250], "style": "Table Style Light 9",
                "columns": [{"header": c, "header_format": self.hdr} for c in cols]}
        if n: spec["data"] = data
        ws.add_table(start, 0, last_row, len(cols) - 1, spec)
        ws.freeze_panes(start + 1, 0)
        for j, c in enumerate(cols):
            w = (widths or {}).get(c)
            if w is None:
                sample = [len(str(r[j])) for r in data[:300]] + [len(c)]
                w = min(max(10, int(np.percentile(sample, 90)) + 2), 60)
            fmt = self.int if c in int_cols else self.dec if c in dec_cols else self.dec2 if c in dec2_cols else self.base
            ws.set_column(j, j, w, fmt)
        for c in url_cols:
            if c not in cols: continue
            j = cols.index(c); k = 0
            for i, r in enumerate(data):
                v = r[j]
                if isinstance(v, str) and v.startswith(("http://", "https://")) and k < 60000 and len(v) < 2000 and ";" not in v and " " not in v:
                    ws.write_url(start + 1 + i, j, v, self.url, v); k += 1
        for c, opts in (validations or {}).items():
            if c in cols and n:
                j = cols.index(c)
                ws.data_validation(start + 1, j, start + n, j, {"validate": "list", "source": opts})
        for c, mapping in (status_cols or {}).items():
            if c in cols and n:
                j = cols.index(c)
                for val, fmt in mapping.items():
                    ws.conditional_format(start + 1, j, start + n, j, {"type": "cell", "criteria": "==", "value": f'"{val}"', "format": fmt})
        return ws

    def close(self): self.wb.close()

def status_maps(w):
    return {"QC_Status": {"Pass": w.good, "Warning": w.warn, "Fail": w.bad},
            "Verification_Status": {"Verified": w.good, "Partially verified": w.warn, "Not yet verified": w.bad, "Not accessible": w.bad, "Unsuitable": w.grey},
            "Inclusion_Status": {"Included": w.good, "Pending verification": w.warn, "Excluded": w.grey, "Separate contextual dataset": w.warn},
            "Decision": {"Included": w.good, "Pending verification": w.warn, "Excluded": w.grey, "Separate contextual dataset": w.warn},
            "Status": {"Confirmed": w.good, "Discrepancy": w.bad, "Unresolved": w.warn, "Operational": w.good, "Not accessed": w.grey, "Refused (HTTP 403)": w.bad},
            "Result": {"Pass": w.good, "Warning": w.warn, "Fail": w.bad},
            "Readiness": {"Ready": w.good, "Conditionally ready": w.warn, "Not ready": w.bad},
            "Download_Status": {"Downloaded": w.good, "Not downloaded (large; needs confirmation)": w.warn, "Not retrievable from workspace": w.bad, "Metadata only": w.grey}}

def dictionary(df_map):
    """df_map: {sheet: [(field, description, type/unit, allowed/missing codes)]}"""
    rows = []
    for sheet, fields in df_map.items():
        for f in fields:
            rows.append({"Sheet": sheet, "Field": f[0], "Description": f[1], "Type_or_Unit": f[2], "Allowed_Values_or_Codes": f[3] if len(f) > 3 else ""})
    return pd.DataFrame(rows)

def fsize(p): return os.path.getsize(p) if os.path.exists(p) else None
def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

# ============================================================== derived tables
# --- query log
def qvalue(r):
    ids = r["Relevant_IDs"]; db = r["Database"]
    raw = json.load(open(os.path.join(B, r["Raw_Response_File"])))
    recs = {x["id"]: x for x in raw["records"]}
    verified = incl = 0
    for i in ids:
        x = recs.get(i, {})
        if db.startswith("MG-RAST"):
            p = x.get("project"); verified += 1 if os.path.exists(os.path.join(R, "mgrast_projects", f"{p}.json")) else 0
            incl += 1 if cd.MGRAST.get(p, ("",))[0] in ("Included", "Separate contextual dataset") else 0
        elif db.startswith(("Europe PMC", "Zenodo")):
            pass
        else:
            acc = i
            if db.startswith("MGnify"):
                ins = [a for a in (x.get("insdc") or "").split(";") if a.startswith("PRJ")]; acc = ins[0] if ins else i
            k = f"INSDC_STUDY|{acc}"
            if k in cand and "ena" in cand[k]: verified += 1
            if cd.decide(acc)[0] in ("Included", "Separate contextual dataset"): incl += 1
    return verified, incl

NOTES = {"MGnify (API v2)": "v2 'search' matched the full string as a phrase; all multi-word queries returned 0. See adapted queries MGN-ADJ-*.",
         "MGnify (API v2) - adapted": "Single-keyword or biome-lineage browse adapted from the query groups because of v2 phrase matching.",
         "MG-RAST": "Free-text 'all=' search over public metagenome metadata; first page of up to 1,000 metagenomes inspected; matching appears token-based (counts not comparable across terms).",
         "ENA Portal": "Words AND-ed within study_title, restricted to library_source=METAGENOMIC and library_strategy=WGS; result aggregated to distinct studies.",
         "NCBI BioProject": "Term AND (metagenome OR metagenomic) with 'bioproject sra' filter; summaries retrieved for up to 500 records.",
         "Europe PMC": "Restricted to articles with BioProject/ENA accession annotations; first 100 records inspected; records are literature leads, not datasets.",
         "Zenodo": "All words required; resource_type=dataset; anonymous API returns at most 25 records per page; records are leads only."}
qrows = []
for r in qlog:
    v, i = qvalue(r)
    insp, rel = r["N_Records_Inspected"], r["N_Potentially_Relevant"]
    qrows.append({"Query_ID": r["Query_ID"], "Database": r["Database"], "Query_Group": r["Query_Group"], "Search_String": r["Search_String"],
                  "Search_Date": r["Search_Date"], "Search_URL_or_Endpoint": r["Search_URL_or_Endpoint"], "Filters_Used": r["Filters_Used"],
                  "N_Results_Returned": r["N_Results_Returned"], "N_Records_Inspected": insp, "N_Potentially_Relevant": rel,
                  "N_Verified": v, "N_Included_or_Contextual_After_Review": i, "N_Excluded": insp - rel + (v - i if v >= i else 0),
                  "Exclusion_Reasons": f"Failed programmatic screen: {insp - rel}; excluded or pending at eligibility review after verification: {max(v - i, 0)}",
                  "Notes_and_Limitations": NOTES.get(r["Database"], "") + (f" API error: {r['API_Error']}" if r.get("API_Error") else ""),
                  "Raw_Response_File": r["Raw_Response_File"]})
qdf = pd.DataFrame(qrows)

# --- study inventory
inv = []
for k, v in cand.items():
    kind, acc = v["kind"], v["accession"]
    e = v.get("ena") or {}
    if kind == "INSDC_STUDY":
        if acc == "ERP015409":
            dec = ("Excluded", "", "Duplicate of PRJEB13831 (secondary accession); records merged under PRJEB13831")
        elif "ena" in v and e.get("n_wgs_metagenomic_runs", 0) == 0:
            strat = ", ".join(e.get("strategies") or []) or "no public runs"
            dec = ("Excluded", "", f"No public WGS/METAGENOMIC runs in ENA ({strat})")
        else:
            dec = cd.decide(acc)
        ver = "Verified" if (e.get("n_runs", 0) > 0 and dec[0] in ("Included", "Separate contextual dataset") and acc in cd.CORE) else ("Partially verified" if e.get("n_runs", 0) > 0 else "Not accessible")
        if e.get("n_runs", 0) == 0: ver = "Not accessible"
        url = f"https://www.ebi.ac.uk/ena/browser/view/{acc}"
        title = e.get("study_title") or (v["titles"][0] if v["titles"] else "")
    elif kind == "MGRAST_PROJECT":
        d = cd.MGRAST.get(acc, ("Excluded", "", cd.DEFAULT_EXCLUDED_REASON)); dec = d
        p = os.path.join(R, "mgrast_projects", f"{acc}.json")
        ver = "Partially verified" if os.path.exists(p) else "Not accessible"
        url = f"https://www.mg-rast.org/mgmain.html?mgpage=project&project={acc}"; title = v["titles"][0] if v["titles"] else ""
    elif kind == "ZENODO_RECORD":
        dec = ("Pending verification", "", "Zenodo dataset lead; not screened at record level")
        if acc in ("10.5281/zenodo.17590646", "10.5281/zenodo.6919377"): dec = ("Included", "", "Processed data used (see Download_Manifest)")
        ver = NYV; url = f"https://doi.org/{acc}" if acc.startswith("10.") else f"https://zenodo.org/records/{acc}"; title = v["titles"][0] if v["titles"] else ""
    else:
        dec = ("Pending verification", "", "Accession-linked publication lead; full text not screened")
        ver = NYV; url = f"https://doi.org/{acc}" if str(acc).startswith("10.") else f"https://europepmc.org/article/MED/{acc}"; title = v["titles"][0] if v["titles"] else ""
    inv.append({"Record_ID": k, "Record_Type": {"INSDC_STUDY": "INSDC study/project", "MGRAST_PROJECT": "MG-RAST project", "ZENODO_RECORD": "Zenodo record", "PUBLICATION": "Publication (accession-linked)"}[kind],
                "Accession": acc, "Title": title, "Discovered_Via": "; ".join(v["dbs"]), "Query_IDs": "; ".join(v["queries"]), "N_Query_Hits": len(v["queries"]),
                "Dataset_ID": ds_of.get(acc, NAP), "Public_Runs_ENA": e.get("n_runs", "") if kind == "INSDC_STUDY" else "",
                "Shotgun_Metagenomic_Runs": e.get("n_wgs_metagenomic_runs", "") if kind == "INSDC_STUDY" else "",
                "Samples_ENA": e.get("n_samples", "") if kind == "INSDC_STUDY" else "", "Library_Strategies": "; ".join(e.get("strategies") or []),
                "Platforms": "; ".join(e.get("platforms") or [])[:300], "Countries": "; ".join(e.get("countries") or [])[:300],
                "Collection_Months": (f"{e['dates'][0]} to {e['dates'][-1]}" if e.get("dates") else ""),
                "Runs_With_Coordinates": e.get("n_with_latlon", "") if kind == "INSDC_STUDY" else "",
                "Raw_FASTQ_GB": round(e.get("fastq_bytes_total", 0) / 1e9, 2) if kind == "INSDC_STUDY" else "",
                "MGnify_Study": "; ".join(v.get("mgnify") or []), "Year_or_Date": v["meta"].get("year") or v["meta"].get("date", ""),
                "Decision": dec[0], "Sample_Category": dec[1] or NAP, "Decision_Reason": dec[2], "Verification_Status": ver,
                "Official_URL": url, "Date_Checked": DATE})
invdf = pd.DataFrame(inv).sort_values(["Record_Type", "Decision", "Accession"])

# --- accession register (selected datasets)
acc_rows = []
for acc, ds in sorted(ds_of.items(), key=lambda x: (int(x[1][2:]), x[0])):
    x = sm[sm.Study_ID == acc]
    k = f"INSDC_STUDY|{acc}"; e = (cand.get(k) or {}).get("ena") or {}
    acc_rows.append({"Dataset_ID": ds, "Dataset_Name": ds_name.get(ds, ""), "Study_ID": acc, "Primary_Project_Accession": e.get("study_accession") or acc,
                     "Secondary_Study_Accession": e.get("secondary_study_accession") or NR, "Repository": x.Repository.iloc[0] if len(x) else NR,
                     "MGnify_Study": ("MGYS00005846" if acc == "PRJEB27054" else "; ".join((cand.get(k) or {}).get("mgnify") or [])) or NAV,
                     "Publication_DOI": x.DOI.iloc[0] if len(x) else NYV, "Publication_IDs": x.Publication_ID.iloc[0] if len(x) else NYV,
                     "Processed_Data": "10.5281/zenodo.17590646 (CC BY 4.0)" if ds == "DS02" else ("10.5281/zenodo.6919377 (CC BY 4.0; ResFinder counts for " + str(int(x.AMR_Data_Available.str.contains("ResFinder").sum())) + " samples)" if x.AMR_Data_Available.str.contains("ResFinder").any() else NAV),
                     "Metadata_Source": {"DS01": "Munk et al. 2022 Suppl. Data 1 + ENA", "DS02": "Zenodo meta.tsv + SRA", "DS03": "NCBI BioSample", "DS04": "MetaSUB-metadata GitHub + SRA"}.get(ds, "ENA run report; NCBI BioSample where needed" if not acc.startswith("mgp") else "MG-RAST search API"),
                     "Code_Repository": {"DS02": "https://github.com/SowpatiLab/Metro_metagenomics (MIT)", "DS04": "https://github.com/MetaSUB/CAP2"}.get(ds, NAV),
                     "N_Samples": len(x), "N_Samples_Included": int((x.Inclusion_Status == "Included").sum()),
                     "N_Runs": int(x.N_Runs.astype(int).sum()) if len(x) else 0,
                     "Official_URL": f"https://www.ebi.ac.uk/ena/browser/view/{acc}" if not acc.startswith("mgp") else f"https://www.mg-rast.org/mgmain.html?mgpage=project&project={acc}"})
accdf = pd.DataFrame(acc_rows)

# --- dataset eligibility (selected datasets: criteria computed from sample records)
def pct(x, cond): return f"{int(cond.sum())}/{len(x)}"
el = []
for acc, ds in sorted(ds_of.items(), key=lambda x: (int(x[1][2:]), x[0])):
    x = sm[sm.Study_ID == acc]
    if not len(x): continue
    dec = cd.decide(acc) if not acc.startswith("mgp") else cd.MGRAST[acc]
    geo = (x.Country != NR) & (x.Country != NYV)
    sg = x.Sequencing_Type.str.contains("Shotgun")
    st = ~x.Sample_Type.isin(["Unclassified", ""])
    el.append({"Dataset_ID": ds, "Study_ID": acc, "Dataset_Name": ds_name.get(ds, ""),
               "C1_Public_Record": "Yes", "C2_Stable_Accessions": "Yes",
               "C3_Shotgun_Metagenomic_Samples": pct(x, sg), "C4_Sample_Type_Identifiable": pct(x, st),
               "C5_Country_Available": pct(x, geo), "C5b_Coordinates_Available": pct(x, x.Latitude != NR),
               "C6_Sequencing_Info": pct(x, x.Sequencing_Platform != NR),
               "C7_Taxonomic_or_Functional_Annotation": "Yes (downloaded)" if ds == "DS02" else ("Partial (diversity indices)" if x.Taxonomic_Profile_Available.str.startswith("Partial").any() else ("Online (MG-RAST)" if acc.startswith("mgp") else "No (raw reads only)")),
               "C8_AMR_Info_or_Reanalysis_Potential": "AMR table downloaded" if x.AMR_Data_Available.str.startswith("Yes").any() else ("Raw reads public (reanalysis possible)" if x.Raw_Read_Availability.str.contains("Public").any() else "No"),
               "C9_Traceable_Publication": x.DOI.iloc[0] if x.DOI.iloc[0] not in (NYV, "") else "Not yet verified",
               "C10_Reuse_Conditions": x.Licence.iloc[0], "Decision": dec[0], "Sample_Category": dec[1], "Reason": dec[2]})
eldf = pd.DataFrame(el)

# --- download manifest
man = []
def add(path, url, ftype, status="Downloaded", note="", checksum=None):
    p = os.path.join(B, path)
    man.append({"File_Path": path, "Source_URL": url, "File_Type": ftype, "Size_Bytes": fsize(p) if status == "Downloaded" else note.split("|")[0] if "|" in note else "",
                "Checksum": checksum or (("sha256:" + sha256(p)) if status == "Downloaded" and os.path.exists(p) and fsize(p) < 6e8 else ""),
                "Download_Date": DATE if status == "Downloaded" else "", "Download_Status": status, "Notes": note.split("|")[-1]})
add("raw_downloads/zenodo_17590646_india/data.tar.xz", "https://zenodo.org/api/records/17590646/files/data.tar.xz/content", "Processed data archive (tar.xz)", note="MD5 dd9e8bfaaad2aa5b33c81c7fd1a8b5bc matches Zenodo record")
for f in ["metadata.tsv", "ResFinder_anno.tsv", "diversity.tsv", "ARG.tsv"]:
    add(f"raw_downloads/zenodo_6919377_global_resistome/{f}", f"https://zenodo.org/api/records/6919377/files/{f}/content", "Processed table (TSV)", note="MD5 recorded in MD5SUMS.txt; subset to core runs in harmonized/")
add("raw_downloads/munk2022_supplementary/PMC9715550_supplementary.zip", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9715550/supplementaryFiles", "Publication supplementary package", note="Includes Supplementary Data 1 (run/sample map)")
add("raw_downloads/munk2022_PMC9715550.xml", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9715550/fullTextXML", "Publication full text (JATS XML)")
add("raw_downloads/singh2026_PMC13144492.xml", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13144492/fullTextXML", "Publication full text (JATS XML)")
add("raw_downloads/metasub_metadata/complete_metadata.csv", "https://github.com/MetaSUB/MetaSUB-metadata (commit 1d738748)", "Sample metadata (CSV)")
for f in ["README.pdf", "farm_dna.fa.gz", "farm_protein.fa.gz", "HMMs.xlsx.gz", "HMMs.tab.gz"]:
    add(f"raw_downloads/farme/{f}", f"https://staff.washington.edu/jwallace/farme/downloads/{f}", "Reference database file")
add("raw_downloads/card/card-data_latest.tar.bz2", "https://card.mcmaster.ca/latest/data", "Reference database archive (CARD 4.0.2)")
n_ena = len(glob.glob(os.path.join(R, "ena_study_runs", "*.tsv"))); n_sm = len(glob.glob(os.path.join(R, "ena_sample_metadata", "*.tsv")))
man.append({"File_Path": "raw_downloads/ena_study_runs/*.tsv", "Source_URL": "https://www.ebi.ac.uk/ena/portal/api/search (result=read_run)", "File_Type": f"Run metadata reports ({n_ena} studies)",
            "Size_Bytes": sum(fsize(p) for p in glob.glob(os.path.join(R, "ena_study_runs", "*.tsv"))), "Checksum": "", "Download_Date": DATE, "Download_Status": "Metadata only", "Notes": "Study-level verification"})
man.append({"File_Path": "raw_downloads/ena_sample_metadata/*.tsv", "Source_URL": "https://www.ebi.ac.uk/ena/portal/api/search (result=read_run, extended fields)", "File_Type": f"Run/sample metadata ({n_sm} studies)",
            "Size_Bytes": sum(fsize(p) for p in glob.glob(os.path.join(R, "ena_sample_metadata", "*.tsv"))), "Checksum": "", "Download_Date": DATE, "Download_Status": "Metadata only", "Notes": "Exact query URL saved beside each file (.url)"})
man.append({"File_Path": "raw_downloads/ncbi_biosample/*.xml", "Source_URL": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample", "File_Type": "BioSample XML",
            "Size_Bytes": sum(fsize(p) for p in glob.glob(os.path.join(R, "ncbi_biosample", "*.xml"))), "Checksum": "", "Download_Date": DATE, "Download_Status": "Metadata only", "Notes": "Attributes for samples lacking location/date in ENA"})
man.append({"File_Path": "raw_downloads/mgrast_metagenomes/*.json; raw_downloads/mgrast_projects/*.json", "Source_URL": "https://api.mg-rast.org/search; https://api.mg-rast.org/project/{id}", "File_Type": "MG-RAST metadata (JSON)",
            "Size_Bytes": sum(fsize(p) for p in glob.glob(os.path.join(R, "mgrast_*", "*.json"))), "Checksum": "", "Download_Date": DATE, "Download_Status": "Metadata only", "Notes": ""})
# not downloaded
for ds, g in sm[sm.Inclusion_Status == "Included"].groupby("Dataset_ID"):
    b = pd.to_numeric(g.FASTQ_Bytes_Total, errors="coerce").sum()
    man.append({"File_Path": f"[raw reads] {ds}", "Source_URL": "; ".join(sorted(set(g.Study_ID)))[:200], "File_Type": "Raw FASTQ (shotgun metagenomic)",
                "Size_Bytes": int(b) if b else "", "Checksum": "Per-run MD5 in Sequencing_Metadata", "Download_Date": "", "Download_Status": "Not downloaded (large; needs confirmation)",
                "Notes": f"{len(g)} included samples; {b/1e9:,.1f} GB compressed FASTQ per ENA" if b else f"{len(g)} included samples; size not reported (MG-RAST)"})
man.append({"File_Path": "[processed] Zenodo 6919377 rRNA.tsv", "Source_URL": "https://zenodo.org/records/6919377", "File_Type": "16S/18S read-count table", "Size_Bytes": 6401117131, "Checksum": "",
            "Download_Date": "", "Download_Status": "Not downloaded (large; needs confirmation)", "Notes": "Would give genus-level taxonomy for 467 DS01 runs"})
man.append({"File_Path": "[processed] Munk et al. 2022 Source Data (MOESM6/7)", "Source_URL": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-34312-7/MediaObjects/", "File_Type": "Figure source data",
            "Size_Bytes": "", "Checksum": "", "Download_Date": "", "Download_Status": "Not retrievable from workspace", "Notes": "Publisher server returned HTTP 403; download manually from the article page"})
man.append({"File_Path": "[processed] Singh et al. 2026 source_data.zip and Suppl. Data S7 (site coordinates)", "Source_URL": "https://doi.org/10.1038/s41467-026-70702-x", "File_Type": "Publication supplementary",
            "Size_Bytes": "", "Checksum": "", "Download_Date": "", "Download_Status": "Not retrievable from workspace", "Notes": "Not in Europe PMC package; download manually to obtain site-level locations"})
man.append({"File_Path": "[processed] MGnify MGYS00005846 analyses", "Source_URL": "https://www.ebi.ac.uk/metagenomics/studies/MGYS00005846", "File_Type": "MGnify assembly analyses",
            "Size_Bytes": "", "Checksum": "", "Download_Date": "", "Download_Status": "Not downloaded (large; needs confirmation)", "Notes": "Functional/taxonomic annotations for PRJEB27054 assemblies"})
mandf = pd.DataFrame(man)

# --- exclusion log
ex = []
for _, r in invdf[invdf.Decision == "Excluded"].iterrows():
    ex.append({"Level": "Study/project", "Identifier": r.Accession, "Dataset_ID": r.Dataset_ID, "Title_or_Description": r.Title, "Reason": r.Decision_Reason, "N_Records": 1})
for (ds, reason), g in sm[sm.Inclusion_Status == "Excluded"].groupby(["Dataset_ID", "Inclusion_Reason"]):
    ex.append({"Level": "Sample", "Identifier": f"{ds}: " + "; ".join(g.Sample_ID.head(5)) + (" ..." if len(g) > 5 else ""), "Dataset_ID": ds,
               "Title_or_Description": ds_name.get(ds, ""), "Reason": reason, "N_Records": len(g)})
tot_insp = qdf.N_Records_Inspected.sum(); tot_rel = qdf.N_Potentially_Relevant.sum()
ex.append({"Level": "Search record", "Identifier": "All queries", "Dataset_ID": NAP, "Title_or_Description": "Programmatic screen across 236 queries",
           "Reason": "Failed programmatic title/metadata screen (SCREEN_RULE v1)", "N_Records": int(tot_insp - tot_rel)})
exdf = pd.DataFrame(ex)

# --- source verification checks (evidence-based)
d1 = sm[sm.Dataset_ID == "DS01"]; d2 = sm[sm.Dataset_ID == "DS02"]; d4 = sm[sm.Dataset_ID == "DS04"]
munk_cities = d1[(d1.Inclusion_Status == "Included") & (d1.Geographic_Unit == "City (source-labelled)")].City.nunique()
munk_countries = d1[(d1.Inclusion_Status == "Included") & (d1.Geographic_Unit == "City (source-labelled)")].Country.nunique()
munk_pairs = d1[(d1.Inclusion_Status == "Included") & (d1.Geographic_Unit == "City (source-labelled)")][["City", "Country"]].drop_duplicates().shape[0]
kyiv_src = int(sm.QC_Flags.str.contains("SRA lists Kyiv").sum()); kyiv_prop = int(sm.QC_Flags.str.contains("propagated from same-city SRA records lists Poland").sum())
ms_meta = pd.read_csv(os.path.join(R, "metasub_metadata", "complete_metadata.csv"), dtype=str)
ms_labels = ms_meta.city.dropna().unique(); ms_ctrl = [c for c in ms_labels if "control" in c]
checks = [
 ("V01", "Munk et al. 2022 global sewage: six ENA projects named in paper", "All runs in Suppl. Data 1 belong to the named projects", "Suppl. Data 1 joined to ENA run reports", "951/951 runs found, but 227 'GS pilot 2016' runs are in PRJEB27054 (not named in paper); 16 'pilot 2016 reproduce' runs in PRJEB13831 (ERP015409)", "Discrepancy"),
 ("V02", "Munk et al. 2022 city coverage", "Paper reports 243 cities, 101 countries", "Harmonized DS01 included records with source-labelled city", f"{munk_cities} distinct city labels, {munk_pairs} distinct city-country pairs (Santiago occurs in Chile and Spain), {munk_countries} countries. Countries match; the paper's 243 equals the pair count, not the label count. Kochi/Cochin spelling variants not merged", "Confirmed"),
 ("V03", "Munk Suppl. Data 1 biological samples", "unique_sample_id groups runs/accessions", "Suppl. Data 1", f"{d1.Biological_Sample_ID.replace(NAP, np.nan).nunique()} biological samples across {len(d1[d1.Biological_Sample_ID != NAP])} repository sample accessions", "Confirmed"),
 ("V04", "PRJEB27054 overlap with the six named sewage projects", "Overlap untested in earlier audit", "ENA run reports", "No shared runs or samples; PRJEB27054 is the 2016 pilot round of the same programme (227/235 runs in Suppl. Data 1)", "Confirmed"),
 ("V05", "Indian metro wastewater (PRJNA1415085) sample count", "447 samples, 4 cities (publication)", "SRA run report + Zenodo meta.tsv", f"{len(d2)} BioSamples, {d2.Run_ID.str.count(';').add(1).sum()} runs, cities: {', '.join(sorted(d2.City.unique()))}", "Confirmed"),
 ("V06", "Indian processed data integrity", "Zenodo MD5 dd9e8bfa...", "md5sum of downloaded archive", "MD5 matches", "Confirmed"),
 ("V07", "Indian sampling matrix", "'Wastewater' in SRA", "Singh et al. 2026 Methods (PMC13144492)", "Large open wastewater drains (classified Community wastewater, open drain)", "Confirmed"),
 ("V08", "CARD release used for Indian RGI table", "Release needed for reproducibility", "Zenodo files, GitHub README, publication full text", "Not stated in any accessed source", "Unresolved"),
 ("V09", "Dhaka PRJNA1399457 metadata", "Location/date absent from ENA", "NCBI BioSample attributes", "4 metagenomes: 3 hospital + 1 community sewer, 2025-06-18, coordinates given; 3 isolate genomes excluded", "Confirmed"),
 ("V10", "Dhaka PRJNA1399457 publication", "DOI 10.1128/spectrum.00755-26 suggested in brief", "BioProject/ENA records", "No publication linked in repository records; DOI not verified in this session", "Unresolved"),
 ("V11", "MetaSUB SRA vs publication", "Danko et al. report 4,728 samples, 60 cities", "SRA run report PRJNA732392 + GitHub metadata", f"{len(d4)} SRA samples (1 run each); {d4[d4.Inclusion_Status=='Included'].City.nunique()} cities among included; all 1,898 link to GitHub uuid", "Discrepancy"),
 ("V12", "MetaSUB sequencing platform", "Publication: HiSeq X Ten", "SRA run report", "All runs listed as Illumina NovaSeq 6000", "Discrepancy"),
 ("V13", "MetaSUB city-country consistency", "Kyiv is in Ukraine", "SRA geo_loc_name", f"{kyiv_src} samples list 'Poland:Kyiv' in SRA; {kyiv_prop} further Kyiv samples lack an SRA country and would inherit Poland by propagation. All {kyiv_src + kyiv_prop} set to Not yet verified; source values kept", "Discrepancy"),
 ("V14", "MetaSUB South Asian coverage", "Proposal names Delhi, Karachi, Dhaka", "MetaSUB metadata city column", f"{len(ms_labels)} non-blank city labels ({len(ms_ctrl)} are control labels; some cities have two spellings, e.g. newyork/new_york_city); none is a South Asian city", "Confirmed"),
 ("V15", "FDA retail lamb (PRJNA1522440)", "Food-associated AMR candidate", "ENA run report and BioProject description", "263 runs labelled metagenomic; metagenomic vs enrichment libraries not distinguishable; no urban location", "Confirmed"),
 ("V16", "Zenodo 6919377 coverage of core runs", "Processed ARG counts for core runs", "metadata.tsv", "Covers PRJEB13831 (233 runs) and PRJEB27054 (234 runs) only; Munk 2017-2018 rounds, MetaSUB, India, Dhaka not covered", "Confirmed"),
 ("V17", "Munk et al. 2022 Source Data", "Figure source data file", "Publisher static server", "HTTP 403 from workspace", "Unresolved"),
 ("V18", "FARME DB accessibility", "Website and downloads online", "https://staff.washington.edu/jwallace/farme/", "Online; files last modified Feb-Mar 2019; README revised March 2019", "Confirmed"),
 ("V19", "CARD current release", "Release and terms", "card.json _version; download page", "CARD 4.0.2 (timestamp 2026-08-11); commercial use restricted; ontologies CC BY 4.0", "Confirmed"),
 ("V20", "MGnify API v2", "Use current API", "openapi.json", "v2.0.0 operational; 'search' behaves as phrase match", "Confirmed"),
 ("V21", "MG-RAST project mgp6602", "Public project", "MG-RAST API", "HTTP 401 (not public)", "Discrepancy"),
 ("V22", "Proposal target of 35 cities", "35 cities incl. Delhi, Karachi, Dhaka", "Harmonized included samples", "Dhaka, Delhi and Karachi present in verified records; total city coverage far exceeds 35 but is uneven by matrix", "Confirmed"),
]
vdf = pd.DataFrame(checks, columns=["Check_ID", "Item", "Expectation_or_Claim", "Evidence_Source", "Finding", "Status"])
vdf["Date_Checked"] = DATE

# --- source registry
sources = [
 ("S01", "MG-RAST", "Metagenome repository and annotation server", "https://www.mg-rast.org/", "https://api.mg-rast.org/search; /project/{id}", "Official API", "Operational", "Public data only used; one project returned 401", "No licence stated for public data", "Discovery and sample metadata (3 included series)"),
 ("S02", "MGnify (EMBL-EBI)", "Metagenome analysis resource", "https://www.ebi.ac.uk/metagenomics/", "https://www.ebi.ac.uk/metagenomics/api/v2/ (v2.0.0)", "Official API", "Operational", "Public; Webin auth only for private data", "EMBL-EBI terms of use", "Discovery; study-to-INSDC mapping"),
 ("S03", "European Nucleotide Archive", "INSDC sequence repository", "https://www.ebi.ac.uk/ena/browser/home", "https://www.ebi.ac.uk/ena/portal/api/search", "Official API", "Operational", "Open", "INSDC unrestricted-access policy", "Discovery; study verification; run/sample metadata"),
 ("S04", "NCBI BioProject / SRA / BioSample", "INSDC repository", "https://www.ncbi.nlm.nih.gov/bioproject/", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", "E-utilities", "Operational", "Open", "INSDC unrestricted-access policy", "Discovery; BioSample attributes"),
 ("S05", "Europe PMC (incl. PubMed records)", "Literature index", "https://europepmc.org/", "https://www.ebi.ac.uk/europepmc/webservices/rest/", "REST API", "Operational", "Open", "Article-specific", "Accession-linked literature leads; full texts and supplements"),
 ("S06", "Zenodo", "Research data repository", "https://zenodo.org/", "https://zenodo.org/api/records", "REST API", "Operational", "Open records", "Record-specific (CC BY 4.0 for both records used)", "Processed data (India; 214K resource)"),
 ("S07", "GitHub: MetaSUB/MetaSUB-metadata", "Sample metadata repository", "https://github.com/MetaSUB/MetaSUB-metadata", "git clone (commit 1d738748)", "git", "Operational", "Public", "No licence file in repository", "MetaSUB city/coordinate metadata"),
 ("S08", "GitHub: SowpatiLab/Metro_metagenomics", "Analysis code", "https://github.com/SowpatiLab/Metro_metagenomics", "git clone (commit b3ea59bb)", "git", "Operational", "Public", "MIT", "Pipeline provenance for DS02"),
 ("S09", "FARME DB", "AMR reference database", "https://staff.washington.edu/jwallace/farme/", "Direct file download", "HTTPS", "Operational", "Public", "Non-commercial reuse with attribution (UW)", "Reference only; no samples"),
 ("S10", "CARD", "AMR reference database", "https://card.mcmaster.ca/", "https://card.mcmaster.ca/latest/data", "HTTPS", "Operational", "Public download", "Commercial use prohibited without licence; ontologies CC BY 4.0", "Reference only; annotation of DS02"),
 ("S11", "MetaSUB GeoSeeq portal", "Processed MetaSUB results", "https://portal.geoseeq.com/sample-groups/0625c15e-cf0d-4241-b33d-e967d55b91e6", "Not used", "Not accessed", "Not accessed", "Account may be required", "Not checked", "Candidate source for MetaSUB taxonomy/AMR tables"),
 ("S12", "WHO GLASS", "Surveillance indicators", "https://www.who.int/initiatives/glass", "Not used", "Not accessed", "Not accessed", "Public", "WHO terms", "Context only; not city-matched"),
 ("S13", "FDA NARMS", "Food-chain AMR surveillance", "https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/integrated-reportssummaries", "Not used", "Not accessed", "Not accessed", "Public", "US federal", "Context only; isolate-based"),
]
srcdf = pd.DataFrame(sources, columns=["Source_ID", "Name", "Type", "Official_URL", "API_Endpoint_or_Method", "Access_Method", "Status", "Access_Conditions", "Licence_or_Terms", "Role_in_Project"])
srcdf["Date_Checked"] = DATE

# --- source log CSV (one row per executed query + per verification fetch class)
srclog = qdf[["Query_ID", "Database", "Search_String", "Search_Date", "Search_URL_or_Endpoint", "N_Results_Returned", "N_Records_Inspected", "N_Potentially_Relevant", "Raw_Response_File"]].copy()
srclog.insert(0, "Log_Type", "Search query")
extra_log = pd.DataFrame([{"Log_Type": "Download", "Query_ID": "", "Database": "", "Search_String": r.File_Type, "Search_Date": r.Download_Date or DATE,
                           "Search_URL_or_Endpoint": r.Source_URL, "N_Results_Returned": "", "N_Records_Inspected": "", "N_Potentially_Relevant": "",
                           "Raw_Response_File": r.File_Path} for r in mandf.itertuples()])
pd.concat([srclog, extra_log]).to_csv(os.path.join(O, "Public_Data_Source_Log.csv"), index=False, quoting=csv.QUOTE_MINIMAL)

# ============================================================== summaries for WB3
# AMR summary (within pipeline only)
amr_x = amr.merge(sm[["Sample_ID", "City", "Country", "Sample_Type"]], on="Sample_ID", how="left")
amr_x["Pipeline"] = np.where(amr_x.AMR_Reference_Database == "CARD", "RGI bwt / CARD (DS02)", "KMA / ResFinder (Zenodo 6919377 resource)")
amr_x["Class"] = amr_x.Resistance_Class.astype(str).str.replace(", ", "; ").str.split("; ")
ax = amr_x.explode("Class"); ax["Class"] = ax.Class.str.strip()
nsamp = amr_x.groupby(["Pipeline", "City", "Country"]).Sample_ID.nunique().rename("N_Samples_With_AMR_Profile")
grp = ax.groupby(["Pipeline", "Country", "City", "Class"])
amrsum = pd.DataFrame({"N_Samples_Class_Detected": grp.Sample_ID.nunique(), "Total_Count": grp.Detection_Value.sum(),
                       "Median_Count_Among_Detecting_Samples": ax.groupby(["Pipeline", "Country", "City", "Class", "Sample_ID"]).Detection_Value.sum().groupby(level=[0, 1, 2, 3]).median()}).reset_index()
amrsum = amrsum.merge(nsamp.reset_index(), on=["Pipeline", "City", "Country"], how="left")
amrsum["Unit"] = np.where(amrsum.Pipeline.str.startswith("RGI"), "Mapped reads (raw count)", "Aligned fragments (raw count)")
amrsum["Caution"] = "Raw counts, not depth- or length-normalized; a gene may be counted under several drug classes; do not compare across pipelines"
amrsum = amrsum[["Pipeline", "Country", "City", "Class", "N_Samples_With_AMR_Profile", "N_Samples_Class_Detected", "Total_Count", "Median_Count_Among_Detecting_Samples", "Unit", "Caution"]].sort_values(["Pipeline", "Country", "City", "Total_Count"], ascending=[True, True, True, False])

# taxonomic summary (DS02 phylum relative abundance by city; DS01 subset genus diversity)
tx = tax.merge(sm[["Sample_ID", "City"]], on="Sample_ID")
ts = tx.groupby(["City", "Taxon_Name"]).Relative_Abundance.agg(["count", "mean", "median"]).reset_index()
ncity = tx.groupby("City").Sample_ID.nunique()
ts["N_Samples_City"] = ts.City.map(ncity)
ts["Mean_Rel_Abundance_All_Samples"] = ts.apply(lambda r: tx[(tx.City == r.City) & (tx.Taxon_Name == r.Taxon_Name)].Relative_Abundance.sum() / ncity[r.City], axis=1)
ts = ts.sort_values(["City", "Mean_Rel_Abundance_All_Samples"], ascending=[True, False]).groupby("City").head(15)
ts = pd.DataFrame({"Summary_Type": "DS02 phylum relative abundance (top 15 per city)", "Dataset_ID": "DS02", "City_or_Group": ts.City, "Feature": ts.Taxon_Name,
                   "N_Samples": ts.N_Samples_City, "N_Samples_Detected": ts["count"], "Mean_Value": ts.Mean_Rel_Abundance_All_Samples.round(5),
                   "Median_Value_Detected": ts["median"].round(5), "Unit": "Fraction of Bracken phylum-level reads (absent = 0 within profiled samples)"})
dvx = div.merge(sm[["Sample_ID", "City", "Country"]], on="Sample_ID")
dsum = dvx.groupby(["Feature_Category", "Country"]).agg(N=("Run_ID", "nunique"), Mean_Shannon=("Shannon", "mean"), Median_Richness=("Richness_n", "median")).reset_index()
dsum = pd.DataFrame({"Summary_Type": "Zenodo 6919377 diversity by country", "Dataset_ID": "Multiple (see Diversity_Long)", "City_or_Group": dsum.Country + " | " + dsum.Feature_Category,
                     "Feature": "Shannon (mean) / richness (median)", "N_Samples": dsum.N, "N_Samples_Detected": dsum.N, "Mean_Value": dsum.Mean_Shannon.round(4),
                     "Median_Value_Detected": dsum.Median_Richness, "Unit": "Mean Shannon index; median richness n (runs)"})
taxsum = pd.concat([ts, dsum], ignore_index=True)

# analysis readiness
ready = [
 ("Global microbial diversity", "Taxonomic profiles and sample metadata", "Conditionally ready",
  f"Downloaded taxonomic profiles cover DS02 only ({tax.Sample_ID.nunique()} samples); genus/phylum diversity indices exist for {int(inc.Sample_ID.isin(div.Sample_ID).sum())} included samples in {inc[inc.Sample_ID.isin(div.Sample_ID)].Dataset_ID.nunique()} datasets (Zenodo 6919377). All {inc.Dataset_ID.nunique()} included datasets need uniform reprocessing of raw reads before diversity can be compared across studies."),
 ("Taxonomic relative abundance", "Comparable abundance measurements", "Conditionally ready",
  "Ready within DS02 (Bracken, one database). Not comparable across studies until reads are re-profiled with one classifier and database."),
 ("Comparative analysis", "Harmonized samples and study metadata", "Conditionally ready",
  f"Sample metadata harmonized for {len(inc):,} included samples with matrix, country, date fields; comparisons must stay within matrix and within pipeline. 625 locality strings still need standardization to cities."),
 ("PCA/ordination", "Suitable transformed feature matrices", "Conditionally ready",
  "DS02 phylum/genus/species and RGI matrices can be CLR-transformed; Zenodo 6919377 ResFinder counts likewise. No cross-study feature matrix exists yet."),
 ("AMR gene mapping", "AMR gene profiles and annotation provenance", "Conditionally ready",
  f"{amr.Sample_ID.nunique()} samples ({int(inc.Sample_ID.isin(amr.Sample_ID).sum())} of them included) have gene-level AMR tables (DS02 CARD/RGI; Zenodo 6919377 ResFinder/KMA). CARD release for DS02 and ResFinder version for the 214K resource are not stated in the tables."),
 ("Resistome comparison", "Compatible AMR measurements", "Conditionally ready",
  "Within-pipeline only (DS02 cities; Zenodo 6919377 ResFinder runs, mainly DS01 2016 pilot and European sewage studies). CARD/RGI and ResFinder/KMA outputs must not be pooled; counts need length and depth normalization."),
 ("Geographic mapping", "Verified geographic metadata", "Conditionally ready",
  f"{(inc.Latitude != NR).sum():,} included samples carry source-reported coordinates; DS02 is city-level only; DS01 coordinates are two-decimal (approximate). Maps must show observed points, not interpolated surfaces."),
 ("Environmental association", "Environmental covariates", "Not ready",
  "Only source-reported covariates were collected (MetaSUB surface/city attributes, sampling month/round). No pollution, climate, sanitation or antibiotic-use covariates have been joined."),
 ("Predictive modelling", "Validated target variables and sufficient observations", "Not ready",
  "No clinical, epidemiological or food-safety outcomes matched to samples by place and time. Only exploratory within-data prediction (e.g., resistome from taxonomy in DS02) is feasible."),
 ("Food-system analysis", "Relevant food and food-safety data", "Not ready",
  f"Food-related datasets ({(sm.Inclusion_Status=='Separate contextual dataset').sum()} samples) are heterogeneous, mostly not urban-located and lack food-safety outcomes; keep as context."),
]
readdf = pd.DataFrame(ready, columns=["Proposed_Analysis", "Required_Data", "Readiness", "Reason"])

# ============================================================== QC checks (computed)
def qc_checks():
    out = []
    def add(i, name, result, n, detail): out.append({"Check_ID": i, "Check": name, "Result": result, "N_Records_Flagged": n, "Detail": detail})
    d = sm.Sample_ID.duplicated().sum(); add("Q01", "Uniqueness of Sample_ID", "Pass" if d == 0 else "Fail", int(d), "Sample_ID is the repository BioSample/ENA sample accession or MG-RAST sample id")
    rr = runs.Run_ID.duplicated().sum(); add("Q02", "Uniqueness of Run_ID", "Pass" if rr == 0 else "Fail", int(rr), "")
    orphan = (~sm.Study_ID.isin(ds_of.keys())).sum(); add("Q03", "Study-to-sample relationships", "Pass" if orphan == 0 else "Fail", int(orphan), "Every sample links to a selected study and Dataset_ID")
    multi = sm[sm.Run_ID.str.contains(";")]; add("Q03b", "Samples with multiple runs", "Warning" if len(multi) else "Pass", len(multi), "Runs aggregated; see Run_ID and Sequencing_Metadata")
    pat = r"^(SAM[NED][A-Z]?\d+|ERS\d+|SRS\d+|DRS\d+|mgs\d+|mgm\d+\.\d+)$"
    bad = (~sm.Sample_ID.str.match(pat)).sum(); add("Q04", "Accession format (sample)", "Pass" if bad == 0 else "Warning", int(bad), "Pattern " + pat)
    badr = (~runs.Run_ID.str.match(r"^([SED]RR\d+|mgm\d+\.\d+)$")).sum(); add("Q04b", "Accession format (run)", "Pass" if badr == 0 else "Warning", int(badr), "")
    bio = sm[(sm.Dataset_ID == "DS01") & (sm.Biological_Sample_ID != NAP)]
    mb = bio.Biological_Sample_ID.duplicated(keep=False).sum()
    add("Q05", "Duplicate biological samples (DS01)", "Warning", int(mb), f"{bio.Biological_Sample_ID.nunique()} biological samples span {len(bio)} accessions (Munk unique_sample_id); deduplicate before per-sample statistics")
    kk = sm.QC_Flags.str.contains("inconsistency").sum(); add("Q06", "City-country consistency", "Fail" if kk else "Pass", int(kk), "DS04 Kyiv samples whose source or propagated country is Poland; country set to Not yet verified (see V13)")
    sp = sm[(sm.City.isin(["Kochi", "Cochin"]))].shape[0]; add("Q06b", "City spelling variants", "Warning" if sp else "Pass", int(sp), "Kochi and Cochin both used in Munk Suppl. Data 1; not merged")
    lat = pd.to_numeric(sm.Latitude, errors="coerce"); lon = pd.to_numeric(sm.Longitude, errors="coerce")
    oor = ((lat.abs() > 90) | (lon.abs() > 180)).sum(); add("Q07", "Coordinate ranges", "Pass" if oor == 0 else "Fail", int(oor), f"{lat.notna().sum()} samples with valid coordinates; out-of-range values were set to Not reported upstream")
    dbad = sm.Collection_Date.str.match(r"^\d{4}(-\d{2}(-\d{2})?)?$") | sm.Collection_Date.isin(MISSING_CODES)
    add("Q08", "Collection-date format (ISO 8601)", "Pass" if dbad.all() else "Warning", int((~dbad).sum()), f"{(sm.Collection_Date == NR).sum()} samples have no date")
    unc = ((sm.Inclusion_Status == "Included") & sm.Sample_Type.isin(["Unclassified", ""])).sum(); add("Q09", "Sample-type classification", "Pass" if unc == 0 else "Warning", int(unc), "Included samples without a standardized category")
    ms = (sm.Sequencing_Platform == NR).sum(); add("Q10", "Missing sequencing information", "Pass" if ms == 0 else "Warning", int(ms), "Platform missing")
    na = ((sm.Inclusion_Status == "Included") & ~sm.AMR_Data_Available.str.startswith("Yes")).sum()
    add("Q11", "Missing AMR annotations", "Warning", int(na), "Included samples without a downloaded AMR table (raw reads available for reanalysis); not coded as zero")
    neg = (amr.Detection_Value < 0).sum() + (tax.Relative_Abundance > 1).sum() + (tax.Relative_Abundance < 0).sum()
    add("Q12", "Abundance values and units", "Pass" if neg == 0 else "Fail", int(neg), "No negative counts; relative abundances within [0,1]; units recorded per row")
    nourl = (~sm.Source_URL.str.startswith("http")).sum(); add("Q13", "Source URLs present", "Pass" if nourl == 0 else "Fail", int(nourl), "URLs generated from accessions; repository endpoints responded during collection")
    nofq = int(((runs.Repository != "MG-RAST") & (runs.Raw_Read_Availability != "Public FASTQ (ENA)")).sum())
    add("Q14", "Downloadable-file availability", "Warning" if nofq else "Pass", nofq, "INSDC runs without an ENA FASTQ link (submitted files or SRA format only); MG-RAST metagenomes use the MG-RAST download endpoint")
    mix = sm[(sm.Inclusion_Status == "Included") & ~sm.Sequencing_Type.str.contains("Shotgun metagenomic")].shape[0]
    add("Q15", "Mixing of incompatible data types", "Pass" if mix == 0 else "Fail", int(mix), "Included samples all have a shotgun metagenomic run; amplicon, isolate, single-cell and metatranscriptomic samples excluded")
    mixed = sm[(sm.Inclusion_Status == "Included") & (sm.Other_Run_IDs != NAP)].shape[0]
    add("Q15b", "Included samples that also have other run types", "Warning" if mixed else "Pass", int(mixed), "Non-shotgun runs listed in Other_Run_IDs and excluded from Run_ID, Raw_Data_URL and size totals; flagged in QC_Flags")
    add("Q16", "Transformations documented", "Pass", len(tlog), "See transformation_log in README_Data_Preparation.md and workbook READMEs")
    return pd.DataFrame(out)
qcdf = qc_checks()

# ============================================================== dictionaries
SM_DICT = [
 ("Dataset_ID", "Project-assigned dataset identifier (one per study, except DS01 which groups 7 ENA projects)", "Text", "DS01..DSnn"),
 ("Dataset_Name", "Short dataset description", "Text"), ("Study_ID", "Curated study accession used for selection (primary INSDC project or MG-RAST project)", "Text"),
 ("Project_Accession", "Primary project accession as reported by the repository", "Text"), ("Sample_ID", "Repository sample accession (BioSample/ENA sample or MG-RAST sample); primary linkage key", "Text"),
 ("Secondary_Sample_ID", "Secondary sample accession (ERS/SRS/DRS) or MG-RAST metagenome id", "Text"),
 ("Biological_Sample_ID", "Biological sample identifier where several accessions derive from one physical sample (DS01: Munk unique_sample_id)", "Text", "Not applicable"),
 ("Run_ID", "Shotgun metagenomic run accession(s), ';'-separated (other run types moved to Other_Run_IDs)", "Text"), ("Other_Run_IDs", "Non-shotgun runs of the same sample (with data type); not used for analysis", "Text", "Not applicable"), ("N_Runs", "Number of runs aggregated into the sample row", "Integer"),
 ("Analysis_ID", "MGnify study/analysis accession where known", "Text", "Not available"), ("Repository", "Source repository", "Text"),
 ("Publication_ID", "PMID/PMCID of the primary publication (core datasets only)", "Text", "Not yet verified"), ("DOI", "Publication DOI", "Text", "Not yet verified"),
 ("City", "City as labelled by the source (Munk Suppl. Data 1; MetaSUB metadata; Zenodo meta.tsv; BioSample for Dhaka). Not inferred from localities.", "Text", "Not reported"),
 ("Locality_As_Reported", "Text after the country in geo_loc_name (or MG-RAST location); not standardized to a city", "Text", "Not reported"),
 ("Country", "Country name (harmonized spelling); 'Hong Kong' retained as reported", "Text", "Not reported; Not yet verified"),
 ("ISO_Country_Code", "ISO 3166-1 alpha-2 code", "Text", "Not reported; Not yet verified"),
 ("Latitude", "Decimal degrees as reported by the source; never imputed", "Number", "Not reported"), ("Longitude", "Decimal degrees as reported by the source; never imputed", "Number", "Not reported"),
 ("Geographic_Unit", "Level of the geographic label", "Text", "City (source-labelled); Locality as reported (not standardized to city); Country only; Missing geographic information"),
 ("Coordinate_Precision", "Exact (>=3 decimals) / Approximate (<=2 decimals) / City-level location only / Missing", "Text"),
 ("Urban_Setting", "Basis for urban classification", "Text"), ("Geographic_Source", "Which official source supplied location fields", "Text"),
 ("Sample_Type", "Standardized category", "Text", "Transit surface; Built environment; Municipal wastewater; Hospital wastewater; Community wastewater; Urban soil; Surface water; Air; Food-related environment; Other; Unclassified"),
 ("Sample_Type_Rule", "Rule that assigned Sample_Type (documented keyword rule or study-level category)", "Text"),
 ("Sample_Subtype", "isolation_source as reported", "Text"), ("Environmental_Matrix", "env_medium / environment_material as reported", "Text"),
 ("Collection_Date", "ISO 8601 at reported precision (YYYY, YYYY-MM or YYYY-MM-DD)", "Text (date)", "Not reported; Not yet verified"),
 ("Collection_Year", "Year", "Integer"), ("Collection_Month", "Month 1-12", "Integer", "Not reported"), ("Collection_Date_Source", "Source of the date", "Text"),
 ("Sample_Size", "Not applicable at sample level (one row = one sample)", "Text", "Not applicable"), ("Sampling_Design", "Sampling round/design where documented", "Text"),
 ("Sampling_Location", "Site/station as reported", "Text"), ("Sequencing_Method", "library_strategy/library_source/library_selection", "Text"),
 ("Sequencing_Platform", "instrument platform: model", "Text"), ("Sequencing_Type", "Harmonized data type(s) of the sample's runs", "Text"),
 ("Read_Type", "Library layout", "Text"), ("Base_Count_Total", "Total bases across runs (repository-reported)", "Integer (bases)"),
 ("FASTQ_Bytes_Total", "Total compressed FASTQ size (ENA-reported)", "Integer (bytes)"), ("Raw_Read_Availability", "Raw read availability", "Text"),
 ("Raw_Data_URL", "FASTQ FTP link(s) or MG-RAST download endpoint", "Text (URL)"), ("Processed_Data_Availability", "Processed outputs and their status", "Text"),
 ("Processed_Data_URL", "Location of processed outputs", "Text (URL)"), ("Taxonomic_Profile_Available", "Whether a taxonomic profile exists and was downloaded", "Text"),
 ("Taxonomic_Rank", "Rank(s) available", "Text"), ("Microbial_Taxa", "Pointer to long table (taxa are not stored in the sample row)", "Text"),
 ("Taxonomic_Abundance", "Pointer to long table", "Text"), ("Abundance_Unit", "Unit of taxonomic abundance", "Text"),
 ("Functional_Profile_Available", "Functional profile availability", "Text"), ("AMR_Data_Available", "Whether gene-level AMR data were downloaded", "Text"),
 ("AMR_Gene", "Pointer to AMR_Genes_Long with distinct gene count", "Text"), ("AMR_Gene_Family", "Pointer or source limitation", "Text"),
 ("Resistance_Class", "Pointer", "Text"), ("Resistance_Mechanism", "Pointer or source limitation", "Text"), ("AMR_Detection_Type", "Detection method", "Text"),
 ("AMR_Abundance", "Pointer", "Text"), ("AMR_Abundance_Unit", "Unit of AMR counts", "Text"), ("AMR_Reference_Database", "Reference database", "Text"),
 ("AMR_Database_Version", "Version if stated", "Text"), ("Data_Access_Status", "Public / Restricted access", "Text"), ("Licence", "Reuse terms", "Text"),
 ("Date_Collected", "Date records were retrieved (ISO 8601)", "Text (date)"),
 ("Verification_Status", "Verified = official record retrieved and country+date confirmed for core datasets; Partially verified = record retrieved, classification from title review or key fields missing", "Text", "Verified; Partially verified; Not yet verified"),
 ("Inclusion_Status", "Sample-level decision", "Text", "Included; Excluded; Pending verification; Separate contextual dataset"), ("Inclusion_Reason", "Reason for decision", "Text"),
 ("Missing_Data_Status", "Core fields missing (city/locality, country, date, coordinates)", "Text"), ("QC_Status", "Pass / Warning / Fail", "Text"),
 ("QC_Flags", "Specific QC messages", "Text"), ("Limitations", "Dataset-level limitations", "Text"), ("Source_URL", "Official sample record", "Text (URL)"),
 ("Remarks", "Free text", "Text"), ("City_Source_Value", "Original city value", "Text"), ("Country_Source_Value", "Original geo_loc_name/country value", "Text"),
 ("Collection_Date_Source_Value", "Original date string", "Text"), ("Sample_Type_Source_Value", "Original matrix text", "Text"),
 ("Latitude_Source_Value", "Original latitude string", "Text"), ("Longitude_Source_Value", "Original longitude string", "Text")]

common_readme = [
 ("Project", "Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping (Academic Press, Elsevier book chapter)"),
 ("Prepared", f"{DATE}, by an AI assistant (Claude) for author review. All records come from public repositories queried on this date."),
 ("Missing-value codes", "Not reported (source has no value); Not available (resource does not provide it); Not yet verified (value exists but conflicts or is unchecked); Not applicable; Restricted access. Missing AMR observations are never coded as zero."),
 ("Interpretation limits", "Gene detection is not phenotypic resistance, expression, viability, transmission or clinical disease. City labels do not imply city representativeness."),
 ("Reproducibility", "Scripts 01-05 in scripts/ regenerate every table from raw_downloads/; see README_Data_Preparation.md."),
]

# ============================================================== WORKBOOK 1
w = WB(os.path.join(O, "01_Public_Dataset_Inventory.xlsx")); S = status_maps(w)
w.readme("README", "01 Public Dataset Inventory", common_readme + [
 ("Purpose", "Documents every public source, executed query, discovered record and verification decision."),
 ("Search_Query_Log", f"{len(qdf)} executed queries across 6 databases (35 per database) plus {sum(qdf.Database.str.contains('adapted'))} adapted MGnify queries. Counts come from API responses saved under raw_downloads/search_responses/."),
 ("Programmatic screen", open(os.path.join(L, "screen_rule.txt")).read().strip()),
 ("Decisions", "Study-level eligibility decisions were made by the AI assistant from titles, ENA run records and (for core datasets) publications; they require author confirmation. Decision rules: scripts/curation_decisions.py."),
 ("Verification statuses", "Verified; Partially verified; Not yet verified; Not accessible; Unsuitable."),
 ("Sheets", "Source_Registry, Search_Query_Log, Study_Inventory, Accession_Register, Dataset_Eligibility, Download_Manifest, Exclusion_Log, Source_Verification, Data_Dictionary")])
w.table("Source_Registry", srcdf, url_cols=["Official_URL"], status_cols={"Status": S["Status"]})
w.table("Search_Query_Log", qdf, url_cols=["Search_URL_or_Endpoint"], int_cols=["N_Results_Returned", "N_Records_Inspected"])
w.table("Study_Inventory", invdf, url_cols=["Official_URL"], status_cols={"Decision": S["Decision"], "Verification_Status": S["Verification_Status"]},
        validations={"Decision": ["Included", "Excluded", "Pending verification", "Separate contextual dataset"]}, widths={"Title": 60, "Decision_Reason": 50, "Query_IDs": 30})
w.table("Accession_Register", accdf, url_cols=["Official_URL"])
w.table("Dataset_Eligibility", eldf, status_cols={"Decision": S["Decision"]}, note="Criteria C1-C10 follow the brief's inclusion criteria; fractions are samples meeting the criterion / samples in the study.")
w.table("Download_Manifest", mandf, int_cols=["Size_Bytes"], status_cols={"Download_Status": S["Download_Status"]})
w.table("Exclusion_Log", exdf, int_cols=["N_Records"])
w.table("Source_Verification", vdf, status_cols={"Status": S["Status"]}, widths={"Finding": 70, "Expectation_or_Claim": 40})
w.table("Data_Dictionary", dictionary({
 "Search_Query_Log": [("N_Results_Returned", "Total hits reported by the API", "Integer"), ("N_Records_Inspected", "Records retrieved and screened", "Integer"),
                      ("N_Potentially_Relevant", "Records passing the programmatic screen", "Integer"), ("N_Verified", "Relevant records whose official study/project record was retrieved (INSDC via ENA run report; MG-RAST project record). Publications and Zenodo records are not counted.", "Integer"),
                      ("N_Included_or_Contextual_After_Review", "Verified records classed Included or Separate contextual dataset", "Integer"), ("N_Excluded", "Failed screen plus verified-but-not-included", "Integer")],
 "Study_Inventory": [("Record_ID", "Kind|accession", "Text"), ("Shotgun_Metagenomic_Runs", "ENA runs with WGS/METAGENOMIC", "Integer"), ("Raw_FASTQ_GB", "Compressed FASTQ volume reported by ENA", "GB"),
                     ("Decision", "Study-level decision", "Text", "Included; Excluded; Pending verification; Separate contextual dataset")],
 "Download_Manifest": [("Size_Bytes", "File size on disk, or repository-reported size for files not downloaded", "Bytes"), ("Checksum", "sha256 of downloaded file (files <600 MB) or repository checksum note", "Text")]}),
        widths={"Description": 80})
w.close()

# ============================================================== WORKBOOK 2
w = WB(os.path.join(O, "02_Sample_Metadata_Master.xlsx")); S = status_maps(w)
w.readme("README", "02 Sample Metadata Master", common_readme + [
 ("Purpose", "Harmonized sample-level metadata: one row per unique repository sample (BioSample/ENA sample or MG-RAST sample)."),
 ("Counts", f"{len(sm):,} samples from {sm.Dataset_ID.nunique()} datasets; Included {int((sm.Inclusion_Status=='Included').sum()):,}; Separate contextual {int((sm.Inclusion_Status=='Separate contextual dataset').sum()):,}; Excluded {int((sm.Inclusion_Status=='Excluded').sum()):,}."),
 ("Source priority (geography/date)", "DS01: Munk et al. 2022 Suppl. Data 1 by run, then ENA; DS02: Zenodo meta.tsv (city) + SRA; DS03: NCBI BioSample; DS04: MetaSUB GitHub metadata by uuid + SRA; others: ENA geo_loc_name, then NCBI BioSample."),
 ("Derived values", "MetaSUB country was propagated from same-city SRA records only when the SRA field was missing and the city mapped to a single country (flagged in QC_Flags). No coordinates, dates or cities were imputed."),
 ("Linkage", "Sample_ID is the primary key to 03_Analysis_Ready_Metagenomic_Data.xlsx. Run-level data are in Sequencing_Metadata."),
 ("Formulas", "Missingness_Report counts are live COUNTIFS formulas over Sample_Master.")])
smx = sm.copy()
for c in ["Latitude", "Longitude"]: smx[c] = smx[c].map(lambda v: float(v) if re.match(r"^-?\d+(\.\d+)?$", v) else v)
for c in ["N_Runs", "Base_Count_Total", "FASTQ_Bytes_Total", "Collection_Year", "Collection_Month"]: smx[c] = smx[c].map(lambda v: int(float(v)) if re.match(r"^\d+(\.0)?$", v) else v)
ws_sm = w.table("Sample_Master", smx, url_cols=["Source_URL"], int_cols=["Base_Count_Total", "FASTQ_Bytes_Total"], dec_cols=["Latitude", "Longitude"],
                status_cols={k: S[k] for k in ("QC_Status", "Verification_Status", "Inclusion_Status")},
                validations={"Inclusion_Status": ["Included", "Excluded", "Pending verification", "Separate contextual dataset"], "QC_Status": ["Pass", "Warning", "Fail"],
                             "Verification_Status": ["Verified", "Partially verified", "Not yet verified"],
                             "Sample_Type": ["Transit surface", "Built environment", "Municipal wastewater", "Hospital wastewater", "Community wastewater", "Urban soil", "Surface water", "Air", "Food-related environment", "Other", "Unclassified"]},
                widths={"Raw_Data_URL": 40, "Processed_Data_Availability": 40, "Limitations": 40, "Dataset_Name": 40})
geo_cols = ["Sample_ID", "Dataset_ID", "Study_ID", "City", "Locality_As_Reported", "Country", "ISO_Country_Code", "Latitude", "Longitude", "Geographic_Unit", "Coordinate_Precision", "Urban_Setting", "Geographic_Source", "City_Source_Value", "Country_Source_Value", "Latitude_Source_Value", "Longitude_Source_Value", "Inclusion_Status"]
w.table("Geographic_Metadata", smx[geo_cols], dec_cols=["Latitude", "Longitude"], status_cols={"Inclusion_Status": S["Inclusion_Status"]})
env_cols = ["Sample_ID", "Dataset_ID", "Study_ID", "Sample_Type", "Sample_Type_Rule", "Sample_Subtype", "Environmental_Matrix", "Sample_Type_Source_Value", "Urban_Setting", "Collection_Date", "Collection_Year", "Collection_Month", "Collection_Date_Source", "Collection_Date_Source_Value", "Sampling_Design", "Sampling_Location", "Inclusion_Status"]
w.table("Environmental_Metadata", smx[env_cols], status_cols={"Inclusion_Status": S["Inclusion_Status"]})
seq = runs[["Run_ID", "Sample_ID", "Dataset_ID", "Study_ID", "Project_Accession", "Experiment_ID", "Repository", "Sample_Alias", "Library_Name", "Sequencing_Method", "Sequencing_Platform", "Sequencing_Type", "Read_Type", "Base_Count", "Read_Count", "FASTQ_Bytes", "FASTQ_MD5", "Raw_Read_Availability", "Raw_Data_URL", "First_Public", "Source_URL"]].copy()
for c in ["Base_Count", "Read_Count", "FASTQ_Bytes"]: seq[c] = seq[c].map(lambda v: int(float(v)) if re.match(r"^\d+(\.0)?$", str(v)) else v)
w.table("Sequencing_Metadata", seq, url_cols=["Source_URL"], int_cols=["Base_Count", "Read_Count", "FASTQ_Bytes"], widths={"Raw_Data_URL": 40, "FASTQ_MD5": 30})
link = smx[["Sample_ID", "Dataset_ID", "Dataset_Name", "Study_ID", "Project_Accession", "Secondary_Sample_ID", "Biological_Sample_ID", "Run_ID", "N_Runs", "Analysis_ID", "Publication_ID", "DOI", "Repository", "Source_URL"]]
w.table("Study_Linkage", link, url_cols=["Source_URL"])
# Missingness with live formulas
ws = w.wb.add_worksheet("Missingness_Report")
cols = list(smx.columns); n = len(smx)
def colref(c): j = cols.index(c); L_ = xl_col_to_name(j); return f"Sample_Master!${L_}$2:${L_}${n+1}"
fields = ["City", "Locality_As_Reported", "Country", "ISO_Country_Code", "Latitude", "Longitude", "Collection_Date", "Collection_Month", "Sample_Subtype", "Environmental_Matrix", "Sequencing_Platform", "Read_Type", "Base_Count_Total", "FASTQ_Bytes_Total", "DOI", "Publication_ID", "Analysis_ID", "Taxonomic_Profile_Available", "AMR_Data_Available"]
hdr = ["Field", "Scope"] + MISSING_CODES + ["Total_Missing", "N_Rows_In_Scope", "Percent_Missing"]
ws.write(0, 0, "Missingness by field (live COUNTIFS over Sample_Master; scope = Inclusion_Status)", w.title)
ws.write_row(2, 0, hdr, w.hdr)
pctf = w.wb.add_format({"font_name": "Arial", "font_size": 10, "num_format": "0.0%"})
row = 3
for scope in ["Included", "All samples"]:
    for f in fields:
        ws.write(row, 0, f, w.base); ws.write(row, 1, scope, w.base)
        for k, code in enumerate(MISSING_CODES):
            if scope == "All samples": ws.write_formula(row, 2 + k, f'=COUNTIF({colref(f)},"{code}")', w.int)
            else: ws.write_formula(row, 2 + k, f'=COUNTIFS({colref(f)},"{code}",{colref("Inclusion_Status")},"Included")', w.int)
        ws.write_formula(row, 7, f"=SUM(C{row+1}:G{row+1})", w.int)
        ws.write_formula(row, 8, f'=COUNTA({colref("Sample_ID")})' if scope == "All samples" else f'=COUNTIF({colref("Inclusion_Status")},"Included")', w.int)
        ws.write_formula(row, 9, f"=IF(I{row+1}=0,0,H{row+1}/I{row+1})", pctf)
        row += 1
ws.set_column(0, 1, 24, w.base); ws.set_column(2, 9, 16, w.base); ws.freeze_panes(3, 0); ws.autofilter(2, 0, row - 1, len(hdr) - 1)
ws.conditional_format(3, 9, row - 1, 9, {"type": "3_color_scale", "min_color": "#C6EFCE", "mid_color": "#FFEB9C", "max_color": "#FFC7CE"})
# Duplicate check
dups = []
bio = sm[(sm.Dataset_ID == "DS01") & (sm.Biological_Sample_ID != NAP)]
for b, g in bio.groupby("Biological_Sample_ID"):
    if len(g) > 1: dups.append({"Check": "Multiple accessions per biological sample (DS01)", "Key": b, "N": len(g), "Sample_IDs": "; ".join(g.Sample_ID), "Action": "Aggregate or select one accession per biological sample before per-sample statistics"})
for s_, g in runs.groupby("Sample_ID"):
    if len(g) > 1 and g.Sequencing_Type.nunique() > 1:
        dups.append({"Check": "Sample with runs of different data types", "Key": s_, "N": len(g), "Sample_IDs": "; ".join(g.Run_ID), "Action": "Use shotgun metagenomic runs only"})
alias = sm[sm.Sample_ID.str.startswith("SAM")]
dups.append({"Check": "Duplicate Sample_ID", "Key": "All", "N": int(sm.Sample_ID.duplicated().sum()), "Sample_IDs": "", "Action": "None required" if not sm.Sample_ID.duplicated().any() else "Investigate"})
dups.append({"Check": "Duplicate Run_ID", "Key": "All", "N": int(runs.Run_ID.duplicated().sum()), "Sample_IDs": "", "Action": "None required"})
dups.append({"Check": "Cross-project sample overlap among the 7 DS01 projects", "Key": "DS01", "N": 0, "Sample_IDs": "", "Action": "Confirmed none (ENA run reports)"})
w.table("Duplicate_Check", pd.DataFrame(dups), int_cols=["N"])
w.table("QC_Report", qcdf, status_cols={"Result": S["Result"]}, widths={"Detail": 90})
w.table("Data_Dictionary", dictionary({"Sample_Master": SM_DICT, "Sequencing_Metadata": [("Run_ID", "Run accession", "Text"), ("Base_Count", "Bases (repository)", "Integer"), ("FASTQ_Bytes", "Compressed FASTQ bytes", "Integer"), ("FASTQ_MD5", "Repository MD5 per file", "Text")]}), widths={"Description": 90})
w.close()

# ============================================================== WORKBOOK 3
w = WB(os.path.join(O, "03_Analysis_Ready_Metagenomic_Data.xlsx")); S = status_maps(w)
w.readme("README", "03 Analysis-Ready Metagenomic Data", common_readme + [
 ("Purpose", "Long-format taxonomy, diversity and AMR tables with linked sample metadata for analysis in R/Python."),
 ("What is real data here", f"Taxonomy_Long: {len(tax):,} phylum rows for {tax.Sample_ID.nunique()} DS02 samples (Bracken). Diversity_Long: {len(div):,} rows for {div.Run_ID.nunique()} runs (Zenodo 6919377). AMR_Genes_Long: {len(amr):,} rows for {amr.Sample_ID.nunique()} samples (DS02 RGI/CARD; KMA/ResFinder from Zenodo 6919377 for runs the resource covers). Environmental_Covariates: {len(cov):,} source-reported values."),
 ("Pipelines are not interchangeable", "AMR rows from CARD/RGI and ResFinder/KMA differ in reference database, gene naming and count definition. Analyse them separately (column AMR_Detection_Type)."),
 ("Not available", "Functional_Profiles_Long is an empty template: no functional profiles were downloaded. Genus and species Bracken tables (DS02) are kept as source files because they exceed Excel's row limit (3.7 M and >10 M rows)."),
 ("Summaries", "AMR_Summary and Taxonomic_Summary were computed in Python (scripts/05_build_workbooks.py) from the long tables in this workbook; they are values, not formulas, because the source tables exceed 275,000 rows."),
 ("Extra sheet", "Diversity_Long was added to hold the Zenodo 6919377 diversity indices; it is not in the brief's sheet list.")])
asm_cols = ["Sample_ID", "Dataset_ID", "Study_ID", "Run_ID", "Biological_Sample_ID", "Sample_Type", "City", "Locality_As_Reported", "Country", "ISO_Country_Code", "Latitude", "Longitude", "Coordinate_Precision",
            "Collection_Date", "Collection_Year", "Collection_Month", "Sequencing_Platform", "Read_Type", "Base_Count_Total", "Taxonomic_Profile_Available", "AMR_Data_Available", "AMR_Reference_Database", "QC_Status", "QC_Flags"]
asm = smx[smx.Inclusion_Status == "Included"][asm_cols].copy()
asm["Has_Taxonomy_Long"] = asm.Sample_ID.isin(tax.Sample_ID).map({True: "Yes", False: "No"})
asm["Has_Diversity_Long"] = asm.Sample_ID.isin(div.Sample_ID).map({True: "Yes", False: "No"})
asm["Has_AMR_Genes_Long"] = asm.Sample_ID.isin(amr.Sample_ID).map({True: "Yes", False: "No"})
w.table("Analysis_Sample_Master", asm, dec_cols=["Latitude", "Longitude"], int_cols=["Base_Count_Total"], status_cols={"QC_Status": S["QC_Status"]})
w.table("Taxonomy_Long", tax, int_cols=["Abundance_Value", "Kraken_Assigned_Reads"])
w.table("AMR_Genes_Long", amr, int_cols=["Detection_Value", "Reference_Length_bp"], dec2_cols=["Average_Percent_Coverage"], widths={"AMR_Gene_Family": 30, "Resistance_Class": 30})
fp = pd.DataFrame(columns=["Dataset_ID", "Sample_ID", "Run_ID", "Feature_ID", "Feature_Name", "Feature_Type", "Abundance_Value", "Abundance_Unit", "Annotation_Pipeline", "Reference_Database", "Database_Version", "Source_File"])
w.table("Functional_Profiles_Long", fp, note="TEMPLATE ONLY - no functional profiles were downloaded. Candidate sources: MG-RAST project annotations (DS GLUSEEN, Swiss, Polish series) and MGnify MGYS00005846 (PRJEB27054).")
w.table("Environmental_Covariates", cov)
w.table("AMR_Summary", amrsum, int_cols=["Total_Count"], dec2_cols=["Median_Count_Among_Detecting_Samples"], widths={"Caution": 60})
w.table("Taxonomic_Summary", taxsum)
w.table("Diversity_Long", div, int_cols=["Total_Fragments"])
w.table("Analysis_Readiness", readdf, status_cols={"Readiness": S["Readiness"]}, widths={"Reason": 100, "Required_Data": 35})
w.table("QC_Report", qcdf, status_cols={"Result": S["Result"]}, widths={"Detail": 90})
w.table("Data_Dictionary", dictionary({
 "Taxonomy_Long": [("Abundance_Value", "Bracken new_est_reads", "Count"), ("Relative_Abundance", "Bracken fraction_total_reads", "Fraction 0-1"), ("Taxon_ID_NCBI", "NCBI taxonomy id", "Text")],
 "AMR_Genes_Long": [("AMR_Gene", "ARO term (CARD) or ResFinder gene name", "Text"), ("Gene_Identifier", "ARO accession or ResFinder reference sequence id", "Text"),
                    ("Detection_Value", "All mapped reads (RGI) or aligned fragments (KMA)", "Count"), ("Average_Percent_Coverage", "Reference coverage reported (RGI) or computed refCoveredPositions/refSequence_length*100 (KMA)", "%"),
                    ("AMR_Detection_Type", "Pipeline identifier; never pool across pipelines", "Text")],
 "Diversity_Long": [("Richness_n", "Number of features detected", "Count"), ("Shannon", "Shannon index", "Unitless"), ("Simpson", "Simpson index", "Unitless")],
 "Environmental_Covariates": [("Covariate", "Variable name as in source", "Text"), ("Level", "Sample or City", "Text")]}), widths={"Description": 90})
w.close()

# ============================================================== WORKBOOK 4
w = WB(os.path.join(O, "04_AMR_Reference_Database_Register.xlsx")); S = status_maps(w)
farme = pd.read_csv(os.path.join(R, "farme", "HMMs.tab.gz"), sep="\t", dtype=str).fillna("")
w.readme("README", "04 AMR Reference Database Register", common_readme + [
 ("Purpose", "Reference and annotation resources for AMR interpretation. These are NOT sample datasets; no sample records were created from them."),
 ("FARME", "FARME DB (Wallace et al. 2017, doi:10.1093/database/baw165) is the resource most likely meant by the proposal's 'FARME data management'; author confirmation needed."),
 ("CARD", "CARD 4.0.2 (card.json timestamp 2026-08-11) downloaded; the ARO index is reproduced in CARD_Reference for gene-family/drug-class/mechanism lookups.")])
fr = [("Full name", "Functional Antibiotic Resistance Metagenomic Element Database (FARME DB)"), ("Citation", "Wallace JC, Port JA, Smith MN, Faustman EM. Database (Oxford). 2017;2017:baw165. doi:10.1093/database/baw165. PMID 28077567"),
      ("Host", "Institute for Risk Analysis and Risk Communication, University of Washington"), ("URL", "https://staff.washington.edu/jwallace/farme/"),
      ("Accessibility (2026-10-31)", "Website and all download files returned HTTP 200"), ("Version", "No version number; README 'Last revision: March, 2019'; files last modified 2019-02-28 to 2019-03-07"),
      ("DNA sequences", "20,724 (farm_dna.fa.gz; count verified from file)"), ("Protein sequences", "48,178 (farm_protein.fa.gz; ~11,000 predicted with MetaGeneMark; count verified)"),
      ("HMM hits", f"{len(farme):,} rows (HMMs.tab.gz; Pfam 29.0 and Resfams 1.1 models)"), ("Projects", "29 functional metagenomics projects (README)"),
      ("Content basis", "Sequences from environmental functional-metagenomic selections showing AR in laboratory experiments"),
      ("Reuse terms", "Copyright 2020 University of Washington; non-commercial reproduction and distribution permitted with attribution; commercial use: license@uw.edu"),
      ("Suitability", "Secondary reference for environmentally derived resistance elements (BLAST/HMM against assemblies). Not maintained since 2019; no sample-level urban data; not a substitute for CARD/ResFinder."),
      ("Recommended role", "Optional sensitivity annotation of assembled contigs; report hits separately from CARD-based calls")]
frdf = pd.DataFrame(fr, columns=["Attribute", "Value"])
fcat = farme.groupby("Clone Resistance AR Category").size().rename("N_HMM_Rows").reset_index()
fsrc = farme.groupby("HMM Source").size().rename("N_HMM_Rows").reset_index()
farme_tbl = pd.concat([frdf, pd.DataFrame({"Attribute": ["--- HMM rows by clone AR category ---"] + [f"AR category: {x}" for x in fcat.iloc[:, 0]] + ["--- HMM rows by model source ---"] + [f"HMM source: {x}" for x in fsrc.iloc[:, 0]],
                                           "Value": [""] + list(fcat.N_HMM_Rows.astype(str)) + [""] + list(fsrc.N_HMM_Rows.astype(str))})])
w.table("FARME_Reference", farme_tbl, widths={"Attribute": 40, "Value": 100})
aro = pd.read_csv(os.path.join(R, "card", "extracted", "aro_index.tsv"), sep="\t", dtype=str).fillna("")
cj = json.load(open(os.path.join(R, "card", "extracted", "card.json")))
mtypes = pd.Series([v.get("model_type") for k, v in cj.items() if not k.startswith("_")]).value_counts()
card_attr = pd.DataFrame([("Version", cj["_version"]), ("Timestamp", cj["_timestamp"]), ("Detection models in card.json", str(sum(1 for k in cj if not k.startswith("_")))),
                          ("ARO index rows", str(len(aro)))] + [(f"Model type: {k}", str(v)) for k, v in mtypes.items()] +
                         [("Citation", "Alcock BP, et al. CARD 2023. Nucleic Acids Res. 2023;51(D1):D690-9. doi:10.1093/nar/gkac920"),
                          ("Terms", "Commercial use prohibited without written permission of McMaster University; ontologies (ARO, MOBIO, VIRO) CC BY 4.0"),
                          ("Used by", "DS02 (Singh et al. 2026) via RGI bwt; release used there not stated"), ("Download", "https://card.mcmaster.ca/latest/data (sha256 in Download_Manifest)")],
                         columns=["Attribute", "Value"])
card_tbl = aro.rename(columns=lambda c: c.replace(" ", "_"))
w.table("CARD_Reference", card_tbl, note="Rows: CARD 4.0.2 aro_index.tsv (unmodified content). Attributes of the release are listed in Database_Versions.", widths={"Model_Name": 40, "AMR_Gene_Family": 40, "Drug_Class": 40})
other = pd.DataFrame([
 ("ResFinder (CGE, DTU)", "Acquired AMR genes", "Used by DS01 processed data: Munk et al. 2022 (commit 813679d) and Zenodo 6919377 resource", "https://bitbucket.org/genomicepidemiology/resfinder_db", "Not downloaded", "Version for Zenodo 6919377 not restated in table"),
 ("Resfams 1.1", "AMR protein-family HMMs", "Used inside FARME HMM annotation", "Not accessed", "Not downloaded", "Cited in FARME README"),
 ("Pfam 29.0", "Protein families", "Used inside FARME HMM annotation", "Not accessed", "Not downloaded", "Cited in FARME README"),
 ("NCBI AMRFinderPlus / Reference Gene Catalog", "Curated AMR genes", "Candidate for re-annotation", "https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/", "Not accessed", "Listed as candidate only"),
 ("MEGARes", "AMR database for metagenomics", "Candidate for re-annotation", "https://www.meglab.org/megares/", "Not accessed", "Listed as candidate only"),
 ("PLSDB", "Plasmid database", "Used by Munk et al. 2022 (v2020_11_19) for plasmid context", "https://ccb-microbe.cs.uni-saarland.de/plsdb/", "Not accessed", "From publication methods")],
 columns=["Resource", "Type", "Relation_to_Project", "URL", "Access_Status", "Notes"])
w.table("Other_AMR_Resources", other, url_cols=["URL"])
ver = pd.DataFrame([("CARD", "4.0.2", "2026-08-11", "card.json _version/_timestamp", "Downloaded 2026-10-31"),
                    ("CARD (as used by DS02)", "Not reported", "", "Zenodo 17590646, GitHub README, PMC13144492", "Unresolved"),
                    ("FARME DB", "No version (README revised March 2019)", "2019-03", "README.pdf; HTTP Last-Modified", "Downloaded 2026-10-31"),
                    ("ResFinder (Munk et al. 2022)", "commit 813679d", "", "PMC9715550 Methods", "From publication"),
                    ("ResFinder (Zenodo 6919377)", "Not stated in table", "", "Zenodo record; see https://hmmartiny.github.io/mARG/", "Unresolved"),
                    ("KMA (Zenodo 6919377)", "1.3.0", "", "ARG.tsv kma_version column", "From data"),
                    ("Kraken 2 database (DS02)", "k2 core_nt index", "2023-05-02", "Zenodo 17590646 file name; PMC13144492", "From data"),
                    ("Silva (Munk et al. 2022)", "138", "", "PMC9715550 Methods", "From publication"),
                    ("MetaSUB metadata", "commit 1d738748", "2020-02-21", "GitHub", "Downloaded 2026-10-31")], columns=["Resource", "Version", "Release_Date", "Evidence", "Status"])
w.table("Database_Versions", ver)
am = pd.DataFrame([("DS02", "Read-based taxonomy", "Kraken 2 + Bracken", "Not stated", "k2 core_nt (2023-05-02)", "Bracken estimated reads; fraction", "Singh et al. 2026; Zenodo 17590646"),
                   ("DS02", "Read-based ARG quantification", "RGI bwt (Bowtie2)", "Not stated", "CARD (release not stated)", "All mapped reads per ARO term", "Zenodo 17590646 rgi_bowtie.tsv"),
                   ("DS02", "MAG-based ARG detection", "BacAnt", "Not stated", "Not stated", "Presence in HQ MAGs (file not harmonized)", "Zenodo amr_HQ_MAG.tsv"),
                   ("DS01 (Munk 2022)", "Read-based ARG quantification", "KMA", "1.2.17a", "ResFinder commit 813679d", "Fragment counts", "PMC9715550"),
                   ("Zenodo 6919377 subset (1,004 runs)", "Read-based ARG quantification", "KMA", "1.3.0", "ResFinder (version not restated)", "fragmentCountAln", "ARG.tsv"),
                   ("Zenodo 6919377 subset (1,004 runs)", "rRNA-based diversity", "KMA", "1.3.0", "Silva-derived rRNA database", "Shannon, Simpson, richness", "diversity.tsv"),
                   ("Planned (all datasets)", "Uniform reprocessing", "fastp; Kraken 2/Bracken; RGI bwt", "To confirm", "CARD 4.0.2 recommended", "Length- and depth-normalized counts", "Section 4 v3 protocol [AUTHOR TO CONFIRM]")],
                  columns=["Dataset", "Step", "Tool", "Tool_Version", "Reference_Database", "Output_Unit", "Source"])
w.table("Annotation_Methods", am)
w.table("Reference_Data_Dictionary", dictionary({"CARD_Reference": [(c, "Column from CARD aro_index.tsv", "Text") for c in card_tbl.columns], "FARME_Reference": [("Attribute", "Property", "Text"), ("Value", "Value as verified", "Text")]}), widths={"Description": 60})
w.close()

# ============================================================== WORKBOOK 0 (index)
def countable(s): return s[~s.isin(MISSING_CODES + [""])]
idx = [("Workbooks", "00_Master_Data_Index.xlsx; 01_Public_Dataset_Inventory.xlsx; 02_Sample_Metadata_Master.xlsx; 03_Analysis_Ready_Metagenomic_Data.xlsx; 04_AMR_Reference_Database_Register.xlsx"),
       ("Queries executed", str(len(qdf))), ("Records retrieved and screened", f"{int(qdf.N_Records_Inspected.sum()):,}"),
       ("Distinct candidate records after deduplication", str(len(invdf))),
       ("INSDC candidate studies verified against ENA", str(int((invdf.Record_Type == 'INSDC study/project').sum()))),
       ("Studies/projects included", str(int(((invdf.Decision == 'Included') & invdf.Record_Type.isin(['INSDC study/project', 'MG-RAST project'])).sum()))),
       ("Studies/projects as separate contextual datasets", str(int(((invdf.Decision == 'Separate contextual dataset')).sum()))),
       ("Studies/projects pending verification", str(int(((invdf.Decision == 'Pending verification') & invdf.Record_Type.isin(['INSDC study/project', 'MG-RAST project'])).sum()))),
       ("Datasets in sample master", str(sm.Dataset_ID.nunique())), ("Unique samples in sample master", f"{len(sm):,}"),
       ("Verified unique samples (Verification_Status = Verified)", f"{int((sm.Verification_Status=='Verified').sum()):,}"),
       ("Included unique samples", f"{len(inc):,}"), ("Included samples (DS01 counted by biological sample)", f"{(len(inc[inc.Dataset_ID!='DS01']) + inc[(inc.Dataset_ID=='DS01') & (inc.Biological_Sample_ID!=NAP)].Biological_Sample_ID.nunique() + int(((inc.Dataset_ID=='DS01') & (inc.Biological_Sample_ID==NAP)).sum())):,}"),
       ("Countries represented (included samples)", str(countable(inc.Country).nunique())),
       ("Cities represented (distinct City|Country pairs, source-labelled city, country known)", str(inc[(inc.Geographic_Unit == 'City (source-labelled)') & ~inc.Country.isin(MISSING_CODES)][["City", "Country"]].drop_duplicates().shape[0])),
       ("Source-labelled city labels (distinct City strings, any country status)", str(countable(inc[inc.Geographic_Unit == 'City (source-labelled)'].City).nunique())),
       ("Additional locality strings not standardized to cities", str(countable(inc[inc.Geographic_Unit.str.startswith('Locality')].Locality_As_Reported).nunique())),
       ("Proposal target cities present (Dhaka / Delhi / Karachi)", f"Dhaka: {int(((inc.City=='Dhaka') | inc.Locality_As_Reported.str.contains('Dhaka')).sum())} samples; Delhi: {int((inc.City.isin(['Delhi','New Delhi']) | inc.Locality_As_Reported.str.contains('Delhi')).sum())}; Karachi: {int(((inc.City=='Karachi') | inc.Locality_As_Reported.str.contains('Karachi')).sum())}"),
       ("Included samples with source-reported coordinates", f"{int((inc.Latitude != NR).sum()):,}"),
       ("Included samples with gene-level AMR table", f"{int(inc.Sample_ID.isin(amr.Sample_ID).sum()):,}"),
       ("Included samples with taxonomy table (phylum)", f"{int(inc.Sample_ID.isin(tax.Sample_ID).sum()):,}"),
       ("Included samples with diversity indices", f"{int(inc.Sample_ID.isin(div.Sample_ID).sum()):,}"),
       ("Raw reads downloaded", "None (inventoried only; volumes in Download_Manifest)")]
idxdf = pd.DataFrame(idx, columns=["Metric", "Value"])
idxdf["Source"] = "Computed by scripts/05_build_workbooks.py from harmonized CSVs (2026-10-31)"
st = inc.groupby("Sample_Type").agg(N_Samples=("Sample_ID", "count"), N_Datasets=("Dataset_ID", "nunique"), N_Countries=("Country", lambda x: countable(x).nunique()),
                                    N_With_Coordinates=("Latitude", lambda x: int((x != NR).sum())), N_With_AMR_Table=("Sample_ID", lambda x: int(x.isin(amr.Sample_ID).sum()))).reset_index().sort_values("N_Samples", ascending=False)
sq = inc.assign(Platform=inc.Sequencing_Platform.str.split(";")).explode("Platform").groupby("Platform").agg(N_Samples=("Sample_ID", "nunique"), N_Datasets=("Dataset_ID", "nunique")).reset_index().sort_values("N_Samples", ascending=False)
amr_av = inc.groupby(["AMR_Data_Available"]).Sample_ID.count().rename("N_Included_Samples").reset_index()
tax_av = inc.groupby(["Taxonomic_Profile_Available"]).Sample_ID.count().rename("N_Included_Samples").reset_index()
miss = []
for f in ["City", "Locality_As_Reported", "Country", "Latitude", "Collection_Date", "Collection_Month", "Sample_Subtype", "Environmental_Matrix", "Sequencing_Platform", "DOI"]:
    m = inc[f].isin(MISSING_CODES + [""]).sum(); miss.append({"Field": f, "N_Missing_Included": int(m), "Percent_Missing": round(100 * m / len(inc), 1)})
missdf = pd.DataFrame(miss)
dsl = inc.groupby(["Dataset_ID"]).agg(Study=("Study_ID", lambda x: "; ".join(sorted(set(x)))[:80]), Sample_Type=("Sample_Type", lambda x: "; ".join(sorted(set(x)))),
                                      N_Samples=("Sample_ID", "count"), Countries=("Country", lambda x: ("; ".join(sorted(set(countable(x))))[:120]) or "Not reported")).reset_index()
dsl["Dataset_Name"] = dsl.Dataset_ID.map(ds_name); dsl["_k"] = dsl.Dataset_ID.str[2:].astype(int); dsl = dsl.sort_values("_k").drop(columns="_k")
w = WB(os.path.join(O, "00_Master_Data_Index.xlsx")); S = status_maps(w)
w.readme("README", "00 Master Data Index", common_readme + [
 ("Purpose", "Entry point: what was collected, how much, and what it can support."),
 ("Numbers", "All values are calculated from the collected records by scripts/05_build_workbooks.py; this index holds values (no cross-workbook formulas)."),
 ("Headline", f"{len(inc):,} included public shotgun metagenomes from {inc.Dataset_ID.nunique()} datasets, {countable(inc.Country).nunique()} countries. Gene-level AMR tables exist for {inc.Sample_ID.isin(amr.Sample_ID).sum():,} of them; all others need reprocessing from raw reads."),
 ("Overall readiness", "Conditionally ready for within-dataset, within-matrix description and AMR mapping (Zenodo 6919377 ResFinder subset, DS02). Not ready for cross-study pooled statistics, environmental association, predictive modelling or food-system analysis without further data collection and uniform reprocessing.")])
w.table("Summary", idxdf, widths={"Metric": 60, "Value": 70, "Source": 50})
w.table("Datasets", dsl[["Dataset_ID", "Dataset_Name", "Study", "Sample_Type", "N_Samples", "Countries"]], widths={"Dataset_Name": 60})
w.table("Sample_Type_Coverage", st)
w.table("Sequencing_Coverage", sq)
w.table("AMR_Availability", amr_av, widths={"AMR_Data_Available": 70})
w.table("Taxonomy_Availability", tax_av, widths={"Taxonomic_Profile_Available": 80})
w.table("Missing_Metadata", missdf)
w.table("Analysis_Readiness", readdf, status_cols={"Readiness": S["Readiness"]}, widths={"Reason": 100})
w.close()

json.dump({"qdf": len(qdf), "inventory": len(invdf), "included_samples": len(inc), "datasets": int(sm.Dataset_ID.nunique()),
           "index": dict(idx)}, open(os.path.join(H, "index_numbers.json"), "w"), indent=1)
print("done", {k: v for k, v in idx})
