"""Step 01 - Recursive inventory of the research directory.

For every file under SOURCE (excluding this output package and the git clone):
size, modification time, SHA-256, file type, integrity check (does the file open?),
duplicate group (identical SHA-256), and a content-based classification of its role
in the chapter. Writes 00_Project_Inventory/Source_File_Inventory.xlsx and a CSV copy.
Original files are opened read-only and never modified.
"""
import csv
import datetime as dt
import hashlib
import json
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from PIL import Image

from amr_paths import SOURCE, INV, EXCLUDE_DIRS, log

Image.MAX_IMAGE_PIXELS = None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def integrity(path: Path) -> str:
    """Light-weight check that the file can be opened by its format reader."""
    ext = path.suffix.lower()
    try:
        if path.stat().st_size == 0:
            return "EMPTY"
        if ext in {".xlsx", ".docx", ".pptx", ".xlsm"}:
            with zipfile.ZipFile(path) as z:
                bad = z.testzip()
            return "OK" if bad is None else f"CORRUPT member {bad}"
        if ext in {".png", ".jpg", ".jpeg", ".tif", ".tiff"}:
            with Image.open(path) as im:
                im.verify()
            return "OK"
        if ext in {".json", ".ipynb"}:
            json.loads(path.read_text(encoding="utf-8"))
            return "OK"
        if ext in {".csv", ".tsv", ".txt", ".md", ".py", ".js", ".log", ".sh"}:
            with open(path, encoding="utf-8") as fh:
                fh.read(65536)
            return "OK"
        if ext in {".tgz", ".gz"}:
            import gzip
            with gzip.open(path) as fh:
                fh.read(65536)
            return "OK"
        if ext == ".pdf":
            return "OK" if path.read_bytes()[:5] == b"%PDF-" else "NOT A PDF"
        return "NOT CHECKED (format)"
    except Exception as exc:  # report, never raise
        return f"ERROR {type(exc).__name__}: {str(exc)[:80]}"


FILETYPE = {
    ".xlsx": "Excel workbook", ".xls": "Excel workbook", ".csv": "CSV table", ".tsv": "TSV table",
    ".docx": "Word document", ".doc": "Word document", ".pdf": "PDF", ".pptx": "PowerPoint",
    ".png": "Image (PNG)", ".tif": "Image (TIFF)", ".tiff": "Image (TIFF)", ".jpg": "Image (JPEG)",
    ".py": "Python script", ".ipynb": "Jupyter notebook", ".js": "JavaScript", ".json": "JSON",
    ".md": "Markdown", ".txt": "Text", ".log": "Log", ".npy": "NumPy array", ".pyc": "Python bytecode",
    ".tgz": "Archive", ".gz": "Archive", ".fasta": "FASTA", ".xml": "XML", ".url": "URL record",
    ".geojson": "GeoJSON", ".sh": "Shell script", ".res": "KMA result", ".mapstat": "KMA mapstat",
}


def classify(rel: str) -> tuple:
    """Return (role, relevance, public_release_note) from the path; refined below by duplicates."""
    r = rel.replace("\\", "/")
    top = r.split("/")[0]
    low = r.lower()
    if low.endswith(".pyc") or "__pycache__" in low:
        return "Python bytecode cache", "Irrelevant (regenerable)", "Exclude"
    if "_to_delete" in low or low.endswith(".writetest"):
        return "Temporary / marked for deletion", "Irrelevant", "Exclude"
    if r == "Supplementary Material.docx":
        return "Supplementary Material (primary reference document)", "Core", "Publish final revised version"
    if top == "08_Final_Submission":
        if "Chapter_Draft" in r and r.endswith(".docx"):
            return "Chapter manuscript draft (Elsevier)", "Related - chapter text", "Withhold (publisher copyright)"
        if "QA_and_Completion" in r:
            return "Chapter QA report", "Related - process record", "Withhold (internal)"
        if "Supplementary_Data_Tables" in r:
            return "Supplementary data workbook (S1-S96, earlier build)", "Core - superseded by revised export", "Superseded"
        if "/Figures/" in r:
            return "Chapter main-text figure", "Related - chapter", "Withhold (main-text figure, publisher copyright)"
        if "Supplementary_Figures" in r:
            return "Early supplementary figure export (pre-v2 numbering)", "Superseded", "Superseded"
        if "Build_Sources" in r:
            if "manuscript.txt" in r or "refs.py" in r or "ref_verification" in r or "content.py" in r:
                return "Chapter build source (manuscript text/references)", "Related - chapter", "Withhold (chapter text)"
            return "Supplement build source (captions, interpretations, scripts)", "Core - provenance", "Internal provenance"
    if top == "Updated_Research_Analysis":
        if "/04_Tables/" in r:
            return "Verified analytical table (supplement source)", "Core", "Publish (as S-table export)"
        if "/05_Figures/" in r:
            return "Verified analytical figure (supplement source)", "Core", "Publish (as S-figure export)"
        if "/01_Data/" in r:
            return "Updated harmonized / reprocessed data", "Core", "Publish cleaned version"
        if "/02_Scripts/" in r or "/03_Notebooks/" in r:
            return "Verified analysis code / notebook", "Core - reproducibility", "Publish (archival copy)"
        if "/06_Logs/" in r:
            return "Processing log", "Supporting", "Internal"
        if "/07_Reports/" in r:
            return "Verification / Dhaka report", "Supporting", "Internal"
        if "00_Backup" in r:
            return "Checksum / backup record", "Supporting", "Internal"
        return "Verified analysis workbook / register", "Core", "Internal"
    if top == "AMR_Chapter_Reproducibility_Bundle":
        if "/raw_downloads/" in r:
            return "Third-party raw download (unmodified)", "Source data", "Withhold (third-party; re-download from source)"
        if "/harmonized/" in r:
            return "Harmonized public metadata / annotation table", "Core", "Publish cleaned version"
        if "/scripts/" in r:
            return "Data-collection script", "Core - reproducibility", "Publish (archival copy)"
        return "Data-collection log / output", "Supporting", "Internal"
    if top == "Analysis_Results":
        return "First analysis run (superseded by audit + update)", "Superseded", "Superseded"
    if top == "Independent_Audit":
        return "Independent audit record / corrected build", "Supporting - audit trail", "Internal"
    if top == "Figure_Table_Interpretations":
        return "Figure/table interpretation drafts", "Supporting", "Internal"
    if top == r:  # root-level file
        name = r
        if name.startswith(("00_", "01_", "02_", "03_", "04_")) and name.endswith(".xlsx"):
            return "Primary data-collection workbook", "Core", "Already public in repository"
        if name.startswith("Public_Data_Source_Log"):
            return "Query / download source log", "Core", "Already public in repository"
        if name.startswith("README_Data_Preparation"):
            return "Data-preparation README", "Core", "Already public (repository README)"
        if "Discussion" in name or "Abstract" in name:
            return "Chapter text draft (discussion/abstract)", "Related - chapter text", "Withhold (chapter text)"
        if name.startswith("Fig1_conceptual"):
            return "Chapter main-text figure (high-resolution)", "Related - chapter", "Withhold (main-text figure)"
        if name.startswith("Data_Collection_and_QC"):
            return "Data collection and QC report", "Supporting", "Internal"
        if name.startswith("AMR_") and name.endswith(".xlsx"):
            return "Analysis results workbook", "Supporting", "Internal"
    return "Other", "Unclassified", "Review"


def main():
    rows = []
    files = [p for p in SOURCE.rglob("*") if p.is_file()
             and p.relative_to(SOURCE).parts[0] not in EXCLUDE_DIRS]
    log(f"Inventory: {len(files)} files found under source directory")
    for i, p in enumerate(sorted(files)):
        rel = str(p.relative_to(SOURCE)).replace("\\", "/")
        st = p.stat()
        role, relevance, release = classify(rel)
        rows.append(dict(
            File_ID=f"F{i+1:05d}", Relative_Path=rel, Top_Folder=rel.split("/")[0] if "/" in rel else "(root)",
            File_Name=p.name, Extension=p.suffix.lower() or "(none)",
            File_Type=FILETYPE.get(p.suffix.lower(), "Other"), Size_Bytes=st.st_size,
            Size_MB=round(st.st_size / 1e6, 3),
            Modified=dt.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M"),
            SHA256=sha256(p), Integrity=integrity(p), Role=role, Relevance=relevance,
            Public_Release=release))
    df = pd.DataFrame(rows)

    # Duplicates: identical content (SHA-256). Canonical copy = path without "_1" suffix, then shortest path.
    df["Duplicate_Group"] = ""
    df["Duplicate_Status"] = "Unique"
    groups = df[df.duplicated("SHA256", keep=False) & (df.Size_Bytes > 0)].groupby("SHA256")
    for k, (h, g) in enumerate(groups, 1):
        order = g.assign(_p=g.File_Name.str.contains(r"_1\.\w+$") * 1, _l=g.Relative_Path.str.len()) \
                 .sort_values(["_p", "_l"])
        df.loc[g.index, "Duplicate_Group"] = f"D{k:04d}"
        df.loc[order.index[0], "Duplicate_Status"] = "Duplicate - canonical copy"
        df.loc[order.index[1:], "Duplicate_Status"] = "Duplicate - redundant copy"
    df["Size_MB"] = df["Size_MB"].astype(float)

    INV.mkdir(parents=True, exist_ok=True)
    df.to_csv(INV / "Source_File_Inventory.csv", index=False, quoting=csv.QUOTE_MINIMAL)

    summ = []
    summ.append(("Source directory", "(research directory; relative paths below)"))
    summ.append(("Inventory date", dt.date.today().isoformat()))
    summ.append(("Files inspected", len(df)))
    summ.append(("Total size (GB)", round(df.Size_Bytes.sum() / 1e9, 3)))
    summ.append(("Files with integrity problems", int((~df.Integrity.isin(["OK", "NOT CHECKED (format)"])).sum())))
    summ.append(("Duplicate groups (identical SHA-256)", df.Duplicate_Group.replace("", pd.NA).nunique()))
    summ.append(("Redundant duplicate copies", int((df.Duplicate_Status == "Duplicate - redundant copy").sum())))
    summary = pd.DataFrame(summ, columns=["Item", "Value"])
    by_folder = df.groupby("Top_Folder").agg(Files=("File_ID", "count"), Size_MB=("Size_MB", "sum")).reset_index()
    by_type = df.groupby("File_Type").agg(Files=("File_ID", "count"), Size_MB=("Size_MB", "sum")).reset_index() \
                .sort_values("Files", ascending=False)
    by_role = df.groupby(["Relevance", "Role", "Public_Release"]).agg(Files=("File_ID", "count"),
                                                                     Size_MB=("Size_MB", "sum")).reset_index()
    dups = df[df.Duplicate_Group != ""].sort_values(["Duplicate_Group", "Duplicate_Status"])
    probs = df[~df.Integrity.isin(["OK", "NOT CHECKED (format)"])]
    core = df[df.Relevance.str.startswith("Core")]

    out = INV / "Source_File_Inventory.xlsx"
    with pd.ExcelWriter(out, engine="openpyxl") as xw:
        summary.to_excel(xw, sheet_name="Summary", index=False)
        by_folder.to_excel(xw, sheet_name="By_Folder", index=False)
        by_type.to_excel(xw, sheet_name="By_File_Type", index=False)
        by_role.to_excel(xw, sheet_name="By_Role_and_Release", index=False)
        core.to_excel(xw, sheet_name="Core_Files", index=False)
        dups.to_excel(xw, sheet_name="Duplicates", index=False)
        probs.to_excel(xw, sheet_name="Integrity_Issues", index=False)
        df.to_excel(xw, sheet_name="All_Files", index=False)
    _format(out)
    log(f"Inventory written: {len(df)} files, {summary.Value[5]} duplicate groups, "
        f"{summary.Value[4]} integrity issues -> {out.name}")
    print(summary.to_string(index=False))
    print(by_role.to_string(index=False))


def _format(path):
    from openpyxl import load_workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    wb = load_workbook(path)
    for ws in wb.worksheets:
        for c in ws[1]:
            c.font = Font(bold=True, color="FFFFFF")
            c.fill = PatternFill("solid", fgColor="1F4E78")
            c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = "A2"
        if ws.max_row > 1:
            ws.auto_filter.ref = ws.dimensions
        for col in ws.columns:
            width = max(len(str(c.value or "")) for c in list(col)[:300])
            ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 70)
    wb.save(path)


if __name__ == "__main__":
    main()
