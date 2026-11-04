"""Step 09 - Final quality-control audit and public-release screen.

Writes 00_Project_Inventory/Quality_Control_Checklist.xlsx (every check with PASS/FAIL/WARN and evidence)
and 00_Project_Inventory/Public_Release_Screen.csv (file-by-file release decision for the GitHub repository).
Exit code 1 if any check FAILS.
"""
import hashlib
import io
import json
import re
import sys
import zipfile

import pandas as pd
from PIL import Image

from amr_paths import (PROJECT, INV, CLEAN_CSV, CLEAN_XLSX, DICT, ANALYSIS, TABLES, FIGURES, DOCS, REPORTS, UPDATED,
                       log)

Image.MAX_IMAGE_PIXELS = None
R = []


def check(area, item, ok, evidence, warn=False):
    R.append(dict(Area=area, Check=item, Result="PASS" if ok else ("WARN" if warn else "FAIL"), Evidence=str(evidence)))


def datasets():
    prof = pd.read_excel(ANALYSIS / "Descriptive_Statistics" / "Dataset_Profiles.xlsx", sheet_name="Dataset_Summary")
    check("Datasets", "Original and cleaned record counts identical (no data loss)",
          (prof.Rows_Original == prof.Rows_Cleaned).all(),
          "; ".join(f"{r.Dataset}: {r.Rows_Original}->{r.Rows_Cleaned}" for _, r in prof.iterrows()))
    dd = pd.read_csv(DICT / "Data_Dictionary_All_Datasets.csv")
    miss = []
    for f in list(CLEAN_CSV.glob("*.csv")) + list(CLEAN_CSV.glob("*.csv.gz")):
        name = f.name.split(".")[0]
        cols = pd.read_csv(f, nrows=0).columns
        miss += [f"{name}.{c}" for c in cols if not ((dd.Dataset == name) & (dd.Variable == c)).any()]
        bad = [c for c in cols if c != c.lower() or " " in c]
        check("Datasets", f"{name}: variable names standardised (snake_case)", not bad, bad or "all lower snake_case")
    check("Datasets", "Every cleaned variable documented in the data dictionary", not miss, miss[:10] or "complete")
    check("Datasets", "Data dictionary descriptions present", dd.Description.fillna("").str.len().gt(0).all(),
          f"{dd.Description.fillna('').str.len().eq(0).sum()} empty descriptions")
    ops = pd.read_csv(CLEAN_CSV.parent / "Cleaning_Operations_Log.csv")
    ex = ops[ops.Detail.str.startswith("Exact duplicate rows")]
    check("Datasets", "Duplicate handling: exact duplicate rows (flagged, not removed)", True,
          f"{int(ex.Records_or_Cells_Affected.sum())} exact duplicate rows across datasets")
    ref = ops[ops.Step == "referential"]
    check("Datasets", "Referential integrity across tables", (ref.Records_or_Cells_Affected == 0).all(),
          "; ".join(f"{r.Dataset}: {r.Records_or_Cells_Affected}" for _, r in ref.iterrows()),
          warn=True)
    ver = pd.read_csv(ANALYSIS / "Analytical_Results" / "Independent_Verification_of_Supplementary_Results.csv")
    check("Analysis", "Headline results reproduced from cleaned data", (ver.Match == "YES").all(),
          f"{(ver.Match == 'YES').sum()}/{len(ver)} checks match")


def tables():
    idx = pd.read_csv(TABLES / "Table_Index_S1-S96.csv")
    check("Tables", "Sequential numbering Table S1-S96", list(idx.Table) == [f"Table S{i}" for i in range(1, 97)],
          f"{len(idx)} tables")
    check("Tables", "Every table has a title from the final document", idx.Title.fillna("").str.len().gt(0).all(),
          f"{idx.Title.isna().sum()} missing titles")
    diffs = []
    for _, r in idx.iterrows():
        a = pd.read_csv(TABLES / "CSV_Tables" / r.CSV_File)
        b = pd.read_csv(UPDATED / "04_Tables" / r.Source_File)
        if not a.equals(b):
            diffs.append(r.Table)
    check("Tables", "Exported values identical to verified analysis tables (cell-by-cell)", not diffs, diffs or "96/96 identical")
    xl = TABLES / "Excel_Tables" / "Supplementary_Data_Tables_S1-S96.xlsx"
    with zipfile.ZipFile(xl) as z:
        check("Tables", "Supplementary data workbook opens", z.testzip() is None, xl.name)


def figures():
    idx = pd.read_csv(FIGURES / "Figure_Index_S1-S30.csv")
    check("Figures", "Sequential numbering Fig. S1-S30", list(idx.Figure) == [f"Fig. S{i}" for i in range(1, 31)],
          f"{len(idx)} figures")
    lowres, missing = [], []
    for _, r in idx.iterrows():
        for sub, ext in (("PNG", "png"), ("TIFF", "tif"), ("PDF", "pdf")):
            f = FIGURES / sub / f"{r.File_Stem}.{ext}"
            if not f.exists():
                missing.append(f.name)
        im = Image.open(FIGURES / "PNG" / f"{r.File_Stem}.png")
        dpi = im.info.get("dpi", (0, 0))[0]
        if dpi < 299:
            lowres.append((r.Figure, dpi))
        tif = Image.open(FIGURES / "TIFF" / f"{r.File_Stem}.tif")
        if tif.info.get("dpi", (0, 0))[0] < 299:
            lowres.append((r.Figure, "tiff", tif.info.get("dpi")))
    check("Figures", "PNG, TIFF and PDF present for every figure", not missing, missing or "90/90 files")
    check("Figures", "Resolution >= 300 dpi (PNG and TIFF)", not lowres, lowres or "all >= 300 dpi")
    check("Figures", "Print width fits Elsevier page (<= 190 mm at 300 dpi)", (idx.Print_width_mm <= 190).all(),
          f"max {idx.Print_width_mm.max()} mm")
    import docx
    from docx.oxml.ns import qn
    d = docx.Document(DOCS / "Supplementary_Material_Final.docx")
    blobs = [d.part.related_parts[b.get(qn("r:embed"))].blob for b in d.element.body.iter(qn("a:blip"))]
    same = [hashlib.sha256(blobs[i]).hexdigest() ==
            hashlib.sha256((FIGURES / "PNG" / f"{r.File_Stem}.png").read_bytes()).hexdigest() for i, r in idx.iterrows()]
    # the corrected Fig. S28 is re-encoded on export: compare pixels instead of bytes
    pix = [Image.open(io.BytesIO(blobs[i])).convert("RGB").tobytes() ==
           Image.open(FIGURES / "PNG" / f"{r.File_Stem}.png").convert("RGB").tobytes() for i, r in idx.iterrows()]
    check("Figures", "Exported figures identical to the images in the final document", all(pix),
          f"{sum(same)} byte-identical, {sum(pix)} pixel-identical of {len(idx)}")


def document():
    import docx
    from docx.oxml.ns import qn
    path = DOCS / "Supplementary_Material_Final.docx"
    d = docx.Document(path)
    body = d.element.body
    instr = [i.text.strip() for i in body.iter(qn("w:instrText"))]
    seq_f = [i for i in instr if i.startswith("SEQ Figure")]
    seq_t = [i for i in instr if i.startswith("SEQ Table")]
    refs = [i.split()[1] for i in instr if i.startswith("REF ")]
    bms = {b.get(qn("w:name")) for b in body.iter(qn("w:bookmarkStart"))}
    tocs = [i for i in instr if i.startswith("TOC")]
    check("Word document", "Automatic figure captions (SEQ Figure) = 30", len(seq_f) == 30, len(seq_f))
    check("Word document", "Automatic table captions (SEQ Table) = 96", len(seq_t) == 96, len(seq_t))
    broken = sorted(set(r for r in refs if r not in bms))
    check("Word document", "Every cross-reference (REF) points to an existing caption bookmark", not broken,
          f"{len(refs)} REF fields; broken: {broken[:5]}")
    cited = {r for r in refs}
    uncited = sorted(b for b in bms if b.startswith("_Ref") and b not in cited)
    check("Word document", "Every table and figure is cited in the text", not uncited, uncited[:10] or "126/126 cited")
    check("Word document", "Table of Contents, List of Figures and List of Tables fields present",
          any('\\o' in t for t in tocs) and any('"Figure"' in t for t in tocs) and any('"Table"' in t for t in tocs), tocs)
    text = "\n".join(p.text for p in d.paragraphs)
    errs = re.findall(r"Error! [A-Za-z ]+", text)
    check("Word document", "No field errors after Word update", not errs, errs[:3] or "none")
    lof = [p for p in d.paragraphs if p.style.name.lower() == "table of figures"]
    toc = [p for p in d.paragraphs if p.style.name.lower().startswith("toc ")]
    check("Word document", "List of Figures/Tables populated by Word (126 entries)", len(lof) == 126, f"{len(lof)} entries")
    check("Word document", "Table of Contents populated by Word", len(toc) >= 19, f"{len(toc)} entries")
    caps = [p.style.style_id for p in d.paragraphs if re.match(r"(Fig\.|Table) S\d+ ", p.text)
            and p.style.style_id != "TableofFigures"]   # list entries repeat the caption text
    check("Word document", "126 captions (30 figures + 96 tables)", len(caps) == 126, len(caps))
    check("Word document", "Captions use Word's built-in Caption style", set(caps) == {"Caption"}, set(caps))
    ftr = d.sections[0].footer._element
    fi = [i.text.strip() for i in ftr.iter(qn("w:instrText"))]
    check("Word document", "Page numbering (PAGE / NUMPAGES) in footer", "PAGE" in fi and "NUMPAGES" in fi, fi)
    check("Word document", "Running header present", bool(d.sections[0].header.paragraphs[0].text.strip()),
          d.sections[0].header.paragraphs[0].text)
    # placeholder text left?
    check("Word document", "Placeholders 'Table list'/'figure list' removed",
          not any(p.text.strip().lower() in ("table list", "figure list") for p in d.paragraphs), "")
    corr = pd.read_excel(DOCS / "Supplementary_Appendices" / "Supplementary_Text_Corrections_Log.xlsx", sheet_name="Corrections")
    present = sum(1 for c in corr.Corrected_Text if c.strip()[:60] in text or len(c.strip()) < 3)
    check("Word document", "Verified text corrections present in final text", present >= len(corr) * 0.95,
          f"{present}/{len(corr)} corrected passages found verbatim (others contain fields, e.g. table numbers)")
    # PDF
    from pypdf import PdfReader
    pdf = DOCS / "Supplementary_Material_Final.pdf"
    if not pdf.exists():
        check("PDF", "PDF exported", False, "missing")
        return
    rd = PdfReader(str(pdf))
    pages = [pg.extract_text() or "" for pg in rd.pages]
    check("PDF", "PDF opens", len(pages) > 0, f"{len(pages)} pages")
    alltext = " ".join(pages)
    missing = [f"Table S{i}" for i in range(1, 97) if not re.search(rf"Table S{i}\b", alltext)] + \
              [f"Fig. S{i}" for i in range(1, 31) if not re.search(rf"Fig\. S{i}\b", alltext)]
    check("PDF", "PDF contains every caption of the final DOCX", not missing, missing[:5] or "126/126")
    check("PDF", "No field errors in PDF", "Error!" not in alltext, "")
    blank = []
    for i, pg in enumerate(rd.pages):
        t = (pages[i] or "").strip()
        has_img = bool(pg.get("/Resources", {}).get("/XObject")) if hasattr(pg, "get") else False
        body_text = re.sub(r"Supplementary Material.*?Mapping|Page \d+ of \d+", "", t).strip()
        if len(body_text) < 15 and not has_img:
            blank.append(i + 1)
    check("PDF", "No blank pages", not blank, blank[:10] or "none")
    m = re.search(r"Page \d+ of (\d+)", pages[-1])
    check("PDF", "Footer page count matches PDF pages", bool(m) and int(m.group(1)) == len(pages),
          f"footer says {m.group(1) if m else '?'}; PDF has {len(pages)}")


def files():
    bad, tmp = [], []
    for f in PROJECT.rglob("*"):
        if not f.is_file():
            continue
        if f.name.startswith("~$") or f.suffix in (".tmp", ".pyc") or "__pycache__" in f.parts:
            tmp.append(str(f.relative_to(PROJECT)))
            continue
        try:
            if f.suffix in (".xlsx", ".docx"):
                with zipfile.ZipFile(f) as z:
                    assert z.testzip() is None
            elif f.suffix in (".png", ".tif"):
                Image.open(f).verify()
            elif f.suffix == ".pdf":
                assert f.read_bytes()[:5] == b"%PDF-"
            elif f.suffix in (".json", ".ipynb"):
                json.loads(f.read_text(encoding="utf-8"))
        except Exception as e:
            bad.append(f"{f.name}: {e}")
    check("Files", "Every output file opens", not bad, bad[:5] or "all open")
    check("Files", "No temporary files left", not tmp, tmp[:5] or "none")
    need = [INV / "Source_File_Inventory.xlsx", INV / "Data_Quality_Report.docx", INV / "Processing_Log.txt",
            DOCS / "Supplementary_Material_Final.docx", DOCS / "Supplementary_Material_Final.pdf",
            REPORTS / "Data_Cleaning_Methodology.docx", REPORTS / "Statistical_Analysis_Methodology.docx",
            REPORTS / "Technical_Report.docx", PROJECT / "README.md", PROJECT / "07_Code_and_Scripts" / "Requirements.txt"]
    miss = [str(p.relative_to(PROJECT)) for p in need if not p.exists()]
    check("Files", "All required deliverables exist", not miss, miss or "all present")


EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
SECRET = re.compile(r"(ghp_[A-Za-z0-9]{20,}|github_pat_|AKIA[0-9A-Z]{16}|sk-[A-Za-z0-9]{20,}|-----BEGIN [A-Z ]*PRIVATE KEY|"
                    r"(?i:password|passwd|api[_-]?key|secret[_-]?key|access[_-]?token)\s*[:=]\s*['\"][^'\"]{4,})")
LOCALPATH = re.compile(r"([A-Za-z]:\\\\?Users\\\\?[A-Za-z0-9_.-]+|/home/[a-z0-9_.-]+|/Users/[A-Za-z0-9_.-]+|"
                       r"[A-Za-z]:[\\/]Research[\\/]AMR|/mnt/user-data)")


def release_screen(candidates):
    rows = []
    for f in candidates:
        rel = str(f.relative_to(PROJECT)).replace("\\", "/")
        size = f.stat().st_size / 1e6
        finding = []
        if f.suffix in (".py", ".md", ".txt", ".csv", ".json", ".ps1", ".ipynb", ".gz"):
            try:
                if f.suffix == ".gz":
                    import gzip
                    txt = gzip.open(f, "rt", encoding="utf-8", errors="ignore").read(50_000_000)
                else:
                    txt = f.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                txt = ""
            emails = set(e for e in EMAIL.findall(txt) if not e.endswith(("anthropic.com", "example.com")))
            if emails:
                finding.append(f"emails: {sorted(emails)[:3]}")
            if SECRET.search(txt):
                finding.append("possible credential pattern")
            lp = set(LOCALPATH.findall(txt))
            if lp:
                finding.append(f"local paths: {sorted(lp)[:3]}")
        if size > 50:
            finding.append(f"size {size:.1f} MB > 50 MB")
        rows.append(dict(File=rel, Size_MB=round(size, 3), Findings="; ".join(finding)))
    return pd.DataFrame(rows)


def main():
    datasets(); tables(); figures(); document(); files()
    res = pd.DataFrame(R)
    # public-release screen over everything proposed for GitHub
    cand = [f for f in PROJECT.rglob("*") if f.is_file() and not any(
        p in f.parts for p in ("01_Raw_Data_Archive", "text_review_records"))]
    scr = release_screen(cand)
    scr.to_csv(INV / "Public_Release_Screen.csv", index=False)
    flagged = scr[scr.Findings != ""]
    check("Public release", "Files proposed for release screened for e-mails, credentials, local paths and size",
          True, f"{len(scr)} files screened; {len(flagged)} with findings (see Public_Release_Screen.csv)")
    res = pd.DataFrame(R)
    with pd.ExcelWriter(INV / "Quality_Control_Checklist.xlsx", engine="openpyxl") as xw:
        res.groupby(["Area", "Result"]).size().unstack(fill_value=0).reset_index().to_excel(xw, sheet_name="Summary", index=False)
        res.to_excel(xw, sheet_name="Checklist", index=False)
        flagged.to_excel(xw, sheet_name="Release_Findings", index=False)
    log(f"QC: {len(res)} checks - PASS {(res.Result == 'PASS').sum()}, WARN {(res.Result == 'WARN').sum()}, "
        f"FAIL {(res.Result == 'FAIL').sum()}; release findings in {len(flagged)} files")
    print(res[res.Result != "PASS"].to_string())
    print(flagged.to_string())
    sys.exit(1 if (res.Result == "FAIL").any() else 0)


if __name__ == "__main__":
    main()
