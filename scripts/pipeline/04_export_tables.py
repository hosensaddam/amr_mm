"""Step 04 - Export the 96 supplementary tables (unabridged) as CSV and as a formatted Excel workbook.

The S-number -> analysis-table mapping is read from the supplement build script (TORDER in
08_Final_Submission/Build_Sources/v2/build_supp.py), titles and notes from the final Word document,
so the exported tables carry exactly the captions used in Supplementary_Material_Final.docx.
Values are copied from the verified CSV files without modification.
"""
import re

import pandas as pd
from docx import Document
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side

from amr_paths import BUILD_V2, UPDATED, TABLES, DOCS, log


def torder():
    src = (BUILD_V2 / "build_supp.py").read_text(encoding="utf-8")
    return eval(re.search(r"TORDER = (\[.*?\n\])", src, re.S).group(1))


def table_map():
    rows, n = [], 0
    files = sorted((UPDATED / "04_Tables").glob("*.csv"))
    for group, ids in torder():
        for tid in ids:
            n += 1
            f = [p for p in files if p.name.startswith(tid + "_")]
            assert len(f) == 1, (tid, f)
            rows.append(dict(Table=f"Table S{n}", S=n, Section=group, Analysis_ID=tid, Source_File=f[0].name,
                             path=f[0]))
    return rows


def captions_from_docx():
    """Title and notes for each table from the final document (after Word updated the fields)."""
    doc = Document(DOCS / "Supplementary_Material_Final.docx")
    P = doc.paragraphs
    out = {}
    for i, p in enumerate(P):
        m = re.match(r"Table S(\d+)\s+(.*)", p.text.strip())
        if m and p.style.style_id == "Caption":
            n = int(m.group(1))
            j = i + 1
            while j < len(P) and (not P[j].text.strip() or re.fullmatch(r"Part \d+ of \d+", P[j].text.strip())):
                j += 1
            out[n] = (m.group(2).strip(), P[j].text.strip() if j < len(P) else "")
    return out


def main():
    csv_dir, xl_dir = TABLES / "CSV_Tables", TABLES / "Excel_Tables"
    csv_dir.mkdir(parents=True, exist_ok=True)
    xl_dir.mkdir(parents=True, exist_ok=True)
    caps = captions_from_docx()
    rows = table_map()
    index = []
    xl = xl_dir / "Supplementary_Data_Tables_S1-S96.xlsx"
    with pd.ExcelWriter(xl, engine="openpyxl") as xw:
        pd.DataFrame().to_excel(xw, sheet_name="Index")
        for r in rows:
            df = pd.read_csv(r["path"])
            title, notes = caps.get(r["S"], ("", ""))
            stem = f"Table_S{r['S']:02d}_{r['Analysis_ID']}"
            df.to_csv(csv_dir / f"{stem}.csv", index=False)
            sheet = f"S{r['S']}"
            pd.DataFrame([[f"{r['Table']}. {title}"], [notes], [f"Source: {r['Source_File']} (verified analysis, 31 Oct 2026)"]]) \
                .to_excel(xw, sheet_name=sheet, index=False, header=False)
            df.to_excel(xw, sheet_name=sheet, index=False, startrow=4)
            index.append(dict(Table=r["Table"], Sheet=sheet, Section=r["Section"], Title=title, Rows=len(df),
                              Columns=df.shape[1], Analysis_ID=r["Analysis_ID"], Source_File=r["Source_File"],
                              CSV_File=f"{stem}.csv"))
        pd.DataFrame(index).to_excel(xw, sheet_name="Index", index=False)
    # formatting
    wb = load_workbook(xl)
    thin = Side(style="thin", color="000000")
    for ws in wb.worksheets:
        if ws.title == "Index":
            hdr_row = 1
            ws.column_dimensions["D"].width = 90
        else:
            ws["A1"].font = Font(bold=True, size=12)
            ws["A2"].font = Font(italic=True, size=9)
            ws["A3"].font = Font(size=9, color="555555")
            hdr_row = 5
        for c in ws[hdr_row]:
            c.font = Font(bold=True, color="FFFFFF")
            c.fill = PatternFill("solid", fgColor="1F4E78")
            c.alignment = Alignment(wrap_text=True, vertical="top")
            c.border = Border(top=thin, bottom=thin)
        ws.freeze_panes = ws.cell(row=hdr_row + 1, column=1)
        for col in ws.iter_cols(min_row=hdr_row, max_row=min(ws.max_row, hdr_row + 200)):
            w = max(len(str(c.value)) if c.value is not None else 0 for c in col)
            letter = col[0].column_letter
            if ws.title != "Index" or letter != "D":
                ws.column_dimensions[letter].width = min(max(9, w + 2), 45)
    wb.save(xl)
    pd.DataFrame(index).to_csv(TABLES / "Table_Index_S1-S96.csv", index=False)
    missing = [i["Table"] for i in index if not i["Title"]]
    log(f"Tables exported: {len(index)} CSV + workbook; tables without caption match: {missing}")


if __name__ == "__main__":
    main()
