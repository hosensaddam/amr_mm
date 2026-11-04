"""Step 05 - Export Figures S1-S30 as PNG (300 dpi), TIFF (LZW) and PDF, with a figure index.

Figures are not re-drawn: the verified analysis figures (Updated_Research_Analysis/05_Figures, PNG 300 dpi
and TIFF 600 dpi) are the ones embedded in the supplement. The image actually embedded in the final
document is used whenever it differs from the analysis PNG (Fig. S28 = Fig24 with corrected title), so the
exported files always match the document. PDFs embed the full-resolution raster (no resampling).
Plotting code: 07_Code_and_Scripts/Original_Analysis_Code (notebooks listed in the index).
"""
import hashlib
import io
import re
import shutil
import sys
import zipfile

import pandas as pd
from PIL import Image

from amr_paths import BUILD_V2, UPDATED, FIGURES, DOCS, log

Image.MAX_IMAGE_PIXELS = None


def sfigs():
    sys.path.insert(0, str(BUILD_V2))
    from supp_figs import SFIGS  # (fid, group, caption, source tables)
    return SFIGS


def docx_images():
    """Embedded images of the final DOCX in document order -> bytes."""
    import docx
    from docx.oxml.ns import qn
    d = docx.Document(DOCS / "Supplementary_Material_Final.docx")
    out = []
    for blip in d.element.body.iter(qn("a:blip")):
        rid = blip.get(qn("r:embed"))
        out.append(d.part.related_parts[rid].blob)
    return out


def main():
    for sub in ("PNG", "TIFF", "PDF", "Editable_Figures"):
        (FIGURES / sub).mkdir(parents=True, exist_ok=True)
    cat = pd.read_csv(UPDATED / "04_Tables" / "R06u_updated_figure_catalog.csv")
    figs = sfigs()
    embedded = docx_images()
    assert len(embedded) == len(figs) == 30, (len(embedded), len(figs))
    fdir = UPDATED / "05_Figures"
    rows = []
    for k, (fid, group, cap, srcs) in enumerate(figs, 1):
        png_src = next(fdir.glob(f"{fid}_*.png"))
        tif_src = next((fdir / "tiff").glob(f"{fid}_*.tif"))
        stem = f"Fig_S{k:02d}_{png_src.stem}"
        emb = embedded[k - 1]
        same = hashlib.sha256(emb).hexdigest() == hashlib.sha256(png_src.read_bytes()).hexdigest()
        im = Image.open(io.BytesIO(emb)).convert("RGB")
        if same:
            shutil.copyfile(png_src, FIGURES / "PNG" / f"{stem}.png")
            shutil.copyfile(tif_src, FIGURES / "TIFF" / f"{stem}.tif")
            tif_dpi = Image.open(tif_src).info.get("dpi", (600, 600))[0]
            note = "Verified analysis figure (identical to document image)"
        else:  # document carries a corrected version: export from it
            im.save(FIGURES / "PNG" / f"{stem}.png", dpi=(300, 300), optimize=True)
            im.save(FIGURES / "TIFF" / f"{stem}.tif", dpi=(300, 300), compression="tiff_lzw")
            tif_dpi = 300
            note = "Corrected version embedded in the supplement (title fixed); exported from the document"
        im.save(FIGURES / "PDF" / f"{stem}.pdf", "PDF", resolution=300.0)
        w, h = im.size
        c = cat[cat.file == png_src.name]
        rows.append(dict(Figure=f"Fig. S{k}", Section=group, Analysis_ID=fid, File_Stem=stem,
                         Width_px=w, Height_px=h, PNG_dpi=300, TIFF_dpi=round(tif_dpi),
                         Print_width_mm=round(w / 300 * 25.4, 1), Print_height_mm=round(h / 300 * 25.4, 1),
                         Min_font_pt=c.min_font_pt.iloc[0] if len(c) else "",
                         Notebook=c.notebook.iloc[0] if len(c) else "", Source_Tables=srcs,
                         Caption_Verified=cap, Note=note))
    idx = pd.DataFrame(rows)
    idx.to_csv(FIGURES / "Figure_Index_S1-S30.csv", index=False)
    idx.to_excel(FIGURES / "Figure_Index_S1-S30.xlsx", index=False)
    # Editable sources: the percent-format notebook sources that draw each figure
    nb_src = UPDATED / "03_Notebooks" / "notebook_sources"
    used = sorted(set(idx.Notebook) - {""})
    for nb in used:
        num = re.sub(r"\D", "", nb)
        for f in nb_src.glob(f"{num}_*.py"):
            shutil.copyfile(f, FIGURES / "Editable_Figures" / f.name)
    (FIGURES / "Editable_Figures" / "README.txt").write_text(
        "Figures S1-S30 are produced by the analysis notebooks listed in ../Figure_Index_S1-S30.csv (column Notebook).\n"
        "The percent-format notebook sources in this folder are the editable plotting code; executed notebooks\n"
        "are in 07_Code_and_Scripts/Original_Analysis_Code. Re-running them regenerates the PNG (300 dpi) and\n"
        "TIFF (600 dpi) files from the verified tables in 04_Tables/CSV_Tables.\n", encoding="utf-8")
    log(f"Figures exported: {len(idx)} x (PNG, TIFF, PDF); from document image: "
        f"{(idx.Note.str.startswith('Corrected')).sum()}; min print width {idx.Print_width_mm.min()} mm")


if __name__ == "__main__":
    main()
