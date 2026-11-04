"""Step 06 - Build Supplementary_Material_Final.docx from the author's Supplementary Material.docx.

What this script does (the source document is opened read-only and never modified):
 1. Applies the verified text corrections listed in
    06_Research_Documents/Supplementary_Appendices/Supplementary_Text_Corrections_Log.xlsx
    (each correction restores the meaning of the verified, author-approved interpretation).
 2. Rebuilds every caption label as an automatic Word caption: "Fig. S{SEQ Figure}" / "Table S{SEQ Table}",
    with a hidden bookmark (_RefFigSn / _RefTabSn) around the number.
 3. Converts every in-text mention of a supplementary table or figure ("Table S12", "Tables S3 and S4",
    "Fig. S5") into a REF cross-reference field (hyperlinked) to that bookmark.
 4. Replaces the "Table list" / "figure list" placeholders with a Table of Contents, List of Figures
    and List of Tables (TOC fields), adds numbered section headings, a running header, a
    "Page X of Y" footer and document properties.
Fields are inserted with placeholder results; 07_update_fields_and_export_pdf.ps1 then lets Microsoft Word
update all fields and export the PDF.
"""
import copy
import json
import re
import shutil
import sys

import pandas as pd
from docx import Document
from docx.enum.text import WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

from amr_paths import SUPP_DOCX_SOURCE, DOCS, BUILD_V2, log

CHAPTER = "Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping"
OUT = DOCS / "Supplementary_Material_Final.docx"
CORR = DOCS / "Supplementary_Appendices" / "Supplementary_Text_Corrections_Log.xlsx"
N_FIG, N_TAB = 30, 96

INTRO = ("This supplement accompanies the chapter “" + CHAPTER + "” (Academic Press, Elsevier). "
         "It reports every accepted statistical figure (Figs. S1–S30) and statistical table (Tables S1–S96) "
         "from the chapter’s verified analysis (31 October 2026), each followed by an interpretation. "
         "Table values were generated programmatically from the verified analysis files without manual transcription; "
         "tables longer than 120 rows are abridged here and given in full in the accompanying Supplementary Data "
         "workbook (one sheet per table). Streams 1–3 use different pipelines, reference databases and abundance "
         "units (RPKM, FPKM, FPKM-16S) and must not be compared in magnitude. Detection of a resistance gene indicates "
         "DNA presence, not phenotypic resistance, expression, viability, exposure or transmission.")

_bm_id = [9000]


def next_bm():
    _bm_id[0] += 1
    return str(_bm_id[0])


# ------------------------------------------------------------------ low-level XML helpers
def make_run(text, rpr=None):
    r = OxmlElement("w:r")
    if rpr is not None:
        r.append(copy.deepcopy(rpr))
    t = OxmlElement("w:t")
    t.set(qn("xml:space"), "preserve")
    t.text = text
    r.append(t)
    return r


def field_runs(instr, result, rpr=None):
    """Complex field: begin / instrText / separate / result / end (list of w:r elements)."""
    out = []
    for kind in ("begin", None, "separate"):
        r = OxmlElement("w:r")
        if rpr is not None:
            r.append(copy.deepcopy(rpr))
        if kind is None:
            it = OxmlElement("w:instrText")
            it.set(qn("xml:space"), "preserve")
            it.text = f" {instr} "
            r.append(it)
        else:
            fc = OxmlElement("w:fldChar")
            fc.set(qn("w:fldCharType"), kind)
            if kind == "begin":
                fc.set(qn("w:dirty"), "true")  # ask Word to refresh on open
            r.append(fc)
        out.append(r)
    out.append(make_run(result, rpr))
    r = OxmlElement("w:r")
    if rpr is not None:
        r.append(copy.deepcopy(rpr))
    fc = OxmlElement("w:fldChar")
    fc.set(qn("w:fldCharType"), "end")
    r.append(fc)
    out.append(r)
    return out


def ptext(p):
    return "".join(t.text or "" for t in p._p.iter(qn("w:t")))


def text_runs(p):
    """Text-carrying runs of the paragraph, including runs nested in hyperlinks (document order)."""
    return [r for r in p._p.iter(qn("w:r")) if r.find(qn("w:t")) is not None]


def replace_in_paragraph(p, old, new):
    """Replace `old` by `new` across runs, keeping the formatting of the run where the match starts."""
    runs = text_runs(p)
    texts = [r.find(qn("w:t")).text or "" for r in runs]
    full = "".join(texts)
    i = full.find(old)
    if i < 0:
        return False
    j = i + len(old)
    pos = 0
    first = True
    for r, t in zip(runs, texts):
        a, b = pos, pos + len(t)
        pos = b
        if b <= i or a >= j:
            continue
        s, e = max(i, a) - a, min(j, b) - a
        tnode = r.find(qn("w:t"))
        tnode.set(qn("xml:space"), "preserve")
        tnode.text = t[:s] + (new if first else "") + t[e:]
        first = False
    return True


def set_paragraph_text(p, new):
    runs = text_runs(p)
    if not runs:
        p.add_run(new)
        return
    runs[0].find(qn("w:t")).text = new
    runs[0].find(qn("w:t")).set(qn("xml:space"), "preserve")
    for r in runs[1:]:
        r.getparent().remove(r)


# ------------------------------------------------------------------ caption and cross-reference conversion
CAP_RE = re.compile(r"^\s*(Fig\.|Table)\s+S(\d+)\b\s*")
REF_GROUP = re.compile(r"\b(Tables?|Figs?\.|Figures?)\s+(S\d+(?:\s*(?:,|and|&|–|-|to)\s*S\d+)*)")
SNUM = re.compile(r"S(\d+)")


def bookmark_name(kind, n):
    return f"_Ref{'Fig' if kind == 'Fig' else 'Tab'}S{n}"


def convert_caption(p, kind, n, doc_caption_style="Caption"):
    """Turn the leading 'Fig. S1 ' / 'Table S1 ' text into label + bookmarked SEQ field."""
    # remove any pre-existing bookmarks in the paragraph
    for el in list(p._p):
        if el.tag in (qn("w:bookmarkStart"), qn("w:bookmarkEnd")):
            p._p.remove(el)
    runs = text_runs(p)
    full = "".join(r.find(qn("w:t")).text or "" for r in runs)
    m = CAP_RE.match(full)
    assert m, full[:60]
    replace_in_paragraph(p, full[:m.end()], "")
    rpr = runs[0].find(qn("w:rPr"))
    lab_rpr = copy.deepcopy(rpr) if rpr is not None else OxmlElement("w:rPr")
    if lab_rpr.find(qn("w:b")) is None:
        lab_rpr.insert(0, OxmlElement("w:b"))
    seq = "Figure" if kind == "Fig" else "Table"
    label = "Fig. " if kind == "Fig" else "Table "
    els = [make_run(label, lab_rpr)]
    bs = OxmlElement("w:bookmarkStart"); bid = next_bm()
    bs.set(qn("w:id"), bid); bs.set(qn("w:name"), bookmark_name(kind, n))
    els.append(bs)
    els.append(make_run("S", lab_rpr))
    els += field_runs(f"SEQ {seq} \\* ARABIC", str(n), lab_rpr)
    be = OxmlElement("w:bookmarkEnd"); be.set(qn("w:id"), bid)
    els.append(be)
    els.append(make_run(" ", lab_rpr))
    ppr = p._p.find(qn("w:pPr"))
    idx = list(p._p).index(ppr) + 1 if ppr is not None else 0
    for k, el in enumerate(els):
        p._p.insert(idx + k, el)
    # Caption style + keep caption with its table (tables) / figure above (figures)
    p.style = doc_caption_style
    if kind == "Table":
        p.paragraph_format.keep_with_next = True


def unwrap_hyperlinks(p):
    """Replace internal w:hyperlink wrappers (static links left by an earlier build) by their runs."""
    for h in p._p.findall(qn("w:hyperlink")):
        parent = h.getparent()
        i = list(parent).index(h)
        for k, child in enumerate(list(h)):
            parent.insert(i + k, child)
        parent.remove(h)


def plain_runs(p):
    """Text runs that are NOT part of a field (instructions or results), in document order. Cross-references
    are only created in plain text, so caption labels (S + SEQ result) and existing fields are never touched."""
    out, depth = [], 0
    for r in p._p.iter(qn("w:r")):
        fc = r.find(qn("w:fldChar"))
        if fc is not None:
            kind = fc.get(qn("w:fldCharType"))
            depth += 1 if kind == "begin" else -1 if kind == "end" else 0
            continue
        if depth == 0 and r.find(qn("w:t")) is not None:
            out.append(r)
    return out


def split_run(r, offset):
    """Split text run r at character offset; returns the right-hand run (inserted after r)."""
    t = r.find(qn("w:t"))
    text = t.text or ""
    right = copy.deepcopy(r)
    t.text = text[:offset]
    t.set(qn("xml:space"), "preserve")
    rt = right.find(qn("w:t"))
    rt.text = text[offset:]
    rt.set(qn("xml:space"), "preserve")
    r.addnext(right)
    return right


def isolate(p, a, b):
    """Split runs so that characters [a, b) of the paragraph text are covered by whole runs; return them."""
    pos = 0
    for r in plain_runs(p):
        n = len(r.find(qn("w:t")).text or "")
        if pos < a < pos + n:
            split_run(r, a - pos)
            break
        pos += n
    pos = 0
    for r in plain_runs(p):
        n = len(r.find(qn("w:t")).text or "")
        if pos < b < pos + n:
            split_run(r, b - pos)
            break
        pos += n
    out, pos = [], 0
    for r in plain_runs(p):
        n = len(r.find(qn("w:t")).text or "")
        if pos >= a and pos + n <= b and n:
            out.append(r)
        pos += n
    return out


def add_crossrefs(p, stats=None):
    """Replace S-number tokens following Table/Fig labels by hyperlinked REF fields (paragraph-level)."""
    unwrap_hyperlinks(p)
    full = "".join(r.find(qn("w:t")).text or "" for r in plain_runs(p))
    spans = []
    for g in REF_GROUP.finditer(full):
        kind = "Table" if g.group(1).startswith("Table") else "Fig"
        for s in SNUM.finditer(g.group(2)):
            n = int(s.group(1))
            if (kind == "Table" and n <= N_TAB) or (kind == "Fig" and n <= N_FIG):
                spans.append((g.start(2) + s.start(), g.start(2) + s.end(), kind, n))
    for a, b, kind, n in reversed(spans):
        runs = isolate(p, a, b)
        if not runs:
            continue
        rpr = runs[0].find(qn("w:rPr"))
        new = field_runs(f"REF {bookmark_name(kind, n)} \\h", f"S{n}", rpr)
        anchor = runs[0]
        for el in new:
            anchor.addprevious(el)
        for r in runs:
            r.getparent().remove(r)
        if stats is not None:
            stats[kind] = stats.get(kind, 0) + 1
            stats.setdefault("targets", set()).add((kind, n))


# ------------------------------------------------------------------ styles, numbering, header/footer
def builtin_style(doc, ui_name, internal_name, base="Normal"):
    """Return a paragraph style that Word recognises as its BUILT-IN style (e.g. internal name 'caption',
    styleId 'Caption'): python-docx marks added styles as custom, which makes Word rename them ('Caption1')."""
    from docx.enum.style import WD_STYLE_TYPE
    for s in doc.styles:
        nm = s.element.find(qn("w:name"))
        if nm is not None and nm.get(qn("w:val")).lower() == internal_name.lower():
            break
    else:
        s = doc.styles.add_style(ui_name, WD_STYLE_TYPE.PARAGRAPH)
        s.base_style = doc.styles[base]
    s.element.set(qn("w:styleId"), ui_name.replace(" ", ""))
    s.element.find(qn("w:name")).set(qn("w:val"), internal_name)
    s.element.attrib.pop(qn("w:customStyle"), None)
    return s


def ensure_styles(doc):
    st = doc.styles
    cap = builtin_style(doc, "Caption", "caption")
    cap.font.size = Pt(10)
    cap.font.color.rgb = None
    cap.font.italic = False
    cap.font.bold = False
    cap.paragraph_format.space_before = Pt(6)
    cap.paragraph_format.space_after = Pt(4)
    cap.paragraph_format.line_spacing = 1.15
    th = builtin_style(doc, "TOC Heading", "TOC Heading")
    th.font.bold = True
    th.font.size = Pt(14)
    th.paragraph_format.space_after = Pt(6)
    for lvl in (1, 2, 3):
        s = builtin_style(doc, f"TOC {lvl}", f"toc {lvl}")
        s.paragraph_format.left_indent = Pt(12 * (lvl - 1))
        s.paragraph_format.space_after = Pt(2)
        s.paragraph_format.line_spacing = 1.0
        s.font.size = Pt(10)
    s = builtin_style(doc, "Table of Figures", "table of figures")
    s.paragraph_format.space_after = Pt(2)
    s.paragraph_format.line_spacing = 1.0
    s.font.size = Pt(9)
    for name, size in (("Heading 1", 16), ("Heading 2", 14), ("Heading 3", 12)):
        h = st[name]
        h.font.name = "Times New Roman"
        h.font.size = Pt(size)
        h.font.bold = True
        h.font.color.rgb = None
        h.paragraph_format.space_before = Pt(12 if name != "Heading 3" else 10)
        h.paragraph_format.space_after = Pt(6)
        h.paragraph_format.keep_with_next = True
        h.paragraph_format.line_spacing = 1.15
        rpr = h.element.get_or_add_rPr()
        fonts = rpr.find(qn("w:rFonts"))
        if fonts is not None:
            for a in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
                fonts.attrib.pop(qn(a), None)


def add_heading_numbering(doc):
    """Multilevel list: Heading 2 -> '1.', Heading 3 -> '1.1' (Heading 1 = document title, unnumbered)."""
    numbering = doc.part.numbering_part.element
    abs_ids = [int(a.get(qn("w:abstractNumId"))) for a in numbering.findall(qn("w:abstractNum"))]
    num_ids = [int(n.get(qn("w:numId"))) for n in numbering.findall(qn("w:num"))]
    aid, nid = max(abs_ids + [0]) + 1, max(num_ids + [0]) + 1
    an = OxmlElement("w:abstractNum"); an.set(qn("w:abstractNumId"), str(aid))
    mlt = OxmlElement("w:multiLevelType"); mlt.set(qn("w:val"), "multilevel"); an.append(mlt)
    for ilvl, (fmt, style) in enumerate([("%1.", "Heading2"), ("%1.%2", "Heading3")]):
        lvl = OxmlElement("w:lvl"); lvl.set(qn("w:ilvl"), str(ilvl))
        for tag, val in (("w:start", "1"), ("w:numFmt", "decimal"), ("w:pStyle", style),
                         ("w:lvlText", fmt), ("w:lvlJc", "left")):
            e = OxmlElement(tag); e.set(qn("w:val"), val); lvl.append(e)
        ppr = OxmlElement("w:pPr"); ind = OxmlElement("w:ind")
        ind.set(qn("w:left"), str(567 + 284 * ilvl)); ind.set(qn("w:hanging"), str(567 + 284 * ilvl)); ppr.append(ind)
        lvl.append(ppr)
        an.append(lvl)
    first_num = numbering.find(qn("w:num"))
    if first_num is not None:
        first_num.addprevious(an)
    else:
        numbering.append(an)
    num = OxmlElement("w:num"); num.set(qn("w:numId"), str(nid))
    ref = OxmlElement("w:abstractNumId"); ref.set(qn("w:val"), str(aid)); num.append(ref)
    numbering.append(num)
    for ilvl, sname in enumerate(("Heading 2", "Heading 3")):
        ppr = doc.styles[sname].element.get_or_add_pPr()
        old = ppr.find(qn("w:numPr"))
        if old is not None:
            ppr.remove(old)
        numpr = OxmlElement("w:numPr")
        il = OxmlElement("w:ilvl"); il.set(qn("w:val"), str(ilvl)); numpr.append(il)
        ni = OxmlElement("w:numId"); ni.set(qn("w:val"), str(nid)); numpr.append(ni)
        ppr.insert(0, numpr)


def header_footer(doc):
    sec = doc.sections[0]
    sec.different_first_page_header_footer = False
    hdr = sec.header
    hdr.is_linked_to_previous = False
    hp = hdr.paragraphs[0]
    set_paragraph_text(hp, "")
    for r in list(hp.runs):
        r._r.getparent().remove(r._r)
    run = hp.add_run(f"Supplementary Material — {CHAPTER}")
    run.italic = True
    run.font.size = Pt(9)
    hp.alignment = 2  # right
    ftr = sec.footer
    fp = ftr.paragraphs[0]
    for el in list(fp._p):
        if el.tag != qn("w:pPr"):
            fp._p.remove(el)
    for extra in ftr.paragraphs[1:]:
        extra._p.getparent().remove(extra._p)
    fp.alignment = 1  # centre
    rpr = OxmlElement("w:rPr"); sz = OxmlElement("w:sz"); sz.set(qn("w:val"), "18"); rpr.append(sz)
    els = [make_run("Page ", rpr)] + field_runs("PAGE", "1", rpr) + [make_run(" of ", rpr)] + \
          field_runs("NUMPAGES", "1", rpr)
    for el in els:
        fp._p.append(el)


def insert_paragraph_after(par, text="", style=None):
    new_p = OxmlElement("w:p")
    par._p.addnext(new_p)
    from docx.text.paragraph import Paragraph
    np_ = Paragraph(new_p, par._parent)
    if text:
        np_.add_run(text)
    if style:
        np_.style = style
    return np_


def toc_paragraph(after, instr, placeholder):
    p = insert_paragraph_after(after)
    for el in field_runs(instr, placeholder):
        p._p.append(el)
    return p


def page_break_after(par):
    p = insert_paragraph_after(par)
    p.add_run().add_break(WD_BREAK.PAGE)
    return p


# ------------------------------------------------------------------ in-text citations of every item
def item_list(kind, nums):
    """'Fig. S4' / 'Figs. S1, S2 and S3' / 'Table S5' / 'Tables S1, S2 and S3'."""
    one, many = ("Fig.", "Figs.") if kind == "Fig" else ("Table", "Tables")
    nums = sorted(nums)
    if len(nums) == 1:
        return f"{one} S{nums[0]}"
    return f"{many} " + ", ".join(f"S{n}" for n in nums[:-1]) + f" and S{nums[-1]}"


def figure_source_map():
    """Table S-number -> figure S-numbers that display it (from the verified figure/source-table mapping)."""
    src = (BUILD_V2 / "build_supp.py").read_text(encoding="utf-8")
    torder = eval(re.search(r"TORDER = (\[.*?\n\])", src, re.S).group(1))
    tab_s, n = {}, 0
    for _, ids in torder:
        for tid in ids:
            n += 1
            tab_s[tid] = n
    sys.path.insert(0, str(BUILD_V2))
    from supp_figs import SFIGS
    out = {}
    for k, (_, _, _, srcs) in enumerate(SFIGS, 1):
        for tid in srcs.split(", "):
            out.setdefault(tab_s[tid], []).append(k)
    return out


def add_section_leadins(doc, caps_style_id="Caption"):
    """Insert 'This section presents ...' after every Heading 3, citing each figure/table in the section."""
    paras = doc.paragraphs
    heads = [i for i, p in enumerate(paras) if p.style.name in ("Heading 2", "Heading 3")]
    n_added = 0
    for k, i in enumerate(heads):
        if paras[i].style.name != "Heading 3":
            continue
        end = heads[k + 1] if k + 1 < len(heads) else len(paras)
        figs, tabs = [], []
        for p in paras[i + 1:end]:
            m = re.match(r"(Fig\.|Table) S(\d+)\b", p.text)
            if m and p.style.style_id == caps_style_id:
                (figs if m.group(1) == "Fig." else tabs).append(int(m.group(2)))
        parts = [item_list("Fig", figs)] if figs else []
        parts += [item_list("Table", tabs)] if tabs else []
        if not parts:
            continue
        lead = insert_paragraph_after(paras[i], "This section presents " + " and ".join(parts) + ".")
        lead.paragraph_format.space_after = Pt(6)
        for r in lead.runs:
            r.italic = True
            r.font.size = Pt(10)
        n_added += 1
    return n_added


def add_table_figure_links(doc, tab2fig, caps_style_id="Caption"):
    """Append 'Displayed in Fig. Sx.' to the notes paragraph of every table that a figure plots."""
    from docx.text.paragraph import Paragraph
    added = 0
    for p in doc.paragraphs:
        m = re.match(r"Table S(\d+)\b", p.text)
        if not (m and p.style.style_id == caps_style_id):
            continue
        n = int(m.group(1))
        if n not in tab2fig:
            continue
        el = p._p.getnext()
        while el is not None:
            if el.tag == qn("w:p"):
                q = Paragraph(el, p._parent)
                t = q.text.strip()
                if t and not re.fullmatch(r"Part \d+ of \d+", t):
                    break
            el = el.getnext()
        if el is None:
            continue
        runs = text_runs(q)
        rpr = runs[-1].find(qn("w:rPr")) if runs else None
        sep = "" if q.text.endswith(" ") else " "
        q._p.append(make_run(f"{sep}Displayed in {item_list('Fig', tab2fig[n])}.", rpr))
        added += 1
    return added


# ------------------------------------------------------------------ main
def locate_fields(paras):
    """Map (item, field) -> paragraph, using the same rules as the review extraction."""
    P = [p.text.strip() for p in paras]
    loc = {}
    for i, t in enumerate(P):
        m = re.match(r"(Table|Fig\.) S(\d+)\b", t)
        if not m:
            continue
        item = f"{'Table' if m.group(1) == 'Table' else 'Fig.'} S{m.group(2)}"
        loc[(item, "caption")] = paras[i]
        if m.group(1) == "Fig.":
            loc[(item, "interpretation")] = paras[i + 1]
            continue
        seg, j = [], i + 1
        while j < len(P) and not re.match(r"Table S\d+\b", P[j]) and not paras[j].style.name.startswith("Heading"):
            if P[j] and not P[j].startswith("Part "):
                seg.append(paras[j])
            j += 1
        if seg:
            loc[(item, "notes")] = seg[0]
        if len(seg) > 1:
            loc[(item, "interpretation")] = seg[-1]
    return loc


def main():
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(SUPP_DOCX_SOURCE, OUT)          # work on a copy only
    doc = Document(OUT)
    ensure_styles(doc)
    cap_style = next(x for x in doc.styles if x.style_id == "Caption")
    th_style = next(x for x in doc.styles if x.style_id == "TOCHeading")
    for p in doc.paragraphs:
        unwrap_hyperlinks(p)
    paras = doc.paragraphs

    # 1. corrections
    corr = pd.read_excel(CORR, sheet_name="Corrections")
    loc = locate_fields(paras)
    applied, failed = 0, []
    for _, c in corr.iterrows():
        p = loc.get((c["Item"], c["Field"]))
        if p is not None and replace_in_paragraph(p, c["Original_Text_In_Document"], c["Corrected_Text"]):
            applied += 1
        else:
            failed.append(f"{c['Item']}/{c['Field']}: {c['Original_Text_In_Document'][:50]}")
    log(f"DOCX: {applied}/{len(corr)} text corrections applied; failed: {len(failed)}")
    for f in failed:
        log(f"  FAILED correction {f}")

    # 2. captions + 3. cross-references
    stats = {}
    caps = {"Fig": [], "Table": []}
    for p in doc.paragraphs:
        m = CAP_RE.match(p.text)
        if m and (p.text.startswith("Fig. S") or p.text.startswith("Table S")) and p.style.name != "Heading 3":
            kind = "Fig" if m.group(1) == "Fig." else "Table"
            n = int(m.group(2))
            # rest of the caption may itself cite other items
            convert_caption(p, kind, n, cap_style)
            caps[kind].append(n)
    for kind, total in (("Fig", N_FIG), ("Table", N_TAB)):
        assert caps[kind] == list(range(1, total + 1)), f"{kind} caption sequence broken: {caps[kind][:10]}"
    n_lead = add_section_leadins(doc)
    n_links = add_table_figure_links(doc, figure_source_map())
    log(f"DOCX: section lead-ins {n_lead}; table->figure links {n_links}")
    for p in doc.paragraphs:
        add_crossrefs(p, stats=stats)
    for t in doc.tables:   # citations inside table cells
        for row in t.rows:
            for cell in row.cells:
                for p in cell.paragraphs:
                    add_crossrefs(p, stats=stats)
    log(f"DOCX: captions Fig {len(caps['Fig'])}, Table {len(caps['Table'])}; "
        f"REF fields: tables {stats.get('Table', 0)}, figures {stats.get('Fig', 0)}")
    cited = stats.get("targets", set())
    not_cited = [f"Table S{n}" for n in range(1, N_TAB + 1) if ("Table", n) not in cited] + \
                [f"Fig. S{n}" for n in range(1, N_FIG + 1) if ("Fig", n) not in cited]

    assert not not_cited, f"items without an in-text citation: {not_cited}"
    # 4. front matter: intro, TOC, List of Figures, List of Tables
    paras = doc.paragraphs
    title = next(p for p in paras if p.style.name == "Heading 1")
    ph = [p for p in paras if p.text.strip().lower() in ("table list", "figure list")]
    intro = insert_paragraph_after(title, INTRO)
    intro.paragraph_format.line_spacing = 1.15
    intro.paragraph_format.space_after = Pt(6)
    for r in intro.runs:
        r.font.size = Pt(11)
    add_crossrefs(intro, stats=stats)
    anchor = intro
    pb0 = page_break_after(anchor)
    h = insert_paragraph_after(pb0, "Contents", th_style)
    toc = toc_paragraph(h, 'TOC \\o "1-3" \\h \\z \\u', "Right-click and choose Update Field to build the contents.")
    pb1 = page_break_after(toc)
    hf = insert_paragraph_after(pb1, "List of Figures", th_style)
    lof = toc_paragraph(hf, 'TOC \\h \\z \\c "Figure"', "List of Figures (update field).")
    pb2 = page_break_after(lof)
    ht = insert_paragraph_after(pb2, "List of Tables", th_style)
    lot = toc_paragraph(ht, 'TOC \\h \\z \\c "Table"', "List of Tables (update field).")
    page_break_after(lot)
    for p in ph:   # placeholders
        p._p.getparent().remove(p._p)
    # Remove empty paragraphs immediately before the first section heading (avoid blank page)
    # 5. headings, figures, part labels
    add_heading_numbering(doc)
    body = list(doc.paragraphs)
    for i, p in enumerate(body):
        if p._p.findall(".//" + qn("w:drawing")):
            p.paragraph_format.keep_with_next = True
            p.alignment = 1
        if re.fullmatch(r"Part \d+ of \d+", p.text.strip()):
            for r in p.runs:
                r.italic = True
                r.font.size = Pt(8)
            p.paragraph_format.keep_with_next = True
    # Supplementary Figures section starts on a new page after the lists; Tables section on a new page
    for p in doc.paragraphs:
        if p.style.name == "Heading 2" and p.text.strip() == "Supplementary Tables":
            p.paragraph_format.page_break_before = True

    # 6. header / footer, properties
    header_footer(doc)
    cp = doc.core_properties
    cp.title = "Supplementary Material"
    cp.subject = CHAPTER
    cp.keywords = "antimicrobial resistance; urban microbiome; metagenomics; supplementary material"
    cp.comments = "Automatic captions (SEQ), cross-references (REF), TOC, List of Figures and List of Tables."
    # Fields are updated by Word in step 07; no 'update fields on open' prompt for readers.
    uf = doc.settings.element.find(qn("w:updateFields"))
    if uf is not None:
        doc.settings.element.remove(uf)
    doc.save(OUT)
    json.dump(dict(corrections_applied=applied, corrections_failed=failed, captions={k: len(v) for k, v in caps.items()},
                   ref_tables=stats.get("Table", 0), ref_figures=stats.get("Fig", 0), not_cited_in_body=not_cited),
              open(DOCS / "Supplementary_Appendices" / "docx_build_report.json", "w"), indent=1)
    log(f"DOCX written: {OUT.name}; every table and figure cited in the text")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
