"""Build supplementary blocks (figures + 96 tables) and the full-data Excel workbook.
Appends 'supp' to build.json written by build_json.py."""
import json, math, re, sys, glob, os
import pandas as pd
sys.path.insert(0, '/home/claude/amr'); sys.path.insert(0, '/home/claude/amr/supp')
from supp_figs import SFIGS, EXTRA_INTERP
from PIL import Image

TDIR = '/mnt/user-data/uploads/Research/AMR/Updated_Research_Analysis/04_Tables/'
FDIR = '/mnt/user-data/uploads/Research/AMR/Updated_Research_Analysis/05_Figures/'
catalog = json.load(open('/home/claude/amr/supp/table_catalog.json'))
interp = {}
for f in sorted(glob.glob('/home/claude/amr/supp/interp_g*.json')):
    interp.update(json.load(open(f)))
figraw = json.load(open('/home/claude/amr/supp/fig_interp_raw.json'))
OVR = {}
try:
    OVR = json.load(open('/home/claude/amr/supp/overrides.json'))
except FileNotFoundError:
    pass
for _k, _v in OVR.items():
    interp.setdefault(_k, {}).update(_v)

TORDER = [
 ("Evidence base and geographic coverage", ["D01","D02","D04","D05","D06","D07","D09","D12","R09u","S01","S03","GI06"]),
 ("Stream 1 (Indian open drains): taxonomic composition, diversity and ordination", ["T00","T02","T03","T04","T06s","T08s","P01","P02","P03b","P03c"]),
 ("Stream 1 (Indian open drains): resistome composition, abundance and determinants", ["A01","A02","A03","A04","A05s","A07s","C01","C03","C04","P04","P05","A09","A10s","A11"]),
 ("Stream 2 (curated global resource): urban matrices and food-related contextual units", ["A12","A13","A14","T11","T12","C08","C08b","C09","C10","C11","C12","C13","C14","C15","P08"]),
 ("Stream 3 (uniform reprocessing): Dhaka and South Asian comparisons", ["DH02","DH03","DH32","DH33","DH06","DH08","DH09","DH10","DH17","DH11","DH12","DH14","DH16","DH23","DH25","DH28","DH29","DH29b","DH30","DH31","GI01","GI02","GI04","GI05"]),
 ("Predictive modelling", ["M01","M02","M03","M04","M05","M06","M07"]),
 ("Spatial analysis (2016 Global Sewage Surveillance pilot)", ["S04","S05","S06","S07"]),
 ("Robustness and sequencing-depth sensitivity", ["SA01","SA02","SA03","SA04","SA05","SD01","SD03","SD04","SD05","SD07"]),
]
assert sum(len(v) for _, v in TORDER) == 96

PCOLS = re.compile(r'(^|_)(p|q|p_bh|p_value|pvalue|p_perm|perm_p|p_adj|q_bh|pval|p_lrt|lrt_p|p_exact|exact_p|p_bh_within_family|q_value)($|_)', re.I)
def is_p(col):
    c = str(col).lower()
    return bool(PCOLS.search(c)) or c in ('p', 'q') or c.startswith('p_') or c.endswith('_p') or 'pvalue' in c or c.endswith('_q') or c.startswith('q_')

def fmt(v, col):
    if v is None or (isinstance(v, float) and math.isnan(v)): return '–'
    if isinstance(v, (bool,)): return str(v)
    if isinstance(v, (int,)) or (isinstance(v, float) and float(v).is_integer() and not is_p(col) and abs(v) < 1e15):
        return f'{int(v):,}' if abs(v) >= 10000 else str(int(v))
    if isinstance(v, float):
        a = abs(v)
        if is_p(col):
            return f'{v:.1e}' if a < 0.001 else f'{v:.3f}'
        if a >= 10000: return f'{v:,.0f}'
        if a >= 100: return f'{v:.1f}'
        if a >= 1: return f'{v:.2f}'
        if a >= 0.01: return f'{v:.3f}'
        if a == 0: return '0'
        return f'{v:.2e}'
    s = str(v)
    return s

def hdr(c):
    s = str(c).replace('_', ' ')
    s = re.sub(r'\bpct\b', '%', s)
    return s

def table_blocks(tid, df, sort_hint):
    nrows, ncols = df.shape
    note_extra = []
    show = df
    if nrows > 120:
        col = sort_hint if sort_hint in df.columns else None
        if col is None:
            num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
            col = num[0] if num else None
        if col is not None:
            show = df.sort_values(col, ascending=False, key=lambda s: s.abs() if pd.api.types.is_numeric_dtype(s) else s).head(30)
            note_extra.append(f'Abridged: the 30 of {nrows} rows with the largest absolute values of "{hdr(col)}" are shown; the full table is in the Supplementary Data workbook.')
        else:
            show = df.head(30)
            note_extra.append(f'Abridged: first 30 of {nrows} rows shown; the full table is in the Supplementary Data workbook.')
    cols = list(show.columns)
    pieces = []
    if len(cols) > 11 and len(show) <= 8:
        # transpose: first column values become headers
        key = cols[0]
        header = [hdr(key)] + [fmt(v, key) for v in show[key]]
        rows = [[hdr(c)] + [fmt(v, c) for v in show[c]] for c in cols[1:]]
        pieces.append((header, rows, ''))
        note_extra.append('Table transposed for layout (statistics as rows).')
    elif len(cols) > 11:
        key = cols[0]; rest = cols[1:]
        chunks = [rest[i:i+9] for i in range(0, len(rest), 9)]
        for k, ch in enumerate(chunks):
            header = [hdr(key)] + [hdr(c) for c in ch]
            rows = [[fmt(r[key], key)] + [fmt(r[c], c) for c in ch] for _, r in show.iterrows()]
            pieces.append((header, rows, f'Part {k+1} of {len(chunks)}'))
        note_extra.append(f'Columns split into {len(chunks)} parts for layout; the first column is repeated.')
    else:
        header = [hdr(c) for c in cols]
        rows = [[fmt(r[c], c) for c in cols] for _, r in show.iterrows()]
        pieces.append((header, rows, ''))
    return pieces, note_extra

def main():
    supp = []
    supp.append(dict(t='h1', text='Supplementary Material'))
    supp.append(dict(t='sp', text='This supplement reports every accepted statistical figure (30, Figs. S1–S30) and statistical table (96, Tables S1–S96) from the chapter’s verified analysis (Updated_Research_Analysis, 31 October 2026), each with an interpretation. Items excluded at screening (superseded outputs, quality-control and process records) are not reproduced. Table values were generated programmatically from the verified CSV files without manual transcription; tables longer than 120 rows are abridged here and given in full in the accompanying Supplementary Data workbook (one sheet per table). Streams 1–3 use different pipelines, reference databases and abundance units (RPKM, FPKM, FPKM-16S) and must not be compared in magnitude. Detection of a resistance gene indicates DNA presence, not phenotypic resistance, expression, viability, exposure or transmission.'))
    supp.append(dict(t='h2', text='Supplementary Figures'))
    tab_map = {}
    _n = 0
    for _g, _ids in TORDER:
        for _t in _ids:
            _n += 1; tab_map[_t] = f'S{_n}'
    fig_map = {}
    cur = None
    for i, (fid, grp, cap, srcs) in enumerate(SFIGS, start=1):
        fig_map[fid] = f'S{i}'
        if grp != cur:
            supp.append(dict(t='h3', text=grp)); cur = grp
        path = FDIR + [f for f in os.listdir(FDIR) if f.startswith(fid)][0]
        if fid == 'Fig24': path = '/home/claude/amr/figs/Fig24_fixed.png'
        W, H = Image.open(path).size
        if fid in EXTRA_INTERP:
            itxt = EXTRA_INTERP[fid]
        else:
            r = figraw[fid]
            itxt = r['interp'] + (' ' + r['limitation'] if r['limitation'] else '')
        supp.append(dict(t='sfig', n=f'S{i}', path=path, iw=W, ih=H, caption=cap, interp=itxt,
                         source='Source data: ' + ', '.join('Table ' + tab_map[x] for x in srcs.split(', ')) + '. Original analysis file: ' + os.path.basename(path).replace('_fixed','') + '.'))
    supp.append(dict(t='h2', text='Supplementary Tables'))
    n = 0
    for grp, ids in TORDER:
        supp.append(dict(t='h3', text=grp))
        for tid in ids:
            n += 1
            meta = catalog[tid]; it = interp.get(tid, {})
            df = pd.read_csv(TDIR + meta['file'])
            pieces, extra = table_blocks(tid, df, it.get('sort_hint', ''))
            title = it.get('title') or meta['title']
            notes = ' '.join(x for x in [it.get('notes', ''), ' '.join(extra)] if x)
            supp.append(dict(t='stable', n=f'S{n}', tid=tid, title=title, pieces=pieces, notes=notes, srcfile=meta['file'],
                             interp=it.get('interpretation', '[INTERPRETATION MISSING]'), nrows=int(df.shape[0]), ncols=int(df.shape[1])))
    ids_sorted = sorted(map(re.escape, tab_map), key=len, reverse=True)
    p1 = re.compile(r'\b(Tables?) (' + '|'.join(ids_sorted) + r')(?![A-Za-z0-9])')
    p2 = re.compile(r'(?<!Table )(?<!Tables )(?<![A-Za-z0-9])(' + '|'.join(ids_sorted) + r')(?![A-Za-z0-9])')
    def fix(s):
        s = p1.sub(lambda m: m.group(1) + ' ' + tab_map[m.group(2)], s)
        s = p2.sub(lambda m: 'Table ' + tab_map[m.group(1)], s)
        return s
    for b in supp:
        for k in ('interp', 'notes'):
            if k in b and isinstance(b[k], str) and not b[k].startswith('Source:'):
                b[k] = fix(b[k])
    data = json.load(open('/home/claude/amr/build.json'))
    data['supp'] = supp
    json.dump(data, open('/home/claude/amr/build.json', 'w'), ensure_ascii=False)
    json.dump(dict(fig_map=fig_map, tab_map=tab_map), open('/home/claude/amr/supp/smaps.json', 'w'), indent=1)
    # Excel workbook with all 96 full tables
    with pd.ExcelWriter('/home/claude/amr/Supplementary_Data_Tables_S1-S96.xlsx', engine='openpyxl') as xw:
        idx = []
        for grp, ids in TORDER:
            for tid in ids:
                meta = catalog[tid]
                idx.append(dict(Table=f'Table {tab_map[tid]}', Section=grp, Title=interp.get(tid, {}).get('title') or meta['title'], Rows=None, Source_CSV=meta['file']))
        for row in idx:
            tid = [k for k, v in tab_map.items() if f'Table {v}' == row['Table']][0]
            df = pd.read_csv(TDIR + catalog[tid]['file']); row['Rows'] = len(df)
        pd.DataFrame(idx).to_excel(xw, sheet_name='Index', index=False)
        for grp, ids in TORDER:
            for tid in ids:
                df = pd.read_csv(TDIR + catalog[tid]['file'])
                df.to_excel(xw, sheet_name=f'Table {tab_map[tid]}', index=False)
    print('supp blocks', len(supp), 'tables', n, 'figs', len(SFIGS))

if __name__ == '__main__':
    main()
