"""Patch amr_common.py of the corrected package for the 31 Oct 2026 update.
Every replacement is asserted to match exactly once so that no silent edit occurs."""
import sys, re
from pathlib import Path
P = Path(sys.argv[1]) / "scripts" / "amr_common.py"
s = P.read_text(encoding="utf-8")

def rep(old, new, count=1):
    global s
    n = s.count(old)
    assert n == count, (n, old[:120])
    s = s.replace(old, new)

# 1. Updated harmonized inputs (originals untouched): AMR_HARM points to the corrected copy
rep('import json, hashlib, platform, warnings, datetime as dt, re\n',
    'import json, hashlib, platform, warnings, datetime as dt, re, os\n')
rep('HARM = BUNDLE / "harmonized"\n',
    'HARM_ORIGINAL = BUNDLE / "harmonized"\n'
    '# UPDATE 2026-10-31: AMR_HARM may point to the corrected copy of the harmonized tables\n'
    '# (sample-level sample-type overrides AUD-04, DS02 site identifiers AUD-11, coordinate provenance AUD-18).\n'
    'HARM = Path(os.environ["AMR_HARM"]) if os.environ.get("AMR_HARM") else HARM_ORIGINAL\n')
rep('ANALYSIS_DATE = "2026-10-31"', 'ANALYSIS_DATE = "2026-10-31"\nUPDATE_LABEL = "Updated analysis (audit corrections + Dhaka reprocessing), 31 October 2026"')

# 2. Publication typography (AUD-15): Elsevier final-size lettering 6.5-9 pt, sans-serif
rep('''plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
})''', '''plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "legend.fontsize": 6.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
# Elsevier artwork sizes (final print size): single column ~85 mm (90 mm max), double column ~175 mm (190 mm max)
SINGLE_COL_IN, DOUBLE_COL_IN = 85 / 25.4, 175 / 25.4''')

# 3. save_fig: RGB PNG (300 dpi) + LZW TIFF (600 dpi), final-size and legibility audit recorded in the catalogue
rep('''def save_fig(fig, num: int, slug: str, title: str, caption: str, nb: str, objective: str, source: str):
    fn = f"Fig{num:02d}_{slug}.png"
    fig.savefig(FIG / fn, dpi=300, facecolor="white")
    plt.close(fig)
    _append("figures", nb, dict(key=fn, number=num, file=fn, title=title, caption=caption,
                                notebook=nb, objective=objective, source_data=source))
    print(f"[figure] Figure {num}: {title} -> figures/{fn}")
    return fn''', '''def figure_layout_audit(fig):
    """Return (min_font_pt, n_overlapping_text_pairs, n_texts) for the visible text of a drawn figure.
    Overlap = two non-empty text boxes intersecting by more than 15% of the smaller box (tick labels,
    titles, legends, annotations). Used to flag, not to fix, layout problems (AUD-15)."""
    import matplotlib.text as mtext
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    boxes, sizes = [], []
    for t in fig.findobj(mtext.Text):
        if not t.get_visible() or not str(t.get_text()).strip() or t.get_alpha() == 0: continue
        try: bb = t.get_window_extent(r)
        except Exception: continue
        if bb.width < 1 or bb.height < 1: continue
        boxes.append(bb); sizes.append(float(t.get_fontsize()))
    ov = 0
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            w = min(a.x1, b.x1) - max(a.x0, b.x0); h = min(a.y1, b.y1) - max(a.y0, b.y0)
            if w > 0 and h > 0 and w * h > 0.15 * min(a.width * a.height, b.width * b.height): ov += 1
    return (min(sizes) if sizes else np.nan), ov, len(boxes)

def save_fig(fig, num: int, slug: str, title: str, caption: str, nb: str, objective: str, source: str):
    """UPDATE 2026-10-31 (AUD-15): figures are laid out at final print width (<= 175 mm) and saved as
    RGB PNG (300 dpi) plus LZW-compressed RGB TIFF (600 dpi, combination artwork) in figures/tiff/."""
    from PIL import Image
    fn = f"Fig{num:02d}_{slug}.png"
    minfont, overlaps, ntext = figure_layout_audit(fig)
    fig.savefig(FIG / fn, dpi=300, facecolor="white")
    (FIG / "tiff").mkdir(exist_ok=True)
    tif = f"Fig{num:02d}_{slug}.tif"
    fig.savefig(FIG / "tiff" / tif, dpi=600, facecolor="white", format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    for pth in (FIG / fn, FIG / "tiff" / tif):  # flatten to RGB (no alpha channel) for production
        im = Image.open(pth)
        if im.mode != "RGB":
            dpi = im.info.get("dpi", (300, 300)); im = im.convert("RGB")
            im.save(pth, dpi=dpi, **({"compression": "tiff_lzw"} if pth.suffix == ".tif" else {}))
    im = Image.open(FIG / fn); wmm, hmm = im.size[0] / 300 * 25.4, im.size[1] / 300 * 25.4
    layout = "single column (<= 90 mm)" if wmm <= 90 else ("double column (<= 190 mm)" if wmm <= 190 else "OVERSIZE (> 190 mm)")
    _append("figures", nb, dict(key=fn, number=num, file=fn, tiff_file=f"tiff/{tif}", title=title, caption=caption,
                                notebook=nb, objective=objective, source_data=source, width_mm=round(wmm, 1),
                                height_mm=round(hmm, 1), print_layout=layout, min_font_pt=minfont,
                                text_overlaps_flagged=overlaps, texts_checked=ntext))
    print(f"[figure] Figure {num}: {title} -> figures/{fn} ({wmm:.0f} x {hmm:.0f} mm; min font {minfont} pt; overlaps {overlaps})")
    return fn''')

# 4. Site-aware (cluster) inference helpers for AUD-11
rep('''def cliffs_delta_from_u(u, n1, n2):''', '''def site_level_kw(df, value, site, group, n_perm=9999, seed=SEED, agg="median"):
    """AUD-11: Kruskal-Wallis on site-level summaries (one value per site), with a permutation p-value in
    which group labels are permuted among whole sites. The site, not the repeated sample, is the unit."""
    d = df[[value, site, group]].dropna()
    s = d.groupby(site, observed=True).agg(v=(value, agg), g=(group, "first"), n_samples=(value, "size"))
    assert (d.groupby(site, observed=True)[group].nunique() == 1).all(), "group must be constant within site"
    k = kruskal_eps2(s.rename(columns={"v": value, "g": group}), value, group)
    rng = np.random.default_rng(seed); vals, labs = s.v.values, s.g.values; cnt = 0
    for _ in range(n_perm):
        lp = rng.permutation(labs)
        H = stats.kruskal(*[vals[lp == u] for u in np.unique(lp)])[0]
        cnt += H >= k["H"] - 1e-12
    k.update(p_perm=(cnt + 1) / (n_perm + 1), n_perm=n_perm, n_sites=len(s), samples=int(s.n_samples.sum()),
             site_summary=agg, sites_per_group="; ".join(f"{g}={c}" for g, c in s.g.value_counts().sort_index().items()))
    return k

def mixed_site_lrt(df, formula, terms, group):
    """AUD-11: linear mixed model with a random intercept for `group` (site). Coefficients from REML;
    each listed fixed term tested by a maximum-likelihood ratio test against the model without it."""
    import statsmodels.formula.api as smf
    lhs, rhs = formula.split("~"); parts = [p.strip() for p in rhs.split("+")]
    d = df.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = smf.mixedlm(formula, d, groups=d[group]).fit(reml=False, method=["lbfgs", "powell"])
        out = {}
        for t in terms:
            red = [p for p in parts if p != t] or ["1"]
            m0 = smf.mixedlm(lhs + "~ " + " + ".join(red), d, groups=d[group]).fit(reml=False, method=["lbfgs", "powell"])
            LR = max(0.0, 2 * (full.llf - m0.llf)); dfl = len(full.fe_params) - len(m0.fe_params)
            out[t] = dict(LR=LR, df=dfl, p=stats.chi2.sf(LR, dfl))
        reml = smf.mixedlm(formula, d, groups=d[group]).fit(reml=True, method=["lbfgs", "powell"])
    vr = float(reml.cov_re.iloc[0, 0]); icc = vr / (vr + reml.scale); ci = reml.conf_int()
    fe = reml.fe_params.index
    coef = pd.DataFrame({"term": fe, "estimate": reml.fe_params.values, "SE": reml.bse_fe.values,
                         "ci_low": ci.loc[fe, 0].values, "ci_high": ci.loc[fe, 1].values, "p_value": reml.pvalues[fe].values})
    coef = pd.concat([coef, pd.DataFrame({"term": ["Site variance (random intercept)", "Residual variance", "ICC (site)"],
                                          "estimate": [vr, reml.scale, icc]})], ignore_index=True)
    return dict(lrt=out, coef=coef, icc=icc, n=int(reml.nobs), n_groups=int(d[group].nunique()), model=reml)

def site_perm_chisq(detected, site, group, n_perm=9999, seed=SEED):
    """AUD-11: chi-square statistic for group x detection, with p-value from permuting group labels among
    whole sites (cluster permutation), so repeated samples of a site are not treated as independent."""
    d = pd.DataFrame({"y": np.asarray(detected, bool), "s": np.asarray(site), "g": np.asarray(group)})
    def chi(g):
        t = pd.crosstab(g, d.y).reindex(columns=[False, True], fill_value=0).values
        if (t.sum(axis=0) == 0).any(): return np.nan
        return stats.chi2_contingency(t, correction=False)[0]
    c0 = chi(d.g)
    if np.isnan(c0): return dict(chi2=np.nan, p_perm=np.nan, n=len(d), n_sites=d.s.nunique())
    sg = d.groupby("s").g.first(); rng = np.random.default_rng(seed); cnt = 0
    for _ in range(n_perm):
        mp = dict(zip(sg.index, rng.permutation(sg.values))); cnt += chi(d.s.map(mp)) >= c0 - 1e-12
    site_prev = d.groupby("s").agg(g=("g", "first"), prev=("y", "mean"))
    return dict(chi2=c0, p_perm=(cnt + 1) / (n_perm + 1), n=len(d), n_sites=d.s.nunique(),
                site_prevalence_by_group="; ".join(f"{g}: median {100*v.median():.0f}% (range {100*v.min():.0f}-{100*v.max():.0f}%)"
                                                   for g, v in site_prev.groupby("g").prev))

def cliffs_delta_from_u(u, n1, n2):''')
P.write_text(s, encoding="utf-8")
print("amr_common.py patched")
