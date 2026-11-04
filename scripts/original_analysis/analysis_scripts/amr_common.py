"""Shared configuration, I/O registry and statistical helpers for the
AMR urban-microbiome analysis (Analysis_Results).

All analyses read the harmonized, unmodified source tables in
AMR_Chapter_Reproducibility_Bundle/.../harmonized and write only inside
Analysis_Results/. Original files are never modified.
"""
from __future__ import annotations
import json, hashlib, platform, warnings, datetime as dt, re, os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

SEED = 20261031
RNG = np.random.default_rng(SEED)
def _find_root():
    """Input root = the research directory holding the harmonized bundle (read only).
    Set AMR_ROOT to override. Outputs go to the Analysis_Results folder that holds this script."""
    import os
    cands = [os.environ.get("AMR_ROOT"), str(Path(__file__).resolve().parents[2]),
             r"C:\Users\USER\Documents\Research\AMR"]
    for c in cands:
        if c and (Path(c) / "AMR_Chapter_Reproducibility_Bundle").exists():
            return Path(c)
    raise FileNotFoundError("Research directory not found; set AMR_ROOT")
ROOT = _find_root()
BUNDLE = ROOT / "AMR_Chapter_Reproducibility_Bundle" / "AMR_Chapter_Public_Data_Collection"
HARM_ORIGINAL = BUNDLE / "harmonized"
# UPDATE 2026-10-31: AMR_HARM may point to the corrected copy of the harmonized tables
# (sample-level sample-type overrides AUD-04, DS02 site identifiers AUD-11, coordinate provenance AUD-18).
HARM = Path(os.environ["AMR_HARM"]) if os.environ.get("AMR_HARM") else HARM_ORIGINAL
RAW = BUNDLE / "raw_downloads"
OUT = Path(__file__).resolve().parents[1]
FIG, TAB, QCD, DER, LOG, CAT, EXT, VER = (OUT / s for s in
    ["figures", "tables", "qc", "derived", "logs", "catalog", "external_reference", "verification"])
for _p in (FIG, TAB, QCD, DER, LOG, CAT, EXT, VER):
    _p.mkdir(parents=True, exist_ok=True)

MISSING_CODES = ["Not reported", "Not available", "Not yet verified", "Not applicable", "Restricted access"]
ANALYSIS_DATE = "2026-10-31"
UPDATE_LABEL = "Updated analysis (audit corrections + Dhaka reprocessing), 31 October 2026"

# ---------------------------------------------------------------- palettes
SAMPLE_TYPE_COLORS = {
    "Municipal wastewater": "#0072B2", "Community wastewater": "#56B4E9",
    "Hospital wastewater": "#D55E00", "Surface water": "#009E73",
    "Transit surface": "#E69F00", "Built environment": "#CC79A7",
    "Air": "#7F7F7F", "Urban soil": "#8C6D31", "Food-related environment": "#882255",
}
CITY_COLORS = {"Delhi": "#0072B2", "Mumbai": "#D55E00", "Chennai": "#009E73", "Kolkata": "#CC79A7"}
SEASON_ORDER = ["Winter (Dec-Feb)", "Pre-monsoon (Mar-May)", "Monsoon (Jun-Sep)", "Post-monsoon (Oct-Nov)"]
SEASON_COLORS = dict(zip(SEASON_ORDER, ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]))
CLASS_ORDER = ["Beta-lactam", "Aminoglycoside", "Tetracycline", "MLS", "Fluoroquinolone/quinolone",
               "Phenicol", "Folate pathway antagonist", "Peptide/polymyxin", "Glycopeptide",
               "Rifamycin", "Fosfomycin", "Oxazolidinone", "Biocide", "Multi-class", "Other", "Unclassified"]
_cmap = plt.get_cmap("tab20")
CLASS_COLORS = {c: matplotlib.colors.to_hex(_cmap(i % 20)) for i, c in enumerate(CLASS_ORDER)}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "legend.fontsize": 6.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
# Elsevier artwork sizes (final print size): single column ~85 mm (90 mm max), double column ~175 mm (190 mm max)
SINGLE_COL_IN, DOUBLE_COL_IN = 85 / 25.4, 175 / 25.4

def season_india(month: int) -> str:
    """India Meteorological Department season for a calendar month."""
    if month in (12, 1, 2): return SEASON_ORDER[0]
    if month in (3, 4, 5): return SEASON_ORDER[1]
    if month in (6, 7, 8, 9): return SEASON_ORDER[2]
    if month in (10, 11): return SEASON_ORDER[3]
    return "Not reported"

# ---------------------------------------------------------------- loaders
def load_csv(name: str, **kw) -> pd.DataFrame:
    return pd.read_csv(HARM / name, low_memory=False, **kw)

def is_missing(s: pd.Series) -> pd.Series:
    return s.isna() | s.astype(str).str.strip().isin(MISSING_CODES + ["", "nan", "None"])

def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.where(~is_missing(s)), errors="coerce")

def load_world():
    import geopandas as gpd
    z = EXT / "ne_50m_admin_0_countries.zip"
    if not z.exists():
        raise FileNotFoundError("Natural Earth file missing; run notebook 01 first")
    w = gpd.read_file(f"zip://{z}")
    w["ISO2"] = w["ISO_A2_EH"].where(w["ISO_A2_EH"] != "-99", w["ISO_A2"])
    return w

# ---------------------------------------------------------------- registry
def _cat_path(kind, nb): return CAT / f"{kind}_{nb}.json"

def reset_catalog(nb: str):
    for kind in ("figures", "tables", "tests", "exclusions", "findings"):
        _cat_path(kind, nb).write_text("[]", encoding="utf-8")

def _append(kind, nb, rec):
    p = _cat_path(kind, nb)
    cur = json.loads(p.read_text(encoding="utf-8")) if p.exists() else []
    cur = [r for r in cur if r.get("key") != rec.get("key")]
    cur.append(rec)
    p.write_text(json.dumps(cur, indent=1, default=str), encoding="utf-8")

def figure_layout_audit(fig):
    """Return (min_font_pt, n_overlapping_text_pairs, n_texts) for the visible text of a drawn figure.
    Overlap = two non-empty text boxes intersecting by more than 15% of the smaller box (tick labels,
    titles, legends, annotations). Used to flag, not to fix, layout problems (AUD-15)."""
    import matplotlib.text as mtext
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    boxes, sizes, rots = [], [], []
    for t in fig.findobj(mtext.Text):
        if not t.get_visible() or not str(t.get_text()).strip() or t.get_alpha() == 0: continue
        try: bb = t.get_window_extent(r)
        except Exception: continue
        if bb.width < 1 or bb.height < 1: continue
        boxes.append(bb); sizes.append(float(t.get_fontsize())); rots.append(float(t.get_rotation()) % 180)
    ov = 0; dpi = fig.dpi
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            if rots[i] == rots[j] and rots[i] % 90 != 0:   # parallel rotated labels: use separation normal to the text direction
                th = np.radians(rots[i]); nx, ny = -np.sin(th), np.cos(th)
                sep = abs((b.x0 + b.x1 - a.x0 - a.x1) / 2 * nx + (b.y0 + b.y1 - a.y0 - a.y1) / 2 * ny)
                if sep >= 0.85 * max(sizes[i], sizes[j]) * dpi / 72: continue
            w = min(a.x1, b.x1) - max(a.x0, b.x0); h = min(a.y1, b.y1) - max(a.y0, b.y0)
            if w > 0 and h > 0 and w * h > 0.15 * min(a.width * a.height, b.width * b.height): ov += 1
    return (min(sizes) if sizes else np.nan), ov, len(boxes)

def save_fig(fig, num: int, slug: str, title: str, caption: str, nb: str, objective: str, source: str):
    """UPDATE 2026-10-31 (AUD-15): figures are laid out at final print width (<= 175 mm) and saved as
    RGB PNG (300 dpi) plus LZW-compressed RGB TIFF (600 dpi, combination artwork) in figures/tiff/."""
    from PIL import Image
    from matplotlib.ticker import NullFormatter, LogLocator, FuncFormatter
    def _lab(v, _):
        if v <= 0: return ""
        e = np.floor(np.log10(v)); m = v / 10 ** e
        return f"{v:g}" if 1e-3 <= v < 1e4 else (f"$10^{{{int(e)}}}$" if abs(m - 1) < 1e-9 else f"${m:g}\\times10^{{{int(e)}}}$")
    for ax_ in fig.axes:   # AUD-15: short-span log axes - label 1-2-5 steps, suppress crowded minor labels
        for axis, lim in ((ax_.xaxis, ax_.get_xlim()), (ax_.yaxis, ax_.get_ylim())):
            if axis.get_scale() == "log" and min(lim) > 0:
                span = np.log10(max(lim) / min(lim)); axis.set_minor_formatter(NullFormatter())
                if span < 2.5: axis.set_major_locator(LogLocator(10, subs=(1.0, 2.0, 5.0))); axis.set_major_formatter(FuncFormatter(_lab))
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
    return fn

def save_table(df: pd.DataFrame, name: str, title: str, nb: str, description: str = "",
               sheet: str | None = None, index=False):
    fn = f"{name}.csv"
    df.to_csv(TAB / fn, index=index, encoding="utf-8")
    _append("tables", nb, dict(key=fn, file=fn, title=title, notebook=nb, description=description,
                               sheet=sheet, rows=int(len(df)), cols=int(df.shape[1])))
    print(f"[table] {title}: {len(df)} rows -> tables/{fn}")
    return fn

def record_test(nb, analysis, test, n, statistic=np.nan, df=None, p=np.nan, effect_name="",
                effect=np.nan, ci_low=np.nan, ci_high=np.nan, family="", note="", objective=""):
    key = f"{analysis}|{test}|{family}"
    _append("tests", nb, dict(key=key, notebook=nb, objective=objective, analysis=analysis, test=test,
                              family=family, n=n, statistic=statistic, df=df, p_value=p,
                              effect_size_name=effect_name, effect_size=effect, ci_low=ci_low,
                              ci_high=ci_high, note=note))

def record_exclusion(nb, analysis, unit, n_excluded, reason, ids_example=""):
    _append("exclusions", nb, dict(key=f"{analysis}|{reason}", notebook=nb, analysis=analysis, unit=unit,
                                   n_excluded=int(n_excluded), reason=reason, example_ids=ids_example))

def record_finding(nb, fid, objective, finding, n, denominator, evidence, ref, interpretation, limitation,
                   kind="Empirical (completed analysis)"):
    _append("findings", nb, dict(key=fid, id=fid, notebook=nb, objective=objective, finding=finding,
                                 sample_size=n, denominator=denominator, statistical_evidence=evidence,
                                 table_or_figure=ref, interpretation=interpretation,
                                 limitations=limitation, evidence_type=kind))

def load_catalog(kind: str) -> pd.DataFrame:
    rows = []
    for p in sorted(CAT.glob(f"{kind}_*.json")):
        rows += json.loads(p.read_text(encoding="utf-8"))
    return pd.DataFrame(rows)

# ---------------------------------------------------------------- statistics
def wilson(k, n, alpha=0.05):
    from statsmodels.stats.proportion import proportion_confint
    if n == 0: return (np.nan, np.nan)
    return proportion_confint(k, n, alpha=alpha, method="wilson")

def bh(p):
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(p, float); out = np.full(p.shape, np.nan); m = ~np.isnan(p)
    if m.any(): out[m] = multipletests(p[m], method="fdr_bh")[1]
    return out

def holm(p):
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(p, float); out = np.full(p.shape, np.nan); m = ~np.isnan(p)
    if m.any(): out[m] = multipletests(p[m], method="holm")[1]
    return out

def kruskal_eps2(df, value, group):
    d = df[[value, group]].dropna()
    groups = [g[value].values for _, g in d.groupby(group, observed=True) if len(g) > 0]
    k, n = len(groups), len(d)
    H, p = stats.kruskal(*groups)
    eps2 = (H - k + 1) / (n - k)
    return dict(H=H, df=k - 1, p=p, epsilon2=eps2, n=n, k=k)

def site_level_kw(df, value, site, group, n_perm=9999, seed=SEED, agg="median"):
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

def ri_lmm(df, formula, group, reml=True):
    """Exact (RE)ML fit of a linear model with one random intercept per `group`, by profiling the
    variance ratio theta = sigma2_group / sigma2 (bounded >= 0). Robust at the boundary (theta -> 0 gives OLS),
    where statsmodels MixedLM can fail with a singular Hessian. Returns dict with fe table, variance components, llf."""
    import statsmodels.formula.api as smf
    from scipy.optimize import minimize_scalar
    mod = smf.ols(formula, df); X = np.asarray(mod.exog, float); y = np.asarray(mod.endog, float)
    rows = mod.data.row_labels if mod.data.row_labels is not None else df.index
    g = df.loc[rows, group].astype(str).values; names = list(mod.exog_names)
    N, p = X.shape; ug, gi = np.unique(g, return_inverse=True); ng = np.bincount(gi)
    def transform(theta):
        c = 1 - 1 / np.sqrt(1 + ng * theta)
        ym = np.bincount(gi, y) / ng; Xm = np.vstack([np.bincount(gi, X[:, j]) / ng for j in range(p)]).T
        return y - c[gi] * ym[gi], X - c[gi][:, None] * Xm[gi]
    def crit(theta, full=False):
        ys, Xs = transform(theta)
        XtX = Xs.T @ Xs; beta = np.linalg.lstsq(Xs, ys, rcond=None)[0]; rss = float(((ys - Xs @ beta) ** 2).sum())
        ld = float(np.log1p(ng * theta).sum())
        if reml:
            s2 = rss / (N - p); sign, ldx = np.linalg.slogdet(XtX)
            ll = -0.5 * ((N - p) * np.log(2 * np.pi * s2) + ld + ldx + (N - p))
        else:
            s2 = rss / N; ll = -0.5 * (N * np.log(2 * np.pi * s2) + ld + N)
        return (ll, beta, s2, XtX) if full else -ll
    best = None
    for hi in (1e-2, 1, 1e2, 1e4):
        r = minimize_scalar(crit, bounds=(0, hi), method="bounded", options={"xatol": 1e-10})
        if best is None or r.fun < best.fun: best = r
    th = float(best.x) if crit(best.x) < crit(0.0) else 0.0
    ll, beta, s2, XtX = crit(th, full=True)
    cov = s2 * np.linalg.pinv(XtX); se = np.sqrt(np.diag(cov)); z = beta / se
    fe = pd.DataFrame({"term": names, "estimate": beta, "SE": se, "ci_low": beta - 1.959964 * se, "ci_high": beta + 1.959964 * se,
                       "p_value": 2 * stats.norm.sf(np.abs(z))})
    return dict(fe=fe, sigma2=s2, sigma2_group=th * s2, icc=th / (1 + th), llf=ll, n=N, n_groups=len(ug), k_fe=p)

def mixed_site_lrt(df, formula, terms, group):
    """AUD-11: linear mixed model with a random intercept for `group` (site). Coefficients from REML;
    each listed fixed term tested by a maximum-likelihood ratio test against the model without it.
    Implemented with ri_lmm (exact profiled likelihood); statsmodels MixedLM is used as a cross-check where it converges."""
    lhs, rhs = formula.split("~"); parts = [p.strip() for p in rhs.split("+")]
    full = ri_lmm(df, formula, group, reml=False); out = {}
    for t in terms:
        red = [p for p in parts if p != t] or ["1"]
        m0 = ri_lmm(df, lhs + "~ " + " + ".join(red), group, reml=False)
        LR = max(0.0, 2 * (full["llf"] - m0["llf"])); dfl = full["k_fe"] - m0["k_fe"]
        out[t] = dict(LR=LR, df=dfl, p=stats.chi2.sf(LR, dfl))
    reml = ri_lmm(df, formula, group, reml=True)
    coef = pd.concat([reml["fe"], pd.DataFrame({"term": ["Site variance (random intercept)", "Residual variance", "ICC (site)"],
                                                "estimate": [reml["sigma2_group"], reml["sigma2"], reml["icc"]]})], ignore_index=True)
    sm_check = np.nan
    try:
        import statsmodels.formula.api as smf
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ms = smf.mixedlm(formula, df, groups=df[group]).fit(reml=True)
        sm_check = float(np.max(np.abs(ms.fe_params.values - reml["fe"].estimate.values)))
    except Exception:
        pass
    return dict(lrt=out, coef=coef, icc=reml["icc"], n=reml["n"], n_groups=reml["n_groups"], statsmodels_max_abs_diff=sm_check)

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

def cliffs_delta_from_u(u, n1, n2):
    return 2 * u / (n1 * n2) - 1

def pairwise_mwu(df, value, group, adjust="holm"):
    d = df[[value, group]].dropna()
    levels = [l for l in d[group].unique()]
    rows = []
    for i in range(len(levels)):
        for j in range(i + 1, len(levels)):
            a = d.loc[d[group] == levels[i], value].values; b = d.loc[d[group] == levels[j], value].values
            if len(a) < 2 or len(b) < 2: continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append(dict(group_1=levels[i], group_2=levels[j], n_1=len(a), n_2=len(b),
                             median_1=np.median(a), median_2=np.median(b), U=u, p_value=p,
                             cliffs_delta=cliffs_delta_from_u(u, len(a), len(b))))
    out = pd.DataFrame(rows)
    if len(out): out["p_adj_" + adjust] = holm(out.p_value) if adjust == "holm" else bh(out.p_value)
    return out

def spearman_ci(x, y, alpha=0.05):
    x = np.asarray(x, float); y = np.asarray(y, float); m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]; n = len(x)
    r, p = stats.spearmanr(x, y)
    se = np.sqrt((1 + r ** 2 / 2) / (n - 3))  # Bonett & Wright (2000)
    z = np.arctanh(r); zc = stats.norm.ppf(1 - alpha / 2)
    return dict(rho=r, p=p, n=n, ci_low=np.tanh(z - zc * se), ci_high=np.tanh(z + zc * se))

def boot_ci(x, func=np.median, n_boot=2000, seed=SEED):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    if len(x) < 3: return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    b = [func(rng.choice(x, len(x), replace=True)) for _ in range(n_boot)]
    return tuple(np.percentile(b, [2.5, 97.5]))

def closure(m):
    m = np.asarray(m, float); s = m.sum(axis=1, keepdims=True)
    return np.divide(m, s, out=np.zeros_like(m), where=s > 0)

def mult_replace(prop, delta=None):
    """Multiplicative zero replacement (Martin-Fernandez et al. 2003) on closed rows."""
    p = closure(prop)
    if delta is None:
        nz = p[p > 0]; delta = 0.65 * nz.min()
    z = p == 0; nz_count = z.sum(axis=1, keepdims=True)
    out = np.where(z, delta, p * (1 - nz_count * delta))
    return closure(out), delta

def clr(prop):
    lp = np.log(prop)
    return lp - lp.mean(axis=1, keepdims=True)

def shannon(counts):
    c = np.asarray(counts, float); c = c[c > 0]
    if c.sum() == 0: return np.nan
    p = c / c.sum(); return float(-(p * np.log(p)).sum())

def gini_simpson(counts):
    c = np.asarray(counts, float); c = c[c > 0]
    if c.sum() == 0: return np.nan
    p = c / c.sum(); return float(1 - (p ** 2).sum())

# ------------------------------------------------- distance-based statistics
def gower_center(D):
    A = -0.5 * np.asarray(D, float) ** 2
    n = A.shape[0]; J = np.eye(n) - np.ones((n, n)) / n
    return J @ A @ J

def _design(df, cols):
    if not cols: return np.ones((len(df), 1))
    X = pd.get_dummies(df[cols].astype(str), drop_first=True).astype(float).values
    return np.column_stack([np.ones(len(df)), X])

def permanova_dbrda(D, meta, term, covariates=(), strata=None, block=None, n_perm=999, seed=SEED):
    """Distance-based PERMANOVA (McArdle & Anderson 2001) for one term, adjusted for
    covariates (sequential). Permutation of the term's rows within `strata`, or of
    whole `block`s when the term is constant within blocks (e.g. study-level factors)."""
    rng = np.random.default_rng(seed)
    meta = meta.reset_index(drop=True)
    G = gower_center(D); N = G.shape[0]; SST = np.trace(G)
    Xc = _design(meta, list(covariates))
    Qc, _ = np.linalg.qr(Xc); Rc = np.eye(N) - Qc @ Qc.T
    Gr = Rc @ G @ Rc
    Xt_full = pd.get_dummies(meta[term].astype(str), drop_first=True).astype(float).values
    dft = Xt_full.shape[1]; dfres = N - Xc.shape[1] - dft

    def F_of(Xt):
        Xr = Rc @ Xt
        Q, R = np.linalg.qr(Xr)
        keep = np.abs(np.diag(R)) > 1e-10
        Q = Q[:, keep]
        ss = float(np.sum((Gr @ Q) * Q))
        ssr = float(np.trace(Gr)) - ss
        return (ss / dft) / (ssr / dfres), ss

    F0, ss0 = F_of(Xt_full)
    labels = meta[term].astype(str).values
    cnt = 0
    for _ in range(n_perm):
        if block is not None:
            b = meta[block].astype(str).values
            ub = pd.unique(b); blab = pd.Series(labels).groupby(b).first()
            perm_lab = dict(zip(ub, rng.permutation(blab.loc[ub].values)))
            lab_p = np.array([perm_lab[x] for x in b])
        elif strata is not None:
            lab_p = labels.copy(); s = meta[strata].astype(str).values
            for u in np.unique(s):
                ix = np.where(s == u)[0]; lab_p[ix] = rng.permutation(labels[ix])
        else:
            lab_p = rng.permutation(labels)
        Xp = pd.get_dummies(pd.Series(lab_p), drop_first=True).astype(float).values
        if Xp.shape[1] != dft: continue
        Fp, _ = F_of(Xp)
        cnt += Fp >= F0
    p = (cnt + 1) / (n_perm + 1)
    return dict(term=term, covariates=";".join(covariates) or "(none)", pseudo_F=F0, df_term=dft,
                df_resid=dfres, R2_partial=ss0 / SST, p_value=p, n=N, permutations=n_perm,
                permutation_scheme=("blocks: " + block) if block else (("within " + strata) if strata else "free"))

def permdisp_test(D, groups, n_perm=999):
    from skbio import DistanceMatrix
    from skbio.stats.distance import permdisp
    ids = [str(i) for i in range(len(groups))]
    dm = DistanceMatrix(np.asarray(D, float), ids)
    r = permdisp(dm, [str(g) for g in groups], permutations=n_perm, seed=SEED)
    return dict(F=float(r["test statistic"]), p_value=float(r["p-value"]), n=int(r["sample size"]),
                groups=int(r["number of groups"]))

def env_info():
    import importlib
    pk = {}
    for m in ["pandas", "numpy", "scipy", "statsmodels", "sklearn", "matplotlib", "seaborn", "openpyxl",
              "skbio", "geopandas", "shapely", "networkx", "nbformat", "nbconvert", "nbclient", "PIL"]:
        try: pk[m] = importlib.import_module(m).__version__
        except Exception as e: pk[m] = f"not available ({e.__class__.__name__})"
    return dict(python=platform.python_version(), executable=__import__("sys").executable,
                platform=platform.platform(), machine=platform.machine(), packages=pk, seed=SEED,
                run_timestamp=dt.datetime.now().isoformat(timespec="seconds"))

def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""): h.update(b)
    return h.hexdigest()

def nb_header(nb):
    np.random.seed(SEED)  # AUDIT FIX AUD-10: seeds legacy global RNG used by seaborn strip-plot jitter
    reset_catalog(nb)
    warnings.filterwarnings("ignore", category=FutureWarning)
    print(f"Notebook {nb} | root = {ROOT}")
    print(json.dumps(env_info()["packages"], indent=0)[:400])

def df_to_md(df: pd.DataFrame, floatfmt="{:.4g}") -> str:
    """Minimal Markdown table writer (no external dependency)."""
    def f(v):
        if isinstance(v, float): return floatfmt.format(v) if not np.isnan(v) else ""
        return str(v).replace("|", "/").replace("\n", " ")
    cols = [str(c) for c in df.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows(): out.append("| " + " | ".join(f(v) for v in r.values) + " |")
    return "\n".join(out)

# ------------------------------------------------- AMR class harmonisation
_CARD_GROUP = {
    "cephalosporin": "Beta-lactam", "penam": "Beta-lactam", "carbapenem": "Beta-lactam", "monobactam": "Beta-lactam",
    "cephamycin": "Beta-lactam", "penem": "Beta-lactam", "oxacephem": "Beta-lactam",
    "aminoglycoside antibiotic": "Aminoglycoside", "tetracycline antibiotic": "Tetracycline", "glycylcycline": "Tetracycline",
    "macrolide antibiotic": "MLS", "lincosamide antibiotic": "MLS", "streptogramin antibiotic": "MLS",
    "streptogramin A antibiotic": "MLS", "streptogramin B antibiotic": "MLS", "pleuromutilin antibiotic": "MLS",
    "fluoroquinolone antibiotic": "Fluoroquinolone/quinolone", "phenicol antibiotic": "Phenicol",
    "diaminopyrimidine antibiotic": "Folate pathway antagonist", "sulfonamide antibiotic": "Folate pathway antagonist",
    "peptide antibiotic": "Peptide/polymyxin", "glycopeptide antibiotic": "Glycopeptide", "rifamycin antibiotic": "Rifamycin",
    "phosphonic acid antibiotic": "Fosfomycin", "oxazolidinone antibiotic": "Oxazolidinone",
    "disinfecting agents and antiseptics": "Biocide",
}
_RF_GROUP = {"Beta_lactam": "Beta-lactam", "Aminoglycoside": "Aminoglycoside", "Tetracycline": "Tetracycline",
             "Macrolide": "MLS", "Lincosamide": "MLS", "Streptogramin_A": "MLS", "Streptogramin_B": "MLS", "Pleuromutilin": "MLS",
             "Quinolone": "Fluoroquinolone/quinolone", "Fluoroquinolone": "Fluoroquinolone/quinolone", "Phenicol": "Phenicol",
             "Folate_pathway_antagonist": "Folate pathway antagonist", "Polymyxin": "Peptide/polymyxin", "Glycopeptide": "Glycopeptide",
             "Rifampicin": "Rifamycin", "Fosfomycin": "Fosfomycin", "Oxazolidinone": "Oxazolidinone", "Steroid_antibacterial": "Other"}

def class_group(cls: str, pipeline: str) -> str:
    """Map a source drug-class string (';'-separated) to one harmonised display group.
    Genes whose classes span more than one group are labelled 'Multi-class'.
    Used for display alongside, never for pooling, the two pipelines."""
    s = str(cls)
    if s in ("", "nan", "Missing") or s in MISSING_CODES: return "Unclassified"
    mp = _CARD_GROUP if pipeline == "CARD" else _RF_GROUP
    parts = [p.strip() for p in s.split(";") if p.strip()]
    groups = {mp.get(p, "Other") for p in parts}
    if pipeline == "CARD" and len(groups) > 1 and "Biocide" in groups: groups.discard("Biocide")
    return groups.pop() if len(groups) == 1 else "Multi-class"

CLINICAL_PANEL = ["NDM", "KPC", "OXA-48-like", "VIM", "IMP", "CTX-M", "MCR", "tet(X)", "qnr", "vanA", "sul"]
# OXA-48-like alleles: union of the earlier list and every allele annotated "OXA-48-like" in ResFinder phenotypes.txt (resfinder-4.7.2),
# excluding blaOXA-54 (annotated as natural/intrinsic in Shewanella oneidensis, not an acquired clinical determinant). Updated 31 Oct 2026.
_OXA48 = r"^blaOXA-(48|162|163|181|199|204|232|244|245|247|252|370|405|416|436|438|439|484|505|514|515|517|519|535|538|546|547|566|567)$"
def clinical_label(gene: str, family: str, pipeline: str):
    g, f = str(gene), str(family)
    if pipeline == "CARD":
        rules = [("NDM", "NDM beta-lactamase" in f), ("KPC", "KPC beta-lactamase" in f), ("OXA-48-like", "OXA-48-like" in f),
                 ("VIM", "VIM beta-lactamase" in f), ("IMP", "IMP beta-lactamase" in f), ("CTX-M", "CTX-M beta-lactamase" in f),
                 ("MCR", "MCR phosphoethanolamine" in f), ("tet(X)", g.startswith("tet(X")), ("qnr", "(qnr)" in f),
                 ("vanA", g == "vanA"), ("sul", "sulfonamide resistant sul" in f)]
    else:
        rules = [("NDM", g.startswith("blaNDM")), ("KPC", g.startswith("blaKPC")), ("OXA-48-like", re.match(_OXA48, g) is not None),
                 ("VIM", g.startswith("blaVIM")), ("IMP", g.startswith("blaIMP")), ("CTX-M", g.startswith("blaCTX-M")),
                 ("MCR", g.startswith("mcr-")), ("tet(X)", g.startswith("tet(X")), ("qnr", g.startswith("qnr")),
                 ("vanA", g.startswith("vanA") or g.startswith("VanHAX")), ("sul", re.match(r"^sul[1-4]", g) is not None)]   # ResFinder names the Tn1546 vanA operon "VanHAX"
    # NOTE (verification 31 Oct 2026): for ResFinder, pass the FULL allele name (template before the first "_", e.g. blaOXA-181, mcr-3.1);
    # the OXA-48-like regex and the MCR rule do not match family-collapsed names such as "blaOXA" or "mcr".
    for lab, ok in rules:
        if ok: return lab
    return None

def rf_gene_name(template: str) -> str:
    allele = str(template).split("_")[0]
    if allele.startswith("bla") or allele.startswith("mcr"):
        return re.sub(r"-\d+(\.\d+)?[a-z]?$", "", allele)
    return allele

_EUK_PHYLA = set("""Arthropoda Chordata Streptophyta Ascomycota Mollusca Platyhelminthes Basidiomycota Nematoda Annelida Cnidaria
Chlorophyta Preaxostyla Echinodermata Rhodophyta Apicomplexa Haptophyta Oomycota Bryozoa Ciliophora Bacillariophyta Mucoromycota
Euglenozoa Porifera Nemertea Blastocladiomycota Cryptomycota Evosea Tubulinea Heterolobosea Rotifera Cercozoa Discosea Perkinsozoa
Prasinodermophyta Sanchytriomycota Olpidiomycota Microsporidia Brachiopoda Chytridiomycota Hemichordata Zoopagomycota Fornicata
Parabasalia Priapulida Foraminifera Gastrotricha Endomyxa Sipuncula Placozoa Xenacoelomorpha Tardigrada Onychophora Kinorhyncha
Ctenophora Chaetognatha Nematomorpha Picozoa Acanthocephala Phoronida Entoprocta Gnathostomulida Dicyemida Cycliophora
Hemimastigophora Orthonectida Loricifera""".split())
_ARC_EXTRA = {"Candidatus Thermoplasmatota", "Candidatus Diapherotrites"}

def phylum_domain(name: str) -> str:
    """Assign a phylum name (NCBI taxonomy as used by Kraken 2) to a domain."""
    n = str(name)
    if n.endswith("viricota"): return "Viruses"
    if "archaeota" in n or n in _ARC_EXTRA: return "Archaea"
    if n in _EUK_PHYLA: return "Eukaryota"
    return "Bacteria"

def term_label(term: str) -> str:
    """'C(Sample_Type, Treatment('X'))[T.Air]' -> 'Air'."""
    m = re.search(r"\[T\.(.+)\]$", str(term))
    return m.group(1) if m else str(term)

def rnd(df: pd.DataFrame, n: int = 4) -> pd.DataFrame:
    """Round float columns to (n + 2) significant digits (never to fixed decimals), so that very
    small p-values and effect sizes are preserved in saved tables."""
    out = df.copy()
    for c in out.columns:
        if pd.api.types.is_float_dtype(out[c]):
            out[c] = out[c].map(lambda x: float(f"{x:.{n + 2}g}") if pd.notna(x) else x)
    return out
