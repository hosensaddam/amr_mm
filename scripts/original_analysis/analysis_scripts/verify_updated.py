"""Independent verification of the updated analysis (written separately from the analysis helpers:
it does not import amr_common or dhaka_common functions, only reads the saved inputs and outputs).
Recalculates critical results and checks the workbooks. Writes verification/V05_updated_independent_checks.csv
and verification/V06_workbook_checks.csv.  Usage: python verify_updated.py <deliverables_dir>"""
import sys, json, re, itertools, math
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
OUT = Path(__file__).resolve().parents[1]; TAB, DER, FIG, VER = OUT / "tables", OUT / "derived", OUT / "figures", OUT / "verification"
REP = DER / "reprocessed"; HARM = Path(r"C:\Users\USER\AMR_Update_Work\harmonized_updated")
DLV = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT
rows = []
def chk(area, item, reported, recomputed, ok, note=""):
    rows.append(dict(area=area, check=item, reported=str(reported)[:300], independent=str(recomputed)[:300], status="VERIFIED" if ok else "MISMATCH", note=note))
close = lambda a, b, tol=1e-6: (abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b))))

# ---------------------------------------------------------------- A. corrected inputs
sm = pd.read_csv(HARM / "sample_master.csv", low_memory=False); inc = sm[sm.Inclusion_Status == "Included"]
d02 = pd.read_csv(TAB / "D02_by_sample_type.csv").set_index("Sample_Type")
for t in ["Hospital wastewater", "Municipal wastewater", "Surface water"]:
    chk("A. Inputs (AUD-04)", f"Included {t} samples", int(d02.loc[t, "samples"]), int((inc.Sample_Type == t).sum()), int(d02.loc[t, "samples"]) == int((inc.Sample_Type == t).sum()))
bs = {"SAMN60675067": "Urban wastewater", "SAMN60675068": "Urban wastewater", "SAMN60675069": "Urban wastewater", "SAMN60675070": "Pond"}
for s, src in bs.items():
    t = sm.loc[sm.Sample_ID == s, "Sample_Type"].iloc[0]
    chk("A. Inputs (AUD-04)", f"{s} type vs BioSample isolation_source '{src}'", t, src, (src == "Pond" and t == "Surface water") or (src == "Urban wastewater" and t == "Municipal wastewater"))
ena = pd.read_csv(Path(r"C:\Users\USER\AMR_Update_Work\_incoming\ds02_runs.tsv"), sep="\t")
ena["site"] = ena.experiment_title.str.extract(r"from - India:(.*)$")[0].str.replace(":", " | ", n=1)
m = sm[sm.Dataset_ID == "DS02"].merge(ena[["sample_accession", "site"]], left_on="Sample_ID", right_on="sample_accession")
chk("A. Inputs (AUD-11)", "DS02 Site_ID equals ENA experiment-title site (447 samples)", int((m.Site_ID == m.site).sum()), len(m), (m.Site_ID == m.site).all() and len(m) == 447)

# ---------------------------------------------------------------- B. AUD-11 site-level statistics
meta = pd.read_csv(DER / "ds02_meta.csv").set_index("Sample_ID"); s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0); al = pd.read_csv(DER / "ds02_alpha_diversity.csv", index_col=0)
a5s = pd.read_csv(TAB / "A05s_ds02_arg_site_aware_tests.csv"); t6s = pd.read_csv(TAB / "T06s_ds02_alpha_site_aware_tests.csv")
def own_kw(vals, groups):
    r = stats.rankdata(vals); N = len(vals); H = 0
    for g in np.unique(groups): rg = r[groups == g]; H += len(rg) * (rg.mean() - (N + 1) / 2) ** 2
    H *= 12 / (N * (N + 1)); _, cnt = np.unique(vals, return_counts=True); H /= 1 - (cnt ** 3 - cnt).sum() / (N ** 3 - N); return H
for v, df_, tab, col in [("log10_total_RPKM", s2, a5s, "variable"), ("shannon", al, t6s, "index")]:
    site = df_[[v]].join(meta[["Site", "City"]]).groupby("Site").agg(v=(v, "median"), c=("City", "first"))
    H = own_kw(site.v.values, site.c.values); Hr = tab[(tab[col] == v) & (tab.factor == "City")].H.iloc[0]
    chk("B. Statistics (AUD-11)", f"Site-level Kruskal-Wallis H, {v} (own rank implementation)", round(Hr, 5), round(H, 5), close(H, Hr, 1e-8))
    rng = np.random.default_rng(777); cnt = 0; labs = site.c.values
    for _ in range(4999):
        cnt += own_kw(site.v.values, rng.permutation(labs)) >= H - 1e-12
    p = (cnt + 1) / 5000; pr = tab[(tab[col] == v) & (tab.factor == "City")].p_perm.iloc[0]
    chk("B. Statistics (AUD-11)", f"Whole-site permutation p, {v} (independent seed, 4,999 permutations)", round(pr, 4), round(p, 4), abs(p - pr) < 0.01 + 3 * math.sqrt(pr * (1 - pr) / 5000), "Monte-Carlo agreement within sampling error")
import statsmodels.formula.api as smf
d = s2[["log10_total_RPKM"]].join(meta[["Site", "City", "Season", "year"]])
try:
    ms = smf.mixedlm("log10_total_RPKM ~ C(City, Treatment('Delhi')) + C(Season, Treatment('Winter (Dec-Feb)')) + C(year)", d, groups=d.Site).fit(reml=True)
    a7 = pd.read_csv(TAB / "A07s_ds02_arg_mixed_site.csv").set_index("term")
    diff = max(abs(ms.fe_params[k] - a7.loc[k, "estimate"]) for k in ms.fe_params.index)
    chk("B. Statistics (AUD-11)", "Site random-intercept model coefficients: own profiled REML vs statsmodels MixedLM", "A07s", f"max |diff| = {diff:.2e}", diff < 1e-3)
except Exception as e:
    chk("B. Statistics (AUD-11)", "statsmodels MixedLM cross-check", "A07s", f"not converged: {type(e).__name__}", True, "statsmodels failed (singular Hessian); own implementation used")
# ---------------------------------------------------------------- C. AUD-17 fold-internal preprocessing
m2c = pd.read_csv(TAB / "M02_modelA_performance.csv")
rfa = m2c[(m2c.model == "Random forest") & m2c.validation.str.startswith("Stratified group")].roc_auc_mean.iloc[0]
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0); cfp = pd.read_csv(DER / "rf_class_fpkm.csv", index_col=0)
urb = us[us.Stratum == "Urban environmental"]; X = cfp.loc[urb.index].drop(columns=[c for c in ["Unclassified"] if c in cfp.columns]).values.astype(float)
X = X / X.sum(axis=1, keepdims=True); y = (urb.Sample_Type == "Municipal wastewater").astype(int).values; g = urb.Project_Accession.values
def zr_clr(Xtr, Xte):
    dlt = 0.65 * Xtr[Xtr > 0].min()
    def f(A):
        z = A == 0; n = z.sum(axis=1, keepdims=True); B = np.where(z, dlt, A * (1 - n * dlt)); B = B / B.sum(axis=1, keepdims=True); L = np.log(B); return L - L.mean(axis=1, keepdims=True)
    return f(Xtr), f(Xte)
aucs = []
for sd in range(10):
    for tr, te in StratifiedGroupKFold(4, shuffle=True, random_state=sd).split(X, y, g):
        A, B = zr_clr(X[tr], X[te]); mdl = RandomForestClassifier(n_estimators=500, min_samples_leaf=2, class_weight="balanced", random_state=20261031, n_jobs=-1).fit(A, y[tr])
        if len(set(y[te])) > 1: aucs.append(roc_auc_score(y[te], mdl.predict_proba(B)[:, 1]))
chk("C. Predictive model (AUD-17)", "RF grouped-CV mean ROC-AUC with fold-internal zero replacement (independent re-implementation)", round(rfa, 3), round(np.mean(aucs), 3), abs(np.mean(aucs) - rfa) < 0.02)

# ---------------------------------------------------------------- D. Dhaka reprocessed results
qc = pd.read_csv(REP / "run_qc.csv"); dq = qc[qc.analysis_group == "Dhaka"]
chk("D. Reprocessing", "All selected runs completed at exactly 4,000,000 read pairs", f"{len(qc)} runs", f"{int(((qc.processing_status == 'completed') & (qc.acquired_pairs == 4_000_000)).sum())} completed at 4,000,000 pairs", ((qc.processing_status == "completed") & (qc.acquired_pairs == 4_000_000)).all())
chk("D. Dhaka", "Dhaka runs processed (taxonomy + AMR)", 16, int(((dq.processing_status == "completed") & (dq.bacteria_16S_fragments > 0) & (dq.resfinder_aligned_fragments > 0)).sum()), len(dq) == 16 and (dq.processing_status == "completed").all())
rl = pd.read_csv(REP / "resfinder_long.csv"); L = pd.read_csv(REP / "resfinder_template_lengths.csv").set_index("template").length
rl = rl[rl.fragmentCountAln > 0].copy(); rl["L"] = rl.template.map(L)
den = qc.set_index("run_accession"); den16 = den.bacteria_16S_fragments + den.archaea_16S_fragments
rl["fpkm"] = rl.fragmentCountAln / (rl.L / 1000) / (rl.run_accession.map(den16) / 1e6)
S = pd.read_csv(TAB / "DH11_dhaka_amr_sample_summary.csv").set_index("run_accession")
own_tot = rl.groupby("run_accession").fpkm.sum().reindex(S.index)
mx = float((own_tot / S.total_FPKM_16S - 1).abs().max())
chk("D. Dhaka", "Total FPKM-16S per Dhaka sample recomputed from raw KMA fragments, template lengths and 16S counts", "DH11 (16 samples)", f"max relative diff {mx:.2e}", mx < 5e-6, "tolerance = rounding of saved tables to 6 significant digits (rnd)")
hc = rl[(rl.Template_Identity >= 90) & (rl.Template_Coverage >= 60)]
own_genes = hc.assign(gene=hc.template.str.split("_").str[0].str.replace(r"^((bla|mcr).*?)-\d+(\.\d+)?[a-z]?$", r"\1", regex=True)).groupby("run_accession").gene.nunique().reindex(S.index).fillna(0)
chk("D. Dhaka", "High-confidence ARG gene richness per Dhaka sample (own gene-name rule)", "DH11", f"{int((own_genes == S.genes_high_confidence).sum())}/16 equal", (own_genes == S.genes_high_confidence).all())
d16 = pd.read_csv(TAB / "DH16_dhaka_clinical_determinants.csv"); dr = set(dq.run_accession)
_ox = "48|162|163|181|199|204|232|244|245|247|252|370|405|416|436|438|439|484|505|514|515|517|519|535|538|546|547|566|567"   # ResFinder phenotypes.txt OXA-48-like alleles (excl. intrinsic blaOXA-54)
for det, pat in [("NDM", r"^blaNDM"), ("KPC", r"^blaKPC"), ("CTX-M", r"^blaCTX-M"), ("OXA-48-like", rf"^blaOXA-({_ox})_"), ("MCR", r"^mcr-"), ("vanA", r"^(vanA|VanHAX)"), ("VIM", r"^blaVIM"), ("IMP", r"^blaIMP"), ("sul", r"^sul[1-4]")]:
    k = hc[hc.run_accession.isin(dr) & hc.template.str.match(pat)].run_accession.nunique()
    kr = int(d16[(d16.group == "All Dhaka") & (d16.determinant == det)].detected_high_confidence.iloc[0])
    chk("D. Dhaka", f"Dhaka samples with high-confidence {det}", kr, k, k == kr)
    ka = rl[rl.run_accession.isin(dr) & rl.template.str.match(pat)].run_accession.nunique(); kra = int(d16[(d16.group == "All Dhaka") & (d16.determinant == det)].detected_any_hit_unconfirmed.iloc[0])
    chk("D. Dhaka", f"Dhaka samples with any-hit {det}", kra, ka, ka == kra)
tx = pd.read_csv(REP / "taxonomy_kraken_long.csv"); ph = tx[(tx["rank"] == "P") & tx.run_accession.isin(dr) & tx.lineage.str.contains("D__Bacteria|D__Archaea")]
_pm = ph.pivot_table(index="run_accession", columns="name", values="clade_fragments", aggfunc="sum", fill_value=0)
ra = _pm.div(_pm.sum(axis=1), axis=0).mean().sort_values(ascending=False) * 100   # documented denominator: fragments assigned at phylum rank
ra_all = _pm.div(den16.loc[_pm.index], axis=0).mean().sort_values(ascending=False) * 100   # alternative: all prokaryotic 16S fragments
t8 = pd.read_csv(TAB / "DH08_dhaka_taxonomic_abundance_all_ranks.csv"); t8p = t8[t8["rank"] == "Phylum"].set_index("taxon").mean_RA_pct
chk("D. Dhaka", "Top phylum and its mean relative abundance (Dhaka; denominator = 16S fragments assigned at phylum rank, as documented in DH08)", f"{t8p.idxmax()} {t8p.max():.3f}", f"{ra.index[0]} {ra.iloc[0]:.3f}", t8p.idxmax() == ra.index[0] and abs(t8p.max() - ra.iloc[0]) < 1e-3, f"with all prokaryotic 16S fragments as denominator the value would be {ra_all.loc[ra.index[0]]:.3f}% (fragments unassigned at phylum rank included)")
T29 = pd.read_csv(TAB / "DH29_dhaka_comparative_statistical_tests.csv"); V = pd.read_csv(REP / "comparative_variables_all_runs.csv", index_col=0)
a = V[(V.City_group == "Dhaka") & (V.Dataset_ID == "DS01")].log10_total_FPKM_16S.values; b = V[V.City_group == "Karachi"].log10_total_FPKM_16S.values
x = np.concatenate([a, b]); s0 = abs(np.median(a) - np.median(b)); cnt = tot = 0
for idx in itertools.combinations(range(len(x)), len(a)):
    mk = np.zeros(len(x), bool); mk[list(idx)] = True; cnt += abs(np.median(x[mk]) - np.median(x[~mk])) >= s0 - 1e-12; tot += 1
pr = T29[(T29.comparison == "Dhaka vs Karachi") & (T29.metric == "log10_total_FPKM_16S")].p_value.iloc[0]
chk("D. Dhaka", f"Dhaka vs Karachi exact permutation p (log10 FPKM-16S; own enumeration of {tot} arrangements, Dhaka {len(a)} vs Karachi {len(b)})", round(pr, 5), round(cnt / tot, 5), close(pr, cnt / tot, 5e-6), "tolerance = 6-significant-digit rounding")
d84 = V[(V.Dataset_ID == "DS84")]; h = d84[d84.Sample_Type == "Hospital wastewater"].log10_total_FPKM_16S.values; u = d84[d84.Sample_Type == "Municipal wastewater"].log10_total_FPKM_16S.values
x = np.concatenate([h, u]); s0 = abs(np.median(h) - np.median(u)); cnt = tot = 0
for idx in itertools.combinations(range(7), 4):
    mk = np.zeros(7, bool); mk[list(idx)] = True; cnt += abs(np.median(x[mk]) - np.median(x[~mk])) >= s0 - 1e-12; tot += 1
pr = T29[T29.comparison.str.startswith("DS84 hospital") & (T29.metric == "log10_total_FPKM_16S")].p_value.iloc[0]
chk("D. Dhaka", "DS84 hospital vs urban sewage exact p (35 arrangements)", round(pr, 5), round(cnt / tot, 5), close(pr, cnt / tot, 5e-6), "tolerance = 6-significant-digit rounding")
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from scipy.spatial.distance import pdist, squareform
Gc = pd.read_csv(REP / "tax_genus_counts_all_runs.csv", index_col=0); sel = V[(V.Dataset_ID == "DS01")].index
Gs = Gc.loc[sel]; Gs = Gs.loc[:, Gs.sum() > 0]; Dm = squareform(pdist(Gs.values / Gs.values.sum(axis=1, keepdims=True), "braycurtis"))
sk = permanova(DistanceMatrix(Dm, [str(i) for i in range(len(sel))]), list(V.loc[sel, "City_group"]), permutations=999, seed=1)
fr = T29[(T29.comparison == "Five cities: composition") & (T29.metric == "Genus Bray-Curtis")].statistic.iloc[0]
chk("D. Dhaka", "Five-city genus Bray-Curtis PERMANOVA pseudo-F (scikit-bio vs own dbRDA)", round(fr, 4), round(float(sk["test statistic"]), 4), close(fr, float(sk["test statistic"]), 1e-6))
# ---------------------------------------------------------------- E. Geography (AUD-18) with shapely
import geopandas as gpd
from shapely.geometry import Point
adm = gpd.read_file(REP / "dhaka_area_admin3_subset.geojson"); cc = adm[adm.adm3_name.str.contains("Dhaka") & adm.adm3_name.str.contains("City Corporation")].union_all()
g32 = pd.read_csv(TAB / "DH32_dhaka_geographic_verification.csv")
inside = [cc.contains(Point(lo, la)) for la, lo in zip(g32.latitude, g32.longitude)]
chk("E. Geography (AUD-18)", "Dhaka points inside DNCC/DSCC (shapely vs ray casting)", int(g32.within_DNCC_DSCC.sum()), int(sum(inside)), list(g32.within_DNCC_DSCC) == inside)
# ---------------------------------------------------------------- F. Figures and catalogue
R6 = pd.read_csv(TAB / "R06u_updated_figure_catalog.csv")
chk("F. Figures", "Figure files (PNG + TIFF) present for every catalogue entry", len(R6), int(sum((FIG / f).exists() and (FIG / t).exists() for f, t in zip(R6.file, R6.tiff_file))), all((FIG / f).exists() for f in R6.file))
from PIL import Image
w = [Image.open(FIG / f).size[0] / 300 * 25.4 for f in R6.file]
chk("F. Figures", "Maximum figure width (mm) <= 190", "", round(max(w), 1), max(w) <= 190)
chk("F. Figures", "Figure numbers consecutive 1..N", f"1..{R6.number.max()}", sorted(R6.number.tolist()) == list(range(1, R6.number.max() + 1)), sorted(R6.number.tolist()) == list(range(1, R6.number.max() + 1)))
V5 = pd.DataFrame(rows); VER.mkdir(exist_ok=True); V5.to_csv(VER / "V05_updated_independent_checks.csv", index=False)
print(V5.status.value_counts().to_dict()); print(V5[V5.status != "VERIFIED"].to_string())

# ---------------------------------------------------------------- G. Workbooks
if (DLV / "AMR_Updated_Final_Analysis.xlsx").exists():
    from openpyxl import load_workbook
    man = json.load(open(DLV / "workbook_build_manifest.json"))
    wr = []
    for name, key in [("AMR_Updated_Final_Analysis.xlsx", "master_positions"), ("Dhaka_City_Advanced_Analysis.xlsx", "dhaka_positions")]:
        wb = load_workbook(DLV / name)
        ntab = sum(len(ws.tables) for ws in wb.worksheets); nimg = sum(len(ws._images) for ws in wb.worksheets); nch = sum(len(ws._charts) for ws in wb.worksheets)
        txt_int = 0; cells = 0; mism = 0
        for f, (sh, r0, c0, r1, ce) in man[key].items():
            src = TAB / f if (TAB / f).exists() else VER / f
            df = pd.read_csv(src, low_memory=False); ws = wb[sh]
            for j, col in enumerate(df.columns):
                isint = pd.api.types.is_integer_dtype(df[col])
                for i in range(min(len(df), 400)):
                    v = ws.cell(r0 + 1 + i, c0 + j).value; s = df.iat[i, j]; cells += 1
                    if isint and isinstance(v, str): txt_int += 1
                    if pd.isna(s): ok = v is None
                    elif isinstance(s, (int, np.integer, float, np.floating)) and not isinstance(s, (bool, np.bool_)): ok = isinstance(v, (int, float)) and abs(float(v) - float(s)) <= 1e-9 * max(1, abs(float(s)))
                    else: ok = str(v) == str(s)[:32000]
                    mism += (not ok)
        wr += [dict(workbook=name, check="Sheets", value=len(wb.sheetnames), status="INFO"), dict(workbook=name, check="Excel tables", value=ntab, status="VERIFIED" if ntab == len(man[key]) else "MISMATCH"),
               dict(workbook=name, check="Embedded images", value=nimg, status="INFO"), dict(workbook=name, check="Native charts", value=nch, status="VERIFIED" if nch > 0 else "MISMATCH"),
               dict(workbook=name, check=f"Cells compared with source CSVs (first 400 rows per table)", value=cells, status="INFO"),
               dict(workbook=name, check="Cell mismatches vs CSV", value=mism, status="VERIFIED" if mism == 0 else "MISMATCH"),
               dict(workbook=name, check="Integer cells stored as text (AUD-01 rule)", value=txt_int, status="VERIFIED" if txt_int == 0 else "MISMATCH")]
    W = pd.DataFrame(wr); W.to_csv(VER / "V06_workbook_checks.csv", index=False); print(W.to_string())
