# %% [markdown]
# # 16 · Dhaka-centred comparative analysis (within Dhaka; Dhaka vs Karachi, Delhi and other South Asian cities)
# **Research objectives:** Objectives 1-3 for Dhaka — differences between Dhaka sample types, datasets and periods, and between Dhaka and comparable South Asian cities.
#
# **Comparability.** All runs (16 Dhaka; 17 Global Sewage Surveillance comparators from Karachi, Kathmandu, Kochi and Colombo; 12 Delhi open-drain samples from DS02) were processed identically (4,000,000 read pairs; fastp; Kraken 2/SILVA 138; KMA/ResFinder), so pipeline and database differences are removed. Study design, matrix, year and instrument still differ and are handled by design:
# * **Matched-protocol comparison (strongest):** Dhaka DS01 (4) versus Karachi (7), Kathmandu (3), Kochi (4) and Colombo (3) — all Global Sewage Surveillance municipal sewage, 2017-2019, Illumina NovaSeq 6000, same sampling and extraction protocol.
# * **Dhaka versus Delhi:** Delhi DS02 open drains (4 sites x 3 time points, 2022-2024) versus Dhaka sewage/community wastewater (8 samples, 5 sampling units). Different studies and matrices; the comparison is matched on wastewater type only, and **the site is the independent unit** (exact permutation on site medians).
# * **Within Dhaka:** the only within-study sample-type contrast is DS84 (hospital 4 vs urban sewage-treatment-plant inlet 3). DS01 contains only municipal sewage and DS03 only hospital/community wastewater, so between-dataset contrasts confound sample type with study, year and instrument and are reported descriptively.
#
# **Statistical approach.** Small groups: exact (or Monte-Carlo) permutation tests of median differences, Mann-Whitney U, Kruskal-Wallis with epsilon-squared, Cliff's delta, ratio of medians with bootstrap 95 % CI; PERMANOVA/PERMDISP for composition; Benjamini-Hochberg within each comparison family. Each table states the independent unit. No result is presented as significant unless supported by the design-appropriate test.

# %%
#BOOTSTRAP
NB = "NB16"
nb_header(NB)
from dhaka_common import *
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
meta = load_run_meta(); meta = meta[(meta.processing_status == "completed") & ~meta.analysis_group.str.startswith("Calibration")].copy()
meta["Unit"] = meta.Sample_ID
meta.loc[meta.Dataset_ID == "DS02", "Unit"] = meta.Site_ID          # Delhi: open-drain site
meta.loc[(meta.City_group == "Dhaka") & (meta.Dataset_ID == "DS01"), "Unit"] = "Dhaka DS01 sewage (site not reported)"
for c in ["Karachi", "Kathmandu", "Kochi", "Colombo"]:
    meta.loc[meta.City_group == c, "Unit"] = meta.loc[meta.City_group == c, "Sample_ID"]
runs = meta.index
Gc = pd.read_csv(REP / "tax_genus_counts_all_runs.csv", index_col=0).loc[runs]; Gc = Gc.loc[:, Gc.sum() > 0]
Pc = pd.read_csv(REP / "tax_phylum_counts_all_runs.csv", index_col=0).loc[runs]
rf = pd.read_csv(REP / "resfinder_annotated_all_runs.csv"); rf = rf[rf.run_accession.isin(runs)]
S = pd.read_csv(REP / "amr_summary_all_runs.csv", index_col=0).loc[runs]
GF = gene_matrix(rf, runs, "FPKM_16S"); GH = gene_matrix(rf, runs, "FPKM_16S", tier="high")
depth = int(Gc.sum(axis=1).min()); Rr = rarefy_counts(Gc.round().astype(int), depth)
V = pd.DataFrame({"log10_total_FPKM_16S": np.log10(S.total_FPKM_16S), "log10_total_RPKM_depth": np.log10(S.total_RPKM_depth), "ARG_genes_high_conf": S.genes_high_confidence,
                  "ARG_gene_shannon": S.gene_shannon_FPKM, "genus_richness_rarefied": (Rr > 0).sum(axis=1), "genus_shannon": [shannon(r) for r in Gc.values]}, index=runs)
V = V.join(meta[["Sample_ID", "City_group", "Dataset_ID", "Sample_Type", "Collection_Date", "Period", "Unit", "instrument_model"]])
V.to_csv(REP / "comparative_variables_all_runs.csv")
METRICS = {"log10_total_FPKM_16S": "log10 total normalized ARG abundance (FPKM-16S)", "ARG_genes_high_conf": "ARG genes detected (high-confidence)",
           "ARG_gene_shannon": "ARG gene Shannon", "genus_richness_rarefied": f"Genus richness (rarefied to {depth:,})", "genus_shannon": "Genus Shannon"}
print("common rarefaction depth (all comparison runs):", depth); print(V.groupby(["City_group", "Dataset_ID", "Sample_Type"]).size())
tests = []
def add_test(level, comparison, metric, unit, n_a, n_b, stat_name, stat, p, p_method, effect_name, effect, extra="", status="Inferential (design-appropriate)"):
    tests.append(dict(level=level, comparison=comparison, metric=metric, independent_unit=unit, n_group_a=n_a, n_group_b=n_b, statistic_name=stat_name, statistic=stat,
                      p_value=p, p_method=p_method, effect_size_name=effect_name, effect_size=effect, detail=extra, status=status))

# %% [markdown]
# ## 6.1 Within-Dhaka comparisons

# %%
Vd = V[V.City_group == "Dhaka"]
# (a) Datasets (descriptive: dataset = study x year x instrument x site set)
for mtr in METRICS:
    k = kruskal_eps2(Vd, mtr, "Dataset_ID")
    add_test("Within Dhaka", "DS01 vs DS03 vs DS84", mtr, "sample (datasets confounded with year, instrument, sites, sample type)", 16, np.nan, "Kruskal-Wallis H", k["H"], k["p"], "asymptotic", "epsilon-squared", k["epsilon2"],
             "; ".join(f"{d}: median {Vd.loc[Vd.Dataset_ID == d, mtr].median():.3g}" for d in DHAKA_DATASETS), status="Descriptive (confounded)")
# (b) Within-study sample-type contrast: DS84 hospital (4) vs urban sewage-treatment-plant inlet (3)
d84 = Vd[Vd.Dataset_ID == "DS84"]; h = d84[d84.Sample_Type == "Hospital wastewater"]; u = d84[d84.Sample_Type == "Municipal wastewater"]
for mtr in METRICS:
    p, how = exact_perm_test(h[mtr].values, u[mtr].values)
    fc = (10 ** (h[mtr].median() - u[mtr].median())) if mtr.startswith("log10") else np.nan
    add_test("Within Dhaka", "DS84 hospital wastewater vs urban sewage (STP inlet), same study", mtr, "sample (distinct sites within one study)", len(h), len(u), "median difference (hospital - urban)", h[mtr].median() - u[mtr].median(), p, f"{how} permutation (35 arrangements)",
             "Cliff's delta", cliffs_delta(h[mtr], u[mtr]), (f"ratio of medians (FPKM) = {fc:.2f}" if mtr.startswith("log10") else ""))
# (c) Composition within DS84
ix84 = [runs.get_loc(r) for r in d84.index]
Dg = squareform(pdist(closure(Gc.values), "braycurtis")); Da = squareform(pdist(closure(GF.values), "braycurtis"))
m84 = d84[["Sample_Type"]].reset_index()
for nm, D in [("Genus Bray-Curtis", Dg), ("ARG Bray-Curtis", Da)]:
    sub = D[np.ix_(ix84, ix84)]
    keep = [i for i, t in enumerate(m84.Sample_Type) if t != "Surface water"]
    sub2 = sub[np.ix_(keep, keep)]; mm2 = m84.iloc[keep].reset_index(drop=True)
    r = permanova_dbrda(sub2, mm2, "Sample_Type", (), n_perm=9999)
    add_test("Within Dhaka", "DS84 hospital vs urban sewage: composition", nm, "sample (DS84, pond excluded)", 4, 3, "pseudo-F", r["pseudo_F"], r["p_value"], "permutation (35 distinct arrangements; p >= 1/35)", "R2", r["R2_partial"])
# (d) Co-located sampling points across datasets
coloc = [("SAMN54495504", "SAMN60675065", "Hospital site A (DS03 2025 vs DS84 2024, 0.02 km)"), ("SAMN54495503", "SAMN60675064", "Hospital site B (DS03 2025 vs DS84 2024, 0.15 km)"),
         ("SAMEA4777503", "SAMN60675067", "Pagla sewage-treatment area (DS01 2017 vs DS84 2024, ~0.13 km)")]
sid = {s: r for r, s in zip(Vd.index, Vd.Sample_ID)}
cl_rows = []
for a, b, lab_ in coloc:
    ia, ib = runs.get_loc(sid[a]), runs.get_loc(sid[b])
    ta = Vd.loc[sid[a], "Sample_Type"]
    others = [runs.get_loc(r) for r in Vd.index if r not in (sid[a], sid[b]) and Vd.loc[r, "Sample_Type"] == ta and Vd.loc[r, "Dataset_ID"] != Vd.loc[sid[a], "Dataset_ID"]]
    for nm, D in [("Genus Bray-Curtis", Dg), ("ARG Bray-Curtis", Da)]:
        cl_rows.append(dict(pair=lab_, sample_a=a, sample_b=b, distance=nm, colocated_dissimilarity=D[ia, ib],
                            median_dissimilarity_same_type_other_dataset=float(np.median([D[ia, j] for j in others])) if others else np.nan, n_reference_pairs=len(others)))
coloc_t = pd.DataFrame(cl_rows)
save_table(coloc_t.pipe(rnd, 4), "DH26_dhaka_colocated_sites", "Dhaka: dissimilarity between samples taken at the same location by different datasets (2017/2024/2025), compared with other cross-dataset pairs of the same sample type (descriptive)", NB, sheet="Dhaka_Comparative_Analysis")
# (e) Temporal (descriptive): DS01 four rounds at Dhaka municipal sewage; DS84 June-November 2024
tmp = Vd.sort_values("Collection_Date")[["Sample_ID", "Dataset_ID", "Sample_Type", "Collection_Date"] + list(METRICS)]
save_table(tmp.reset_index().pipe(rnd, 4), "DH27_dhaka_temporal_series", "Dhaka: sample metrics in collection-date order (DS01 2017-2018 rounds; DS84 Jun-Nov 2024; DS03 Jun 2025) - descriptive only; different sites between datasets", NB, sheet="Dhaka_Comparative_Analysis")
coloc_t

# %% [markdown]
# ## 6.2 Dhaka versus Karachi and other Global Sewage Surveillance cities (matched protocol, 2017-2019)

# %%
import math
gss = V[(V.Dataset_ID == "DS01")].copy()
gss["City_group"] = pd.Categorical(gss.City_group, ["Dhaka", "Karachi", "Kathmandu", "Kochi", "Colombo"])
desc_gss = gss.groupby("City_group", observed=True)[list(METRICS)].agg(["median", "min", "max", "count"]); desc_gss.columns = ["_".join(c) for c in desc_gss.columns]
save_table(desc_gss.reset_index().pipe(rnd, 4), "DH28_south_asia_sewage_city_summary", "Matched-protocol comparison (Global Sewage Surveillance, municipal sewage, 2017-2019, NovaSeq): per-city summaries of normalized ARG abundance, ARG richness/diversity and genus diversity (uniform reprocessing)", NB, sheet="Dhaka_Comparative_Analysis")
for mtr in METRICS:
    k = kruskal_eps2(gss.assign(City_group=gss.City_group.astype(str)), mtr, "City_group")
    add_test("Dhaka vs South Asian sewage (DS01)", "Five cities", mtr, "sewage sample (one per city-round; Colombo: 3 sites, 1 round; same protocol)", len(gss), np.nan, "Kruskal-Wallis H", k["H"], k["p"], "asymptotic", "epsilon-squared", k["epsilon2"])
    a = gss.loc[gss.City_group == "Dhaka", mtr].values
    for c in ["Karachi", "Kathmandu", "Kochi", "Colombo"]:
        b = gss.loc[gss.City_group == c, mtr].values
        p, how = exact_perm_test(a, b)
        if mtr.startswith("log10"):
            fc, lo, hi = fold_change_ci(10 ** a, 10 ** b); extra = f"ratio of medians Dhaka/{c} = {fc:.2f} (bootstrap 95% CI {lo:.2f}-{hi:.2f})"
        else: extra = ""
        add_test("Dhaka vs South Asian sewage (DS01)", f"Dhaka vs {c}", mtr, "sewage sample (one per city-round; Colombo: 3 sites, 1 round; same protocol)", len(a), len(b), "median difference (Dhaka - comparator)", np.median(a) - np.median(b), p, f"{how} permutation", "Cliff's delta", cliffs_delta(a, b), extra)
ixg = [runs.get_loc(r) for r in gss.index]; mg = gss[["City_group"]].astype(str).reset_index()
for nm, D in [("Genus Bray-Curtis", Dg), ("ARG Bray-Curtis", Da)]:
    sub = D[np.ix_(ixg, ixg)]
    r = permanova_dbrda(sub, mg, "City_group", (), n_perm=9999); d = permdisp_test(sub, mg.City_group.values, 9999)
    add_test("Dhaka vs South Asian sewage (DS01)", "Five cities: composition", nm, "sewage sample (one per city-round; Colombo: 3 sites, 1 round)", len(gss), np.nan, "pseudo-F", r["pseudo_F"], r["p_value"], "permutation (9,999)", "R2", r["R2_partial"], f"PERMDISP F = {d['F']:.2f}, p = {d['p_value']:.3f}")
    dk = mg.index[mg.City_group.isin(["Dhaka", "Karachi"])].values
    r2 = permanova_dbrda(sub[np.ix_(dk, dk)], mg.loc[dk].reset_index(drop=True), "City_group", (), n_perm=9999)
    n_dh, n_ka = int((mg.City_group == "Dhaka").sum()), int((mg.City_group == "Karachi").sum())
    add_test("Dhaka vs Karachi (DS01)", "Dhaka vs Karachi: composition", nm, "sewage sample (one per city-round; Colombo: 3 sites, 1 round)", n_dh, n_ka, "pseudo-F", r2["pseudo_F"], r2["p_value"], f"permutation ({math.comb(n_dh + n_ka, n_dh):,} distinct arrangements)", "R2", r2["R2_partial"])
desc_gss.iloc[:, ::4]

# %% [markdown]
# ## 6.3 Dhaka versus Delhi (uniformly reprocessed DS02 subset; site-aware)

# %%
dd = V[(V.City_group == "Dhaka") & V.Sample_Type.isin(["Municipal wastewater", "Community wastewater"])]
de = V[V.City_group == "Delhi"]
U_dh = dd.groupby("Unit")[list(METRICS)].median(); U_de = de.groupby("Unit")[list(METRICS)].median()
print("Dhaka sewage units:", list(U_dh.index)); print("Delhi site units:", list(U_de.index))
for mtr in METRICS:
    p, how = exact_perm_test(U_dh[mtr].values, U_de[mtr].values)
    extra = ""
    if mtr.startswith("log10"):
        fc, lo, hi = fold_change_ci(10 ** dd[mtr].values, 10 ** de[mtr].values); extra = f"sample-level ratio of medians Dhaka/Delhi = {fc:.2f} (bootstrap 95% CI {lo:.2f}-{hi:.2f}; descriptive)"
    add_test("Dhaka vs Delhi", "Dhaka sewage/community wastewater vs Delhi open drains", mtr, f"sampling unit (Dhaka {len(U_dh)} units; Delhi {len(U_de)} sites)", len(U_dh), len(U_de), "median difference of unit medians (Dhaka - Delhi)", U_dh[mtr].median() - U_de[mtr].median(), p, f"{how} permutation of units", "Cliff's delta (units)", cliffs_delta(U_dh[mtr], U_de[mtr]), extra)
    a_all = V.loc[V.City_group == "Dhaka", mtr]; u_, p_ = stats.mannwhitneyu(a_all, de[mtr])
    add_test("Dhaka vs Delhi", "All Dhaka samples vs Delhi open drains (sample level)", mtr, "sample (Delhi samples repeated within 4 sites; descriptive)", len(a_all), len(de), "Mann-Whitney U", u_, p_, "asymptotic", "Cliff's delta", cliffs_delta(a_all, de[mtr]), status="Descriptive (pseudo-replicated; mixed matrices)")
ixdd = [runs.get_loc(r) for r in list(dd.index) + list(de.index)]
mdd = pd.concat([dd, de])[["City_group", "Unit"]].reset_index()
for nm, D in [("Genus Bray-Curtis", Dg), ("ARG Bray-Curtis", Da)]:
    sub = D[np.ix_(ixdd, ixdd)]
    r = permanova_dbrda(sub, mdd, "City_group", (), block="Unit", n_perm=9999); d = permdisp_test(sub, mdd.City_group.values, 9999)
    add_test("Dhaka vs Delhi", "Dhaka sewage vs Delhi open drains: composition", nm, "sampling unit (whole-unit permutation)", len(dd), len(de), "pseudo-F", r["pseudo_F"], r["p_value"], f"permutation of whole units ({math.comb(len(U_dh) + len(U_de), len(U_dh)):,} arrangements)", "R2", r["R2_partial"], f"PERMDISP F = {d['F']:.2f}, p = {d['p_value']:.3f} (sample level)")
# Unit-definition sensitivity for Dhaka vs Delhi (added at verification, 31 Oct 2026): the DS01 sewage rounds have no reported site,
# and DS84 UWW1 (SAMN60675067) lies ~0.13 km from the ENA coordinate of DS01 SAMEA4777503 (Pagla sewage-treatment area, see notebook 17).
_alt = {"Primary: DS01 rounds pooled as one unit; DS84 UWW1 separate": dd.Unit.copy()}
_u = dd.Unit.copy(); _u[dd.Sample_ID == "SAMN60675067"] = "Dhaka DS01 sewage (site not reported)"
_alt["Pagla merged: DS84 UWW1 pooled with the DS01 unit"] = _u
_u = dd.Unit.copy(); _ix = dd.Dataset_ID == "DS01"; _u[_ix] = list(dd.index[_ix])
_alt["DS01 rounds as separate units (treats repeated rounds as independent)"] = _u
_u = dd.Unit.copy(); _u[dd.Sample_ID.isin(["SAMEA4777503", "SAMN60675067"])] = "Pagla STP area (DS01 SAMEA4777503 + DS84 UWW1)"; _u[(dd.Dataset_ID == "DS01") & (dd.Sample_ID != "SAMEA4777503")] = "Dhaka DS01 other rounds (site not reported)"
_alt["Pagla split: SAMEA4777503 + UWW1 as one unit; other DS01 rounds as a second unit"] = _u
_rows = []
for nm_, uu in _alt.items():
    Ud = dd.assign(Unit=uu).groupby("Unit")[list(METRICS)].median()
    for mtr in METRICS:
        p, how = exact_perm_test(Ud[mtr].values, U_de[mtr].values)
        _rows.append(dict(unit_definition=nm_, metric=mtr, n_dhaka_units=len(Ud), n_delhi_sites=len(U_de), median_difference_Dhaka_minus_Delhi=Ud[mtr].median() - U_de[mtr].median(), p_exact=p, p_method=f"{how} permutation of units"))
SENS = pd.DataFrame(_rows)
save_table(SENS.pipe(rnd, 5), "DH29b_dhaka_delhi_unit_definition_sensitivity", f"Dhaka vs Delhi: exact permutation p-values under {len(_alt)} definitions of the Dhaka sampling unit (DS01 sites unreported; DS84 UWW1 co-located with DS01 at the Pagla sewage-treatment area)", NB, sheet="Dhaka_Statistical_Tests")
_sab = SENS[SENS.metric == "log10_total_FPKM_16S"]
T = pd.DataFrame(tests)
T["p_BH_within_level"] = np.nan
for lv, ix in T[T.status.str.startswith("Inferential")].groupby("level").groups.items(): T.loc[ix, "p_BH_within_level"] = bh(T.loc[ix, "p_value"].values)
for r_ in T.itertuples():
    record_test(NB, f"{r_.level}: {r_.comparison} [{r_.metric}]", f"{r_.statistic_name} ({r_.p_method})", int(r_.n_group_a + (0 if pd.isna(r_.n_group_b) else r_.n_group_b)), r_.statistic, None, r_.p_value, r_.effect_size_name, r_.effect_size,
                family=f"Dhaka comparative: {r_.level}" + ("" if r_.status.startswith("Inferential") else " [descriptive]"), note=f"unit: {r_.independent_unit}; {r_.detail}", objective="Dhaka comparative")
save_table(T.pipe(rnd, 5), "DH29_dhaka_comparative_statistical_tests", "Dhaka comparative analysis: all tests with the independent unit, group sizes, statistic, permutation/asymptotic p-value, BH adjustment within comparison level, effect size and inferential status", NB, sheet="Dhaka_Statistical_Tests")
T[T.status.str.startswith("Inferential")][["level", "comparison", "metric", "n_group_a", "n_group_b", "p_value", "p_BH_within_level", "effect_size"]].round(3)

# %% [markdown]
# ## 6.4 Class composition and clinically important determinants by city group

# %%
grp = V.City_group.where(V.City_group != "Dhaka", "Dhaka " + V.Dataset_ID)
CLm = rf.pivot_table(index="run_accession", columns="Class_Group", values="FPKM_16S", aggfunc="sum", fill_value=0).reindex(runs).fillna(0)
CLs = 100 * CLm.div(CLm.sum(axis=1), axis=0)
cls_city = CLs.groupby(grp).mean().T
save_table(cls_city.reset_index().rename(columns={"index": "Class_Group"}).pipe(rnd, 3), "DH30_class_share_by_city_group", "Mean share (%) of total FPKM-16S by antibiotic-class group for Dhaka datasets and comparator cities (uniform reprocessing)", NB, sheet="Dhaka_Comparative_Analysis")
crow = []
for labc in CLINICAL_PANEL:
    hc_runs = set(rf[(rf.Clinical == labc) & rf.high_confidence].run_accession); any_runs = set(rf[rf.Clinical == labc].run_accession)
    for g_, ix in V.groupby(grp).groups.items():
        k1 = len(set(ix) & hc_runs); k2 = len(set(ix) & any_runs); lo, hi = wilson(k1, len(ix))
        crow.append(dict(determinant=labc, group=g_, n=len(ix), detected_high_confidence=k1, pct_high_confidence=100 * k1 / len(ix), ci_low=100 * lo, ci_high=100 * hi, detected_any_hit_unconfirmed=k2))
clin_c = pd.DataFrame(crow)
save_table(clin_c.pipe(rnd, 2), "DH31_clinical_determinants_by_city_group", "High-confidence (primary) and any-hit (unconfirmed) detection of clinically important determinants by Dhaka dataset and comparator city (uniform reprocessing; Wilson 95% CI)", NB, sheet="Dhaka_Comparative_Analysis")
cls_city.round(1)

# %%
fig, axes = plt.subplots(1, 4, figsize=(6.9, 2.5))
for ax, mtr in zip(axes, ["log10_total_FPKM_16S", "ARG_genes_high_conf", "genus_shannon", "genus_richness_rarefied"]):
    for i, (t_, col) in enumerate([(t, c) for t, c in TYPE_COL.items() if t in set(Vd.Sample_Type)]):
        g = Vd[Vd.Sample_Type == t_]
        for d_, mkr in [("DS01", "o"), ("DS03", "s"), ("DS84", "^")]:
            gg = g[g.Dataset_ID == d_]
            ax.scatter(np.full(len(gg), i) + np.random.default_rng(1).uniform(-.15, .15, len(gg)), gg[mtr], marker=mkr, s=16, color=col, edgecolors="k", linewidths=.3)
        ax.hlines(g[mtr].median(), i - .3, i + .3, color="k", lw=.8)
    ax.set_xticks(range(len([t for t in TYPE_COL if t in set(Vd.Sample_Type)]))); ax.set_xticklabels([t.split()[0][:4] + "." for t in TYPE_COL if t in set(Vd.Sample_Type)], fontsize=6)
    ax.set_ylabel(METRICS[mtr].replace(" (", "\n("), fontsize=6)
fig.legend([plt.Line2D([], [], marker=m, ls="", mfc="w", mec="k", ms=4) for m in "os^"], ["DS01", "DS03", "DS84"], frameon=False, ncol=3, fontsize=6, loc="lower center", bbox_to_anchor=(.5, -.08))
fig.tight_layout(rect=(0, .06, 1, 1))
dh84 = T[(T.comparison.str.startswith("DS84 hospital")) & (T.metric == "log10_total_FPKM_16S")].iloc[0]
save_fig(fig, 39, "dhaka_within_city_comparison", "Within-Dhaka comparison of resistome and microbiome metrics by verified sample type and dataset",
         f"Each point is one Dhaka sample (marker = dataset; colour = verified sample type: Muni. municipal/urban sewage, Hosp. hospital, Comm. community, Surf. pond); bars are medians. The only within-study sample-type contrast (DS84: 4 hospital vs 3 urban sewage-treatment-plant inlets) gave an exact permutation p = {dh84.p_value:.3f} for log10 FPKM-16S (Cliff's delta {dh84.effect_size:.2f}; Table DH29). Other contrasts confound sample type with study, year and instrument.",
         NB, "Dhaka comparative (6.1)", "tables/DH29; derived/reprocessed/comparative_variables_all_runs.csv")

# %%
order_g = ["Dhaka DS01", "Dhaka DS03", "Dhaka DS84", "Karachi", "Kathmandu", "Kochi", "Colombo", "Delhi"]
gcol = {"Dhaka DS01": "#D55E00", "Dhaka DS03": "#E8A07A", "Dhaka DS84": "#A23B00", **{k: CITY_COL[k] for k in ["Karachi", "Kathmandu", "Kochi", "Colombo", "Delhi"]}}
RAp = 100 * Pc.div(Pc.sum(axis=1), axis=0); topp = RAp.mean().sort_values(ascending=False).index[:8]
ph_g = RAp[topp].groupby(grp).mean().reindex(order_g); ph_g["Other"] = 100 - ph_g.sum(axis=1)
pc_all = pcoa(DistanceMatrix(Dg, list(runs))); sc = pc_all.samples.iloc[:, :2].values; pe = pc_all.proportion_explained.values
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.3), gridspec_kw={"width_ratios": [1.1, 1]})
left = np.zeros(len(ph_g)); pal = dict(zip(list(topp) + ["Other"], sns.color_palette("tab10", 8) + [(.82, .82, .82)]))
for c in ph_g.columns:
    axes[0].barh(range(len(ph_g)), ph_g[c], left=left, color=pal[c], label=c, edgecolor="white", linewidth=.2); left += ph_g[c].values
axes[0].set_yticks(range(len(ph_g))); axes[0].set_yticklabels([f"{g} (n={int((grp == g).sum())})" for g in ph_g.index], fontsize=6); axes[0].invert_yaxis()
axes[0].set_xlabel("Mean relative abundance\n(% of phylum-assigned prokaryotic 16S)"); axes[0].set_title("A  Phylum composition"); axes[0].legend(frameon=False, fontsize=6, ncol=3, loc="upper center", bbox_to_anchor=(.45, -.17))
for g_ in order_g:
    ix = (grp == g_).values
    axes[1].scatter(sc[ix, 0], sc[ix, 1], s=18, color=gcol[g_], edgecolors="k", linewidths=.3, label=g_, marker="^" if g_ == "Delhi" else "o")
axes[1].set_xlabel(f"PCoA1 ({100*pe[0]:.0f}%)"); axes[1].set_ylabel(f"PCoA2 ({100*pe[1]:.0f}%)"); axes[1].set_title("B  Genus Bray-Curtis PCoA"); axes[1].legend(frameon=False, fontsize=6, ncol=2, loc="upper center", bbox_to_anchor=(.5, -.17))
fig.tight_layout()
save_fig(fig, 40, "dhaka_vs_delhi_south_asia_microbiome", "Microbiome composition of Dhaka compared with Delhi and other South Asian cities (uniform reprocessing)",
         f"(A) Mean phylum relative abundance (% of prokaryotic 16S fragments assigned at phylum rank; Kraken 2/SILVA 138) by Dhaka dataset and comparator city; (B) PCoA of genus-level Bray-Curtis dissimilarity (Kraken 2 clade fragments) for all {len(runs)} comparison runs (circles: Global Sewage Surveillance municipal sewage and Dhaka datasets; triangles: Delhi DS02 open drains). All runs analysed at 4,000,000 read pairs with one pipeline; matrices, years and studies differ (see Table DH29 for design-appropriate tests).",
         NB, "Dhaka comparative (6.2)", "tables/DH28-DH30")

# %%
pa = pcoa(DistanceMatrix(Da, list(runs))); sa_ = pa.samples.iloc[:, :2].values; pea = pa.proportion_explained.values
fig, axes = plt.subplots(1, 3, figsize=(6.9, 3.9), gridspec_kw={"width_ratios": [1, 1.1, 1]})
for i, g_ in enumerate(order_g):
    v = V.loc[grp == g_, "log10_total_FPKM_16S"]
    axes[0].scatter(np.full(len(v), i) + np.random.default_rng(2).uniform(-.15, .15, len(v)), v, s=12, color=gcol[g_], edgecolors="k", linewidths=.3)
    axes[0].hlines(v.median(), i - .3, i + .3, color="k", lw=.8)
axes[0].set_xticks(range(len(order_g))); axes[0].set_xticklabels([g.replace("Dhaka ", "Dhk ") for g in order_g], rotation=60, ha="right", fontsize=6)
axes[0].set_ylabel("log10 total normalized ARG\nabundance (FPKM-16S)"); axes[0].set_title("A  ARG abundance")
cc = cls_city[order_g]; cc = cc.loc[[c for c in CLASS_ORDER if c in cc.index]]; bottom = np.zeros(len(order_g))
for c in cc.index:
    axes[1].bar(range(len(order_g)), cc.loc[c], bottom=bottom, color=CLASS_COLORS[c], label=c, edgecolor="white", linewidth=.2); bottom += cc.loc[c].values
axes[1].set_xticks(range(len(order_g))); axes[1].set_xticklabels([g.replace("Dhaka ", "Dhk ") for g in order_g], rotation=60, ha="right", fontsize=6); axes[1].set_ylim(0, 100)
axes[1].set_ylabel("Mean share of FPKM-16S (%)"); axes[1].set_title("B  Class composition")
h_, l_ = axes[1].get_legend_handles_labels()
fig.legend(h_, l_, frameon=False, fontsize=6, ncol=5, loc="lower center", bbox_to_anchor=(.5, 0.0), title="Panel B antibiotic class", title_fontsize=6.5)
for g_ in order_g:
    ix = (grp == g_).values
    axes[2].scatter(sa_[ix, 0], sa_[ix, 1], s=14, color=gcol[g_], edgecolors="k", linewidths=.3, marker="^" if g_ == "Delhi" else "o")
axes[2].set_xlabel(f"PCoA1 ({100*pea[0]:.0f}%)"); axes[2].set_ylabel(f"PCoA2 ({100*pea[1]:.0f}%)"); axes[2].set_title("C  ARG Bray-Curtis PCoA")
fig.subplots_adjust(left=.09, right=.99, top=.92, bottom=.40, wspace=.42)
dl_ = T[(T.level == "Dhaka vs Delhi") & (T.metric == "log10_total_FPKM_16S") & T.status.str.startswith("Inferential")].iloc[0]
save_fig(fig, 41, "dhaka_vs_delhi_south_asia_resistome", "Resistome of Dhaka compared with Delhi and other South Asian cities (uniform ResFinder/KMA reprocessing)",
         f"(A) log10 total normalized ARG abundance (FPKM-16S) per run with group medians; (B) mean antibiotic-class composition (share of FPKM-16S); (C) PCoA of Bray-Curtis dissimilarity of relative ARG gene abundance (triangles: Delhi DS02 open drains). Dhaka sewage versus Delhi open drains, unit-level exact permutation: p = {dl_.p_value:.3f} ({dl_.independent_unit}). Normalized relative abundance, not absolute ARG quantity; ARG detection is not phenotypic resistance.",
         NB, "Dhaka comparative (6.2)", "tables/DH29-DH31")

# %%
GSS_N = ", ".join(f"{c} {int((gss.City_group == c).sum())}" for c in ["Dhaka", "Karachi", "Kathmandu", "Kochi", "Colombo"])
fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.6))
cities5 = ["Dhaka", "Karachi", "Kathmandu", "Kochi", "Colombo"]
for ax, mtr in zip(axes, ["log10_total_FPKM_16S", "ARG_genes_high_conf", "genus_shannon"]):
    for i, c in enumerate(cities5):
        v = gss.loc[gss.City_group == c, mtr]
        ax.scatter(np.full(len(v), i) + np.random.default_rng(3).uniform(-.15, .15, len(v)), v, s=16, color=CITY_COL[c], edgecolors="k", linewidths=.3)
        ax.hlines(v.median(), i - .3, i + .3, color="k", lw=.8)
    t_ = T[(T.comparison == "Dhaka vs Karachi") & (T.metric == mtr)].iloc[0]
    ax.set_xticks(range(5)); ax.set_xticklabels([f"{c}\n({int((gss.City_group == c).sum())})" for c in cities5], fontsize=6); ax.set_ylabel(METRICS[mtr].replace(" (", "\n("), fontsize=6)
    ax.set_title(f"Dhaka vs Karachi: exact p = {t_.p_value:.3f}", fontsize=6.3)
fig.tight_layout()
save_fig(fig, 42, "dhaka_vs_karachi_matched_sewage", "Dhaka compared with Karachi and other South Asian cities under a matched sewage-surveillance protocol",
         f"Global Sewage Surveillance municipal sewage samples collected 2017-2019 and sequenced on Illumina NovaSeq 6000 ({GSS_N}), reprocessed with one pipeline at 4,000,000 read pairs. Points are sewage samples (city-rounds; the three Colombo samples are three sites sampled on one date) with medians; panel titles give the exact permutation p-value for Dhaka vs Karachi (all tests and effect sizes in Table DH29). Few samples per city; results describe these sampling rounds, not city-wide AMR burden.",
         NB, "Dhaka comparative (6.2)", "tables/DH28-DH29")

# %%
iv = T[T.status.str.startswith("Inferential")]
def pick(level, comp, mtr):
    r = T[(T.level == level) & (T.comparison == comp) & (T.metric == mtr)]
    return r.iloc[0] if len(r) else None
k5 = pick("Dhaka vs South Asian sewage (DS01)", "Five cities", "log10_total_FPKM_16S"); kk = T[(T.comparison == "Dhaka vs Karachi: composition") & (T.metric == "ARG Bray-Curtis")].iloc[0]
dkf = T[(T.comparison == "Dhaka vs Karachi") & (T.metric == "log10_total_FPKM_16S")].iloc[0]
dkg = T[(T.comparison == "Dhaka vs Karachi") & (T.metric == "ARG_genes_high_conf")].iloc[0]
record_finding(NB, "DF06", "Dhaka comparative (within Dhaka)", f"Within the one study that sampled both, DS84 hospital wastewater (n = 4) and urban sewage-treatment-plant inlets (n = 3) showed a median difference in log10 normalized ARG abundance of {dh84.statistic:.2f} log10 units, which is not statistically significant (exact permutation p = {dh84.p_value:.3f}; Cliff's delta {dh84.effect_size:.2f}); between-dataset differences in Dhaka confound sample type with study, year, instrument and sites. Samples taken at the same hospital or sewage-treatment locations by different studies (2017-2025) are listed with their dissimilarities in Table DH26.",
               16, "Dhaka samples", "Exact permutation test; Cliff's delta; descriptive Kruskal-Wallis for datasets", "Tables DH26-DH27, DH29; Figure 39",
               f"With 4 versus 3 samples there are only 35 label arrangements, so p-values move in steps of 1/35; the groups were completely separated (Cliff's delta {dh84.effect_size:.2f}) and {int(round(dh84.p_value * 35))} of the 35 arrangements were at least as extreme (p = {dh84.p_value:.3f}). The contrast indicates direction and magnitude in one study but cannot establish a general hospital-versus-sewage difference.",
               "Very small groups; one pond and one community sample; no site replication within DS84.", kind="Empirical (new reprocessing)")
record_finding(NB, "DF07", "Dhaka comparative (Dhaka vs Karachi / South Asia)", f"Under a matched Global Sewage Surveillance protocol (2017-2019), Dhaka municipal sewage had a ratio of median normalized ARG abundance relative to Karachi of {dkf.detail.split('= ')[1] if '=' in dkf.detail else 'n/a'} (exact permutation p = {dkf.p_value:.3f}; Cliff's delta {dkf.effect_size:.2f}); across five cities the Kruskal-Wallis epsilon-squared was {k5.effect_size:.2f} (p = {k5.p_value:.3f}); Dhaka-Karachi resistome composition: PERMANOVA R2 = {kk.effect_size:.2f}, p = {kk.p_value:.3f} (BH q = {kk.p_BH_within_level:.3f}; dispersion was not tested for this pair; five-city PERMDISP: {T[(T.comparison == 'Five cities: composition') & (T.metric == 'ARG Bray-Curtis')].detail.iloc[0]}). Dhaka sewage carried fewer high-confidence ARG genes than Karachi (median difference {dkg.statistic:.1f} genes; exact p = {dkg.p_value:.3f}, BH q = {dkg.p_BH_within_level:.3f}; Cliff's delta {dkg.effect_size:.2f}).",
               len(gss), f"Global Sewage Surveillance samples ({GSS_N})", "Exact permutation; Kruskal-Wallis; PERMANOVA/PERMDISP; bootstrap CI of ratio of medians", "Tables DH28-DH31; Figures 40-42",
               "This is the strongest available Dhaka-versus-city comparison because sampling programme, matrix, period, instrument and bioinformatics are matched (4 of the 7 Karachi runs come from a separate Global Sewage Surveillance monthly series, PRJEB51229; see the limitations register); it remains limited to 3-7 samples per city.",
               "Few samples per city; one or few sewage sites per city; 2017-2019 only.", kind="Empirical (new reprocessing)")
dlf = T[(T.level == "Dhaka vs Delhi") & (T.metric == "log10_total_FPKM_16S") & T.status.str.startswith("Inferential")].iloc[0]
dlc = T[(T.level == "Dhaka vs Delhi") & (T.metric == "ARG Bray-Curtis")].iloc[0]
import re as _re
_pd = float(_re.search(r"p = ([0-9.eE-]+)", str(dlc.detail)).group(1)) if _re.search(r"p = ([0-9.eE-]+)", str(dlc.detail)) else np.nan
DISP_TXT = (f"dispersion also differed ({dlc.detail}), so the composition difference may partly reflect differences in spread rather than location" if _pd < 0.05 else f"dispersion did not differ detectably ({dlc.detail})")
record_finding(NB, "DF08", "Dhaka comparative (Dhaka vs Delhi)", f"Dhaka sewage/community wastewater ({len(U_dh)} sampling units, {len(dd)} samples) and Delhi open drains ({len(U_de)} sites, {len(de)} samples) reprocessed identically differed in log10 normalized ARG abundance by a median of {dlf.statistic:.2f} log10 units between unit medians (exact permutation p = {dlf.p_value:.3f}; range {_sab.p_exact.min():.3f}-{_sab.p_exact.max():.3f} across {len(_alt)} definitions of the Dhaka sampling unit, Table DH29b; BH-adjusted within the Dhaka-Delhi family q = {dlf.p_BH_within_level:.3f}; {dlf.detail}); resistome composition PERMANOVA (whole-unit permutation) R2 = {dlc.effect_size:.2f}, p = {dlc.p_value:.3f} (q = {dlc.p_BH_within_level:.3f}); {DISP_TXT}.",
               len(dd) + len(de), "Dhaka sewage + Delhi open-drain samples", "Exact permutation on sampling-unit medians; PERMANOVA with whole-unit permutation", "Tables DH29-DH31; Figures 40-41",
               "Pipeline differences are removed, but study, matrix (sewer/STP inlet vs open drain), years (2017-2025 vs 2022-2024) and instruments differ, and the abundance p-value depends on how Dhaka sampling units are defined and does not survive BH adjustment; the contrast is a matched-matrix descriptive comparison, not evidence of a city difference.",
               f"{len(U_dh)} (Dhaka) and {len(U_de)} (Delhi) independent units; different studies and periods; Delhi subset ({len(de)} of 100 DS02 Delhi samples).", kind="Empirical (new reprocessing)")
