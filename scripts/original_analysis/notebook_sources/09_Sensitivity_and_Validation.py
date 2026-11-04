# %% [markdown]
# # 09 · Sensitivity and Robustness Analyses
# **Research objectives:** robustness of the key results for Objectives 1, 2 and 4.
#
# Each key estimate is recomputed under alternative, defensible analytic choices: detection thresholds, normalisation denominators, aggregation level (to address pseudo-replication), choice of duplicate project, exclusion of an older KMA version, and random seeds. A result is described as robust when its direction and approximate magnitude persist across specifications.

# %%
#BOOTSTRAP
NB = "NB09"
nb_header(NB)
import seaborn as sns
import statsmodels.formula.api as smf
from scipy.spatial.distance import pdist, squareform
m2 = pd.read_csv(DER / "ds02_meta.csv").set_index("Sample_ID")
al = pd.read_csv(DER / "ds02_card_allele_long.csv")
s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0).loc[m2.index]
ph = pd.read_csv(DER / "ds02_prokaryote_phylum_counts.csv", index_col=0).loc[m2.index]

# %% [markdown]
# ## 1. DS02: city contrasts in total normalized ARG abundance under alternative specifications

# %%
al["kb"] = al.Reference_Length_bp / 1000
prok_M = ph.sum(axis=1) / 1e6
variants = {
    "Primary: all RGI hits, per million sequenced reads": s2.total_RPKM,
    "Allele coverage >= 50%": s2.total_RPKM_cov50,
    "Allele coverage >= 80%": al[al.Average_Percent_Coverage >= 80].groupby("Sample_ID").RPKM.sum().reindex(m2.index),
    "Denominator: prokaryote-classified reads (Bracken)": (al.reads / al.kb).groupby(al.Sample_ID).sum().reindex(m2.index) / prok_M,
}
dq = pd.read_csv(TAB / "QC16_depth_outliers.csv")
rows, kw = [], []
def fit_city(y, data, label):
    d = data.assign(y=np.log10(y)).dropna(subset=["y"])
    f = smf.ols("y ~ C(City, Treatment('Delhi')) + C(Season) + C(year)", d).fit(cov_type="HC3")
    for c in ["Chennai", "Kolkata", "Mumbai"]:
        t = f"C(City, Treatment('Delhi'))[T.{c}]"
        rows.append(dict(specification=label, contrast=f"{c} vs Delhi", n=int(f.nobs), fold_change=10 ** f.params[t], ci_low=10 ** f.conf_int().loc[t, 0], ci_high=10 ** f.conf_int().loc[t, 1], p_value=f.pvalues[t]))
    k = kruskal_eps2(d, "y", "City"); kw.append(dict(specification=label, n=k["n"], H=k["H"], p=k["p"], epsilon2=k["epsilon2"]))
base = m2[["City", "Season", "year", "Collection_Date"]]
for lab, y in variants.items(): fit_city(y, base, lab)
fit_city(s2.total_RPKM, base[~base.index.isin(dq.Sample_ID)], "Primary, excluding depth-flagged samples")
tclus = pd.read_csv(DER / "ds02_taxonomy_clusters.csv", index_col=0).taxonomic_cluster
d_c = base.assign(y=np.log10(s2.total_RPKM), cl=tclus)
f_c = smf.ols("y ~ C(City, Treatment('Delhi')) + C(Season) + C(year) + C(cl)", d_c).fit(cov_type="HC3")
for c in ["Chennai", "Kolkata", "Mumbai"]:
    t = f"C(City, Treatment('Delhi'))[T.{c}]"
    rows.append(dict(specification="Primary, adjusted for DS02 taxonomic cluster", contrast=f"{c} vs Delhi", n=int(f_c.nobs), fold_change=10 ** f_c.params[t], ci_low=10 ** f_c.conf_int().loc[t, 0], ci_high=10 ** f_c.conf_int().loc[t, 1], p_value=f_c.pvalues[t]))
cm = base.assign(v=s2.total_RPKM).groupby(["City", "Collection_Date", "Season", "year"], as_index=False).v.apply(lambda s: np.exp(np.log(s).mean()))
fit_city(cm.v.values, cm.drop(columns="v"), "City-month geometric means (n = city-months)")
sens1 = pd.DataFrame(rows); kw = pd.DataFrame(kw)
save_table(sens1.pipe(rnd, 5), "SA01_ds02_city_contrast_sensitivity", "Sensitivity of DS02 city fold changes in total normalized ARG abundance (OLS with season and year; HC3)", NB, sheet="Sensitivity_Analysis")
save_table(kw.pipe(rnd, 6), "SA02_ds02_city_kruskal_sensitivity", "Sensitivity of the DS02 city effect (Kruskal-Wallis epsilon-squared) across specifications", NB, sheet="Sensitivity_Analysis")
for _, r in kw.iterrows():
    record_test(NB, f"Sensitivity DS02 normalized ARG abundance by city: {r.specification}", "Kruskal-Wallis", int(r.n), r.H, 3, r.p, "epsilon-squared", r.epsilon2, family="Sensitivity", objective="Objective 2")
kw.round(4)

# %% [markdown]
# ## 2. ResFinder: matrix contrasts and food-related contrast under alternative specifications

# %%
rfl = pd.read_csv(DER / "rf_template_long_all_units.csv")
ru = pd.read_csv(DER / "rf_units_meta.csv", index_col=0)
def unit_summary(df, value="FPKM", units=None):
    t = df.groupby("unit_id")[value].sum().rename("total"); u = ru.loc[t.index].join(t)
    u["y"] = np.log10(u.total); return u if units is None else u.loc[u.index.isin(units)]
prim = ru.index[ru.primary & ru.has_ARG]
alt13 = ru.index[(ru.Project_Accession == "PRJEB13831") & ru.has_ARG]
p27 = ru.index[(ru.Project_Accession == "PRJEB27054")]
specs = {"Primary (PRJEB27054 pilot; all hits; FPKM per bacterial fragments)": unit_summary(rfl, units=prim),
         "PRJEB13831 replaces PRJEB27054": unit_summary(rfl, units=prim.difference(p27).union(alt13)),
         "Both pilot projects included (double counts pilot samples)": unit_summary(rfl, units=prim.union(alt13)),
         "Exclude studies run with KMA 1.2.21": unit_summary(rfl, units=prim[~ru.loc[prim, "kma_versions"].str.contains("1.2.21")]),
         "FPKM per million trimmed fragments": unit_summary(rfl, value="FPKM_trimmed", units=prim),
         "Templates with run coverage >= 50% only": unit_summary(rfl[rfl.max_run_cov >= 50], units=prim)}
rows = []
for lab, u in specs.items():
    urb = u[u.Stratum == "Urban environmental"]
    f = smf.mixedlm("y ~ C(Sample_Type, Treatment('Municipal wastewater'))", urb, groups=urb.Project_Accession).fit(reml=True)
    for t in f.fe_params.index:
        if "Sample_Type" not in t: continue
        rows.append(dict(specification=lab, contrast=f"{term_label(t)} vs municipal wastewater", n_units=len(urb), n_studies=urb.Project_Accession.nunique(),
                         fold_change=10 ** f.fe_params[t], ci_low=10 ** f.conf_int().loc[t, 0], ci_high=10 ** f.conf_int().loc[t, 1], p_value=f.pvalues[t]))
    fm = u[u.Sample_Type.isin(["Municipal wastewater", "Food-related environment"])]
    g = smf.mixedlm("y ~ C(Stratum, Treatment('Urban environmental'))", fm, groups=fm.Project_Accession).fit(reml=True)
    t = [x for x in g.fe_params.index if "Stratum" in x][0]
    rows.append(dict(specification=lab, contrast="Food-related (contextual) vs municipal wastewater", n_units=len(fm), n_studies=fm.Project_Accession.nunique(),
                     fold_change=10 ** g.fe_params[t], ci_low=10 ** g.conf_int().loc[t, 0], ci_high=10 ** g.conf_int().loc[t, 1], p_value=g.pvalues[t]))
sens2 = pd.DataFrame(rows)
stl = specs[list(specs)[0]]; stl = stl[stl.Stratum == "Urban environmental"].groupby(["Project_Accession", "Sample_Type"], as_index=False).y.median()
k = kruskal_eps2(stl, "y", "Sample_Type")
record_test(NB, "Sensitivity ResFinder: study-level median log10 FPKM by matrix", "Kruskal-Wallis", k["n"], k["H"], k["df"], k["p"], "epsilon-squared", k["epsilon2"], family="Sensitivity", note="unit = study x matrix", objective="Objective 2")
save_table(sens2.pipe(rnd, 5), "SA03_rf_matrix_contrast_sensitivity", "Sensitivity of ResFinder matrix and food-related contrasts (mixed models with study random intercept)", NB, sheet="Sensitivity_Analysis")
save_table(pd.DataFrame([dict(specification="Study-level medians (unit = study x matrix)", **k)]).pipe(rnd, 6), "SA04_rf_study_level_kruskal", "Sensitivity: matrix differences using study-level medians", NB, sheet="Sensitivity_Analysis")
sens2.round(3)

# %%
fig, axes = plt.subplots(1, 2, figsize=(6.9, 5.6))
for ax, df, ttl in [(axes[0], sens1, "A  DS02 city contrasts (CARD/RGI)"), (axes[1], sens2, "B  ResFinder matrix contrasts")]:
    specs_ = df.specification.unique().tolist(); con = df.contrast.unique().tolist(); y0 = 0; ticks, labels = [], []
    pal = dict(zip(specs_, sns.color_palette("tab10", len(specs_))))
    for c in con:
        d = df[df.contrast == c]
        for j, (_, r) in enumerate(d.iterrows()):
            yy = y0 + j * .12
            ax.errorbar(r.fold_change, yy, xerr=[[r.fold_change - r.ci_low], [r.ci_high - r.fold_change]], fmt="o", ms=3, lw=.8, color=pal[r.specification],
                        label=r.specification if c == con[0] else None)
        ticks.append(y0 + .12 * (len(d) - 1) / 2); labels.append(c); y0 += .12 * len(d) + .35
    ax.axvline(1, color="k", lw=.6, ls="--"); ax.set_xscale("log"); ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=6.8); ax.invert_yaxis()
    ax.set_xlabel("Fold change in total normalized ARG abundance (95% CI, log scale)"); ax.set_title(ttl, fontsize=7.5)
    ax.legend(frameon=False, fontsize=6.0, loc="upper center", bbox_to_anchor=(.5, -.1), ncol=1)
fig.tight_layout()
save_fig(fig, 29, "sensitivity_forest", "Robustness of key normalized-ARG-abundance contrasts to analytic choices",
         "Fold changes (95% CI) under alternative specifications. (A) DS02 city contrasts versus Delhi from OLS adjusted for season and year (HC3), varying detection threshold, denominator, depth-outlier exclusion and aggregation to city-months. (B) ResFinder matrix contrasts versus municipal wastewater from mixed models with a study random intercept, varying the pilot project used, KMA version, denominator and coverage filter; the food-related contrast uses contextual data.",
         NB, "Objectives 2 / 5 (robustness)", "tables/SA01, SA03")

# %% [markdown]
# ## 3. Other robustness checks

# %%
from skbio import DistanceMatrix
from skbio.stats.distance import permanova as skperm
oth = []
alpha = pd.read_csv(DER / "ds02_alpha_diversity.csv", index_col=0)
for v in ["observed_richness", "rarefied_richness"]:
    k = kruskal_eps2(alpha, v, "City"); oth.append(dict(check=f"DS02 phylum {v} by city (Kruskal-Wallis)", n=k["n"], statistic=k["H"], p_value=k["p"], effect=k["epsilon2"], effect_name="epsilon-squared"))
fam = pd.read_csv(DER / "ds02_card_family_rpkm.csv", index_col=0).loc[m2.index]
Pz, _ = mult_replace(fam.loc[:, (fam > 0).mean() >= .5].values)
D_ait = squareform(pdist(clr(Pz), "euclidean"))
r = permanova_dbrda(D_ait, m2[["City"]].reset_index(), "City", n_perm=999)
oth.append(dict(check="DS02 resistome Aitchison distance (families in >= 50% of samples) ~ City (PERMANOVA)", n=r["n"], statistic=r["pseudo_F"], p_value=r["p_value"], effect=r["R2_partial"], effect_name="R2"))
us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0); cfp = pd.read_csv(DER / "rf_class_fpkm.csv", index_col=0)
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
from sklearn.metrics import roc_auc_score
urb = us[us.Stratum == "Urban environmental"]; Xc = cfp.loc[urb.index].drop(columns=[c for c in ["Unclassified"] if c in cfp.columns])
Pz, _ = mult_replace(Xc.values); X = clr(Pz); y = (urb.Sample_Type == "Municipal wastewater").astype(int).values
aucs = []
for sd in [1, 2, 3, 4, 5]:
    p = cross_val_predict(RandomForestClassifier(n_estimators=500, min_samples_leaf=2, class_weight="balanced", random_state=sd, n_jobs=-1), X, y, groups=urb.Project_Accession.values, cv=LeaveOneGroupOut(), method="predict_proba")[:, 1]
    aucs.append(roc_auc_score(y, p))
oth.append(dict(check="Model A random forest ROC-AUC across 5 random seeds (LOGO)", n=len(y), statistic=np.mean(aucs), p_value=np.nan, effect=np.std(aucs), effect_name=f"SD across seeds; range {min(aucs):.3f}-{max(aucs):.3f}"))
oth = pd.DataFrame(oth)
save_table(oth.pipe(rnd, 5), "SA05_other_robustness_checks", "Other robustness checks: richness rarefaction, alternative distance, model seeds", NB, sheet="Sensitivity_Analysis")
oth

# %%
pr = sens1[sens1.specification.str.startswith("Primary:")]
same_dir = sens1.groupby("contrast").apply(lambda d: bool((np.sign(np.log(d.fold_change)) == np.sign(np.log(d.fold_change.iloc[0]))).all()), include_groups=False)
pr2 = sens2[sens2.specification.str.startswith("Primary")]
same_dir2 = sens2.groupby("contrast").apply(lambda d: bool((np.sign(np.log(d.fold_change)) == np.sign(np.log(d.fold_change.iloc[0]))).all()), include_groups=False)
record_finding(NB, "F19", "Robustness (Objectives 2/5)", f"DS02 city contrasts kept the same direction in {int(same_dir.sum())} of {len(same_dir)} contrasts across {sens1.specification.nunique()} specifications (Kruskal-Wallis epsilon-squared range {kw.epsilon2.min():.2f}-{kw.epsilon2.max():.2f}); ResFinder contrasts kept their direction in {int(same_dir2.sum())} of {len(same_dir2)} contrasts across {sens2.specification.nunique()} specifications.",
               len(m2), "DS02 samples; ResFinder units", "OLS (HC3) and mixed models under alternative specifications", "Tables SA01-SA05; Figure 29",
               ("Directions of the main contrasts are not artefacts of a single threshold, denominator or duplicate-project choice" if (same_dir.all() and same_dir2.all()) else "Some contrasts change direction under alternative specifications (" + ", ".join(list(same_dir.index[~same_dir]) + list(same_dir2.index[~same_dir2])) + "); these should be reported as not robust")
               + "; CI widths change, so magnitudes should be quoted with their specification.",
               "Sensitivity analyses cannot remove confounding by study, site or sequencing depth that is shared by all specifications.")
