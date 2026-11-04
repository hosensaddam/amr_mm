# %% [markdown]
# # 05 · Comparative and Multivariate Analysis
# **Research objectives:** Objective 1 (community structure across cities/matrices), Objective 2 (resistome structure and its association with microbial community characteristics), Objective 3 (environmental/seasonal factors).
#
# **Design principles.** Comparisons are made *within* a pipeline and, wherever possible, within a study. Distance-based PERMANOVA (McArdle & Anderson 2001) is implemented with sequential adjustment for covariates and restricted permutations: within strata for nested terms, and at the level of whole study blocks when the tested factor is constant within studies (sample type). Multivariate dispersion (PERMDISP) is tested alongside each PERMANOVA so that dispersion differences are not read as location differences. Linear mixed models treat study as a random effect. Associations are not interpreted causally.

# %%
#BOOTSTRAP
NB = "NB05"
nb_header(NB)
import seaborn as sns, networkx as nx
import statsmodels.formula.api as smf
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import procrustes
from skbio import DistanceMatrix
from skbio.stats.distance import permanova as sk_permanova, mantel
from skbio.stats.ordination import pcoa
m2 = pd.read_csv(DER / "ds02_meta.csv").set_index("Sample_ID")
ph = pd.read_csv(DER / "ds02_prokaryote_phylum_counts.csv", index_col=0).loc[m2.index]
fam = pd.read_csv(DER / "ds02_card_family_rpkm.csv", index_col=0).loc[m2.index]
fam50 = pd.read_csv(DER / "ds02_card_family_rpkm_cov50.csv", index_col=0).loc[m2.index]
cls = pd.read_csv(DER / "ds02_card_class_rpkm.csv", index_col=0).loc[m2.index]
alpha = pd.read_csv(DER / "ds02_alpha_diversity.csv", index_col=0).loc[m2.index]
s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0).loc[m2.index]
meta2 = m2[["City", "Season", "year", "Site"]].reset_index()
N_PERM = 999

# %% [markdown]
# ## A1. Distance matrices (DS02)

# %%
ra = closure(ph.values)
D_tax_bc = squareform(pdist(ra, "braycurtis"))
clr_tax = clr(closure(ph.values + 0.5))
D_tax_ait = squareform(pdist(clr_tax, "euclidean"))
fam_rel = closure(fam.values)
D_arg_bc = squareform(pdist(fam_rel, "braycurtis"))
D_arg_jac = squareform(pdist((fam50.values > 0), "jaccard"))
np.save(DER / "ds02_D_tax_bc.npy", D_tax_bc); np.save(DER / "ds02_D_tax_aitchison.npy", D_tax_ait); np.save(DER / "ds02_D_arg_bc.npy", D_arg_bc)
print({k: (round(float(v[np.triu_indices_from(v, 1)].mean()), 3)) for k, v in {"tax_BC": D_tax_bc, "tax_Aitchison": D_tax_ait, "ARG_BC": D_arg_bc, "ARG_Jaccard": D_arg_jac}.items()})

# %% [markdown]
# ## A2. PERMANOVA and PERMDISP (DS02)

# %%
res = []
specs = [("Prokaryotic phyla", "Bray-Curtis", D_tax_bc), ("Prokaryotic phyla", "Aitchison (CLR)", D_tax_ait),
         ("ARG gene families", "Bray-Curtis (relative RPKM)", D_arg_bc), ("ARG gene families", "Jaccard (presence, cov>=50%)", D_arg_jac)]
for feat, dist, D in specs:
    # UPDATE (AUD-11): city is constant within site, so city labels are permuted among whole sites (block = Site);
    # season and year vary within site and are permuted within site. A free-permutation city row is kept only as a descriptive/cross-check value.
    rf_ = permanova_dbrda(D, meta2, "City", (), n_perm=N_PERM); rf_.update(features=feat, distance=dist, covariates="(none; free permutation, descriptive)"); res.append(rf_)
    for term, cov, strata, block in [("City", (), None, "Site"), ("Season", ("City",), "Site", None), ("year", ("City", "Season"), "Site", None)]:
        r = permanova_dbrda(D, meta2, term, cov, strata=strata, block=block, n_perm=N_PERM)
        r.update(features=feat, distance=dist); res.append(r)
        record_test(NB, f"DS02 {feat} ({dist}) ~ {term}" + (f" | {'+'.join(cov)}" if cov else ""), "PERMANOVA (dbRDA pseudo-F)", r["n"], r["pseudo_F"], f"{r['df_term']}, {r['df_resid']}", r["p_value"], "partial R2", r["R2_partial"], family="DS02 PERMANOVA (site-aware, primary)", note=r["permutation_scheme"], objective="Objective 1/2/3")
    for fac in ["City", "Season"]:
        d = permdisp_test(D, meta2[fac].astype(str).values, N_PERM)
        res.append(dict(term=fac, covariates="(PERMDISP)", pseudo_F=d["F"], df_term=d["groups"] - 1, df_resid=d["n"] - d["groups"], R2_partial=np.nan, p_value=d["p_value"], n=d["n"], permutations=N_PERM, permutation_scheme="PERMDISP (distance to centroid)", features=feat, distance=dist))
        record_test(NB, f"DS02 {feat} ({dist}) dispersion by {fac}", "PERMDISP", d["n"], d["F"], f"{d['groups']-1}, {d['n']-d['groups']}", d["p_value"], family="DS02 PERMDISP [sample-level, descriptive only]", objective="Objective 1/2")
perm = pd.DataFrame(res)[["features", "distance", "term", "covariates", "n", "df_term", "df_resid", "pseudo_F", "R2_partial", "p_value", "permutations", "permutation_scheme"]]
save_table(perm.pipe(rnd, 5), "C01_ds02_permanova_permdisp", "DS02 PERMANOVA (sequential): city tested by permuting whole sites (19), season/year permuted within site (AUD-11, primary); free-permutation city rows and PERMDISP are descriptive", NB, sheet="Beta_Diversity")
ids = [str(i) for i in range(len(meta2))]
skr = sk_permanova(DistanceMatrix(D_tax_bc, ids), list(meta2.City.astype(str).values), permutations=N_PERM, seed=SEED)
own = perm[(perm.distance == "Bray-Curtis") & (perm.term == "City") & (perm.covariates == "(none; free permutation, descriptive)")].iloc[0]
xval = pd.DataFrame([dict(check="PERMANOVA implementation vs scikit-bio (prokaryotic BC ~ City, free permutation)", own_pseudo_F=own.pseudo_F, skbio_pseudo_F=float(skr["test statistic"]),
                          abs_difference=abs(own.pseudo_F - float(skr["test statistic"])), own_p=own.p_value, skbio_p=float(skr["p-value"]))])
save_table(xval, "C02_permanova_crosscheck", "Cross-check of the custom PERMANOVA against scikit-bio", NB, sheet="Verification_Report")
print(xval.T); perm.round(4)

# %% [markdown]
# ## A3. Concordance between community structure and resistome (DS02)

# %%
mt = mantel(DistanceMatrix(D_tax_ait, ids), DistanceMatrix(D_arg_bc, ids), method="spearman", permutations=N_PERM, seed=SEED)
mt2 = mantel(DistanceMatrix(D_tax_bc, ids), DistanceMatrix(D_arg_bc, ids), method="spearman", permutations=N_PERM, seed=SEED)
pa = pcoa(DistanceMatrix(D_tax_ait, ids), number_of_dimensions=5).samples.values; pb = pcoa(DistanceMatrix(D_arg_bc, ids), number_of_dimensions=5).samples.values
_, _, m2_obs = procrustes(pa, pb)
rng = np.random.default_rng(SEED); null = [procrustes(pa, pb[rng.permutation(len(pb))])[2] for _ in range(N_PERM)]
p_pro = (np.sum(np.array(null) <= m2_obs) + 1) / (N_PERM + 1)
conc = pd.DataFrame([dict(test="Mantel (Spearman): taxonomy Aitchison vs ARG Bray-Curtis", statistic=mt[0], p_value=mt[1], n=mt[2], permutations=N_PERM),
                     dict(test="Mantel (Spearman): taxonomy Bray-Curtis vs ARG Bray-Curtis", statistic=mt2[0], p_value=mt2[1], n=mt2[2], permutations=N_PERM),
                     dict(test="Procrustes m2 (PCoA 5 axes; lower = more concordant)", statistic=m2_obs, p_value=p_pro, n=len(pa), permutations=N_PERM)])
conc["procrustes_correlation_r"] = [np.nan, np.nan, np.sqrt(1 - m2_obs)]
for _, r in conc.iterrows():
    record_test(NB, "DS02 taxonomy-resistome concordance", r.test, int(r.n), r.statistic, None, r.p_value, "Mantel r / Procrustes m2", r.statistic, family="DS02 concordance", objective="Objective 2")
save_table(conc.pipe(rnd, 6), "C03_ds02_taxonomy_resistome_concordance", "DS02: concordance of prokaryotic community and resistome structure (Mantel, Procrustes)", NB, sheet="Correlation_Analysis")
conc

# %%
topph = ph.sum().sort_values(ascending=False).index[:15]
clr_df = pd.DataFrame(clr_tax, index=ph.index, columns=ph.columns)[topph]
cls_keep = cls.sum().sort_values(ascending=False).index[:10]
lcls = np.log10(cls[cls_keep].replace(0, np.nan)).fillna(np.log10(cls[cls_keep][cls[cls_keep] > 0].min().min() / 2))
rows = []
for p_ in topph:
    for c_ in cls_keep:
        r = spearman_ci(clr_df[p_], lcls[c_]); rows.append(dict(phylum=p_, class_group=c_, **r))
cor = pd.DataFrame(rows); cor["p_BH"] = bh(cor.p)
save_table(cor.pipe(rnd, 6), "C04_ds02_phylum_class_spearman", "DS02: Spearman correlations between prokaryotic phylum CLR values and log10 class-group RPKM (BH-adjusted, 150 tests)", NB, sheet="Correlation_Analysis")
for _, r in cor.iterrows():
    record_test(NB, f"DS02 {r.phylum} CLR vs {r.class_group} log10 RPKM", "Spearman", int(r.n), r.rho, None, r.p, "rho", r.rho, r.ci_low, r.ci_high, family="DS02 phylum-class correlations (BH within family)", objective="Objective 2")
R = cor.pivot(index="phylum", columns="class_group", values="rho").loc[topph, cls_keep]; Q = cor.pivot(index="phylum", columns="class_group", values="p_BH").loc[topph, cls_keep]
fig, ax = plt.subplots(figsize=(6.9, 5.4))
ann = R.round(2).astype(str) + np.where(Q < .001, "***", np.where(Q < .01, "**", np.where(Q < .05, "*", "")))
sns.heatmap(R, cmap="RdBu_r", center=0, vmin=-.8, vmax=.8, annot=ann, fmt="", annot_kws={"fontsize": 6}, linewidths=.3, cbar_kws={"label": "Spearman rho"}, ax=ax)
ax.set_xlabel("Resistance-class group (log10 RPKM)"); ax.set_ylabel("Prokaryotic phylum (CLR)"); ax.set_title("Phylum-resistome associations in DS02 wastewater (n = 447)")
plt.setp(ax.get_xticklabels(), rotation=40, ha="right")
save_fig(fig, 15, "ds02_phylum_class_correlation", "Correlation between prokaryotic phylum composition and resistance-class abundance (DS02)",
         "Spearman correlations between centred log-ratio (CLR) abundances of the 15 most abundant prokaryotic phyla and log10 RPKM of the 10 most abundant resistance-class groups. * q < 0.05, ** q < 0.01, *** q < 0.001 (Benjamini-Hochberg over 150 tests). Correlations are compositional associations across samples and do not identify ARG hosts.",
         NB, "Objective 2", "tables/C04_ds02_phylum_class_spearman.csv")

# %%
dfm = s2.join(alpha[["shannon"]]).assign(clr_proteo=clr_df["Proteobacteria"], City=m2.City.astype(str), Season=m2.Season.astype(str))
ols = smf.ols("log10_total_RPKM ~ shannon + clr_proteo + C(City, Treatment('Delhi')) + C(Season, Treatment('Winter (Dec-Feb)'))", dfm).fit(cov_type="HC3")
o = pd.DataFrame({"estimate": ols.params, "SE_HC3": ols.bse, "ci_low": ols.conf_int()[0], "ci_high": ols.conf_int()[1], "p_value": ols.pvalues}).reset_index().rename(columns={"index": "term"})
save_table(o.pipe(rnd, 5), "C05_ds02_arg_load_vs_community_ols", f"DS02: log10 total ARG RPKM ~ phylum Shannon + Proteobacteria CLR + city + season (OLS, HC3; R2 = {ols.rsquared:.3f})", NB, sheet="Correlation_Analysis")
for t in ["shannon", "clr_proteo"]:
    record_test(NB, "DS02 log10 total ARG RPKM ~ community covariates + city + season", f"OLS coefficient: {t}", int(ols.nobs), ols.tvalues[t], int(ols.df_resid), ols.pvalues[t], "coefficient (log10 RPKM per unit)", ols.params[t], ols.conf_int().loc[t, 0], ols.conf_int().loc[t, 1], family="DS02 ARG-community models", objective="Objective 2")
r_sh = spearman_ci(dfm.shannon, dfm.log10_total_RPKM); r_pr = spearman_ci(dfm.clr_proteo, dfm.log10_total_RPKM)
for nm, r in [("phylum Shannon", r_sh), ("Proteobacteria CLR", r_pr)]:
    record_test(NB, f"DS02 {nm} vs log10 total ARG RPKM (unadjusted)", "Spearman", r["n"], r["rho"], None, r["p"], "rho", r["rho"], r["ci_low"], r["ci_high"], family="DS02 ARG-community models", objective="Objective 2")
print(r_sh, r_pr); o.round(4)

# %% [markdown]
# ### A4. ARG co-occurrence network (DS02, exploratory)
# Families detected (coverage ≥ 50 %) in ≥ 50 % of samples; Spearman correlation of log10 RPKM; edges kept for ρ ≥ 0.75 and BH q < 0.001. Co-occurrence can reflect shared hosts, shared mobile elements, or shared sources of faecal input; it is not evidence of physical linkage.

# %%
keep = fam50.columns[(fam50 > 0).mean() >= .5]
keep = fam[keep].mean().sort_values(ascending=False).index[:60]
X = np.log10(fam[keep].where(fam[keep] > 0)).fillna(np.log10(fam[keep][fam[keep] > 0].min().min() / 2))
rho, pv = stats.spearmanr(X.values)
iu = np.triu_indices(len(keep), 1)
edges = pd.DataFrame({"a": keep[iu[0]], "b": keep[iu[1]], "rho": rho[iu], "p": pv[iu]}); edges["q"] = bh(edges.p)
E = edges[(edges.rho >= .75) & (edges.q < .001)].copy()
grp = pd.read_csv(DER / "ds02_card_allele_long.csv", usecols=["AMR_Gene_Family", "Class_Group"]).drop_duplicates("AMR_Gene_Family").set_index("AMR_Gene_Family").Class_Group
G = nx.Graph(); G.add_nodes_from(keep); G.add_weighted_edges_from(E[["a", "b", "rho"]].values)
deg = pd.DataFrame({"AMR_Gene_Family": list(G.nodes), "degree": [G.degree(n) for n in G.nodes], "Class_Group": [grp.get(n) for n in G.nodes],
                    "mean_RPKM": [fam[n].mean() for n in G.nodes]}).sort_values("degree", ascending=False)
comps = sorted(nx.connected_components(G), key=len, reverse=True)
save_table(E.pipe(rnd, 6), "C06_ds02_arg_network_edges", f"DS02 ARG co-occurrence network edges (rho >= 0.75, q < 0.001); {len(E)} edges among {len(keep)} families", NB, sheet="Correlation_Analysis")
save_table(deg.pipe(rnd, 4), "C07_ds02_arg_network_nodes", "DS02 ARG co-occurrence network: node degree and class", NB, sheet="Correlation_Analysis")
fig, ax = plt.subplots(figsize=(6.9, 5.6))
H = G.subgraph([n for n in G.nodes if G.degree(n) > 0])
pos = nx.spring_layout(H, seed=SEED, k=.45)
nx.draw_networkx_edges(H, pos, ax=ax, width=[(d["weight"] - .7) * 6 for _, _, d in H.edges(data=True)], alpha=.35, edge_color="grey")
nx.draw_networkx_nodes(H, pos, ax=ax, node_color=[CLASS_COLORS.get(grp.get(n), "#999") for n in H.nodes], node_size=[30 + 25 * np.sqrt(fam[n].mean()) for n in H.nodes], edgecolors="k", linewidths=.3)
top_lab = deg[deg.degree > 0].head(30).AMR_Gene_Family
_ids = {n: i + 1 for i, n in enumerate(sorted(H.nodes, key=lambda n: -fam[n].mean()))}  # AUD-15: numbered nodes + full-name key (no truncation/overlap)
nx.draw_networkx_labels(H, pos, labels={n: str(_ids[n]) for n in H.nodes}, font_size=5.8, font_weight="bold", ax=ax)
import textwrap as _tw
_keys = [f"{_ids[n]} {n}" for n in sorted(H.nodes, key=lambda n: _ids[n])]; _h = (len(_keys) + 1) // 2
ax.text(0.0, -0.01, "\n".join(_keys[:_h]), transform=ax.transAxes, fontsize=5.5, va="top", ha="left"); ax.text(0.5, -0.01, "\n".join(_keys[_h:]), transform=ax.transAxes, fontsize=5.5, va="top", ha="left")
for c in [c for c in CLASS_ORDER if c in set(grp.get(n) for n in H.nodes)]:
    ax.scatter([], [], color=CLASS_COLORS[c], label=c, s=25)
ax.legend(frameon=False, fontsize=6.2, loc="lower left"); ax.axis("off")
ax.set_title(f"ARG family co-occurrence network, DS02 ({H.number_of_nodes()} connected families, {H.number_of_edges()} edges; largest component {len(comps[0])})", fontsize=7.5)
save_fig(fig, 16, "ds02_arg_cooccurrence_network", "Co-occurrence network of ARG families in Indian open-drain wastewater (DS02)",
         f"Nodes: CARD gene families detected (coverage ≥ 50%) in ≥ 50% of the 447 samples (top 60 by mean RPKM; isolated nodes hidden); node size ∝ √(mean RPKM); colour = class group. Edges: Spearman ρ ≥ 0.75 between log10 RPKM with Benjamini-Hochberg q < 0.001. Exploratory; co-occurrence does not demonstrate co-location on mobile elements.",
         NB, "Objective 2", "tables/C06_ds02_arg_network_edges.csv; tables/C07_ds02_arg_network_nodes.csv")

# %% [markdown]
# ## B. ResFinder/KMA resource: matrix effects, study heterogeneity and diversity associations

# %%
us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0)
gfp = pd.read_csv(DER / "rf_gene_fpkm.csv", index_col=0).loc[us.index]
urb = us[us.Stratum == "Urban environmental"].copy()
REF = "Municipal wastewater"
def mixed(formula, data, label):
    mdl = smf.mixedlm(formula, data, groups=data.Project_Accession).fit(reml=True)
    ml = smf.mixedlm(formula, data, groups=data.Project_Accession).fit(reml=False)
    null_f = formula.split("~")[0] + "~ 1" + "".join(" + " + t for t in formula.split("~")[1].split("+") if "Sample_Type" not in t and "Stratum" not in t and t.strip() != "1")
    ml0 = smf.mixedlm(null_f, data, groups=data.Project_Accession).fit(reml=False)
    lr = 2 * (ml.llf - ml0.llf); dfl = len(ml.fe_params) - len(ml0.fe_params); p = stats.chi2.sf(lr, dfl)
    t = pd.DataFrame({"estimate": mdl.fe_params, "SE": mdl.bse_fe, "ci_low": mdl.conf_int().loc[mdl.fe_params.index, 0], "ci_high": mdl.conf_int().loc[mdl.fe_params.index, 1], "p_value": mdl.pvalues[mdl.fe_params.index]})
    vr = float(mdl.cov_re.iloc[0, 0]); icc = vr / (vr + mdl.scale)
    extra = pd.DataFrame({"estimate": [vr, mdl.scale, icc, lr, p]}, index=["Study variance", "Residual variance", "ICC (study)", "LRT chi2 (matrix term, ML)", "LRT p-value"])
    t = pd.concat([t, extra]); t.index.name = "term"
    record_test(NB, label, "Linear mixed model, LRT for matrix term", int(mdl.nobs), lr, dfl, p, "ICC (study)", icc, family="ResFinder mixed models", note=formula, objective="Objective 2/5")
    return t.reset_index(), mdl
t1, mdl1 = mixed(f"log10_total_FPKM ~ C(Sample_Type, Treatment('{REF}'))", urb, "ResFinder urban: log10 total FPKM ~ sample type + (1|study)")
t1["fold_change_vs_MW"] = np.where(t1.term.str.contains("Sample_Type"), 10 ** t1.estimate, np.nan)
t2, _ = mixed(f"gene_richness ~ C(Sample_Type, Treatment('{REF}')) + log10_bacterial_fragments", urb, "ResFinder urban: gene richness ~ sample type + log10 bacterial fragments + (1|study)")
fm = us[us.Sample_Type.isin([REF, "Food-related environment"])].copy()
t3, _ = mixed("log10_total_FPKM ~ C(Stratum, Treatment('Urban environmental'))", fm, "ResFinder: food-related (contextual) vs municipal wastewater, log10 total FPKM + (1|study)")
t3["fold_change"] = np.where(t3.term.str.contains("Stratum"), 10 ** t3.estimate, np.nan)
# UPDATE (AUD-13): identifiability of matrix effects given the study structure
sm_struct = urb.groupby("Sample_Type").agg(units=("Project_Accession", "size"), studies=("Project_Accession", "nunique"),
                                           study_list=("Project_Accession", lambda s: "; ".join(f"{k} ({v})" for k, v in s.value_counts().items()))).reset_index()
sm_struct["largest_study_share_pct"] = [100 * urb[urb.Sample_Type == t].Project_Accession.value_counts().iloc[0] / n for t, n in zip(sm_struct.Sample_Type, sm_struct.units)]
sm_struct["matrix_effect_identifiable_from_study_effect"] = np.where(sm_struct.studies >= 3, "Partly (>= 3 studies; still study-limited)", np.where(sm_struct.studies == 2, "Weakly (2 studies)", "No - single-study matrix: matrix and study effects are not separable"))
t1["single_study_matrix"] = [("Yes" if sm_struct.set_index("Sample_Type").studies.get(term_label(x), 99) == 1 else "No") if "Sample_Type" in str(x) else "" for x in t1.term]
t1["studies_for_matrix"] = [sm_struct.set_index("Sample_Type").studies.get(term_label(x), np.nan) if "Sample_Type" in str(x) else np.nan for x in t1.term]
save_table(sm_struct.pipe(rnd, 3), "C08b_rf_matrix_study_structure", "AUD-13: study structure of each urban matrix in the ResFinder resource and identifiability of matrix effects", NB, sheet="Comparative_Analysis")
# Within-study matrix contrasts: only studies containing >= 2 matrices with >= 3 units each
ws_rows = []
for prj, d in urb.groupby("Project_Accession"):
    cnt = d.Sample_Type.value_counts(); ok = cnt[cnt >= 3].index.tolist()
    if len(ok) < 2: continue
    ref = REF if REF in ok else ok[0]
    for t_ in [x for x in ok if x != ref]:
        a_ = d.loc[d.Sample_Type == t_, "log10_total_FPKM"]; b_ = d.loc[d.Sample_Type == ref, "log10_total_FPKM"]
        u_, p_ = stats.mannwhitneyu(a_, b_, alternative="two-sided")
        ws_rows.append(dict(study=prj, matrix=t_, reference=ref, n_matrix=len(a_), n_reference=len(b_), median_log10_FPKM_matrix=a_.median(), median_log10_FPKM_reference=b_.median(),
                            fold_change_of_medians=10 ** (a_.median() - b_.median()), mann_whitney_U=u_, p_value=p_, cliffs_delta=cliffs_delta_from_u(u_, len(a_), len(b_))))
ws = pd.DataFrame(ws_rows, columns=["study", "matrix", "reference", "n_matrix", "n_reference", "median_log10_FPKM_matrix", "median_log10_FPKM_reference", "fold_change_of_medians", "mann_whitney_U", "p_value", "cliffs_delta"])
if len(ws): ws["p_holm"] = holm(ws.p_value)
save_table(ws.pipe(rnd, 5), "C08c_rf_within_study_matrix_contrasts", "AUD-13: within-study matrix contrasts of log10 total normalized ARG abundance (FPKM), studies with >= 2 matrices of >= 3 units (no between-study confounding)", NB, sheet="Comparative_Analysis",
           description="Empty if no study sampled two or more matrices with at least three units each.")
print(sm_struct); print(ws)
save_table(t1.pipe(rnd, 5), "C08_rf_mixed_arg_load", f"Mixed model: log10 total ARG FPKM by urban matrix (reference municipal wastewater), study random intercept; n = {len(urb)} units, {urb.Project_Accession.nunique()} studies", NB, sheet="Comparative_Analysis")
save_table(t2.pipe(rnd, 5), "C09_rf_mixed_gene_richness", "Mixed model: ARG gene richness by urban matrix adjusted for log10 bacterial fragments, study random intercept", NB, sheet="Comparative_Analysis")
save_table(t3.pipe(rnd, 5), "C10_rf_mixed_food_vs_wastewater", f"Mixed model: food-related environments (contextual) vs municipal wastewater, log10 total FPKM; n = {len(fm)}", NB, sheet="Comparative_Analysis")
t1.round(3)

# %%
st = us.groupby(["Project_Accession", "Dataset_ID", "Sample_Type", "Stratum"]).log10_total_FPKM.agg(["mean", "std", "count"]).reset_index()
st["se"] = st["std"] / np.sqrt(st["count"]); tcrit = stats.t.ppf(.975, (st["count"] - 1).clip(lower=1))
st["ci_low"] = st["mean"] - tcrit * st.se; st["ci_high"] = st["mean"] + tcrit * st.se
def dersimonian(m, se):
    w = 1 / se ** 2; mf = np.sum(w * m) / w.sum(); Qs = np.sum(w * (m - mf) ** 2); k = len(m)
    tau2 = max(0, (Qs - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())); ws = 1 / (se ** 2 + tau2)
    mr = np.sum(ws * m) / ws.sum(); ser = np.sqrt(1 / ws.sum()); I2 = max(0, (Qs - (k - 1)) / Qs) if Qs > 0 else 0
    return dict(k=k, pooled=mr, ci_low=mr - 1.96 * ser, ci_high=mr + 1.96 * ser, tau2=tau2, Q=Qs, Q_df=k - 1, Q_p=stats.chi2.sf(Qs, k - 1), I2=I2)
meta_rows = []
for t_ in st.Sample_Type.unique():
    d = st[(st.Sample_Type == t_) & (st["count"] >= 2)]
    if len(d) >= 2:
        r = dersimonian(d["mean"].values, d.se.values); r["Sample_Type"] = t_; meta_rows.append(r)
        record_test(NB, f"Between-study heterogeneity of mean log10 FPKM: {t_}", "Cochran Q (DerSimonian-Laird)", int(d["count"].sum()), r["Q"], r["Q_df"], r["Q_p"], "I-squared", r["I2"], family="ResFinder heterogeneity", note=f"k = {r['k']} studies; tau2 = {r['tau2']:.3f}", objective="Objective 2")
mh = pd.DataFrame(meta_rows)
save_table(st.pipe(rnd, 4), "C11_rf_study_level_means", "ResFinder: study-level mean log10 total FPKM by matrix (95% t-intervals)", NB, sheet="Comparative_Analysis")
save_table(mh.pipe(rnd, 5), "C12_rf_between_study_heterogeneity", "Random-effects (DerSimonian-Laird) pooled mean log10 FPKM and between-study heterogeneity by matrix (studies with >= 2 units)", NB, sheet="Comparative_Analysis")
names = pd.read_json(HARM / "dataset_names.json", typ="series")
fig, ax = plt.subplots(figsize=(6.0, 7.2))
sto = st.sort_values(["Stratum", "Sample_Type", "mean"], ascending=[False, True, True]).reset_index(drop=True)
for i, r in sto.iterrows():
    ax.errorbar(r["mean"], i, xerr=[[r["mean"] - r.ci_low], [r.ci_high - r["mean"]]] if r["count"] > 1 else None, fmt="s", ms=4, color=SAMPLE_TYPE_COLORS.get(r.Sample_Type), lw=.9, capsize=2)
ax.set_yticks(range(len(sto))); ax.set_yticklabels([f"{r.Dataset_ID} {r.Project_Accession} · {r.Sample_Type} (n={int(r['count'])})" for _, r in sto.iterrows()], fontsize=6.2)
for _, r in mh.iterrows():
    ys = sto.index[sto.Sample_Type == r.Sample_Type]
    ax.plot([r.pooled] * 2, [ys.min() - .4, ys.max() + .4], color=SAMPLE_TYPE_COLORS.get(r.Sample_Type), lw=1.2, ls="--")
ax.set_xlabel("Mean log10 total ARG FPKM (95% CI)"); ax.set_title("Study-level normalized ARG abundance by matrix (ResFinder/KMA)\ndashed line = random-effects pooled mean", fontsize=7.5)
save_fig(fig, 17, "rf_study_forest", "Between-study heterogeneity of normalized ARG abundance (forest plot, ResFinder/KMA)",
         "Mean log10 total ARG FPKM per study and matrix with 95% t-intervals (squares without whiskers: single unit). Dashed lines: DerSimonian-Laird random-effects pooled means for matrices with ≥ 2 studies; I² and Cochran's Q are in Table C12. The spread between studies within a matrix is comparable to the spread between matrices.",
         NB, "Objective 2", "tables/C11_rf_study_level_means.csv; tables/C12_rf_between_study_heterogeneity.csv")
mh.round(3)

# %%
rel = closure(gfp.values)
D_rf = squareform(pdist(rel, "braycurtis"))
np.save(DER / "rf_D_gene_bc.npy", D_rf)
mrf = us[["Sample_Type", "Stratum", "Project_Accession"]].reset_index()
mrf["block"] = mrf.Project_Accession + "|" + mrf.Sample_Type
iu_ = mrf.index[mrf.Stratum == "Urban environmental"].values
Du = D_rf[np.ix_(iu_, iu_)]; mu = mrf.loc[iu_].reset_index(drop=True)
rr = []
r = permanova_dbrda(Du, mu, "Sample_Type", (), block="block", n_perm=N_PERM); r.update(set="Urban units"); rr.append(r)
r = permanova_dbrda(Du, mu, "Project_Accession", ("Sample_Type",), strata="Sample_Type", n_perm=N_PERM); r.update(set="Urban units"); rr.append(r)
r = permanova_dbrda(D_rf, mrf, "Sample_Type", (), block="block", n_perm=N_PERM); r.update(set="Urban + food-related contextual"); rr.append(r)
d = permdisp_test(Du, mu.Sample_Type.values, N_PERM)
rr.append(dict(set="Urban units", term="Sample_Type", covariates="(PERMDISP)", pseudo_F=d["F"], df_term=d["groups"] - 1, df_resid=d["n"] - d["groups"], R2_partial=np.nan, p_value=d["p_value"], n=d["n"], permutations=N_PERM, permutation_scheme="PERMDISP"))
rp = pd.DataFrame(rr)[["set", "term", "covariates", "n", "df_term", "df_resid", "pseudo_F", "R2_partial", "p_value", "permutation_scheme"]]
for _, r in rp.iterrows():
    record_test(NB, f"ResFinder gene composition (Bray-Curtis) ~ {r.term} | {r.covariates} [{r.set}]", "PERMDISP" if "PERMDISP" in r.covariates else "PERMANOVA (dbRDA pseudo-F)", int(r.n), r.pseudo_F, f"{r.df_term}, {r.df_resid}", r.p_value, "partial R2", r.R2_partial, family="ResFinder PERMANOVA", note=r.permutation_scheme, objective="Objective 2")
save_table(rp.pipe(rnd, 5), "C13_rf_permanova", "ResFinder gene composition (Bray-Curtis on relative FPKM): PERMANOVA with study-block permutation, study-within-matrix test, and PERMDISP", NB, sheet="Beta_Diversity")
rp

# %%
dd = us.dropna(subset=["genus_shannon"]).copy()
rows = []
for grp_, d in [("All urban units", dd[dd.Stratum == "Urban environmental"]), ("Food-related (contextual)", dd[dd.Stratum != "Urban environmental"])] + [(t_, d) for t_, d in dd[dd.Stratum == "Urban environmental"].groupby("Sample_Type")]:
    for x, y in [("genus_shannon", "gene_shannon_fragments"), ("genus_richness", "gene_richness"), ("genus_shannon", "log10_total_FPKM")]:
        if len(d) < 8: continue
        r = spearman_ci(d[x], d[y]); rows.append(dict(group=grp_, x=x, y=y, **r))
dc = pd.DataFrame(rows); dc["p_BH"] = bh(dc.p)
for _, r in dc.iterrows():
    record_test(NB, f"ResFinder {r.x} vs {r.y} [{r.group}]", "Spearman", int(r.n), r.rho, None, r.p, "rho", r.rho, r.ci_low, r.ci_high, family="ResFinder diversity associations", objective="Objective 2")
save_table(dc.pipe(rnd, 6), "C14_rf_diversity_associations", "ResFinder: Spearman associations between bacterial (genus) diversity and ARG diversity/normalized abundance (BH-adjusted)", NB, sheet="Correlation_Analysis")
du = dd[dd.Stratum == "Urban environmental"]
mm = smf.mixedlm(f"gene_shannon_fragments ~ genus_shannon + C(Sample_Type, Treatment('{REF}')) + log10_bacterial_fragments", du, groups=du.Project_Accession).fit(reml=True)
mmt = pd.DataFrame({"estimate": mm.fe_params, "SE": mm.bse_fe, "ci_low": mm.conf_int().loc[mm.fe_params.index, 0], "ci_high": mm.conf_int().loc[mm.fe_params.index, 1], "p_value": mm.pvalues[mm.fe_params.index]}).reset_index().rename(columns={"index": "term"})
save_table(mmt.pipe(rnd, 5), "C15_rf_arg_diversity_mixed", "Mixed model: ARG gene Shannon ~ genus Shannon + matrix + log10 bacterial fragments + (1|study), urban units", NB, sheet="Correlation_Analysis")
record_test(NB, "ResFinder ARG Shannon ~ genus Shannon (+ matrix, depth, study RE)", "Mixed-model coefficient (Wald z)", int(mm.nobs), mm.tvalues["genus_shannon"], None, mm.pvalues["genus_shannon"], "coefficient", mm.params["genus_shannon"], mm.conf_int().loc["genus_shannon", 0], mm.conf_int().loc["genus_shannon", 1], family="ResFinder diversity associations", objective="Objective 2")
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.2))
for ax, (x, y, xl, yl) in zip(axes, [("genus_shannon", "gene_shannon_fragments", "Genus Shannon (SILVA rRNA)", "ARG gene Shannon"), ("genus_shannon", "log10_total_FPKM", "Genus Shannon (SILVA rRNA)", "log10 total ARG FPKM")]):
    for t_, d in dd.groupby("Sample_Type"):
        ax.scatter(d[x], d[y], s=9, alpha=.6, color=SAMPLE_TYPE_COLORS.get(t_), label=f"{t_} (n={len(d)})", edgecolors="none")
    rr_ = dc[(dc.group == "All urban units") & (dc.x == x) & (dc.y == y)].iloc[0]
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(f"Urban units: Spearman ρ = {rr_.rho:.2f} (95% CI {rr_.ci_low:.2f} to {rr_.ci_high:.2f})", fontsize=7.2)
axes[1].legend(frameon=False, fontsize=6.2, loc="lower left")
fig.tight_layout()
save_fig(fig, 18, "rf_diversity_vs_resistome", "Bacterial diversity versus ARG diversity and normalized ARG abundance across matrices (ResFinder/KMA)",
         "Each point is an analytic unit. Spearman coefficients (Bonett-Wright 95% CI) are for urban units pooled; within-matrix correlations and the study-adjusted mixed model are in Tables C14-C15. Pooled correlations can be driven by between-matrix and between-study differences (ecological fallacy).",
         NB, "Objective 2", "tables/C14_rf_diversity_associations.csv; tables/C15_rf_arg_diversity_mixed.csv")

# %%
pc = perm.set_index(["distance", "term", "covariates"])
c_bc = pc.loc[("Bray-Curtis", "City", "(none)")]; c_arg = pc.loc[("Bray-Curtis (relative RPKM)", "City", "(none)")]; s_arg = pc.loc[("Bray-Curtis (relative RPKM)", "Season", "City")]
record_finding(NB, "F10", "Objective 1/2/3", f"Within DS02, city explained {100*c_bc.R2_partial:.1f}% of variation in prokaryotic community composition (Bray-Curtis PERMANOVA pseudo-F = {c_bc.pseudo_F:.1f}, p = {c_bc.p_value:.3f}) and {100*c_arg.R2_partial:.1f}% of resistome composition (pseudo-F = {c_arg.pseudo_F:.1f}, p = {c_arg.p_value:.3f}); season, after adjusting for city, explained {100*s_arg.R2_partial:.1f}% of resistome variation (p = {s_arg.p_value:.3f}).",
               len(meta2), "DS02 samples", "PERMANOVA (999 permutations): city labels permuted among the 19 whole sites; season permuted within site (AUD-11); PERMDISP descriptive", "Table C01",
               f"City structures both profiles more than season does. PERMDISP by city: taxonomy p = {pc.loc[('Bray-Curtis','City','(PERMDISP)')].p_value:.3f}, resistome p = {pc.loc[('Bray-Curtis (relative RPKM)','City','(PERMDISP)')].p_value:.3f}; where PERMDISP is significant, part of the PERMANOVA signal reflects unequal within-city dispersion rather than centroid shifts.",
               "Inference on city rests on 19 sites; season and year are partly collinear within a 25-month series.")
record_finding(NB, "F11", "Objective 2", f"Community structure and resistome structure were concordant in DS02 (Mantel r = {mt[0]:.2f}, p = {mt[1]:.3f}; Procrustes correlation = {np.sqrt(1-m2_obs):.2f}, p = {p_pro:.3f}).",
               len(meta2), "DS02 samples", "Mantel (Spearman, 999 permutations); Procrustes on 5 PCoA axes", "Table C03; Figure 15",
               "Samples with more similar prokaryotic communities tend to carry more similar resistomes, consistent with (but not proof of) community composition shaping the resistome.",
               "Phylum resolution; both profiles derive from the same reads and share depth and extraction biases.")
fc = t1[t1.term.str.contains("Sample_Type")].copy()
lower = [term_label(x) for x, h in zip(fc.term, fc.ci_high) if 10 ** h < 1]
higher = [term_label(x) for x, l in zip(fc.term, fc.ci_low) if 10 ** l > 1]
ss_ = sm_struct.set_index("Sample_Type").studies
record_finding(NB, "F12", "Objective 2", "Study-adjusted (random-intercept) fold changes in total normalized ARG abundance (FPKM) relative to municipal wastewater were " + "; ".join(f"{term_label(x)} {10**e:.2f} (95% CI {10**l:.2f}-{10**h:.2f}; {int(ss_.get(term_label(x), 0))} {'study - single-study matrix, not separable from study effect' if ss_.get(term_label(x), 0) == 1 else 'studies'})" for x, e, l, h in zip(fc.term, fc.estimate, fc.ci_low, fc.ci_high)) + f"; the study ICC was {t1.set_index('term').loc['ICC (study)','estimate']:.2f}. Within-study contrasts were possible for {ws.study.nunique() if len(ws) else 0} stud{'y' if (ws.study.nunique() if len(ws) else 0) == 1 else 'ies'} (Table C08c).",
               len(urb), f"Urban units ({urb.Project_Accession.nunique()} studies)", f"Mixed model LRT p = {t1.set_index('term').loc['LRT p-value','estimate']:.2e}", "Tables C08, C11-C12; Figure 17",
               f"Single-study matrices ({', '.join(ss_[ss_ == 1].index) or 'none'}) cannot be interpreted as independent environmental-matrix effects; their fold changes describe one study and depend on the random-effect assumption. Estimates for the remaining matrices are study-limited and not generalizable (AUD-13).",
               "Matrix is confounded with study, KMA version and sequencing depth; few studies per matrix.")
