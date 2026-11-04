# %% [markdown]
# # 06 · Principal Component Analysis and Ordination
# **Research objectives:** Objective 1 (community structure) and Objective 2 (resistome structure).
#
# **Rationale and data representation.** Relative abundances are compositional (rows sum to a constant), so PCA on raw proportions is inappropriate. Each matrix is (i) zero-replaced by multiplicative replacement (δ = 0.65 × smallest non-zero proportion), (ii) centred log-ratio (CLR) transformed and (iii) analysed by PCA on the covariance matrix (CLR values are not rescaled, so loadings reflect log-ratio variance; this is PCA of Aitchison geometry). Features with low prevalence are removed first to limit the influence of zero replacement. For the ResFinder gene matrix, which is sparse (many genes absent from most units), PCA is replaced by principal coordinates analysis (PCoA) and non-metric multidimensional scaling (NMDS) of Bray–Curtis dissimilarities.
#
# PCA axes summarise covariation; they do not identify causes.

# %%
#BOOTSTRAP
NB = "NB06"
nb_header(NB)
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
m2 = pd.read_csv(DER / "ds02_meta.csv").set_index("Sample_ID")
ph = pd.read_csv(DER / "ds02_prokaryote_phylum_counts.csv", index_col=0).loc[m2.index]
cls = pd.read_csv(DER / "ds02_card_class_rpkm.csv", index_col=0).loc[m2.index]

def clr_pca(M, prevalence=0.9, n_comp=10, min_mean=1e-4):
    P = closure(M.values)
    keep = ((M > 0).mean() >= prevalence) & (P.mean(axis=0) >= min_mean)
    Mk = M.loc[:, keep]
    Pz, delta = mult_replace(Mk.values)
    Z = clr(Pz)
    pca = PCA(n_components=min(n_comp, Z.shape[1]), random_state=SEED).fit(Z)
    scores = pd.DataFrame(pca.transform(Z), index=M.index, columns=[f"PC{i+1}" for i in range(pca.n_components_)])
    load = pd.DataFrame(pca.components_.T, index=Mk.columns, columns=scores.columns)
    ev = pd.DataFrame({"PC": scores.columns, "eigenvalue": pca.explained_variance_, "explained_variance_pct": 100 * pca.explained_variance_ratio_,
                       "cumulative_pct": 100 * np.cumsum(pca.explained_variance_ratio_)})
    contrib = 100 * load ** 2 / (load ** 2).sum()
    total_var = Z.var(axis=0, ddof=1).sum()
    ev["pct_of_total_CLR_variance"] = 100 * ev.eigenvalue / total_var
    return scores, load, ev, contrib, dict(features_in=M.shape[1], features_kept=int(keep.sum()), delta=delta)

def scree(ax, ev, title):
    ax.bar(range(1, len(ev) + 1), ev.pct_of_total_CLR_variance, color="#4C72B0"); ax2 = ax.twinx()
    ax2.plot(range(1, len(ev) + 1), ev.pct_of_total_CLR_variance.cumsum(), "o-", color="#C44E52", ms=3); ax2.set_ylim(0, 100); ax2.set_ylabel("Cumulative %", fontsize=7)
    ax.set_xlabel("Principal component"); ax.set_ylabel("% of total CLR variance"); ax.set_title(title, fontsize=7.5); ax.set_xticks(range(1, len(ev) + 1))

def biplot(ax, scores, load, groups, colors, title, ev, n_arrows=8, markers=None):
    for g in [g for g in colors if g in set(groups)]:
        ix = groups == g
        ax.scatter(scores.loc[ix, "PC1"], scores.loc[ix, "PC2"], s=9, alpha=.6, color=colors[g], label=f"{g} (n={int(ix.sum())})", edgecolors="none")
    L = load[["PC1", "PC2"]]; top = (L ** 2).sum(axis=1).sort_values(ascending=False).index[:n_arrows]
    sc = 0.8 * np.abs(scores[["PC1", "PC2"]]).max().min() / np.abs(L.loc[top]).max().max()
    for i_, f in enumerate(top):  # AUD-15: numbered arrow tips + key below the panel (no overlapping labels)
        ax.arrow(0, 0, L.loc[f, "PC1"] * sc, L.loc[f, "PC2"] * sc, color="k", lw=.6, head_width=.08 * sc / 10, alpha=.8, length_includes_head=True)
        _r = np.hypot(L.loc[f, "PC1"], L.loc[f, "PC2"]) * sc; _rmax = np.hypot(L.loc[top, "PC1"], L.loc[top, "PC2"]).max() * sc
        _k = max(1.1, (0.18 + 0.06 * i_) * _rmax / max(_r, 1e-9))   # push short-vector labels outwards (staggered radius)
        ax.text(L.loc[f, "PC1"] * sc * _k, L.loc[f, "PC2"] * sc * _k, str(i_ + 1), fontsize=6.0, weight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=.7))
    import textwrap as _tw
    ax.text(0, -0.42, _tw.fill("Loadings: " + "; ".join(f"{i_ + 1} {f}" for i_, f in enumerate(top)), 70), transform=ax.transAxes, fontsize=5.5, va="top")
    ax.axhline(0, color="grey", lw=.4); ax.axvline(0, color="grey", lw=.4)
    ax.set_xlabel(f"PC1 ({ev.pct_of_total_CLR_variance.iloc[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({ev.pct_of_total_CLR_variance.iloc[1]:.1f}%)")
    ax.set_title(title, fontsize=7.5); ax.legend(frameon=False, fontsize=5.8, markerscale=1.3, loc="upper center", bbox_to_anchor=(.5, -.17), ncol=2, handletextpad=.2, columnspacing=.6)

# %% [markdown]
# ## 1. PCA of prokaryotic phylum composition (DS02)

# %%
sc_t, ld_t, ev_t, ct_t, info_t = clr_pca(ph, prevalence=0.9, min_mean=1e-3)
print(info_t)
save_table(ev_t.pipe(rnd, 4), "P01_ds02_taxonomy_pca_variance", f"DS02 CLR-PCA of prokaryotic phyla: explained variance ({info_t['features_kept']} of {info_t['features_in']} phyla with prevalence >= 90% and mean RA >= 0.1%)", NB, sheet="PCA_Results")
lt = ld_t[["PC1", "PC2", "PC3"]].join(ct_t[["PC1", "PC2", "PC3"]].add_suffix("_contribution_pct")).sort_values("PC1_contribution_pct", ascending=False)
save_table(lt.pipe(rnd, 5).reset_index().rename(columns={"index": "Phylum"}), "P02_ds02_taxonomy_pca_loadings", "DS02 CLR-PCA loadings and % contribution of phyla to PC1-PC3", NB, sheet="PCA_Results")
sct = sc_t[["PC1", "PC2", "PC3"]].join(m2[["City", "Season", "Collection_Date"]])
sct.to_csv(DER / "ds02_taxonomy_pca_scores.csv")
pc_city = []
for pc in ["PC1", "PC2", "PC3"]:
    for f in ["City", "Season"]:
        k = kruskal_eps2(sct, pc, f); pc_city.append(dict(component=pc, factor=f, **k))
        record_test(NB, f"DS02 taxonomy {pc} score by {f}", "Kruskal-Wallis", k["n"], k["H"], k["df"], k["p"], "epsilon-squared", k["epsilon2"], family="DS02 PCA scores", objective="Objective 1")
pc_city = pd.DataFrame(pc_city); pc_city["p_BH"] = bh(pc_city.p)
save_table(pc_city.pipe(rnd, 6), "P03_ds02_taxonomy_pc_scores_by_group", "DS02 taxonomy PC scores by city and season (Kruskal-Wallis, BH)", NB, sheet="PCA_Results")
ev_t.head(5).round(2)

# %%
fig = plt.figure(figsize=(6.9, 3.9)); gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.3, .8])
a0, a1, a2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
biplot(a0, sc_t, ld_t, m2.City.values, CITY_COLORS, "A  Scores by city with top loadings", ev_t)
for g in SEASON_ORDER:
    ix = m2.Season.values == g; a1.scatter(sc_t.loc[ix, "PC1"], sc_t.loc[ix, "PC2"], s=9, alpha=.6, color=SEASON_COLORS[g], label=f"{g} (n={int(ix.sum())})", edgecolors="none")
a1.set_xlabel(a0.get_xlabel()); a1.set_ylabel(a0.get_ylabel()); a1.set_title("B  Scores by season", fontsize=7.5); a1.legend(frameon=False, fontsize=5.8, loc="upper center", bbox_to_anchor=(.5, -.17), ncol=1, handletextpad=.2)
scree(a2, ev_t, "C  Scree plot")
fig.tight_layout()
save_fig(fig, 19, "ds02_taxonomy_pca", "Principal component analysis of prokaryotic phylum composition (DS02, CLR)",
         f"PCA of centred log-ratio transformed prokaryotic phylum composition ({info_t['features_kept']} phyla present in ≥ 90% of samples with mean relative abundance ≥ 0.1%; multiplicative zero replacement) for 447 samples. (A) Scores coloured by city with the eight highest-magnitude loading vectors (arrows scaled for display). (B) Scores coloured by season. (C) Percentage of total CLR variance per component (bars) and cumulative (line). Loadings and contributions: Tables P01-P02.",
         NB, "Objective 1", "tables/P01-P03; derived/ds02_taxonomy_pca_scores.csv")

# %% [markdown]
# ### 1b. Discrete structure in the taxonomic ordination (DS02)
# The score plot shows discrete clusters that do not align with city or season. Clusters are identified by k-means on PC1–PC3 (k chosen by silhouette, 2–6) and cross-tabulated against every available covariate. Because site identifiers and processing batches are not in the deposited metadata, the cause (persistent site characteristics vs a technical batch) cannot be determined here; the cluster label is therefore used as an adjustment variable in sensitivity analyses (notebook 09), not as a finding about ecology.

# %%
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
Z3 = sc_t[["PC1", "PC2", "PC3"]].values
sil = {k: silhouette_score(Z3, KMeans(k, n_init=20, random_state=SEED).fit_predict(Z3)) for k in range(2, 7)}
kbest = max(sil, key=sil.get)
lab = KMeans(kbest, n_init=20, random_state=SEED).fit_predict(Z3)
order_ = pd.Series(sc_t.PC1.values).groupby(lab).mean().sort_values().index
lab = pd.Series(lab, index=sc_t.index).map({o: f"T{i+1}" for i, o in enumerate(order_)})
s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0).loc[m2.index]
clus = pd.DataFrame({"taxonomic_cluster": lab}).join(m2[["City", "Season", "year", "reads_total"]])
clus.to_csv(DER / "ds02_taxonomy_clusters.csv")
from scipy.stats import chi2_contingency
cx = []
for f in ["City", "Season", "year"]:
    t = pd.crosstab(clus.taxonomic_cluster, clus[f]); c2, p, dfree, _ = chi2_contingency(t)
    v = np.sqrt(c2 / (t.values.sum() * (min(t.shape) - 1)))
    cx.append(dict(covariate=f, chi2=c2, df=dfree, p_value=p, cramers_v=v))
    record_test(NB, f"DS02 taxonomic cluster vs {f}", "Chi-square test of independence", len(clus), c2, dfree, p, "Cramer's V", v, family="DS02 taxonomic clusters", objective="QC / Objective 1")
k_d = kruskal_eps2(clus.assign(ld=np.log10(clus.reads_total)), "ld", "taxonomic_cluster"); cx.append(dict(covariate="log10 sequenced reads (Kruskal-Wallis)", chi2=k_d["H"], df=k_d["df"], p_value=k_d["p"], cramers_v=k_d["epsilon2"]))
k_a = kruskal_eps2(clus.join(s2.log10_total_RPKM), "log10_total_RPKM", "taxonomic_cluster"); cx.append(dict(covariate="log10 total ARG RPKM (Kruskal-Wallis)", chi2=k_a["H"], df=k_a["df"], p_value=k_a["p"], cramers_v=k_a["epsilon2"]))
record_test(NB, "DS02 log10 total ARG RPKM by taxonomic cluster", "Kruskal-Wallis", k_a["n"], k_a["H"], k_a["df"], k_a["p"], "epsilon-squared", k_a["epsilon2"], family="DS02 taxonomic clusters", objective="QC / Objective 2")
cx = pd.DataFrame(cx).rename(columns={"chi2": "statistic", "cramers_v": "effect (Cramer's V or epsilon-squared)"})
ra_t = closure(ph.values); ra_t = pd.DataFrame(100 * ra_t, index=ph.index, columns=ph.columns)
prof = ra_t[lt.index[:6]].groupby(lab).median().round(3)
tabc = pd.crosstab(clus.taxonomic_cluster, clus.City).join(prof)
tabc.insert(0, "n", clus.taxonomic_cluster.value_counts())
save_table(tabc.reset_index(), "P03b_ds02_taxonomic_clusters", f"DS02 discrete taxonomic clusters (k-means on PC1-PC3, k = {kbest}; silhouette {sil[kbest]:.2f}): size, city mix and median relative abundance (%) of the top-loading phyla", NB, sheet="PCA_Results")
save_table(cx.pipe(rnd, 6), "P03c_ds02_cluster_associations", "DS02 taxonomic clusters: association with city, season, year, sequencing depth and normalized ARG abundance", NB, sheet="PCA_Results")
record_finding(NB, "F13b", "Data quality / Objective 1", f"DS02 prokaryotic profiles fall into {kbest} discrete clusters (silhouette {sil[kbest]:.2f}) that are only weakly related to city (Cramer's V = {cx.iloc[0,4]:.2f}) and season (V = {cx.iloc[1,4]:.2f}) and are defined mainly by large differences in low-abundance phyla such as {', '.join(lt.index[:2])}.",
               len(clus), "DS02 samples", "k-means + silhouette; chi-square tests", "Tables P03b-P03c; Figure 19",
               "The clusters may reflect persistent site differences or a processing batch; neither can be verified without site identifiers. City contrasts are re-estimated with cluster adjustment in notebook 09.",
               "Site IDs and sequencing/classification batches are not in the deposited metadata.", "Data-quality finding")
tabc

# %% [markdown]
# ## 2. PCA of resistance-class composition (DS02)

# %%
cls_k = cls.drop(columns=[c for c in ["Unclassified"] if c in cls.columns])
sc_r, ld_r, ev_r, ct_r, info_r = clr_pca(cls_k, prevalence=0.9, min_mean=0)
save_table(ev_r.pipe(rnd, 4), "P04_ds02_resistome_pca_variance", f"DS02 CLR-PCA of resistance-class composition: explained variance ({info_r['features_kept']} class groups)", NB, sheet="PCA_Results")
save_table(ld_r[["PC1", "PC2", "PC3"]].join(ct_r[["PC1", "PC2", "PC3"]].add_suffix("_contribution_pct")).pipe(rnd, 5).reset_index().rename(columns={"index": "Class_Group"}), "P05_ds02_resistome_pca_loadings", "DS02 CLR-PCA loadings of resistance-class groups", NB, sheet="PCA_Results")
sc_r.join(m2[["City", "Season"]]).to_csv(DER / "ds02_resistome_pca_scores.csv")
k1 = kruskal_eps2(sc_r.join(m2.City), "PC1", "City")
record_test(NB, "DS02 resistome-class PC1 score by City", "Kruskal-Wallis", k1["n"], k1["H"], k1["df"], k1["p"], "epsilon-squared", k1["epsilon2"], family="DS02 PCA scores", objective="Objective 2")
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.9), gridspec_kw={"width_ratios": [1.4, .8]})
biplot(axes[0], sc_r, ld_r, m2.City.values, CITY_COLORS, "A  Resistance-class composition: scores and loadings", ev_r, n_arrows=10)
scree(axes[1], ev_r, "B  Scree plot"); fig.tight_layout()
save_fig(fig, 20, "ds02_resistome_pca", "Principal component analysis of resistance-class composition (DS02, CLR)",
         f"PCA of CLR-transformed class-group composition of total ARG RPKM ({info_r['features_kept']} class groups present in ≥ 90% of samples) for 447 samples, coloured by city, with loading vectors. PC1 by city: Kruskal-Wallis ε² = {k1['epsilon2']:.2f}.",
         NB, "Objective 2", "tables/P04-P05; derived/ds02_resistome_pca_scores.csv")

# %% [markdown]
# ## 3. Resistome ordination across matrices (ResFinder/KMA): CLR-PCA of classes, PCoA and NMDS of genes

# %%
us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0)
cfp = pd.read_csv(DER / "rf_class_fpkm.csv", index_col=0).loc[us.index]
D = np.load(DER / "rf_D_gene_bc.npy")
cf = cfp.drop(columns=[c for c in ["Unclassified", "Other"] if c in cfp.columns])
sc_c, ld_c, ev_c, ct_c, info_c = clr_pca(cf, prevalence=0.5, min_mean=0)
save_table(ev_c.pipe(rnd, 4), "P06_rf_class_pca_variance", f"ResFinder CLR-PCA of class composition: explained variance ({info_c['features_kept']} class groups present in >= 50% of units)", NB, sheet="PCA_Results")
save_table(ld_c[["PC1", "PC2", "PC3"]].join(ct_c[["PC1", "PC2", "PC3"]].add_suffix("_contribution_pct")).pipe(rnd, 5).reset_index().rename(columns={"index": "Class_Group"}), "P07_rf_class_pca_loadings", "ResFinder CLR-PCA loadings of class groups", NB, sheet="PCA_Results")
ids = [str(i) for i in range(len(us))]
po = pcoa(DistanceMatrix(D, ids))
pe = po.proportion_explained.iloc[:10]
pcs = po.samples.iloc[:, :3].set_axis(["PCoA1", "PCoA2", "PCoA3"], axis=1); pcs.index = us.index
neg = float((po.eigvals < 0).sum())
nmds = MDS(n_components=2, metric_mds=False, metric="precomputed", init="random", random_state=SEED, n_init=4, max_iter=300, normalized_stress="auto").fit(D)
nm = pd.DataFrame(nmds.embedding_, index=us.index, columns=["NMDS1", "NMDS2"])
ordv = pd.DataFrame({"axis": [f"PCoA{i+1}" for i in range(len(pe))], "proportion_explained_pct": 100 * pe.values})
ordv.loc[len(ordv)] = dict(axis="NMDS (2D) stress (normalised Stress-1)", proportion_explained_pct=float(nmds.stress_))
ordv.loc[len(ordv)] = dict(axis="Negative PCoA eigenvalues (count)", proportion_explained_pct=neg)
save_table(ordv.pipe(rnd, 4), "P08_rf_gene_pcoa_nmds", "ResFinder gene composition (Bray-Curtis): PCoA proportion explained and NMDS stress", NB, sheet="PCA_Results")
pcs.join(nm).join(sc_c[["PC1", "PC2"]]).join(us[["Sample_Type", "Stratum", "Project_Accession"]]).to_csv(DER / "rf_ordination_scores.csv")
fig, axes = plt.subplots(1, 3, figsize=(6.9, 3.6))
st = us.Sample_Type.values
biplot(axes[0], sc_c, ld_c, st, SAMPLE_TYPE_COLORS, "A  CLR-PCA of class composition", ev_c, n_arrows=8)
for ax, (x, y), ttl, xl, yl in [(axes[1], ("PCoA1", "PCoA2"), "B  PCoA of gene composition (Bray-Curtis)", f"PCoA1 ({100*pe.iloc[0]:.1f}%)", f"PCoA2 ({100*pe.iloc[1]:.1f}%)"),
                                (axes[2], ("NMDS1", "NMDS2"), f"C  NMDS (stress = {nmds.stress_:.3f})", "NMDS1", "NMDS2")]:
    dfp = pcs.join(nm)
    for g in [g for g in SAMPLE_TYPE_COLORS if g in set(st)]:
        ix = st == g; ax.scatter(dfp.loc[ix, x], dfp.loc[ix, y], s=8, alpha=.6, color=SAMPLE_TYPE_COLORS[g], label=g, edgecolors="none")
    ax.set_title(ttl, fontsize=7.5); ax.set_xlabel(xl); ax.set_ylabel(yl)
axes[2].legend(frameon=False, fontsize=5.6, loc="upper center", bbox_to_anchor=(.5, -.25), ncol=1)
fig.tight_layout()
save_fig(fig, 21, "rf_resistome_ordination", "Ordination of resistome composition across urban matrices and food-related environments (ResFinder/KMA)",
         f"(A) PCA of CLR-transformed class composition ({info_c['features_kept']} class groups present in ≥ 50% of units), with loadings. (B) PCoA and (C) NMDS of Bray-Curtis dissimilarities between relative gene-level FPKM profiles. n = {len(us)} units from {us.Project_Accession.nunique()} studies; food-related environments are contextual. Separation by matrix is confounded with study (Table C13).",
         NB, "Objective 2 / 5", "tables/P06-P08; derived/rf_ordination_scores.csv")

# %%
top_t = lt.index[:3].tolist()
e_pc1 = pc_city[(pc_city.component == "PC1") & (pc_city.factor == "City")].epsilon2.iloc[0]
record_finding(NB, "F13", "Objective 1", f"The first two CLR principal components of prokaryotic phylum composition explained {ev_t.pct_of_total_CLR_variance.iloc[0]:.1f}% and {ev_t.pct_of_total_CLR_variance.iloc[1]:.1f}% of total CLR variance; PC1 was driven mainly by {', '.join(top_t)}, and PC1 scores differed by city (epsilon-squared = {e_pc1:.2f}).",
               len(sc_t), "DS02 samples", "CLR-PCA; Kruskal-Wallis on scores", "Tables P01-P03; Figure 19",
               ("PC1 is only weakly related to city" if e_pc1 < .06 else "PC1 partly separates cities") + "; the dominant ordination structure is the discrete clustering described in F13b, not geography.", "Phylum resolution; PCA axes are descriptive.")
record_finding(NB, "F14", "Objective 2/5", f"PCoA of ResFinder gene composition explained {100*pe.iloc[0]:.1f}% and {100*pe.iloc[1]:.1f}% of Bray-Curtis variation on the first two axes (NMDS stress = {nmds.stress_:.3f}); matrices occupied partly distinct regions of ordination space.",
               len(us), "ResFinder analytic units", "PCoA, NMDS; PERMANOVA in Table C13", "Tables P06-P08; Figure 21",
               "Resistome composition differs between matrices, but visual separation largely coincides with study membership.", "Confounding of matrix with study; sparse gene matrix.")
