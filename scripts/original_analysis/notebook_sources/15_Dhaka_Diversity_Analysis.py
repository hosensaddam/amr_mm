# %% [markdown]
# # 15 · Dhaka City microbial and resistome diversity (alpha and beta)
# **Research objective:** Objectives 1-2 for Dhaka — taxonomic richness, alpha diversity and between-sample (beta) structure of the Dhaka microbiome and resistome.
#
# **Input.** Genus-level Kraken 2 clade fragment counts (prokaryotic 16S; no Bracken redistribution) and ResFinder gene-level FPKM-16S for the 16 Dhaka runs (notebooks 13-14), all at 4,000,000 read pairs.
#
# **Methods.** Alpha: observed genus richness after rarefaction (random subsampling without replacement, fixed seed) to the smallest Dhaka genus-assigned 16S count; Shannon (natural log), Gini-Simpson and Pielou evenness; ARG gene richness (high-confidence) and Shannon. Beta: Bray-Curtis on relative abundance and Aitchison (Euclidean distance of CLR after a 0.5 pseudo-count) for genera; Bray-Curtis on relative FPKM-16S for ARG genes; PCoA and CLR-PCA (SVD). PERMANOVA (dbRDA pseudo-F) and PERMDISP for dataset and verified sample type are **descriptive**: n = 16, unbalanced, with sample type, dataset, year and instrument confounded (see notebook 16 for the within-study contrast).
#
# **Assumptions and limitations.** Compositional data are analysed with ratio-based methods; diversity estimates from ~2,000-30,000 16S fragments per sample are approximate; richness is depth-dependent and therefore rarefied.

# %%
#BOOTSTRAP
NB = "NB15"
nb_header(NB)
from dhaka_common import *
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
meta = load_run_meta(); meta = meta[meta.processing_status == "completed"]
dh = meta[meta.analysis_group == "Dhaka"].sort_values(["Dataset_ID", "Sample_Type", "Collection_Date"]); runs = dh.index
Gc = pd.read_csv(REP / "tax_genus_counts_all_runs.csv", index_col=0).loc[runs]
Gc = Gc.loc[:, Gc.sum() > 0]
rf = pd.read_csv(REP / "resfinder_annotated_all_runs.csv"); S = pd.read_csv(REP / "amr_summary_all_runs.csv", index_col=0)
GF = gene_matrix(rf, runs, "FPKM_16S")
depth = int(Gc.sum(axis=1).min()); print("genus-assigned 16S fragments:", Gc.sum(axis=1).describe().round(0).to_dict(), "rarefaction depth:", depth)

# %% [markdown]
# ## 1. Alpha diversity

# %%
R = rarefy_counts(Gc.round().astype(int), depth)
A = pd.DataFrame({"genus_fragments": Gc.sum(axis=1), "observed_genera": (Gc > 0).sum(axis=1), "rarefied_genus_richness": (R > 0).sum(axis=1),
                  "genus_shannon": [shannon(r) for r in Gc.values], "genus_gini_simpson": [gini_simpson(r) for r in Gc.values],
                  "rarefied_genus_shannon": [shannon(r) for r in R.values]}, index=runs)
A["genus_pielou"] = A.genus_shannon / np.log(A.observed_genera)
A = A.join(S.loc[runs, ["genes_high_confidence", "gene_shannon_FPKM", "gene_simpson_FPKM", "total_FPKM_16S"]]).join(dh[["Sample_ID", "Dataset_ID", "Sample_Type", "Collection_Date", "Period"]])
A.to_csv(REP / "dhaka_alpha_diversity.csv")
save_table(A.reset_index().pipe(rnd, 4), "DH17_dhaka_alpha_diversity", f"Dhaka: genus alpha diversity (rarefied to {depth:,} genus-assigned 16S fragments; Shannon, Gini-Simpson, Pielou) and ARG diversity (high-confidence gene richness, Shannon of gene FPKM-16S) per sample", NB, sheet="Dhaka_Alpha_Diversity")
rs = spearman_ci(np.log10(A.genus_fragments), A.observed_genera); rs2 = spearman_ci(np.log10(A.genus_fragments), A.rarefied_genus_richness)
dchk = pd.DataFrame([dict(check="Observed genus richness vs log10 genus fragments (depth dependence)", **rs), dict(check="Rarefied genus richness vs log10 genus fragments", **rs2)])
save_table(dchk.pipe(rnd, 4), "DH18_dhaka_alpha_depth_check", "Dhaka: depth dependence of genus richness before and after rarefaction (Spearman, n = 16)", NB, sheet="Dhaka_Alpha_Diversity")
by = A.groupby(["Dataset_ID", "Sample_Type"])[["rarefied_genus_richness", "genus_shannon", "genes_high_confidence", "gene_shannon_FPKM"]].agg(["median", "min", "max", "count"])
by.columns = ["_".join(c) for c in by.columns]
save_table(by.reset_index().pipe(rnd, 3), "DH19_dhaka_alpha_by_dataset_type", "Dhaka: alpha-diversity summaries by dataset and verified sample type (descriptive)", NB, sheet="Dhaka_Alpha_Diversity")
by

# %% [markdown]
# ## 2. Beta diversity, PCoA and CLR-PCA

# %%
RA = closure(Gc.values); D_bc = squareform(pdist(RA, "braycurtis"))
CLRg = clr(closure(Gc.values + 0.5)); D_ait = squareform(pdist(CLRg, "euclidean"))
D_arg = squareform(pdist(closure(GF.values), "braycurtis"))
ids = list(runs)
pc_bc = pcoa(DistanceMatrix(D_bc, ids)); pc_arg = pcoa(DistanceMatrix(D_arg, ids))
U, sv, Vt = np.linalg.svd(CLRg - CLRg.mean(axis=0), full_matrices=False); ev = sv ** 2 / (sv ** 2).sum()
pcs = pd.DataFrame({"PCoA1_genus_BC": pc_bc.samples.iloc[:, 0].values, "PCoA2_genus_BC": pc_bc.samples.iloc[:, 1].values,
                    "PCoA1_ARG_BC": pc_arg.samples.iloc[:, 0].values, "PCoA2_ARG_BC": pc_arg.samples.iloc[:, 1].values,
                    "CLR_PC1": U[:, 0] * sv[0], "CLR_PC2": U[:, 1] * sv[1]}, index=runs).join(dh[["Sample_ID", "Dataset_ID", "Sample_Type"]])
var = pd.DataFrame({"axis": ["1", "2", "3"], "PCoA_genus_BC_pct": 100 * pc_bc.proportion_explained.values[:3], "PCoA_ARG_BC_pct": 100 * pc_arg.proportion_explained.values[:3], "CLR_PCA_genus_pct": 100 * ev[:3]})
load = pd.DataFrame({"genus": Gc.columns, "CLR_PC1_loading": Vt[0], "CLR_PC2_loading": Vt[1]}).assign(abs1=lambda d: d.CLR_PC1_loading.abs()).sort_values("abs1", ascending=False).drop(columns="abs1").head(20)
save_table(pcs.reset_index().pipe(rnd, 4), "DH20_dhaka_ordination_scores", "Dhaka: PCoA (genus Bray-Curtis; ARG Bray-Curtis) and CLR-PCA (genus) sample scores", NB, sheet="Dhaka_PCA_Results")
save_table(var.pipe(rnd, 3), "DH21_dhaka_ordination_variance", "Dhaka: variance explained by the first three ordination axes", NB, sheet="Dhaka_PCA_Results")
save_table(load.pipe(rnd, 4), "DH22_dhaka_clr_pca_loadings", "Dhaka: 20 genera with the largest absolute CLR-PC1 loadings", NB, sheet="Dhaka_PCA_Results")
mdh = dh[["Dataset_ID", "Sample_Type", "Period"]].reset_index()
prow = []
for nm, D in [("Genus Bray-Curtis", D_bc), ("Genus Aitchison (CLR)", D_ait), ("ARG gene Bray-Curtis (relative FPKM-16S)", D_arg)]:
    for term in ["Dataset_ID", "Sample_Type"]:
        r = permanova_dbrda(D, mdh, term, (), n_perm=9999); r.update(distance=nm, status="Descriptive (n = 16; dataset, sample type, year and platform confounded)"); prow.append(r)
        record_test(NB, f"Dhaka {nm} ~ {term}", "PERMANOVA (dbRDA pseudo-F; free permutation)", r["n"], r["pseudo_F"], f"{r['df_term']}, {r['df_resid']}", r["p_value"], "R2", r["R2_partial"], family="Dhaka beta diversity [descriptive]", note="confounded design", objective="Objective 1-2 (Dhaka)")
        d = permdisp_test(D, mdh[term].astype(str).values, 9999)
        prow.append(dict(term=term, covariates="(PERMDISP)", pseudo_F=d["F"], df_term=d["groups"] - 1, df_resid=d["n"] - d["groups"], R2_partial=np.nan, p_value=d["p_value"], n=d["n"], permutations=9999, permutation_scheme="PERMDISP", distance=nm, status="Descriptive"))
perm = pd.DataFrame(prow)[["distance", "term", "covariates", "n", "df_term", "df_resid", "pseudo_F", "R2_partial", "p_value", "permutation_scheme", "status"]]
save_table(perm.pipe(rnd, 5), "DH23_dhaka_beta_permanova", "Dhaka: PERMANOVA and PERMDISP of genus and ARG composition by dataset and sample type (descriptive; confounded design, n = 16)", NB, sheet="Dhaka_Beta_Diversity")
dl = []
iu = np.triu_indices(len(runs), 1)
for nm, D in [("Genus Bray-Curtis", D_bc), ("ARG Bray-Curtis", D_arg)]:
    for i, j in zip(*iu):
        a, b = runs[i], runs[j]
        dl.append(dict(distance=nm, run_a=a, run_b=b, dataset_pair="-".join(sorted([dh.loc[a, "Dataset_ID"], dh.loc[b, "Dataset_ID"]])), within_dataset=dh.loc[a, "Dataset_ID"] == dh.loc[b, "Dataset_ID"], value=D[i, j]))
dl = pd.DataFrame(dl)
dsum = dl.groupby(["distance", "dataset_pair"]).value.agg(["median", "min", "max", "count"]).reset_index()
save_table(dsum.pipe(rnd, 4), "DH24_dhaka_beta_within_between_datasets", "Dhaka: Bray-Curtis dissimilarity within and between datasets (median, range of pairwise values)", NB, sheet="Dhaka_Beta_Diversity")
mt = mantel_r = None
from skbio.stats.distance import mantel
mt = mantel(DistanceMatrix(D_ait, ids), DistanceMatrix(D_arg, ids), method="spearman", permutations=9999, seed=SEED)
save_table(pd.DataFrame([dict(test="Mantel (Spearman): genus Aitchison vs ARG Bray-Curtis (Dhaka)", r=mt[0], p_value=mt[1], n=mt[2], permutations=9999)]).pipe(rnd, 4), "DH25_dhaka_microbiome_resistome_mantel", "Dhaka: concordance of microbiome and resistome structure (Mantel)", NB, sheet="Dhaka_Beta_Diversity")
record_test(NB, "Dhaka genus (Aitchison) vs ARG (Bray-Curtis) structure", "Mantel (Spearman)", int(mt[2]), mt[0], None, mt[1], "Mantel r", mt[0], family="Dhaka microbiome-resistome association", objective="Objective 2 (Dhaka)")
perm

# %%
fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.6))
mk = {"DS01": "o", "DS03": "s", "DS84": "^"}
for ax, (x, y, xl, yl, ttl) in zip(axes, [("PCoA1_genus_BC", "PCoA2_genus_BC", f"PCoA1 ({var.PCoA_genus_BC_pct[0]:.0f}%)", f"PCoA2 ({var.PCoA_genus_BC_pct[1]:.0f}%)", "A  Genus, Bray-Curtis"),
                                          ("CLR_PC1", "CLR_PC2", f"PC1 ({var.CLR_PCA_genus_pct[0]:.0f}%)", f"PC2 ({var.CLR_PCA_genus_pct[1]:.0f}%)", "B  Genus, CLR-PCA"),
                                          ("PCoA1_ARG_BC", "PCoA2_ARG_BC", f"PCoA1 ({var.PCoA_ARG_BC_pct[0]:.0f}%)", f"PCoA2 ({var.PCoA_ARG_BC_pct[1]:.0f}%)", "C  ARG genes, Bray-Curtis")]):
    for (d_, t_), g in pcs.groupby(["Dataset_ID", "Sample_Type"]):
        ax.scatter(g[x], g[y], marker=mk[d_], s=24, color=TYPE_COL.get(t_, "k"), edgecolors="k", linewidths=.3, label=f"{d_} · {t_.replace(' wastewater', ' WW')} ({len(g)})")
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(ttl); ax.axhline(0, color=".85", lw=.5, zorder=0); ax.axvline(0, color=".85", lw=.5, zorder=0)
h, l = axes[0].get_legend_handles_labels(); fig.legend(h, l, frameon=False, fontsize=5.8, ncol=3, loc="lower center", bbox_to_anchor=(.5, -.12))
fig.tight_layout(rect=(0, .08, 1, 1))
save_fig(fig, 37, "dhaka_beta_diversity", "Between-sample structure of the Dhaka microbiome and resistome",
         "Ordination of the 16 Dhaka samples: (A) PCoA of genus Bray-Curtis dissimilarity; (B) PCA of CLR-transformed genus counts (0.5 pseudo-count); (C) PCoA of Bray-Curtis dissimilarity of relative ARG gene abundance (FPKM-16S). Marker = dataset (circle DS01, square DS03, triangle DS84); colour = verified sample type (blue municipal/urban sewage, vermilion hospital, light blue community, green pond). Dataset, year, instrument and sample type are partly confounded; PERMANOVA in Table DH23 is descriptive.",
         NB, "Objectives 1-2 (Dhaka)", "tables/DH20-DH23")

# %%
fig, axes = plt.subplots(1, 4, figsize=(6.9, 2.4))
for ax, (c, yl) in zip(axes, [("rarefied_genus_richness", f"Genus richness\n(rarefied, {depth:,} fragments)"), ("genus_shannon", "Genus Shannon"), ("genes_high_confidence", "ARG genes\n(high-confidence)"), ("gene_shannon_FPKM", "ARG Shannon")]):
    for i, d_ in enumerate(DHAKA_DATASETS):
        g = A[A.Dataset_ID == d_]
        ax.scatter(np.full(len(g), i) + np.linspace(-.18, .18, len(g)), g[c], c=[TYPE_COL.get(t, "k") for t in g.Sample_Type], s=16, edgecolors="k", linewidths=.3)
        ax.hlines(g[c].median(), i - .28, i + .28, color="k", lw=.9)
    ax.set_xticks(range(3)); ax.set_xticklabels(DHAKA_DATASETS); ax.set_ylabel(yl, fontsize=6.3)
axes[0].set_title("A", loc="left"); axes[1].set_title("B", loc="left"); axes[2].set_title("C", loc="left"); axes[3].set_title("D", loc="left")
from matplotlib.lines import Line2D
fig.legend([Line2D([], [], marker="o", ls="", mfc=v, mec="k", mew=.3, ms=4) for v in TYPE_COL.values()], list(TYPE_COL), frameon=False, ncol=4, fontsize=5.8, loc="lower center", bbox_to_anchor=(.5, -.1))
fig.tight_layout(rect=(0, .07, 1, 1))
save_fig(fig, 38, "dhaka_alpha_diversity", "Alpha diversity of the Dhaka microbiome and resistome by dataset",
         f"(A) Rarefied genus richness ({depth:,} genus-assigned 16S fragments per sample); (B) genus Shannon diversity; (C) number of ARG genes detected at high confidence; (D) Shannon diversity of gene-level FPKM-16S. Points are samples coloured by verified sample type; bars are dataset medians (DS01 n = 4, DS03 n = 4, DS84 n = 8). Descriptive; datasets differ in year, instrument and sampling sites.",
         NB, "Objectives 1-2 (Dhaka)", "tables/DH17, DH19")

# %%
record_finding(NB, "DF05", "Objective 1 (Dhaka)", f"Rarefied genus richness in Dhaka samples ranged {int(A.rarefied_genus_richness.min())}-{int(A.rarefied_genus_richness.max())} (median {A.rarefied_genus_richness.median():.0f}; {depth:,} fragments) and genus Shannon diversity {A.genus_shannon.min():.2f}-{A.genus_shannon.max():.2f}; dataset explained {100*perm[(perm.distance=='Genus Bray-Curtis')&(perm.term=='Dataset_ID')&(perm.covariates=='(none)')].R2_partial.iloc[0]:.0f}% of genus Bray-Curtis variation (descriptive PERMANOVA) and microbiome and resistome structure were {'concordant' if mt[1] < .05 else 'not significantly concordant'} (Mantel r = {mt[0]:.2f}, permutation p = {mt[1]:.4f}).",
               16, "Dhaka samples", "Rarefied richness; PERMANOVA/PERMDISP (descriptive); Mantel", "Tables DH17-DH25; Figures 37-38",
               "Between-dataset differences are visible, but they combine study, year, instrument, site and sample-type effects and cannot be attributed to any single factor.",
               "n = 16; confounded design; 16S-fragment based diversity; relative composition only.", kind="Empirical (new reprocessing)")
