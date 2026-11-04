# %% [markdown]
# # 13 · Dhaka City taxonomic profiling (16S rRNA reads from shotgun metagenomes)
# **Research objective:** Objective 1 for Dhaka — composition and dominant taxa of the Dhaka wastewater microbiome at phylum, class, order, family and genus level.
#
# **Input.** Kraken 2 (SILVA 138 16S index) clade counts and Bracken genus re-estimates for the 16 Dhaka runs, all processed at 4,000,000 read pairs (notebook 12).
#
# **Method.** Relative abundance at each rank = Kraken 2 clade fragments of a taxon / all prokaryotic (Bacteria + Archaea) 16S fragments assigned at that rank (the share of fragments that could be assigned at each rank is reported in DH10; 16S fragments of ~150 nt often stop at domain level). Kraken 2 clade counts (no redistribution) are the primary estimates at every rank, including genus. Bracken genus re-estimates are reported only as a sensitivity table (DH10b): with the 16S index, Bracken redistributes the many domain-level fragments to genera and can strongly inflate genera with distinctive k-mer profiles (e.g. Klebsiella), so it is not used for inference. Enterobacteriaceae genera (Klebsiella, Escherichia-Shigella, Enterobacter) are poorly separable by 16S fragments. **Species are not reported**: the SILVA 138 Kraken 2 index has no species-level taxonomy, and 16S fragments of ~150 nt cannot resolve species reliably.
#
# **Assumptions and limitations.** Relative abundance is compositional and is not microbial load; 16S copy-number variation is not corrected; read-based 16S classification of shotgun data captures only rRNA-containing fragments (0.1-0.9 % of read pairs here).

# %%
#BOOTSTRAP
NB = "NB13"
nb_header(NB)
from dhaka_common import *
import seaborn as sns
meta = load_run_meta(); meta = meta[meta.processing_status == "completed"]
dh = meta[meta.analysis_group == "Dhaka"].sort_values(["Dataset_ID", "Sample_Type", "Collection_Date"])
runs = dh.index
lab = {r: f"{m.Dataset_ID} {m.Sample_ID[-4:]} ({m.Sample_Type.split()[0][:4]}.)" for r, m in dh.iterrows()}
RANKS = {"P": "Phylum", "C": "Class", "O": "Order", "F": "Family"}
mats = {RANKS[k]: taxonomy_matrix(k, meta.index) for k in RANKS}
mats["Genus"] = taxonomy_matrix("G", meta.index)
bracken_G = taxonomy_matrix("G", meta.index, source="bracken")
for k, M in mats.items(): M.to_csv(REP / f"tax_{k.lower()}_counts_all_runs.csv")
den = meta.bacteria_16S_fragments + meta.archaea_16S_fragments
print({k: M.shape for k, M in mats.items()})

# %% [markdown]
# ## 1. Relative abundance at each rank (Dhaka)

# %%
rows = []
for rk, M in mats.items():
    Md = M.loc[runs]; tot = Md.sum(axis=1)
    RA = 100 * Md.div(tot, axis=0)
    s = pd.DataFrame({"rank": rk, "taxon": RA.columns, "mean_RA_pct": RA.mean().values, "median_RA_pct": RA.median().values, "max_RA_pct": RA.max().values,
                      "prevalence_n": (Md > 0).sum().values, "prevalence_pct": 100 * (Md > 0).mean().values})
    for d_ in DHAKA_DATASETS: s[f"mean_RA_pct_{d_}"] = RA.loc[dh.Dataset_ID == d_].mean().values
    rows.append(s.sort_values("mean_RA_pct", ascending=False))
ra_tab = pd.concat(rows, ignore_index=True)
save_table(ra_tab.pipe(rnd, 4), "DH08_dhaka_taxonomic_abundance_all_ranks", "Dhaka: mean, median and maximum relative abundance (% of prokaryotic 16S fragments assigned at that rank) and prevalence of every phylum, class, order, family and genus (Kraken 2 clade fragments / SILVA 138), with dataset-specific means; 16 samples", NB, sheet="Dhaka_Taxonomic_Abundance")
dom = ra_tab.groupby("rank", sort=False).head(5)[["rank", "taxon", "mean_RA_pct", "prevalence_n"]]
dom

# %%
RAp = 100 * mats["Phylum"].loc[runs].div(mats["Phylum"].loc[runs].sum(axis=1), axis=0)
RAg = 100 * mats["Genus"].loc[runs].div(mats["Genus"].loc[runs].sum(axis=1), axis=0)
per_sample = []
for r in runs:
    for rk, RA in [("Phylum", RAp), ("Genus", RAg)]:
        top = RA.loc[r].sort_values(ascending=False).head(5)
        per_sample.append(dict(run_accession=r, Sample_ID=dh.loc[r, "Sample_ID"], Dataset_ID=dh.loc[r, "Dataset_ID"], Sample_Type=dh.loc[r, "Sample_Type"], rank=rk,
                               dominant_taxa="; ".join(f"{t} {v:.1f}%" for t, v in top.items())))
save_table(pd.DataFrame(per_sample), "DH09_dhaka_dominant_taxa_by_sample", "Dhaka: five dominant phyla and genera in each sample (% of prokaryotic 16S fragments assigned at that rank)", NB, sheet="Dhaka_Taxonomic_Abundance")
unc = pd.DataFrame({"run_accession": runs, "Sample_ID": dh.Sample_ID.values,
                    "pct_16S_unassigned_below_phylum": 100 * (1 - mats["Phylum"].loc[runs].sum(axis=1) / den.loc[runs]).values,
                    "pct_16S_unassigned_below_family": 100 * (1 - mats["Family"].loc[runs].sum(axis=1) / den.loc[runs]).values,
                    "pct_16S_unassigned_below_genus": 100 * (1 - mats["Genus"].loc[runs].sum(axis=1) / den.loc[runs]).values,
                    "genus_assigned_fragments": mats["Genus"].loc[runs].sum(axis=1).values, "prokaryotic_16S_fragments": den.loc[runs].values})
save_table(unc.pipe(rnd, 3), "DH10_dhaka_taxonomic_resolution", "Dhaka: share of prokaryotic 16S fragments not resolved to phylum, family or genus (Kraken 2 LCA), and fragments available for genus-level estimates", NB, sheet="Dhaka_Taxonomic_Abundance")
unc.describe().round(2)

# %%
# Sensitivity: Bracken re-estimation vs Kraken 2 clade counts (genus)
RBk = 100 * mats["Genus"].loc[runs].div(mats["Genus"].loc[runs].sum(axis=1), axis=0); RBb = 100 * bracken_G.loc[runs].div(bracken_G.loc[runs].sum(axis=1), axis=0)
topk = RBk.mean().sort_values(ascending=False).index[:15]
sens = pd.DataFrame({"genus": topk, "mean_RA_pct_kraken_clade_primary": RBk[topk].mean().values, "mean_RA_pct_bracken": RBb.reindex(columns=topk).fillna(0).mean().values})
sens["ratio_bracken_to_kraken"] = sens.mean_RA_pct_bracken / sens.mean_RA_pct_kraken_clade_primary
save_table(sens.pipe(rnd, 3), "DH10b_dhaka_genus_bracken_sensitivity", "Dhaka genus sensitivity: mean relative abundance from Kraken 2 clade fragments (primary) vs Bracken re-estimation (redistributes domain-level 16S fragments; not used for inference)", NB, sheet="Dhaka_Taxonomic_Abundance")
sens

# %%
top_p = RAp.mean().sort_values(ascending=False).index[:10].tolist()
top_g = RAg.mean().sort_values(ascending=False).index[:12].tolist()
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.8), sharey=True)
pal_p = dict(zip(top_p + ["Other phyla"], sns.color_palette("tab10", 10) + [(.82, .82, .82)]))
pal_g = dict(zip(top_g + ["Other genera"], sns.color_palette("tab20", 12) + [(.82, .82, .82)]))
for ax, RA, top, pal, ttl, other in [(axes[0], RAp, top_p, pal_p, "A  Phylum", "Other phyla"), (axes[1], RAg, top_g, pal_g, "B  Genus", "Other genera")]:
    T = RA[top].copy(); T[other] = 100 - T.sum(axis=1); left = np.zeros(len(T))
    for c in T.columns:
        ax.barh(range(len(T)), T[c].values, left=left, color=pal[c], label=c, edgecolor="white", linewidth=.2); left += T[c].values
    ax.set_xlim(0, 100); ax.set_xlabel("Relative abundance (% of assigned 16S fragments)"); ax.set_title(ttl)
    ax.legend(frameon=False, fontsize=5.6, ncol=2, loc="upper center", bbox_to_anchor=(.5, -.16))
axes[0].set_yticks(range(len(runs))); axes[0].set_yticklabels([lab[r] for r in runs], fontsize=6); axes[0].invert_yaxis()
fig.tight_layout()
save_fig(fig, 32, "dhaka_taxonomic_composition", "Taxonomic composition of Dhaka wastewater and pond metagenomes (phylum and genus)",
         f"Relative abundance (% of prokaryotic 16S fragments assigned at the rank shown; Kraken 2 clade fragments, SILVA 138 16S index) of the ten most abundant phyla (A) and twelve most abundant genera (B) in each of the 16 Dhaka samples, each analysed at 4,000,000 read pairs; {unc.pct_16S_unassigned_below_family.min():.0f}-{unc.pct_16S_unassigned_below_family.max():.0f} % of 16S fragments could not be assigned below family (Table DH10), and Enterobacteriaceae genera are poorly separable by 16S fragments. Labels: dataset, last four characters of the BioSample accession and sample type (Muni. municipal/urban sewage, Hosp. hospital, Comm. community, Surf. pond). Relative abundance is a proportion of 16S fragments, not microbial load.",
         NB, "Objective 1 (Dhaka)", "tables/DH08, DH09")

# %%
G = RAg[RAg.mean().sort_values(ascending=False).index[:25]]
fig, ax = plt.subplots(figsize=(6.9, 4.8))
sns.heatmap(np.log10(G.T + 0.01), cmap="mako_r", cbar_kws={"label": "log10 relative abundance (%)", "shrink": .6}, linewidths=.2, ax=ax,
            xticklabels=[lab[r] for r in G.index], yticklabels=True)
ax.tick_params(axis="x", labelsize=5.8, rotation=90); ax.tick_params(axis="y", labelsize=6)
for i, d_ in enumerate(dh.Dataset_ID): ax.add_patch(plt.Rectangle((i, -0.9), 1, .7, color=DATASET_COL[d_], clip_on=False))
ax.set_title("Twenty-five most abundant genera in Dhaka samples", pad=16)
save_fig(fig, 33, "dhaka_genus_heatmap", "Genus-level relative abundance across the 16 Dhaka samples",
         "Heat map of log10 relative abundance (% of genus-assigned prokaryotic 16S fragments; Kraken 2 clade counts) of the 25 genera with the highest mean abundance; 0.01 % added before log transformation for display only. Colour strip above columns: dataset (DS01 blue, DS03 green, DS84 vermilion).",
         NB, "Objective 1 (Dhaka)", "tables/DH08_dhaka_taxonomic_abundance_all_ranks.csv")

# %%
tp = ra_tab[ra_tab["rank"] == "Phylum"].head(3); tg = ra_tab[ra_tab["rank"] == "Genus"].head(5)
record_finding(NB, "DF02", "Objective 1 (Dhaka)", "Across the 16 Dhaka samples (relative abundance as % of 16S fragments assigned at each rank) the prokaryotic community was dominated by " + ", ".join(f"{r.taxon} (mean {r.mean_RA_pct:.1f}%)" for r in tp.itertuples()) + "; the most abundant genera were " + ", ".join(f"{r.taxon} ({r.mean_RA_pct:.1f}%)" for r in tg.itertuples()) + ".",
               16, "Dhaka samples", "Descriptive (mean relative abundance, prevalence)", "Tables DH08-DH10; Figures 32-33",
               "Dhaka wastewater carries a typical sewage-associated community; composition differs visibly between datasets and sample types (notebook 16).",
               f"16S fragments from shotgun reads; only {100 - unc.pct_16S_unassigned_below_phylum.median():.0f}% (median) of prokaryotic 16S fragments were assigned at phylum rank and {100 - unc.pct_16S_unassigned_below_genus.median():.0f}% at genus rank (percentages use these assigned fragments as denominator); species not resolvable; Enterobacteriaceae genera poorly separable; relative, not absolute abundance; one pond sample only.", kind="Empirical (new reprocessing)")
