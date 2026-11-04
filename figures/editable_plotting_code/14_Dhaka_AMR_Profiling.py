# %% [markdown]
# # 14 · Dhaka City antimicrobial-resistance gene profiling (ResFinder/KMA)
# **Research objective:** Objective 2 for Dhaka — occurrence, normalized abundance, diversity, antibiotic classes and mechanisms of acquired ARGs in the 16 Dhaka metagenomes.
#
# **One primary pipeline and database for all samples:** KMA 1.7.1 (`-1t1`, one template per read pair) against the ResFinder database (git tag resfinder-4.7.2, 3,180 templates); class and mechanism from the ResFinder phenotype table of the same release. No CARD/RGI result is combined with these values.
#
# **Quantities.** Aligned fragments (`fragmentCountAln`, as in the Martiny et al. resource) per template; **FPKM-16S** = fragments / template kb / million prokaryotic 16S fragments (relative to bacterial abundance); **RPKM-depth** = fragments / template kb / million quality-filtered read pairs. Both are *normalized relative abundances*, not absolute ARG quantities per volume or per cell.
#
# **Detection tiers.** *High-confidence (primary):* KMA template identity ≥ 90 % and template coverage ≥ 60 % (ResFinder default thresholds). *Lenient (sensitivity):* any aligned fragment — an unconfirmed upper bound. Detection is DNA presence; it does not establish phenotypic resistance, expression, viability, host, transmission or clinical infection.

# %%
#BOOTSTRAP
NB = "NB14"
nb_header(NB)
from dhaka_common import *
import seaborn as sns, textwrap
meta = load_run_meta(); meta = meta[meta.processing_status == "completed"]
rf = resfinder_table(meta)
rf.to_csv(REP / "resfinder_annotated_all_runs.csv", index=False)
dh = meta[meta.analysis_group == "Dhaka"].sort_values(["Dataset_ID", "Sample_Type", "Collection_Date"]); runs = dh.index
lab = {r: f"{m.Dataset_ID} {m.Sample_ID[-4:]} ({m.Sample_Type.split()[0][:4]}.)" for r, m in dh.iterrows()}
rfd = rf[rf.run_accession.isin(runs)]
print(len(rf), "template rows;", rfd.run_accession.nunique(), "Dhaka runs with ARG alignments;", rfd.high_confidence.sum(), "high-confidence template detections in Dhaka")

# %% [markdown]
# ## 1. Sample-level AMR summaries

# %%
def summarise(rfx, runs_):
    out = []
    for r in runs_:
        d = rfx[rfx.run_accession == r]; h = d[d.high_confidence]
        gfp = d.groupby("gene").FPKM_16S.sum()
        out.append(dict(run_accession=r, total_FPKM_16S=d.FPKM_16S.sum(), total_RPKM_depth=d.RPKM_depth.sum(), total_FPKM_16S_high_conf=h.FPKM_16S.sum(),
                        aligned_fragments=int(d.fragmentCountAln.sum()), genes_any_hit=d.gene.nunique(), genes_high_confidence=h.gene.nunique(),
                        alleles_high_confidence=h.template.nunique(), classes_high_confidence=h.Class_Group.nunique(),
                        gene_shannon_FPKM=shannon(gfp.values), gene_simpson_FPKM=gini_simpson(gfp.values)))
    return pd.DataFrame(out).set_index("run_accession")
S = summarise(rf, meta.index); S.to_csv(REP / "amr_summary_all_runs.csv")
Sd = S.loc[runs].join(dh[["Sample_ID", "Dataset_ID", "Sample_Type", "Collection_Date"]])
save_table(Sd.reset_index().pipe(rnd, 4), "DH11_dhaka_amr_sample_summary", "Dhaka: per-sample normalized ARG abundance (FPKM-16S, RPKM-depth), gene/allele/class richness (high-confidence) and ARG diversity (ResFinder/KMA, 4 M read pairs)", NB, sheet="Dhaka_AMR_Genes")
Sd.describe().T.round(2)

# %% [markdown]
# ## 2. Gene-level occurrence and normalized abundance

# %%
GF = gene_matrix(rf, runs, "FPKM_16S"); GH = gene_matrix(rf, runs, "FPKM_16S", tier="high")
cls_of = rf.drop_duplicates("gene").set_index("gene").Class_Group; mech_of = rf.drop_duplicates("gene").set_index("gene").Mechanism
gt = pd.DataFrame({"Class_Group": GF.columns.map(cls_of), "Mechanism": GF.columns.map(mech_of), "mean_FPKM_16S": GF.mean().values, "median_FPKM_16S": GF.median().values,
                   "share_of_total_pct": 100 * GF.sum().values / GF.values.sum(), "detected_any_hit_n": (GF > 0).sum().values,
                   "detected_high_conf_n": (GH.reindex(columns=GF.columns, fill_value=0) > 0).sum().values}, index=GF.columns)
for d_ in DHAKA_DATASETS: gt[f"detected_high_conf_{d_}"] = (GH.reindex(columns=GF.columns, fill_value=0).loc[dh.Dataset_ID == d_] > 0).sum().values
ci = [wilson(k, len(runs)) for k in gt.detected_high_conf_n]
gt["prevalence_high_conf_pct"] = 100 * gt.detected_high_conf_n / len(runs); gt["prev_ci_low"] = [100 * a for a, b in ci]; gt["prev_ci_high"] = [100 * b for a, b in ci]
gt = gt.sort_values("mean_FPKM_16S", ascending=False); gt.index.name = "gene"
save_table(gt.reset_index().pipe(rnd, 4), "DH12_dhaka_amr_gene_table", f"Dhaka: all {len(gt)} ResFinder genes with aligned fragments - class, mechanism, mean/median FPKM-16S, share of total, any-hit and high-confidence detection (Wilson 95% CI) overall and by dataset", NB, sheet="Dhaka_AMR_Genes")
long_dh = rfd[["run_accession", "template", "gene", "Class", "Class_Group", "Mechanism", "Phenotype", "fragmentCountAln", "length_bp", "Template_Identity", "Template_Coverage", "Depth", "high_confidence", "FPKM_16S", "RPKM_depth", "Clinical"]].merge(dh[["Sample_ID", "Dataset_ID", "Sample_Type"]], left_on="run_accession", right_index=True)
save_table(long_dh.pipe(rnd, 4), "DH13_dhaka_amr_template_long", "Dhaka: template-level ResFinder/KMA results (aligned fragments, identity, coverage, depth, detection tier, FPKM-16S, RPKM-depth) for every sample", NB, sheet="Dhaka_AMR_Genes")
gt.head(15)

# %% [markdown]
# ## 3. Antibiotic classes and resistance mechanisms

# %%
CL = rf[rf.run_accession.isin(runs)].pivot_table(index="run_accession", columns="Class_Group", values="FPKM_16S", aggfunc="sum", fill_value=0).reindex(runs).fillna(0)
ME = rf[rf.run_accession.isin(runs)].pivot_table(index="run_accession", columns="Mechanism", values="FPKM_16S", aggfunc="sum", fill_value=0).reindex(runs).fillna(0)
CLs = 100 * CL.div(CL.sum(axis=1), axis=0); MEs = 100 * ME.div(ME.sum(axis=1), axis=0)
ct = pd.DataFrame({"share_of_total_FPKM_pct": 100 * CL.sum() / CL.values.sum(), "mean_sample_share_pct": CLs.mean(), "median_FPKM_16S": CL.median()})
for d_ in DHAKA_DATASETS: ct[f"mean_share_pct_{d_}"] = CLs.loc[dh.Dataset_ID == d_].mean()
for t_ in sorted(dh.Sample_Type.unique()): ct[f"mean_share_pct_{t_}"] = CLs.loc[dh.Sample_Type == t_].mean()
ct = ct.sort_values("share_of_total_FPKM_pct", ascending=False)
save_table(ct.reset_index().rename(columns={"index": "Class_Group"}).pipe(rnd, 3), "DH14_dhaka_resistance_classes", "Dhaka: antibiotic-class groups - share of total FPKM-16S, mean per-sample share, by dataset and by sample type", NB, sheet="Dhaka_Resistance_Classes")
mt = pd.DataFrame({"share_of_total_FPKM_pct": 100 * ME.sum() / ME.values.sum(), "mean_sample_share_pct": MEs.mean()}).sort_values("share_of_total_FPKM_pct", ascending=False)
save_table(mt.reset_index().rename(columns={"index": "Mechanism"}).pipe(rnd, 3), "DH15_dhaka_resistance_mechanisms", "Dhaka: resistance mechanisms (ResFinder phenotype table) - share of total FPKM-16S", NB, sheet="Dhaka_Resistance_Classes")
CLs.round(1).join(dh.Dataset_ID)

# %% [markdown]
# ## 4. Clinically important determinants (high-confidence vs lenient)

# %%
rows = []
for labc in CLINICAL_PANEL:
    d = rfd[rfd.Clinical == labc]
    for g in ["All Dhaka"] + DHAKA_DATASETS + sorted(dh.Sample_Type.unique()):
        u = runs if g == "All Dhaka" else (dh.index[dh.Dataset_ID == g] if g in DHAKA_DATASETS else dh.index[dh.Sample_Type == g])
        k1 = d[d.high_confidence & d.run_accession.isin(u)].run_accession.nunique(); k2 = d[d.run_accession.isin(u)].run_accession.nunique(); lo, hi = wilson(k1, len(u))
        rows.append(dict(determinant=labc, group=g, n=len(u), detected_high_confidence=k1, pct_high_confidence=100 * k1 / len(u), ci_low=100 * lo, ci_high=100 * hi,
                         detected_any_hit_unconfirmed=k2, pct_any_hit_unconfirmed=100 * k2 / len(u),
                         genes_high_confidence="; ".join(sorted(d[d.high_confidence & d.run_accession.isin(u)].allele.unique())) or "-"))
clin = pd.DataFrame(rows)
save_table(clin.pipe(rnd, 2), "DH16_dhaka_clinical_determinants", "Dhaka: clinically important determinants - high-confidence detection (primary; identity >= 90%, coverage >= 60%; Wilson 95% CI) and any-hit detection (unconfirmed upper bound), overall, by dataset and by sample type, with detected alleles", NB, sheet="Dhaka_AMR_Genes")
clin[clin.group == "All Dhaka"][["determinant", "detected_high_confidence", "pct_high_confidence", "detected_any_hit_unconfirmed", "genes_high_confidence"]]

# %%
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.9), gridspec_kw={"width_ratios": [1.25, 1]})
order_c = [c for c in CLASS_ORDER if c in CLs.columns]; left = np.zeros(len(CLs))
for c in order_c:
    axes[0].barh(range(len(CLs)), CLs[c].values, left=left, color=CLASS_COLORS[c], label=c, edgecolor="white", linewidth=.2); left += CLs[c].values
axes[0].set_yticks(range(len(CLs))); axes[0].set_yticklabels([lab[r] for r in CLs.index], fontsize=6); axes[0].invert_yaxis(); axes[0].set_xlim(0, 100)
axes[0].set_xlabel("Share of total FPKM-16S (%)"); axes[0].set_title("A  Antibiotic-class composition"); axes[0].legend(frameon=False, fontsize=5.6, ncol=2, loc="upper center", bbox_to_anchor=(.45, -.15))
for t_, col in TYPE_COL.items():
    ix = Sd.Sample_Type == t_
    axes[1].scatter(Sd.loc[ix, "total_FPKM_16S"], Sd.loc[ix, "genes_high_confidence"], s=22, color=col, label=f"{t_} (n={int(ix.sum())})", edgecolors="k", linewidths=.3)
axes[1].set_xscale("log"); from matplotlib.ticker import NullFormatter; axes[1].xaxis.set_minor_formatter(NullFormatter()); axes[1].set_xlabel("Total normalized ARG abundance\n(FPKM-16S, log scale)"); axes[1].set_ylabel("ARG genes detected (high-confidence)")
axes[1].set_title("B  Abundance vs gene richness"); axes[1].legend(frameon=False, fontsize=5.8, loc="upper center", bbox_to_anchor=(.5, -.25))
fig.tight_layout()
save_fig(fig, 34, "dhaka_resistance_classes", "Antibiotic-class composition and normalized ARG abundance in Dhaka samples (ResFinder/KMA)",
         "(A) Share of each sample's total normalized ARG abundance (FPKM-16S: aligned fragments per kb of ResFinder template per million prokaryotic 16S fragments) by harmonised antibiotic-class group; genes conferring resistance to several groups = Multi-class. (B) Total FPKM-16S versus the number of ARG genes detected at high confidence (identity >= 90 %, coverage >= 60 %), by verified sample type. 4,000,000 read pairs per sample. ARG detection is not phenotypic resistance.",
         NB, "Objective 2 (Dhaka)", "tables/DH11, DH14")

# %%
topg = gt.index[:35]
H = np.log10(GF[topg] + 1)
fig, ax = plt.subplots(figsize=(6.9, 6.6))
sns.heatmap(H.T, cmap="rocket_r", cbar_kws={"label": "log10 (FPKM-16S + 1)", "shrink": .5}, linewidths=.2, ax=ax, xticklabels=[lab[r] for r in H.index])
hc = (GH.reindex(columns=topg, fill_value=0) > 0)
for j, g in enumerate(topg):
    for i, r in enumerate(H.index):
        if hc.loc[r, g]: ax.text(i + .5, j + .5, "•", ha="center", va="center", fontsize=6, color="white")
ax.set_yticklabels([f"{g} [{cls_of.get(g)}]" for g in topg], fontsize=5.8); ax.tick_params(axis="x", labelsize=5.8, rotation=90); ax.set_ylabel(""); ax.set_xlabel("")
for i, d_ in enumerate(H.index.map(dh.Dataset_ID)): ax.add_patch(plt.Rectangle((i, -1.1), 1, .8, color=DATASET_COL[d_], clip_on=False))
ax.set_title("Thirty-five most abundant ARGs in Dhaka samples (• = high-confidence detection)", pad=16)
save_fig(fig, 35, "dhaka_amr_heatmap", "Abundance of the 35 most abundant ARGs across the 16 Dhaka samples (ResFinder/KMA)",
         "log10 (FPKM-16S + 1) per gene (alleles summed) and sample; white dots mark high-confidence detection (at least one allele with identity >= 90 % and coverage >= 60 %); cells without dots are any-hit alignments only (unconfirmed). Brackets: class group. Colour strip: dataset (DS01 blue, DS03 green, DS84 vermilion).",
         NB, "Objective 2 (Dhaka)", "tables/DH12, DH13")

# %%
fig, ax = plt.subplots(figsize=(3.35, 3.4))
t = clin[clin.group == "All Dhaka"].set_index("determinant").reindex(CLINICAL_PANEL); y = np.arange(len(t))
ax.hlines(y, t.pct_high_confidence, t.pct_any_hit_unconfirmed, color="#999999", lw=.8)
ax.errorbar(t.pct_high_confidence, y, xerr=[(t.pct_high_confidence - t.ci_low).clip(lower=0), (t.ci_high - t.pct_high_confidence).clip(lower=0)], fmt="o", ms=3.5, color="#D55E00", lw=.8, capsize=1.5, label="High-confidence (95% CI)")
ax.plot(t.pct_any_hit_unconfirmed, y, "o", ms=3.5, mfc="white", mec="#D55E00", label="Any hit (unconfirmed)")
ax.set_yticks(y); ax.set_yticklabels(CLINICAL_PANEL); ax.invert_yaxis(); ax.set_xlim(-3, 103); ax.set_xlabel("Dhaka samples with determinant (%), n = 16")
ax.legend(frameon=False, fontsize=5.8, loc="upper center", bbox_to_anchor=(.45, -.18), ncol=1); ax.set_title("Clinically important determinants, Dhaka")
save_fig(fig, 36, "dhaka_clinical_determinants", "Detection of clinically important ARG determinants in Dhaka samples",
         "Percentage of the 16 Dhaka samples in which each determinant was detected at high confidence (filled; ResFinder identity >= 90 % and coverage >= 60 %; Wilson 95 % CI) and with any aligned fragment (open; unconfirmed upper bound). OXA-48-like = blaOXA-48 family alleles; sul = sul1-sul4. DNA detection in wastewater does not indicate carriage, infection or phenotypic resistance.",
         NB, "Objective 2/3 (Dhaka)", "tables/DH16_dhaka_clinical_determinants.csv")

# %%
allc = clin[clin.group == "All Dhaka"].set_index("determinant")
top_cls = ct.index[:3].tolist()
cls_txt = ", ".join(f"{c} ({ct.loc[c, 'share_of_total_FPKM_pct']:.0f}%)" for c in top_cls)
_ww = Sd[Sd.Sample_Type != "Surface water"]
record_finding(NB, "DF03", "Objective 2 (Dhaka)", f"Dhaka samples carried {int(Sd.genes_high_confidence.min())}-{int(Sd.genes_high_confidence.max())} ARG genes at high confidence (median {Sd.genes_high_confidence.median():.1f}; genes = ResFinder gene names, with bla and mcr alleles collapsed to families, e.g. blaOXA-1/-10/-181 = blaOXA); total normalized ARG abundance ranged {Sd.total_FPKM_16S.min():,.0f}-{Sd.total_FPKM_16S.max():,.0f} FPKM-16S (median {Sd.total_FPKM_16S.median():,.0f}). The lower bounds come from the single pond-water sample; the {len(_ww)} wastewater samples carried {int(_ww.genes_high_confidence.min())}-{int(_ww.genes_high_confidence.max())} genes and {_ww.total_FPKM_16S.min():,.0f}-{_ww.total_FPKM_16S.max():,.0f} FPKM-16S. The leading class groups by share of FPKM-16S pooled across the 16 samples were {cls_txt}.",
               16, "Dhaka samples", "Descriptive (ResFinder/KMA; high-confidence tier)", "Tables DH11-DH15; Figures 34-35",
               "Dhaka wastewater carries a diverse acquired resistome; abundance is normalized to prokaryotic 16S fragments. Pooled class shares are weighted towards the samples with the highest abundance; per-sample mean shares are given in DH14/DH15.",
               "Normalized relative abundance, not absolute load; read mapping without assembly or host attribution; one database release.", kind="Empirical (new reprocessing)")
record_finding(NB, "DF04", "Objective 2/3 (Dhaka)", "High-confidence detection of clinically important determinants in Dhaka (n = 16): " + "; ".join(f"{d} {int(allc.loc[d,'detected_high_confidence'])}/16" for d in CLINICAL_PANEL) + ". Any-hit (unconfirmed) detection: " + "; ".join(f"{d} {int(allc.loc[d,'detected_any_hit_unconfirmed'])}/16" for d in CLINICAL_PANEL) + ".",
               16, "Dhaka samples", "Wilson 95% CIs (high-confidence tier)", "Table DH16; Figure 36",
               "Carbapenemase and ESBL gene DNA is present in Dhaka wastewater; the high-confidence counts are the evidence-supported values, the any-hit counts an unconfirmed upper bound.",
               "Gene DNA presence only; not phenotypic resistance, carriage or infection; 4 M read pairs limits sensitivity for rare genes.", kind="Empirical (new reprocessing)")
