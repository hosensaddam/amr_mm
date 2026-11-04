# %% [markdown]
# # 19 · Sensitivity analyses: sequencing depth (AUD-16) and analysed depth of the reprocessed set
# **Research objective:** establish whether shallow libraries change any descriptive denominator or inferential result, without automatically removing them.
#
# **Thresholds.** 10 Mbp, 100 Mbp, 500 Mbp and 1 Gbp of sequenced bases per sample (ENA base counts; units of the ResFinder resource use the sum over their runs). For each threshold we report the samples retained, and re-run: (i) sample-type denominators of the evidence base; (ii) DS02 site-aware city tests for normalized ARG abundance and Shannon diversity; (iii) the ResFinder study-random-intercept model of total normalized ARG abundance by matrix, genus Shannon and gene richness by matrix, and clinically important determinant detection; (iv) for the reprocessed set, ARG gene richness after binomial thinning of aligned fragments to 50 % and 25 % of the analysed depth.
#
# **Interpretation rule.** A result is *robust* if its direction is unchanged and its 95 % CI/p-value conclusion is unchanged across thresholds; otherwise it is flagged as depth-sensitive.

# %%
#BOOTSTRAP
NB = "NB19"
nb_header(NB)
import statsmodels.formula.api as smf
md = pd.read_csv(DER / "metadata_set.csv", low_memory=False)
inc = md[md.Inclusion_Status == "Included"].copy(); inc["bases"] = pd.to_numeric(inc.bases, errors="coerce")
TH = [0, 1e7, 1e8, 5e8, 1e9]; THL = {0: "none", 1e7: ">= 10 Mbp", 1e8: ">= 100 Mbp", 5e8: ">= 500 Mbp", 1e9: ">= 1 Gbp"}
print("included:", len(inc), "| <10 Mbp:", int((inc.bases < 1e7).sum()), "| <100 Mbp:", int((inc.bases < 1e8).sum()))

# %% [markdown]
# ## 1. Evidence-base denominators by depth threshold

# %%
rows = []
for t in TH:
    k = inc[inc.bases >= t] if t else inc
    for st_, n in k.Sample_Type.value_counts().items(): rows.append(dict(threshold=THL[t], Sample_Type=st_, samples=int(n)))
    rows.append(dict(threshold=THL[t], Sample_Type="ALL", samples=len(k)))
den = pd.DataFrame(rows).pivot(index="Sample_Type", columns="threshold", values="samples").fillna(0).astype(int)[[THL[t] for t in TH]]
den["lost_at_100Mbp"] = den["none"] - den[">= 100 Mbp"]
save_table(den.reset_index(), "SD01_depth_threshold_denominators", "AUD-16: included samples by sample type before and after each sequencing-depth threshold (samples are not removed from the primary analysis)", NB, sheet="Sensitivity_Analysis")
low = inc[inc.bases < 1e8].groupby(["Dataset_ID", "Platform_Family", "Sample_Type"]).size().rename("samples_lt_100Mbp").reset_index().sort_values("samples_lt_100Mbp", ascending=False)
save_table(low, "SD02_low_depth_samples_by_dataset", "AUD-16: samples with < 100 Mbp by dataset, platform and sample type", NB, sheet="Sensitivity_Analysis")
den

# %% [markdown]
# ## 2. DS02 (CARD/RGI): site-aware city tests by depth threshold

# %%
m2 = pd.read_csv(DER / "ds02_meta.csv").set_index("Sample_ID"); s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0); al = pd.read_csv(DER / "ds02_alpha_diversity.csv", index_col=0)
d = s2[["log10_total_RPKM", "family_richness_cov50"]].join(al[["shannon"]]).join(m2[["City", "Season", "year", "Site", "Base_Count"]])
d["Base_Count"] = pd.to_numeric(d.Base_Count, errors="coerce")
rows = []
for t in TH:
    k = d[d.Base_Count >= t] if t else d
    for v in ["log10_total_RPKM", "shannon", "family_richness_cov50"]:
        ks = site_level_kw(k, v, "Site", "City", n_perm=4999)
        mm = mixed_site_lrt(k, f"{v} ~ C(City) + C(Season) + C(year)", ["C(City)"], "Site")
        rows.append(dict(threshold=THL[t], variable=v, samples=len(k), sites=ks["n_sites"], site_KW_H=ks["H"], site_perm_p=ks["p_perm"], mixed_LRT_p_city=mm["lrt"]["C(City)"]["p"], ICC_site=mm["icc"]))
sd3 = pd.DataFrame(rows)
save_table(sd3.pipe(rnd, 5), "SD03_ds02_depth_sensitivity_site_aware", "AUD-16: DS02 site-aware city tests (site-level Kruskal-Wallis, mixed-model LRT) at each depth threshold", NB, sheet="Sensitivity_Analysis")
sd3

# %% [markdown]
# ## 3. ResFinder resource: matrix model, diversity and determinant detection by depth threshold

# %%
us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0); r2u = pd.read_csv(DER / "rf_run_to_unit.csv")
rm = load_csv("run_master.csv")[["Run_ID", "Base_Count"]]; rm["Base_Count"] = pd.to_numeric(rm.Base_Count, errors="coerce")
ub = r2u.merge(rm, on="Run_ID", how="left").groupby("unit_id").Base_Count.sum(min_count=1)
us["unit_bases"] = us.index.map(ub)
urb = us[(us.Stratum == "Urban environmental")].copy()
print("ResFinder urban units with base counts:", urb.unit_bases.notna().sum(), "of", len(urb))
REF = "Municipal wastewater"; rows = []; rich = []
for t in TH:
    k = urb[urb.unit_bases >= t] if t else urb
    mdl = smf.mixedlm(f"log10_total_FPKM ~ C(Sample_Type, Treatment('{REF}'))", k, groups=k.Project_Accession).fit(reml=True)
    for term in mdl.fe_params.index:
        if "Sample_Type" not in term: continue
        lo, hi = mdl.conf_int().loc[term]
        rows.append(dict(threshold=THL[t], units=len(k), studies=k.Project_Accession.nunique(), matrix=term_label(term), fold_change_vs_MW=10 ** mdl.fe_params[term], ci_low=10 ** lo, ci_high=10 ** hi,
                         units_in_matrix=int((k.Sample_Type == term_label(term)).sum())))
    for st_, g in k.groupby("Sample_Type"):
        rich.append(dict(threshold=THL[t], Sample_Type=st_, units=len(g), median_genus_shannon=g.genus_shannon.median(), median_gene_richness=g.gene_richness.median(), median_total_FPKM=g.total_FPKM.median()))
sd4 = pd.DataFrame(rows); sd5 = pd.DataFrame(rich)
base = sd4[sd4.threshold == "none"].set_index("matrix")
sd4["same_direction_as_primary"] = [np.sign(np.log(r.fold_change_vs_MW)) == np.sign(np.log(base.loc[r.matrix, "fold_change_vs_MW"])) for r in sd4.itertuples()]
sd4["CI_excludes_1"] = (sd4.ci_low > 1) | (sd4.ci_high < 1)
save_table(sd4.pipe(rnd, 4), "SD04_rf_matrix_model_depth_sensitivity", "AUD-16: ResFinder study-random-intercept model of total normalized ARG abundance by matrix (fold change vs municipal wastewater) at each depth threshold", NB, sheet="Sensitivity_Analysis")
save_table(sd5.pipe(rnd, 4), "SD05_rf_diversity_depth_sensitivity", "AUD-16: ResFinder genus Shannon, ARG gene richness and total FPKM by matrix at each depth threshold", NB, sheet="Sensitivity_Analysis")
sd4

# %% [markdown]
# ## 4. Reprocessed set: ARG gene richness after thinning the analysed depth

# %%
REP = DER / "reprocessed"
rf = pd.read_csv(REP / "resfinder_annotated_all_runs.csv"); qc = pd.read_csv(REP / "run_qc.csv").set_index("run_accession")
rng = np.random.default_rng(SEED); rows = []
for r, g in rf.groupby("run_accession"):
    for frac in [1.0, 0.5, 0.25]:
        k = rng.binomial(g.fragmentCountAln.values.astype(int), frac) if frac < 1 else g.fragmentCountAln.values
        kk = pd.Series(k, index=g.index).groupby(g.gene).sum()
        rows.append(dict(run_accession=r, analysis_group=qc.loc[r, "analysis_group"], depth_fraction=frac, analysed_pairs=int(4e6 * frac), genes_any_hit=int((kk > 0).sum()), genes_ge2_fragments=int((kk >= 2).sum())))
th = pd.DataFrame(rows)
sm_ = th.groupby(["analysis_group", "depth_fraction"])[["genes_any_hit", "genes_ge2_fragments"]].median().reset_index()
save_table(th, "SD06_reprocessed_depth_thinning", "Reprocessed set: ARG gene richness (any aligned fragment; >= 2 fragments) after binomial thinning of aligned fragments to 50% and 25% of the 4 M-pair depth", NB, sheet="Sensitivity_Analysis")
save_table(sm_.pipe(rnd, 2), "SD07_reprocessed_depth_thinning_summary", "Reprocessed set: median ARG gene richness by group and analysed-depth fraction", NB, sheet="Sensitivity_Analysis")
sm_

# %%
fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.7))
for v, c in [("log10_total_RPKM", "#0072B2"), ("shannon", "#009E73"), ("family_richness_cov50", "#CC79A7")]:
    t = sd3[sd3.variable == v]; axes[0].plot(range(len(t)), t.site_perm_p, "-o", ms=3, color=c, label=v.replace("_", " "))
axes[0].axhline(.05, color="k", ls="--", lw=.6); axes[0].set_xticks(range(len(TH))); axes[0].set_xticklabels([THL[t] for t in TH], rotation=45, ha="right", fontsize=6)
axes[0].set_ylabel("Site-level permutation p (city)"); axes[0].set_title("A  DS02 city tests"); axes[0].legend(frameon=False, fontsize=6)
for m_, g in sd4.groupby("matrix"):
    axes[1].errorbar(range(len(g)), g.fold_change_vs_MW, yerr=[g.fold_change_vs_MW - g.ci_low, g.ci_high - g.fold_change_vs_MW], fmt="-o", ms=3, lw=.8, capsize=1.5, label=m_, color=SAMPLE_TYPE_COLORS.get(m_))
axes[1].set_yscale("log"); axes[1].axhline(1, color="k", lw=.6); axes[1].set_xticks(range(len(TH))); axes[1].set_xticklabels([THL[t] for t in TH], rotation=45, ha="right", fontsize=6)
axes[1].set_ylabel("Fold change vs municipal WW (95% CI)"); axes[1].set_title("B  ResFinder matrix model"); axes[1].legend(frameon=False, fontsize=6)
for gname, g in sm_.groupby("analysis_group"):
    axes[2].plot(g.depth_fraction * 4, g.genes_any_hit, "-o", ms=3, label=gname.replace(" comparator", "").replace(" (2016 pilot unit in ResFinder resource)", ""))
axes[2].set_xlabel("Analysed depth (million read pairs)"); axes[2].set_ylabel("Median ARG genes (any hit)"); axes[2].set_title("C  Reprocessed set thinning"); axes[2].legend(frameon=False, fontsize=6)
fig.tight_layout()
save_fig(fig, 45, "depth_sensitivity", "Sensitivity of key results to sequencing depth (AUD-16)",
         "(A) Site-level whole-site permutation p-values for DS02 city contrasts after excluding samples below each depth threshold. (B) Study-adjusted fold changes of total normalized ARG abundance by matrix (ResFinder resource) at each threshold. (C) Median ARG gene richness in the reprocessed set after binomial thinning of aligned fragments to 1 and 2 million read pairs. Thresholds in sequenced bases (ENA); samples are not removed from the primary analysis.",
         NB, "Robustness (AUD-16)", "tables/SD01-SD07")

# %%
dep_sens = sd4[~sd4.same_direction_as_primary].matrix.unique().tolist()
r_ = sd3.groupby("variable").site_perm_p.agg(lambda s: "robust" if ((s < .05).all() or (s >= .05).all()) else "depth-sensitive")
record_finding(NB, "F20", "Robustness (AUD-16)", f"Of {len(inc):,} included samples, {int(den.loc['ALL','none'] - den.loc['ALL','>= 10 Mbp'])} have < 10 Mbp and {int(den.loc['ALL','lost_at_100Mbp'])} < 100 Mbp; they affect descriptive denominators only. DS02 site-aware city conclusions were " + "; ".join(f"{k}: {v}" for k, v in r_.items()) + f" across thresholds of 10 Mbp-1 Gbp; ResFinder matrix fold changes kept their direction at every threshold{'' if not dep_sens else ' except for ' + ', '.join(dep_sens)}.",
               len(inc), "Included samples (all); DS02 samples; ResFinder urban units", "Re-analysis at 10 Mbp, 100 Mbp, 500 Mbp and 1 Gbp thresholds", "Tables SD01-SD07; Figure 45",
               "Low-depth samples do not drive the inferential results; they should nonetheless be excluded from any future reprocessing that estimates richness or detection.",
               "Base counts are ENA totals (before QC); ResFinder unit depth summed over runs.", kind="Sensitivity analysis")
