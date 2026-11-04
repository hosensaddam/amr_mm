# %% [markdown]
# # 18 · Integration of the reprocessed set with the global analysis (calibration and cross-pipeline validation)
# **Research objective:** decide, with evidence, how the uniformly reprocessed Dhaka/comparator data relate to the two producer-processed collections, so that no incompatible values are pooled.
#
# **Validation design.**
# 1. *ResFinder resource calibration.* Six runs of five 2016 Global Sewage Surveillance pilot units that are already in the Martiny et al. ResFinder/KMA resource (two Indian units, Colombo, Kathmandu and Karachi) were selected for reprocessing with our pipeline. Gene-level FPKM from the resource and from our pipeline are compared for the same units (Spearman correlation, log-ratio of totals, presence/absence agreement). **This section runs only for calibration runs whose processing_status is 'completed'. If none are complete, the calibration is reported as not performed, Figure 44 is not produced, and Dhaka is not placed on the resource scale.**
# 2. *CARD/RGI cross-check (AUD-12).* The Delhi DS02 samples reprocessed with ResFinder/KMA (all completed runs of the 12 selected) are compared with the producer CARD/RGI calls for the same samples: clinically important determinants detected at high confidence by ResFinder are compared with the CARD high-confidence and any-hit tiers, indicating which CARD tier is closer to an independent call.
#
# **Decision rule.** The reprocessed set is treated as a **third, separate analytical group**. It is pooled with neither the ResFinder resource nor DS02 CARD/RGI; Dhaka is placed in the global context only descriptively (percentile of calibrated values) if, and only if, calibration shows a stable proportional relation.

# %%
#BOOTSTRAP
NB = "NB18"
nb_header(NB)
from dhaka_common import *
import seaborn as sns
meta = load_run_meta(); meta = meta[meta.processing_status == "completed"]
rf = pd.read_csv(REP / "resfinder_annotated_all_runs.csv")
r2u = pd.read_csv(DER / "rf_run_to_unit.csv"); tl = pd.read_csv(DER / "rf_template_long_all_units.csv"); us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0)
meta_all = load_run_meta()
cal_sel = meta_all[meta_all.analysis_group.str.startswith("Calibration")]
cal = meta[meta.analysis_group.str.startswith("Calibration")]
cal_units = r2u[r2u.Run_ID.isin(cal.index)][["Run_ID", "unit_id", "primary", "Sample_ID"]]
CAL_DONE = len(cal_units) > 0
print(f"Calibration runs selected: {len(cal_sel)}; completed: {len(cal)}; units available: {cal_units.unit_id.nunique()}")
print(cal_units)

# %% [markdown]
# ## 1. Calibration against the ResFinder resource (same physical samples)

# %%
rows, conc = [], []
for unit, g in (cal_units.groupby("unit_id") if CAL_DONE else []):
    rr = list(g.Run_ID)
    ours = rf[rf.run_accession.isin(rr)].copy()
    # unit-level: sum fragments over the unit's runs and renormalise by summed 16S fragments
    den = (meta.loc[rr, "bacteria_16S_fragments"] + meta.loc[rr, "archaea_16S_fragments"]).sum()
    og = ours.groupby("gene").agg(frag=("fragmentCountAln", "sum"), length=("length_bp", "mean"))
    og["FPKM_ours"] = og.frag / (og.length / 1000) / (den / 1e6)
    theirs = tl[tl.unit_id == unit].groupby("gene").FPKM.sum().rename("FPKM_resource")
    j = og[["FPKM_ours"]].join(theirs, how="outer").fillna(0)
    both = j[(j.FPKM_ours > 0) & (j.FPKM_resource > 0)]
    rho = stats.spearmanr(np.log10(both.FPKM_ours), np.log10(both.FPKM_resource))[0] if len(both) > 5 else np.nan
    top = j.sort_values("FPKM_resource", ascending=False).head(50)
    rho50 = stats.spearmanr(top.FPKM_ours, top.FPKM_resource)[0]
    rows.append(dict(unit_id=unit, runs=";".join(rr), City=us.loc[unit, "City"] if unit in us.index else "", genes_resource=int((j.FPKM_resource > 0).sum()), genes_ours=int((j.FPKM_ours > 0).sum()),
                     genes_both=len(both), jaccard_presence=len(both) / int(((j.FPKM_ours > 0) | (j.FPKM_resource > 0)).sum()),
                     spearman_log_FPKM_shared_genes=rho, spearman_top50_resource_genes=rho50, total_FPKM_resource=j.FPKM_resource.sum(), total_FPKM_ours=j.FPKM_ours.sum(),
                     log10_ratio_total_ours_vs_resource=np.log10(j.FPKM_ours.sum() / j.FPKM_resource.sum()),
                     share_of_resource_FPKM_in_shared_genes=100 * both.FPKM_resource.sum() / j.FPKM_resource.sum()))
    for gname, v in j.iterrows(): conc.append(dict(unit_id=unit, gene=gname, FPKM_ours=v.FPKM_ours, FPKM_resource=v.FPKM_resource))
calib = pd.DataFrame(rows); conc = pd.DataFrame(conc)
if not CAL_DONE:
    calib = pd.DataFrame([dict(unit_id="(none)", runs="; ".join(cal_sel.index), City="", status=f"NOT PERFORMED: 0 of {len(cal_sel)} selected calibration runs were processed (raw-read processing not completed in this update)")])
    conc = pd.DataFrame(columns=["unit_id", "gene", "FPKM_ours", "FPKM_resource"])
save_table(calib.pipe(rnd, 4), "GI01_resfinder_resource_calibration", "Calibration: our uniform pipeline vs the Martiny et al. ResFinder resource for the same 2016 pilot sewage units (gene-level agreement and ratio of total normalized ARG abundance)", NB, sheet="Global_Integration")
save_table(conc.pipe(rnd, 4), "GI02_resfinder_calibration_gene_pairs", "Calibration: gene-level FPKM in our pipeline and in the resource for each calibration unit", NB, sheet="Global_Integration")
lr = calib.log10_ratio_total_ours_vs_resource if CAL_DONE else pd.Series(dtype=float)
cal_stable = bool(CAL_DONE and len(lr) >= 3 and lr.std() < 0.15 and calib.spearman_log_FPKM_shared_genes.min() > 0.6)
print("calibration performed:", CAL_DONE, "| log10 ratio mean/SD:", lr.mean(), lr.std(), "| stable:", cal_stable)
calib

# %%
_lab = lambda c, u: f"{u} (city not reported)" if (not isinstance(c, str) or c in ("", "Not reported")) else c
if CAL_DONE:
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.1))
    for unit, g in conc.groupby("unit_id"):
        gg = g[(g.FPKM_ours > 0) & (g.FPKM_resource > 0)]
        axes[0].scatter(gg.FPKM_resource, gg.FPKM_ours, s=6, alpha=.6, label=f"{_lab(calib.set_index('unit_id').loc[unit, 'City'], unit)} (ρ={calib.set_index('unit_id').loc[unit, 'spearman_log_FPKM_shared_genes']:.2f})")
    pos = conc[(conc.FPKM_ours > 0) & (conc.FPKM_resource > 0)]
    lim = [min(pos.FPKM_resource.min(), pos.FPKM_ours.min()) / 2, max(pos.FPKM_resource.max(), pos.FPKM_ours.max()) * 2]
    axes[0].plot(lim, lim, "k--", lw=.6); axes[0].set_xscale("log"); axes[0].set_yscale("log"); axes[0].set_xlabel("Gene FPKM, ResFinder resource (Martiny et al.)"); axes[0].set_ylabel("Gene FPKM-16S, this pipeline")
    axes[0].legend(frameon=False, fontsize=5.6); axes[0].set_title("A  Same units, two pipelines")
    axes[1].bar(range(len(calib)), calib.log10_ratio_total_ours_vs_resource, color="#7F7F7F"); axes[1].axhline(0, color="k", lw=.6)
    axes[1].set_xticks(range(len(calib))); axes[1].set_xticklabels([_lab(c, u) for c, u in zip(calib.City, calib.unit_id)], rotation=45, ha="right", fontsize=6)
    axes[1].set_ylabel("log10 (total FPKM this pipeline /\ntotal FPKM resource)"); axes[1].set_title("B  Ratio of total normalized abundance")
    fig.tight_layout()
    save_fig(fig, 44, "resfinder_calibration", "Calibration of the uniform reprocessing pipeline against the ResFinder resource for the same sewage samples",
             f"(A) Gene-level normalized ARG abundance for genes detected by both pipelines in {len(calib)} 2016 Global Sewage Surveillance pilot units (axes log scale; dashed line = equality; ρ = Spearman correlation of log values). (B) log10 ratio of total normalized abundance (this pipeline / resource). The pipelines differ in KMA version (1.7.1 vs 1.3.x), ResFinder release, analysed depth (4 M pairs vs full library) and the 16S denominator (Kraken 2 SILVA 16S fragments vs KMA rRNA alignment), so values are not interchangeable without calibration.",
             NB, "Global integration", "tables/GI01-GI02")
else:
    print("Figure 44 not produced: no calibration run was processed.")

# %% [markdown]
# ## 2. Delhi: independent ResFinder/KMA calls versus producer CARD/RGI calls (AUD-12 cross-check)

# %%
de = meta[meta.analysis_group == "DS02 Delhi comparator"]
sid = dict(zip(de.index, de.Sample_ID))
amr = load_csv("amr_genes_long.csv"); amr = amr[(amr.AMR_Reference_Database == "CARD") & amr.Sample_ID.isin(de.Sample_ID)]
amr["Clinical"] = [clinical_label(g, f, "CARD") for g, f in zip(amr.AMR_Gene, amr.AMR_Gene_Family)]
amr["hc"] = (amr.Average_Percent_Coverage >= 50) & (pd.to_numeric(amr.Completely_Mapped_Reads, errors="coerce").fillna(0) >= 2)
rfd = rf[rf.run_accession.isin(de.index)]
rows = []
for labc in CLINICAL_PANEL:
    for r_, s_ in sid.items():
        rf_hc = bool(((rfd.run_accession == r_) & (rfd.Clinical == labc) & rfd.high_confidence).any())
        rf_any = bool(((rfd.run_accession == r_) & (rfd.Clinical == labc)).any())
        c_hc = bool(((amr.Sample_ID == s_) & (amr.Clinical == labc) & amr.hc).any()); c_any = bool(((amr.Sample_ID == s_) & (amr.Clinical == labc)).any())
        rows.append(dict(determinant=labc, Sample_ID=s_, ResFinder_high_conf=rf_hc, ResFinder_any_hit=rf_any, CARD_high_conf=c_hc, CARD_any_hit=c_any))
x = pd.DataFrame(rows)
agg = x.groupby("determinant").agg(n=("Sample_ID", "size"), ResFinder_high_conf=("ResFinder_high_conf", "sum"), ResFinder_any_hit=("ResFinder_any_hit", "sum"), CARD_high_conf=("CARD_high_conf", "sum"), CARD_any_hit=("CARD_any_hit", "sum"),
                                   agree_RFhc_CARDhc=("ResFinder_high_conf", lambda s: (s == x.loc[s.index, "CARD_high_conf"]).mean()),
                                   agree_RFhc_CARDany=("ResFinder_high_conf", lambda s: (s == x.loc[s.index, "CARD_any_hit"]).mean())).reindex(CLINICAL_PANEL).reset_index()
save_table(x, "GI03_delhi_card_vs_resfinder_calls", f"AUD-12 cross-check: per-sample detection of clinically important determinants in {len(de)} Delhi DS02 samples by producer CARD/RGI (high-confidence and any-hit tiers) and by independent ResFinder/KMA reprocessing", NB, sheet="Global_Integration")
save_table(agg.pipe(rnd, 3), "GI04_delhi_card_vs_resfinder_agreement", f"AUD-12 cross-check: number of Delhi samples (of {len(de)}) with each determinant by tier and pipeline, and agreement of the independent ResFinder high-confidence call with each CARD tier", NB, sheet="Global_Integration")
s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0)
S = pd.read_csv(REP / "amr_summary_all_runs.csv", index_col=0)
tot = pd.DataFrame({"Sample_ID": [sid[r] for r in de.index], "CARD_total_RPKM": [s2.loc[sid[r], "total_RPKM"] for r in de.index], "ResFinder_total_RPKM_depth": S.loc[de.index, "total_RPKM_depth"].values,
                    "ResFinder_total_FPKM_16S": S.loc[de.index, "total_FPKM_16S"].values})
rt = spearman_ci(tot.CARD_total_RPKM, tot.ResFinder_total_RPKM_depth)
save_table(tot.pipe(rnd, 4), "GI05_delhi_card_vs_resfinder_totals", f"Delhi ({len(de)} samples): producer CARD/RGI total RPKM vs independent ResFinder/KMA totals; Spearman rho = {rt['rho']:.2f} (95% CI {rt['ci_low']:.2f}-{rt['ci_high']:.2f})", NB, sheet="Global_Integration")
record_test(NB, "Delhi: CARD/RGI total RPKM vs ResFinder/KMA total RPKM-depth (same samples)", "Spearman", len(de), rt["rho"], None, rt["p"], "rho", rt["rho"], rt["ci_low"], rt["ci_high"], family="Cross-pipeline validation", objective="Global integration")
agg

# %% [markdown]
# ## 3. Updated target-city table and placement of Dhaka in the global context

# %%
d08 = pd.read_csv(TAB / "D08_target_cities.csv")
rp = meta[~meta.analysis_group.str.startswith("Calibration")].groupby(["City_group", "Dataset_ID"]).size()
d08["reprocessed_uniform_pipeline_runs"] = [int(rp.get((c, d), 0)) for c, d in zip(d08.target_city, d08.Dataset_ID)]
d08["note"] = np.where(d08.reprocessed_uniform_pipeline_runs > 0, "Uniformly reprocessed in this update (notebooks 12-16)", "")
save_table(d08, "GI06_target_cities_updated", "Proposal target cities after the update: producer-processed coverage (D08) plus runs uniformly reprocessed in this update", NB, sheet="Global_Integration")
mw = us[(us.Sample_Type == "Municipal wastewater") & us.primary.astype(bool)] if "primary" in us.columns else us[us.Sample_Type == "Municipal wastewater"]
V = pd.read_csv(REP / "comparative_variables_all_runs.csv", index_col=0)
place = []
if cal_stable:
    adj = 10 ** lr.mean()
    for c, g in V[V.Sample_Type == "Municipal wastewater"].groupby("City_group"):
        cal_v = 10 ** g.log10_total_FPKM_16S / adj
        place.append(dict(group=c, n=len(g), median_calibrated_FPKM=np.median(cal_v), percentile_in_resource_municipal_wastewater=100 * (mw.total_FPKM.values[:, None] < cal_v.values[None, :]).mean(),
                          note="Exploratory: calibrated by the mean log-ratio of 5 pilot units; not pooled"))
else:
    place.append(dict(group="(not computed)", n=np.nan, median_calibrated_FPKM=np.nan, percentile_in_resource_municipal_wastewater=np.nan,
                      note=(f"Calibration not stable enough to place reprocessed values on the resource scale (SD of log10 ratio = {lr.std():.2f}; minimum gene-level rho = {calib.spearman_log_FPKM_shared_genes.min():.2f})" if CAL_DONE else "Calibration not performed (no calibration run processed); reprocessed values are NOT placed on the ResFinder-resource scale")))
place = pd.DataFrame(place)
save_table(place.pipe(rnd, 3), "GI07_dhaka_global_placement", "Placement of reprocessed municipal-wastewater values within the ResFinder-resource municipal wastewater distribution (only if calibration is stable; descriptive, never pooled)", NB, sheet="Global_Integration")
place

# %%
cal_txt = (f"For {len(calib)} sewage units present in the ResFinder resource, gene-level normalized abundance from our pipeline correlated with the resource values (Spearman rho {calib.spearman_log_FPKM_shared_genes.min():.2f}-{calib.spearman_log_FPKM_shared_genes.max():.2f} on shared genes; presence Jaccard {calib.jaccard_presence.min():.2f}-{calib.jaccard_presence.max():.2f}), with a log10 ratio of totals of {lr.mean():.2f} ± {lr.std():.2f} (mean ± SD); calibration was judged {'stable' if cal_stable else 'not stable'} for descriptive placement."
           if CAL_DONE else f"Calibration against the ResFinder resource was NOT performed: none of the {len(cal_sel)} selected 2016-pilot calibration runs was processed, so reprocessed Dhaka values are not placed on the resource scale and remain a separate analytical group.")
better_hc = agg.agree_RFhc_CARDhc.mean() >= agg.agree_RFhc_CARDany.mean()
record_finding(NB, "DF10", "Global integration", cal_txt + f" In {len(de)} Delhi samples, independent ResFinder high-confidence calls agreed with the CARD high-confidence tier in {100*agg.agree_RFhc_CARDhc.mean():.0f}% and with the CARD any-hit tier in {100*agg.agree_RFhc_CARDany.mean():.0f}% of determinant-sample pairs; CARD and ResFinder sample totals correlated with rho = {rt['rho']:.2f} (95% CI {rt['ci_low']:.2f} to {rt['ci_high']:.2f}).",
               (len(calib) if CAL_DONE else 0) + len(de), "calibration units + Delhi samples", "Spearman; Jaccard; agreement proportions", "Tables GI01-GI07" + ("; Figure 44" if CAL_DONE else ""),
               ("The reprocessed set is kept as a separate analytical group; the Delhi cross-check supports citing the CARD high-confidence tier (not the any-hit tier) as the evidence-based detection frequency (AUD-12)." if better_hc else "The reprocessed set is kept as a separate analytical group; agreement with the CARD any-hit tier exceeded that with the high-confidence tier, so the CARD high-confidence frequencies may underestimate detection (AUD-12)."),
               ("Five calibration units; database releases and KMA versions differ; " if CAL_DONE else "No calibration performed; ") + f"{len(de)} Delhi samples only; agreement is between two annotation pipelines, not against phenotypic resistance.", kind="Empirical (validation)")
