# %% [markdown]
# # 11 · Audit issue corrections: before-after verification (AUD-04, AUD-11 to AUD-20)
# **Research objective:** document, for every unresolved issue of the Independent Research Audit, what was changed, how the affected outputs moved, and whether the correction is verified.
#
# **Inputs.** Outputs of notebooks 01-10 re-executed from fresh kernels on the updated inputs (this folder), and the previous corrected outputs of the independent audit (`_previous_corrected_outputs/`, copied unchanged from `Independent_Audit/03_Corrected_Analysis_Results`).
#
# **Method.** Table-by-table comparison of the affected quantities; independent re-computation of the AUD-11 site-level statistic and of the AUD-17 fold-internal pre-processing; text search for the superseded terminology (AUD-19); figure-catalogue checks (AUD-15, AUD-20).

# %%
#BOOTSTRAP
NB = "NB11"
nb_header(NB)
PREV = OUT / "_previous_corrected_outputs" / "tables"
def prev(name): return pd.read_csv(PREV / f"{name}.csv")
def cur(name): return pd.read_csv(TAB / f"{name}.csv")
rows = []
def issue(i, status, change, before, after, evidence): rows.append(dict(issue=i, status=status, correction_implemented=change, before=before, after=after, verification_evidence=evidence))

# %% [markdown]
# ## AUD-04 — Dhaka DS84 sample types

# %%
d2p, d2c = prev("D02_by_sample_type").set_index("Sample_Type"), cur("D02_by_sample_type").set_index("Sample_Type")
sm = load_csv("sample_master.csv"); fixed = sm[sm.Sample_Type_Correction.astype(str).str.startswith("AUD-04")][["Sample_ID", "Sample_Type_Original_Harmonized", "Sample_Type", "Sample_Subtype", "Sample_Type_Rule"]]
save_table(fixed, "AC01_aud04_sample_type_corrections", "AUD-04: sample-level sample-type corrections (DS84) with source evidence", NB, sheet="Audit_Corrections")
assert len(fixed) == 4 and set(fixed.Sample_ID) == {"SAMN60675067", "SAMN60675068", "SAMN60675069", "SAMN60675070"}
issue("AUD-04", "RESOLVED", "Sample-level overrides in the updated sample master: UWW1-3 -> Municipal wastewater (BioSample isolation_source 'Urban wastewater'; BioProject PRJNA1475805: 'three urban sewage treatment plant inlets'); Control -> Surface water (isolation_source 'Pond'). Notebooks 01-10 re-run.",
      f"Hospital wastewater {int(d2p.loc['Hospital wastewater','samples'])}; Municipal {int(d2p.loc['Municipal wastewater','samples'])}; Surface water {int(d2p.loc['Surface water','samples'])}",
      f"Hospital wastewater {int(d2c.loc['Hospital wastewater','samples'])}; Municipal {int(d2c.loc['Municipal wastewater','samples'])}; Surface water {int(d2c.loc['Surface water','samples'])}",
      "AC01; D02, D05, D08 re-generated; BioSample/BioProject records archived in 01_Data/ena_records")
fixed

# %% [markdown]
# ## AUD-11 — DS02 pseudo-replication: sample-level vs site-aware inference

# %%
cmp = []
t6p = prev("T06_ds02_alpha_kruskal"); t6s = cur("T06s_ds02_alpha_site_aware_tests")
for idx in ["shannon", "rarefied_richness"]:
    a = t6p[(t6p["index"] == idx) & (t6p.factor == "City")].iloc[0]; b = t6s[(t6s["index"] == idx) & (t6s.factor == "City")].iloc[0]; c = t6s[(t6s["index"] == idx) & (t6s.factor == "C(City)")].iloc[0]
    cmp.append(dict(analysis=f"DS02 {idx} by city", previous_unit="447 samples", previous_test="Kruskal-Wallis", previous_p=a.p, updated_unit=f"{int(b.n_units)} sites", updated_test="KW on site medians (whole-site permutation) | mixed-model LRT", updated_p_site_KW=b.p_perm, updated_p_mixed_LRT=c.p))
a5p = prev("A05_ds02_arg_kruskal"); a5s = cur("A05s_ds02_arg_site_aware_tests")
for v in ["log10_total_RPKM", "family_richness_cov50", "family_shannon"]:
    a = a5p[(a5p.variable == v) & (a5p.factor == "City")].iloc[0]; b = a5s[(a5s.variable == v) & (a5s.factor == "City")].iloc[0]; c = a5s[(a5s.variable == v) & (a5s.factor == "C(City)")].iloc[0]
    cmp.append(dict(analysis=f"DS02 {v} by city", previous_unit="447 samples", previous_test="Kruskal-Wallis", previous_p=a.p, updated_unit=f"{int(b.n_units)} sites", updated_test="KW on site medians (whole-site permutation) | mixed-model LRT", updated_p_site_KW=b.p_perm, updated_p_mixed_LRT=c.p))
c1p = prev("C01_ds02_permanova_permdisp"); c1c = cur("C01_ds02_permanova_permdisp")
for dist in ["Bray-Curtis", "Bray-Curtis (relative RPKM)"]:
    a = c1p[(c1p.distance == dist) & (c1p.term == "City") & (c1p.covariates == "(none)")].iloc[0]; b = c1c[(c1c.distance == dist) & (c1c.term == "City") & (c1c.covariates == "(none)")].iloc[0]
    cmp.append(dict(analysis=f"DS02 PERMANOVA {dist} ~ City", previous_unit="447 samples (free permutation)", previous_test="PERMANOVA", previous_p=a.p_value, updated_unit="19 sites (whole-site permutation)", updated_test="PERMANOVA, block = Site", updated_p_site_KW=b.p_value, updated_p_mixed_LRT=np.nan))
a10p = prev("A10_ds02_clinical_by_city_tests"); a10s = cur("A10s_ds02_clinical_by_city_site_aware")
for r in a10s.itertuples():
    pp = a10p[a10p.determinant == r.determinant]
    cmp.append(dict(analysis=f"DS02 {r.determinant} detection by city", previous_unit="447 samples", previous_test="Chi-square (cov >= 50%)", previous_p=pp.p.iloc[0] if len(pp) else np.nan, updated_unit=f"{r.sites} sites", updated_test="Chi-square, whole-site permutation (high-confidence tier)", updated_p_site_KW=r.p_perm_whole_sites, updated_p_mixed_LRT=np.nan))
cmp = pd.DataFrame(cmp)
cmp["conclusion_changed_at_0.05"] = (cmp.previous_p < .05) != (cmp.updated_p_site_KW < .05)
save_table(cmp.pipe(rnd, 5), "AC02_aud11_pseudoreplication_before_after", "AUD-11: p-values for DS02 city contrasts before (447 samples treated as independent) and after correction (19 identified sites as the unit)", NB, sheet="Audit_Corrections")
# independent recomputation of one site-level statistic (own code, no helper)
m2 = pd.read_csv(DER / "ds02_meta.csv").set_index("Sample_ID"); s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0)
sm_ = s2.join(m2.Site).groupby("Site").agg(v=("log10_total_RPKM", "median"), c=("City", "first"))
H_own = stats.kruskal(*[g.v.values for _, g in sm_.groupby("c")])[0]
H_rep = a5s[(a5s.variable == "log10_total_RPKM") & (a5s.factor == "City")].H.iloc[0]
assert abs(H_own - H_rep) < 1e-8, (H_own, H_rep)
issue("AUD-11", "RESOLVED", "19 DS02 site identities recovered from ENA/SRA experiment titles (PRJNA1415085; 'wastewater sample from - India:<City>:<Site>'); Site added to ds02_meta; site-level Kruskal-Wallis with whole-site permutation, site random-intercept mixed models (ML LRT), whole-site-permutation PERMANOVA and chi-square implemented in NB03-NB05; sample-level tests relabelled descriptive; findings F05, F07, F08, F10 rewritten.",
      "City p-values from 447 pseudo-replicated samples (e.g., ARG abundance KW p = %.1e)" % cmp.loc[cmp.analysis == "DS02 log10_total_RPKM by city", "previous_p"].iloc[0],
      "Site-level whole-site permutation p = %.3f; mixed-model LRT p = %.2g (19 sites)" % tuple(cmp.loc[cmp.analysis == "DS02 log10_total_RPKM by city", ["updated_p_site_KW", "updated_p_mixed_LRT"]].iloc[0]),
      f"AC02; T06s, T08s, A05s, A07s, A10s, C01; independent site-level H recomputed = {H_own:.4f} (matches)")
cmp

# %% [markdown]
# ## AUD-12 — detection tiers; AUD-13 — study-matrix confounding

# %%
a9 = cur("A09_ds02_clinical_determinants"); a9p = prev("A09_ds02_clinical_determinants")
t12 = a9[a9.group == "All"][["determinant", "pct_high_confidence", "hc_ci_low", "hc_ci_high", "pct_cov80", "pct_any_hit_unconfirmed", "uncertainty_range_pp"]].merge(
      a9p[a9p.group == "All"][["determinant", "pct_cov50", "pct_any_hit"]].rename(columns={"pct_cov50": "previous_strict_cov50_pct", "pct_any_hit": "previous_lenient_pct"}), on="determinant")
save_table(t12.pipe(rnd, 2), "AC03_aud12_detection_tiers", "AUD-12: DS02 detection of clinically important determinants - previous strict/lenient bounds vs updated high-confidence (primary), >= 80% coverage and any-hit (unconfirmed) tiers", NB, sheet="Audit_Corrections")
issue("AUD-12", "RESOLVED (within available metadata)", "High-confidence tier (allele coverage >= 50% and >= 2 completely mapped reads) is primary; >= 80% coverage tier added; any-hit reported only as an unconfirmed upper bound; captions/findings rewritten; independent ResFinder/KMA cross-check of 12 Delhi samples (NB18). Identity, MAPQ and CARD release are not in the deposited RGI output (documented).",
      "Ranges cited as prevalence (e.g., KPC %.0f-%.0f%%)" % tuple(a9p[(a9p.group == "All") & (a9p.determinant == "KPC")][["pct_cov50", "pct_any_hit"]].iloc[0]),
      "KPC high-confidence %.0f%% (primary); any hit %.0f%% (unconfirmed)" % tuple(a9[(a9.group == "All") & (a9.determinant == "KPC")][["pct_high_confidence", "pct_any_hit_unconfirmed"]].iloc[0]), "AC03; A09, A10s, F08, Figure 11; GI03-GI04")
c8b = cur("C08b_rf_matrix_study_structure"); c8c = cur("C08c_rf_within_study_matrix_contrasts")
issue("AUD-13", "RESOLVED (reporting and design)", "Matrix-study structure table (C08b) flags single-study matrices; within-study matrix contrasts computed where any study sampled >= 2 matrices with >= 3 units (C08c); F12 rewritten to state single-study matrices beside each estimate and not to generalise.",
      "Fold changes for single-study matrices presented as matrix effects",
      f"Single-study matrices flagged: {', '.join(c8b[c8b.studies == 1].Sample_Type) or 'none'}; within-study contrasts available: {len(c8c)}", "C08, C08b, C08c, F12")
t12

# %% [markdown]
# ## AUD-17 — predictive-model leakage

# %%
m2p, m2c = prev("M02_modelA_performance"), cur("M02_modelA_performance")
lk = m2p[["model", "validation", "roc_auc_mean", "pr_auc_mean", "balanced_accuracy_mean"]].merge(m2c[["model", "validation", "roc_auc_mean", "pr_auc_mean", "balanced_accuracy_mean"]], on=["model", "validation"], suffixes=("_previous", "_fold_internal"))
lk["delta_roc_auc"] = lk.roc_auc_mean_fold_internal - lk.roc_auc_mean_previous
save_table(lk.pipe(rnd, 4), "AC04_aud17_model_leakage_before_after", "AUD-17: Model A performance with transductive zero replacement (previous) vs fold-internal zero replacement and CLR (updated)", NB, sheet="Audit_Corrections")
rfg = lk[(lk.model == "Random forest") & lk.validation.str.startswith("Stratified group")].iloc[0]
issue("AUD-17", "RESOLVED", "Zero-replacement constant and CLR moved inside every model pipeline (FoldwiseZeroReplaceCLR fitted on training folds only); complete Model A workflow re-executed (M01-M04, Figure 27).",
      f"RF grouped-CV ROC-AUC {rfg.roc_auc_mean_previous:.3f}", f"RF grouped-CV ROC-AUC {rfg.roc_auc_mean_fold_internal:.3f} (change {rfg.delta_roc_auc:+.3f})", "AC04; M01-M04; Figure 27")
lk

# %% [markdown]
# ## AUD-15, AUD-19, AUD-20, AUD-16, AUD-18, AUD-14

# %%
figs = load_catalog("figures"); f30 = figs[figs.number <= 30].sort_values("number")
fa = f30[["number", "file", "width_mm", "height_mm", "print_layout", "min_font_pt", "text_overlaps_flagged"]].copy()
from PIL import Image
pp = OUT / "_previous_corrected_outputs" / "figures"
fa["previous_width_mm"] = [round(Image.open(pp / f).size[0] / 300 * 25.4, 1) if (pp / f).exists() else np.nan for f in fa.file]
fa["tiff_exists"] = [(FIG / "tiff" / f.replace(".png", ".tif")).exists() for f in fa.file]
save_table(fa, "AC05_aud15_figure_audit", "AUD-15: final print width, minimum font size, flagged text overlaps and TIFF export for Figures 1-30 (previous widths for comparison)", NB, sheet="Audit_Corrections")
issue("AUD-15", "RESOLVED" if (fa.width_mm <= 190).all() and fa.tiff_exists.all() else "PARTIALLY RESOLVED",
      f"All figures re-laid out at Elsevier double-column width (<= 175-190 mm), {fa.min_font_pt.min():.1f}-9 pt Arial/sans-serif lettering, Figure 27 re-arranged 2 x 2, Figure 10 labels wrapped; RGB PNG 300 dpi + LZW TIFF 600 dpi exported; legibility audit recorded per figure.",
      f"{int((fa.previous_width_mm > 250).sum())} figures wider than 25 cm (max {fa.previous_width_mm.max():.0f} mm)",
      f"max width {fa.width_mm.max():.0f} mm; min font {fa.min_font_pt.min():.1f} pt; figures with flagged overlaps: {int((fa.text_overlaps_flagged > 0).sum())}", "AC05; figures/ and figures/tiff/")
import re as _re
SCAN_TITLE = "AUD-19: remaining occurrences of 'ARG load', 'AMR load' or 'resistance-gene burden' in notebook sources, tables and catalogues (expected: none)"
hits = []
for p in sorted((OUT / "notebook_sources").glob("*.py")) + sorted(TAB.glob("*.csv")) + sorted(CAT.glob("*.json")):
    if p.name.startswith("AC0") or p.name.startswith("11_Audit"): continue  # the scan pattern itself lives in this notebook
    txt = p.read_text(encoding="utf-8", errors="ignore").replace(SCAN_TITLE, "")   # exclude this scan's own title where it is catalogued (self-match)
    n = len(_re.findall(r"\bARG load\b|\btotal ARG load\b|\bAMR load\b|resistance-gene burden", txt, flags=_re.I))
    if n: hits.append(dict(file=p.name, occurrences=n))
hits = pd.DataFrame(hits, columns=["file", "occurrences"])
n_hits = int(hits.occurrences.sum()) if len(hits) else 0
if not len(hits): hits = pd.DataFrame([dict(file="(no file with occurrences)", occurrences=0)])   # explicit zero row keeps an Excel table in the workbook
save_table(hits, "AC06_aud19_terminology_scan", SCAN_TITLE, NB, sheet="Audit_Corrections")
issue("AUD-19", "RESOLVED" if n_hits == 0 else "PARTIALLY RESOLVED", "Terminology replaced throughout ('normalized ARG abundance'; RPKM/FPKM defined as depth- or 16S-normalized relative abundance) in titles, captions, axis labels, findings and table titles.",
      "'total ARG load' in findings, tables and figures", f"{n_hits} remaining occurrences", "AC06")
r06 = cur("R06_figure_catalog")
issue("AUD-20", "RESOLVED" if len(r06) == 30 and set(r06.number) == set(range(1, 31)) else "OPEN", "R06 written after Figure 30 (NB10) with an assertion of 30 figures; Dhaka figures (31-45) catalogued in the updated register (NB20).",
      "29 figures in R06", f"{len(r06)} figures in R06; all PNG and TIFF files present: {bool(fa.tiff_exists.all())}", "R06; AC05")
issue("AUD-16", "RESOLVED (sensitivity analysis; no automatic exclusion)", "Depth-threshold sensitivity (10 Mbp, 100 Mbp, 500 Mbp, 1 Gbp) for denominators, DS02 site-aware tests and ResFinder matrix models; thinning analysis for the reprocessed set (NB19).", "No sensitivity analysis", "See SD01-SD07 and F20", "NB19")
issue("AUD-18", "RESOLVED (provenance and boundary verification)", "Both DS01 coordinate sources retained (ENA sample-level primary for SAMEA4777503, publication city-level as alternative; not averaged); all 16 Dhaka points tested against official DNCC/DSCC polygons (COD-AB adm3 v03); SAMN60675067 inside DSCC near Pagla STP.", "One city-level coordinate for 4 DS01 samples; SAMN60675067 unconfirmed", "See DH32-DH34, Figure 43", "NB17")
issue("AUD-14", "RESOLVED (at workbook build)", "Verification_Report worksheet regenerated from the final verification run of the updated workbook (Updated_Statistical_Verification and Updated_Final_Verification_Report); table counts computed from the final workbook.", "Pre-final V02 embedded", "Final checks embedded after the last build", "Workbook Verification_Report")
AC = pd.DataFrame(rows)
save_table(AC, "AC00_audit_issue_resolution_register", "Resolution of every open audit issue (AUD-04, AUD-11 to AUD-20): correction implemented, before/after values and verification evidence", NB, sheet="Audit_Corrections")
AC[["issue", "status"]]
