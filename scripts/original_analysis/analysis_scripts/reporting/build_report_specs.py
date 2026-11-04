"""Build JSON specs for the two Word reports from the final tables (data-driven; no hand-typed results).
Usage: python build_report_specs.py <tables_dir> <figures_dir> <out_dir>"""
import sys, json
from pathlib import Path
import pandas as pd, numpy as np
from PIL import Image
T, FG, O = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); O.mkdir(exist_ok=True, parents=True)
rd = lambda n: pd.read_csv(T / f"{n}.csv", low_memory=False)
def f(x, d=2):
    if pd.isna(x): return "–"
    if isinstance(x, (int, np.integer)): return f"{x:,}"
    x = float(x)
    if x != 0 and (abs(x) < 10 ** -d or abs(x) >= 1e6): return f"{x:.2e}"
    return f"{x:,.{d}f}"
def tab(df, cols, names=None, widths=None, size=17):
    rows = [[(f(v, 3) if isinstance(v, (float, np.floating)) else ("–" if pd.isna(v) else str(v))) for v in r] for r in df[cols].itertuples(index=False)]
    return {"columns": names or cols, "rows": rows, "widths": widths, "size": size}
def fig(num, width_px=600):
    cat = rd("R06u_updated_figure_catalog").set_index("number").loc[num]
    p = FG / cat.file; w, h = Image.open(p).size
    small = O / f"_fig{num:02d}.png"
    im = Image.open(p); im.thumbnail((2000, 2000)); im.save(small)
    return {"image": str(small), "w": w, "h": h, "width_px": width_px, "caption": f"**Figure {num}.** {cat.title}. {cat.caption}"}
F = rd("R03u_updated_findings_register").set_index("id")
_q = rd("DH05_sequencing_qc_all_runs")
N_RUNS = len(_q)
_c = _q[_q.analysis_group == "DS01 comparator"]
_cc = _c.groupby("City_group").size() if "City_group" in _c.columns else _c.groupby("City").size()
_yrs = pd.to_datetime(_c.Collection_Date.astype(str).str[:10], errors="coerce", format="mixed").dt.year
GSS_TXT = f"{len(_c)} Global Sewage Surveillance municipal sewage samples from " + ", ".join(f"{k} ({v})" for k, v in _cc.sort_values(ascending=False).items()) + f" collected {int(_yrs.min())}-{int(_yrs.max())}"
N_DELHI = int((_q.analysis_group == "DS02 Delhi comparator").sum())
_g1 = rd("GI01_resfinder_resource_calibration")
CAL_DONE = "spearman_log_FPKM_shared_genes" in _g1.columns
CAL_TXT = (f"{int((_q.analysis_group.str.startswith('Calibration')).sum())} runs of {len(_g1)} 2016 pilot sewage units already present in the ResFinder resource ({', '.join(('one Indian unit without a reported city' if str(c) in ('Not reported', 'nan', '') else str(c)) for c in _g1.City)}) were processed as a calibration set." if CAL_DONE else "The calibration set was not processed.")
AC = rd("AC00_audit_issue_resolution_register")
qc = rd("DH06_dhaka_sequencing_qc"); inv = rd("DH01_dhaka_run_inventory"); MIN_PAIRS = float(pd.to_numeric(inv.ENA_read_pairs, errors="coerce").min()) / 1e6; S = rd("DH11_dhaka_amr_sample_summary")
cl = rd("DH16_dhaka_clinical_determinants"); ct = rd("DH14_dhaka_resistance_classes"); ta = rd("DH08_dhaka_taxonomic_abundance_all_ranks")
T29 = rd("DH29_dhaka_comparative_statistical_tests"); geo = rd("DH32_dhaka_geographic_verification"); gi = rd("GI01_resfinder_resource_calibration"); g4 = rd("GI04_delhi_card_vs_resfinder_agreement")
al = rd("DH17_dhaka_alpha_diversity"); c28 = rd("DH28_south_asia_sewage_city_summary"); V4 = rd("V04_updated_analysis_checks")
lim = rd("R04u_updated_limitations")

# ------------------------------------------------------------------ Dhaka comparative report
B = []
B += [{"h1": "1. Purpose and scope"},
      {"p": "This report presents the first sample-level microbiome and antimicrobial-resistance (AMR) profiles for the 16 Dhaka City metagenomes included in the chapter evidence base, and compares Dhaka with Karachi, Delhi and other South Asian cities using one uniform bioinformatic workflow. It extends, and does not replace, the corrected global analysis. All values are traceable to the tables of AMR_Updated_Final_Analysis.xlsx and Dhaka_City_Advanced_Analysis.xlsx; gene detection is reported as DNA presence and never as phenotypic resistance, expression, transmission or clinical infection."},
      {"h1": "2. Data and methods"},
      {"h2": "2.1 Samples"},
      {"p": f"The {len(inv)} Dhaka samples come from three public datasets: DS01 (Global Sewage Surveillance, 4 municipal sewage samples, 2017-2018, Illumina NovaSeq 6000), DS03 (BCSIR hospital and community wastewater, 4 samples, June 2025, NextSeq 2000) and DS84 (PRJNA1475805, 8 samples, June-November 2024, NextSeq 500). After the AUD-04 correction, DS84 comprises four hospital wastewater samples, three urban sewage-treatment-plant inlets (classified as municipal wastewater on the basis of the BioSample isolation source and the BioProject description) and one pond-water comparator. Comparators were {GSS_TXT}, and {N_DELHI} Delhi open-drain samples from DS02 (three time points at each of the four Delhi sites, 2022-2024). {CAL_TXT} All {N_RUNS} selected runs were processed."},
      {"table": tab(inv, ["Sample_ID", "run_accession", "Dataset_ID", "Sample_Type", "Collection_Date", "instrument_model", "ENA_Gbp"], ["BioSample", "Run", "Dataset", "Verified sample type", "Date", "Instrument", "ENA Gbp"], [.15, .13, .09, .22, .12, .17, .12]),
       "title": "**Table 1.** Dhaka samples, runs and verified sample types (source: DH01)."},
      {"h2": "2.2 Uniform processing workflow"},
      {"p": f"For every run, the first 4,000,000 read pairs were streamed from the European Nucleotide Archive, so that all samples were analysed at the same depth (the smallest Dhaka library contains {MIN_PAIRS:.2f} million read pairs). Reads were quality-controlled with fastp 1.4.0 (adapter auto-detection, poly-G trimming, sliding-window Q20 trimming, minimum length 50 nt). Prokaryotic 16S rRNA fragments were classified with Kraken 2.17.1 against the SILVA 138 16S index (confidence 0.1), with Bracken re-estimation at genus level. Acquired ARGs were identified with KMA 1.7.1 against the ResFinder database (release tag resfinder-4.7.2), assigning each read pair to one template. Normalized ARG abundance is expressed as FPKM-16S (aligned fragments per kilobase of template per million prokaryotic 16S fragments) and, as a sensitivity measure, per million quality-filtered read pairs. High-confidence detection requires template identity of at least 90 % and template coverage of at least 60 % (ResFinder defaults); detection from any aligned fragment is reported separately as an unconfirmed upper bound."},
      {"h2": "2.3 Statistical design"},
      {"p": "The independent unit is stated for every test (Table DH29). Groups of three to seven samples were compared with exact permutation tests of median differences, Mann-Whitney U tests, Cliff's delta and ratios of medians with bootstrap confidence intervals; composition was compared with PERMANOVA and PERMDISP on Bray-Curtis dissimilarities. For Dhaka versus Delhi, the Delhi site and the Dhaka sampling unit are the independent units, and PERMANOVA permutes whole units. Between-dataset contrasts within Dhaka are descriptive because dataset, year, instrument, site and sample type are confounded; the only within-study sample-type contrast is DS84 (four hospital and three urban sewage samples), for which only 35 label arrangements exist, so exact p-values move in steps of 1/35 (0.029)."},
      {"h1": "3. Results"},
      {"h2": "3.1 Data usability and quality"},
      {"p": F.loc["DF01", "finding"]}, fig(31),
      {"h2": "3.2 Dhaka microbial composition and diversity"},
      {"p": F.loc["DF02", "finding"]}, {"p": F.loc["DF05", "finding"]},
      {"table": tab(ta[ta["rank"].isin(["Phylum", "Genus"])].groupby("rank", sort=False).head(8), ["rank", "taxon", "mean_RA_pct", "median_RA_pct", "prevalence_n", "mean_RA_pct_DS01", "mean_RA_pct_DS03", "mean_RA_pct_DS84"],
                    ["Rank", "Taxon", "Mean %", "Median %", "Prevalence (of 16)", "Mean % DS01", "Mean % DS03", "Mean % DS84"], [.09, .25, .1, .1, .12, .11, .11, .12]),
       "title": "**Table 2.** Dominant phyla and genera in Dhaka samples (% of prokaryotic 16S fragments; source: DH08)."},
      fig(32), fig(37), fig(38),
      {"h2": "3.3 Dhaka resistome"},
      {"p": F.loc["DF03", "finding"]}, {"p": F.loc["DF04", "finding"]},
      {"table": tab(cl[cl.group == "All Dhaka"], ["determinant", "detected_high_confidence", "pct_high_confidence", "ci_low", "ci_high", "detected_any_hit_unconfirmed", "genes_high_confidence"],
                    ["Determinant", "High-confidence (n of 16)", "%", "95% CI low", "95% CI high", "Any hit (unconfirmed)", "Alleles detected (high-confidence)"], [.12, .13, .08, .1, .1, .13, .34]),
       "title": "**Table 3.** Clinically important determinants in Dhaka (source: DH16)."},
      fig(34), fig(35), fig(36),
      {"h2": "3.4 Within-Dhaka comparisons"},
      {"p": F.loc["DF06", "finding"]}, fig(39),
      {"h2": "3.5 Dhaka versus Karachi and other South Asian sewage (matched protocol)"},
      {"p": F.loc["DF07", "finding"]},
      {"table": tab(T29[T29.level.str.startswith("Dhaka vs South") & T29.comparison.str.startswith("Dhaka vs ")], ["comparison", "metric", "n_group_a", "n_group_b", "statistic", "p_value", "p_BH_within_level", "effect_size", "detail"],
                    ["Comparison", "Metric", "n Dhaka", "n comparator", "Median difference", "Exact p", "BH q", "Cliff's delta", "Ratio of medians (FPKM)"], [.11, .18, .07, .08, .1, .08, .08, .09, .21], size=15),
       "title": "**Table 4.** Dhaka versus each Global Sewage Surveillance comparator city (source: DH29)."},
      fig(42),
      {"h2": "3.6 Dhaka versus Delhi"},
      {"p": F.loc["DF08", "finding"]},
      {"table": tab(T29[T29.level == "Dhaka vs Delhi"], ["comparison", "metric", "independent_unit", "n_group_a", "n_group_b", "p_value", "effect_size", "status"],
                    ["Comparison", "Metric", "Independent unit", "n a", "n b", "p", "Effect size", "Status"], [.16, .16, .2, .06, .06, .08, .1, .18], size=15),
       "title": "**Table 5.** Dhaka versus Delhi tests (source: DH29)."},
      {"table": tab(rd("DH29b_dhaka_delhi_unit_definition_sensitivity"), ["unit_definition", "metric", "n_dhaka_units", "n_delhi_sites", "median_difference_Dhaka_minus_Delhi", "p_exact"], ["Dhaka unit definition", "Metric", "Dhaka units", "Delhi sites", "Median difference", "Exact p"], [.34, .2, .09, .09, .14, .14], size=15),
       "title": "**Table 5b.** Sensitivity of the Dhaka-versus-Delhi tests to the definition of the Dhaka sampling unit (source: DH29b)."},
      fig(40), fig(41),
      {"h2": "3.7 Geography of the Dhaka samples"},
      {"p": F.loc["DF09", "finding"]}, fig(43),
      {"h2": "3.8 Integration with the global analysis"},
      {"p": F.loc["DF10", "finding"]}] + ([fig(44),
      {"p": "Because the calibration was stable, the reprocessed municipal-wastewater values were placed descriptively within the distribution of municipal wastewater in the ResFinder resource after dividing by the mean calibration ratio (Table 6). This placement is exploratory: it uses 3-7 samples per city, a single calibration factor estimated from five units, and different sampling years from most of the resource. It must not be read as a ranking of cities by AMR risk."},
      {"table": tab(rd("GI07_dhaka_global_placement"), ["group", "n", "median_calibrated_FPKM", "percentile_in_resource_municipal_wastewater"], ["Group", "n samples", "Median calibrated FPKM", "Percentile within resource municipal wastewater"], [.25, .15, .3, .3]),
       "title": "**Table 6.** Exploratory placement of reprocessed municipal-wastewater samples within the ResFinder-resource municipal-wastewater distribution (source: GI07)."}] if CAL_DONE else []) + [
      {"h1": "4. Interpretation"},
      {"bullets": [F.loc[i, "interpretation"] for i in ["DF03", "DF04", "DF06", "DF07", "DF08", "DF10"]]},
      {"h1": "5. Limitations"},
      {"bullets": [f"**{r.domain}.** {r.limitation} {r.consequence_for_inference}" for r in lim[lim.domain.str.contains("Dhaka|reprocess|Cross-city|Pipeline integration|Taxonomic resolution", case=False)].itertuples()] +
                  ["Detection of ARG DNA in wastewater does not establish phenotypic resistance, expression, viability, bacterial host, transmission to people or clinical risk; no hotspot, risk ranking or food-chain inference is made."]},
      {"h1": "6. Using these results in the chapter"},
      {"p": "The following wording rules follow from the design and should be applied when the Dhaka results are written into Sections 6.1-6.3 of the chapter."},
      {"bullets": ["**Supported:** 'ARG DNA of clinically important carbapenemase families (e.g., NDM, KPC) was detected at high confidence in n of 16 Dhaka wastewater metagenomes' (with the any-hit count only as an unconfirmed upper bound).",
                   "**Supported:** 'Under a matched sewage-surveillance protocol, Dhaka and Karachi sewage did not differ detectably in normalized ARG abundance, but their resistome composition differed' (state n per city and the exact p-values; dispersion was not tested for this pair).",
                   "**Supported with qualification:** Dhaka versus Delhi differences, described as a comparison of different studies, matrices and periods analysed with one pipeline, not as a city effect; the abundance contrast depends on how Dhaka sampling units are defined (Table DH29b) and does not survive multiple-testing adjustment.",
                   "**Not supported:** statements that Dhaka has a higher AMR 'burden', 'load' or 'risk' than other cities; that genes are expressed, functional or carried by viable pathogens; that wastewater AMR is transmitted to people or food; or that any location is an AMR hotspot.",
                   "**Not supported:** generalizing hospital-versus-sewage differences from DS84 (4 vs 3 samples in one study; 35 possible arrangements) to Dhaka as a whole."]},
      {"h1": "7. Traceability"},
      {"p": "Processing logs (fastp reports, Kraken 2 reports, KMA result and mapstat files) for every run, the ENA manifest with URLs, file sizes and MD5 checksums, all notebooks (executed from fresh kernels) and all tables are stored in Updated_Research_Analysis. Prepared with AI assistance (Claude) for author review; the authors should verify interpretations before use in the chapter."}]
json.dump({"title": "Dhaka City Metagenomic Microbiome and AMR Analysis: Comparative Report", "subtitle": "Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping (Academic Press, Elsevier)",
           "meta": ["Report date: 31 October 2026", "Deliverables: Dhaka_City_Advanced_Analysis.xlsx; AMR_Updated_Final_Analysis.xlsx; notebooks 12-18", "Status: analysis complete for all 16 Dhaka samples; findings require author verification before publication"],
           "blocks": B}, open(O / "dhaka_report.json", "w"), indent=0)

# ------------------------------------------------------------------ Final verification report
V5 = pd.read_csv(T.parent / "verification" / "V05_updated_independent_checks.csv") if (T.parent / "verification" / "V05_updated_independent_checks.csv").exists() else pd.DataFrame()
V6 = pd.read_csv(T.parent / "verification" / "V06_workbook_checks.csv") if (T.parent / "verification" / "V06_workbook_checks.csv").exists() else pd.DataFrame()
nb = pd.read_csv(O / "notebook_runs.csv") if (O / "notebook_runs.csv").exists() else pd.DataFrame()
V7 = pd.read_csv(T.parent / "verification" / "V07_raw_read_reproducibility_check.csv") if (T.parent / "verification" / "V07_raw_read_reproducibility_check.csv").exists() else pd.DataFrame()
V8 = pd.read_csv(T.parent / "verification" / "V08_independent_agent_checks.csv") if (T.parent / "verification" / "V08_independent_agent_checks.csv").exists() else pd.DataFrame()
V8_TXT = (f"A separate verification agent, given only the raw per-run tool outputs, the sample metadata, the stated metric definitions and the final tables (not the analysis code), recomputed key Dhaka results with its own code. Its first pass (53 checks) found 10 disagreements, including the clinical-determinant matching error described in Section 3.3; after correction a second pass (79 checks) found 2 remaining minor issues (a missing OXA-48-like allele and one wording problem), which were also corrected. The final pass on the delivered tables gives: " + ", ".join(f"{k} {v}" for k, v in V8.verdict.value_counts().items()) + " (Table 4e)." ) if len(V8) else "The separate-agent recalculation was not available when this report was generated."
INTEG = json.load(open(O / "integrity.json"))
CORRECTIONS = json.load(open(O / "corrections.json")) if (O / "corrections.json").exists() else []
AUTHOR_DECISIONS = json.load(open(O / "author_decisions.json")) if (O / "author_decisions.json").exists() else []
R6_ = rd("R06u_updated_figure_catalog"); a2 = rd("AC02_aud11_pseudoreplication_before_after"); a4 = rd("AC04_aud17_model_leakage_before_after"); a5 = rd("AC05_aud15_figure_audit")
B = [{"h1": "1. Verdict"},
     {"p": f"All open issues of the Independent Research Audit were addressed: {int(AC.status.str.startswith('RESOLVED').sum())} of {len(AC)} are resolved and {int((~AC.status.str.startswith('RESOLVED')).sum())} partially resolved (see Table 1). The 16 Dhaka samples now have verified taxonomic and AMR profiles. {int((V5.status == 'VERIFIED').sum()) if len(V5) else 0} of {len(V5)} independent recalculations agreed with the reported values, {int((V4.result == 'PASS').sum())} of {len(V4)} internal consistency checks passed, and the workbook checks are summarised in Table 4. Integrity of the original data: {INTEG['files_checked']} original files were re-verified against their pre-update SHA-256 checksums at the end of the update; {INTEG['changed_or_missing']} had changed or were missing."},
     {"h1": "2. Audit issues: correction and re-verification"},
     {"table": tab(AC, ["issue", "status", "before", "after", "verification_evidence"], ["Issue", "Status", "Before", "After", "Evidence"], [.08, .16, .26, .28, .22], size=15), "title": "**Table 1.** Resolution of the open audit issues (source: AC00)."},
     {"h2": "2.1 Pseudo-replication (AUD-11)"},
     {"p": "The 19 DS02 open-drain sites were identified from the ENA/SRA experiment titles of PRJNA1415085. City contrasts were re-estimated with the site as the unit: Kruskal-Wallis on site medians with a whole-site permutation p-value, linear mixed models with a site random intercept (maximum-likelihood ratio tests), PERMANOVA permuting whole sites, and chi-square statistics with whole-site permutation. Sample-level p-values are retained only as descriptive values."},
     {"table": tab(a2, ["analysis", "previous_p", "updated_p_site_KW", "updated_p_mixed_LRT", "conclusion_changed_at_0.05"], ["Analysis", "Previous p (447 samples)", "Site-level p", "Mixed-model LRT p", "Conclusion changed"], [.36, .16, .16, .16, .16], size=16),
      "title": "**Table 2.** DS02 city contrasts before and after the AUD-11 correction (source: AC02)."},
     {"h2": "2.2 Predictive-model leakage (AUD-17)"},
     {"table": tab(a4, ["model", "validation", "roc_auc_mean_previous", "roc_auc_mean_fold_internal", "delta_roc_auc"], ["Model", "Validation", "ROC-AUC previous", "ROC-AUC fold-internal", "Change"], [.2, .4, .13, .14, .13], size=16),
      "title": "**Table 3.** Model A performance with transductive versus fold-internal zero replacement (source: AC04)."},
     {"h2": "2.3 Figures (AUD-15, AUD-20)"},
     {"p": f"All {len(a5)} global figures were re-laid out at Elsevier double-column width (maximum {a5.width_mm.max():.0f} mm), with sans-serif lettering of at least {R6_.min_font_pt.min():.1f} pt at print size, and exported as RGB PNG (300 dpi) and LZW TIFF (600 dpi). The updated catalogue lists {int((R6_.file != '(not produced)').sum())} figures (1-30 global; 31-{int(R6_.number.max())} Dhaka, integration and sensitivity). Four Dhaka-stage figures (31, 41, 43, 45) were found below the 5.5 pt floor in this verification and were corrected (lettering raised to 6 pt; overlapping legend in Figure 41B moved)."},
     {"h1": "3. Independent recalculation"},
     {"table": tab(V5, ["area", "check", "reported", "independent", "status"], ["Area", "Check", "Reported", "Independent", "Status"], [.14, .38, .17, .19, .12], size=15) if len(V5) else {"columns": ["note"], "rows": [["not available"]]},
      "title": "**Table 4a.** Independent recalculation of critical results (verify_updated.py; source: V05)."},
     {"table": tab(V6, ["workbook", "check", "value", "status"], ["Workbook", "Check", "Value", "Status"], [.32, .4, .14, .14], size=16) if len(V6) else {"columns": ["note"], "rows": [["not available"]]},
      "title": "**Table 4b.** Workbook verification (source: V06)."},
     {"table": tab(V4, ["check", "result", "detail"], ["Check", "Result", "Detail"], [.6, .12, .28], size=16), "title": "**Table 4c.** Internal consistency checks computed in notebook 20 (source: V04)."},
     {"h2": "3.1 Raw-read completeness and end-to-end reproducibility"},
     {"p": f"All {N_RUNS} selected runs (16 Dhaka, {len(_c)} Global Sewage Surveillance comparators, {N_DELHI} Delhi DS02 samples and the calibration runs) were processed at exactly 4,000,000 read pairs. Ten runs that were incomplete at the start of this verification (all calibration runs and one run each from Karachi, Kochi, Colombo and Delhi) were completed with the same pinned scripts; the second pipeline version differs from the first only by a download-retry loop. Re-collecting all per-run outputs reproduced the previously collected 41 runs exactly (taxonomy, Bracken and ResFinder tables identical row for row). One Dhaka run (SRR39037538) was then re-downloaded and re-processed end to end in a separate folder; Table 4d compares the outputs."},
     {"table": tab(V7, ["file", "result"], ["Output file", "Re-run versus original"], [.4, .6], size=17) if len(V7) else {"columns": ["note"], "rows": [["not available"]]},
      "title": "**Table 4d.** End-to-end reproducibility of one Dhaka run (source: V07)."},
     {"h2": "3.2 Recalculation by a separate verification agent"},
     {"p": V8_TXT},
     {"table": tab(V8.assign(after=V8.status_after_final_correction.fillna("").str.slice(0, 160)) if "status_after_final_correction" in V8.columns else V8.assign(after=""), ["check_id", "description", "reported_value", "independent_value", "verdict", "after"], ["ID", "Check", "Reported", "Independent", "Agent verdict", "After final correction"], [.07, .27, .17, .17, .1, .22], size=14) if len(V8) else {"columns": ["note"], "rows": [["not available"]]},
      "title": "**Table 4e.** Recalculation from raw per-run outputs by a separate agent without access to the analysis code (source: V08)."},
     {"h2": "3.3 Problems found and corrected during this verification"},
     {"bullets": CORRECTIONS},
     {"h1": "4. Notebook execution"},
     {"table": tab(nb, ["notebook", "status", "seconds", "finished"], ["Notebook", "Status", "Seconds", "Finished"], [.5, .15, .15, .2], size=16) if len(nb) else {"columns": ["note"], "rows": [["see logs/notebook_run_log.jsonl"]]},
      "title": "**Table 5.** Final fresh-kernel execution of notebooks 01-20 on the local Python 3.12.10 environment."},
     {"h1": "5. Remaining limitations"},
     {"bullets": [f"**{r.domain}.** {r.limitation}" for r in lim.itertuples()]},
     {"h1": "6. Decisions required from the authors"},
     {"bullets": AUTHOR_DECISIONS},
     {"note": "Prepared with AI assistance (Claude) as a verification record for author review. Original data and previous outputs were not modified; SHA-256 checksums of the original files are stored in 00_Backup_and_Checksums."}]
json.dump({"title": "Updated Final Verification Report", "subtitle": "Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping",
           "meta": ["Verification date: 31 October 2026", "Scope: audit corrections (AUD-04, AUD-11 to AUD-20), Dhaka reprocessing (16 samples), comparative analyses, workbooks, notebooks and figures"],
           "blocks": B}, open(O / "verification_report.json", "w"), indent=0)
print("specs written")
