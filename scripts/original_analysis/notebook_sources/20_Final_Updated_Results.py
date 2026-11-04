# %% [markdown]
# # 20 · Final updated results: findings register, catalogues and verification checks
# **Research objective:** compile the corrected global analysis (notebooks 01-10), the audit corrections (11), the Dhaka reprocessing and comparisons (12-17), the integration and sensitivity analyses (18-19) into one traceable results base for the chapter, organised by the chapter's results sections 6.1-6.3 and the Dhaka comparative subsection.
#
# Every finding below is traceable to a table or figure produced in this run; the register states its denominator, statistical evidence, interpretation and limitation. Food-system implications (6.3) are limited to contextual evidence: no food-chain transmission or business impact is inferred from environmental metagenomes.

# %%
#BOOTSTRAP
NB = "NB20"
nb_header(NB)
import textwrap
figs = load_catalog("figures").sort_values("number"); tabs = load_catalog("tables"); tests = load_catalog("tests"); finds = load_catalog("findings")
print(len(figs), "figures;", len(tabs), "tables;", len(tests), "tests;", len(finds), "findings")

# %%
SECTION = {"F01": "Evidence base", "F02": "Evidence base", "F03": "Evidence base", "F04": "6.1 Microbiome environmental hazards", "F05": "6.1 Microbiome environmental hazards", "F06": "6.1 Microbiome environmental hazards",
           "F10": "6.1 Microbiome environmental hazards", "F11": "6.2 AMR and public health in global cities", "F13": "6.1 Microbiome environmental hazards", "F13b": "Data quality",
           "F07": "6.2 AMR and public health in global cities", "F08": "6.2 AMR and public health in global cities", "F09": "6.3 AMR and food systems (contextual)", "F12": "6.2 AMR and public health in global cities",
           "F14": "6.2 AMR and public health in global cities", "F15": "Spatial and predictive assessment", "F16": "Spatial and predictive assessment", "F17": "Spatial and predictive assessment",
           "F18": "Spatial and predictive assessment", "F19": "Robustness", "F20": "Robustness"}
F = finds.copy()
F["results_section"] = [SECTION.get(i, "Dhaka comparative findings" if str(i).startswith("DF") else "Other") for i in F.id]
F.loc[F.id.isin(["DF02", "DF05"]), "results_section"] = "6.1 Microbiome environmental hazards (Dhaka)"
F.loc[F.id.isin(["DF03", "DF04"]), "results_section"] = "6.2 AMR and public health (Dhaka)"
F.loc[F.id.isin(["DF06", "DF07", "DF08"]), "results_section"] = "Dhaka comparative findings"
F.loc[F.id.isin(["DF01", "DF09", "DF10"]), "results_section"] = "Dhaka data usability, geography and integration"
key = lambda s: s.str.extract(r"(\d+)")[0].astype(int) + s.str.startswith("DF").astype(int) * 100 + s.str.contains("b").astype(float) / 10
F = F.sort_values("id", key=key)
F = F[["id", "results_section", "objective", "evidence_type", "finding", "sample_size", "denominator", "statistical_evidence", "table_or_figure", "interpretation", "limitations", "notebook"]]
save_table(F, "R03u_updated_findings_register", f"UPDATED findings register ({len(F)} findings): corrected global findings F01-F20 and new Dhaka findings DF01-DF10, organised by results section", NB, sheet="Results_Interpretation")
assert F.id.is_unique and F.finding.str.len().min() > 40
F.groupby("results_section").id.apply(list)

# %%
T = tests.copy(); T["p_value"] = pd.to_numeric(T.p_value, errors="coerce"); T["p_adj_BH_within_family"] = np.nan
for fam_, ix in T.groupby("family").groups.items(): T.loc[ix, "p_adj_BH_within_family"] = bh(T.loc[ix, "p_value"].values)
T["inferential_status"] = np.where(T.family.str.contains(r"descriptive", case=False), "Descriptive only (pseudo-replicated or confounded)", "Inferential")
cols = ["notebook", "objective", "family", "inferential_status", "analysis", "test", "n", "statistic", "df", "p_value", "p_adj_BH_within_family", "effect_size_name", "effect_size", "ci_low", "ci_high", "note"]
save_table(T[cols].sort_values(["notebook", "family", "analysis"]), "R02u_updated_statistical_tests_register", f"UPDATED register of all {len(T)} statistical tests and model metrics (BH within family; inferential status flagged)", NB, sheet="Statistical_Tests")
fc = [c for c in ["number", "file", "tiff_file", "title", "caption", "objective", "notebook", "source_data", "width_mm", "height_mm", "print_layout", "min_font_pt", "text_overlaps_flagged"] if c in figs.columns]
# Figure numbers are analysis-internal. A number is RESERVED (not produced) only when its input data do not exist; it is listed explicitly, never silently skipped.
RESERVED = {44: "Calibration of the uniform pipeline against the ResFinder resource (NB18) - NOT PRODUCED: none of the 6 calibration runs was processed"}
missing = sorted(set(range(1, int(figs.number.max()) + 1)) - set(figs.number.astype(int)))
reserved_rows = pd.DataFrame([dict(number=n, file="(not produced)", tiff_file="(not produced)", title=RESERVED.get(n, "UNDOCUMENTED GAP"), caption="", notebook="NB18" if n == 44 else "") for n in missing])
figcat = pd.concat([figs[fc], reserved_rows[[c for c in reserved_rows.columns if c in fc]]], ignore_index=True).sort_values("number")
save_table(figcat, "R06u_updated_figure_catalog", f"UPDATED figure catalogue: {len(figs)} figures produced (1-30 global; 31-{int(figs.number.max())} Dhaka, integration and sensitivity); reserved/not produced: {', '.join(str(n) for n in missing) or 'none'}", NB, sheet="Figure_Gallery")
assert figs.number.is_unique and all(n in RESERVED for n in missing), f"undocumented figure-number gaps: {missing}"
assert all((FIG / f).exists() for f in figs.file) and all((FIG / f).exists() for f in figs.tiff_file)
print(figs[["number", "width_mm", "min_font_pt", "text_overlaps_flagged"]].describe().round(1))

# %%
lim = pd.read_csv(TAB / "R04_limitations.csv")
_q = pd.read_csv(DER / "reprocessed" / "run_qc.csv")
_nd = _q[_q.processing_status != "completed"]
N_SEL, N_NOT_DONE = len(_q), len(_nd)
NOT_DONE_TXT = "; ".join(f"{r.run_accession} ({r.analysis_group}, {r.country}, {r.sample_alias})" for r in _nd.itertuples())
_cal_done = int(((_q.processing_status == "completed") & _q.analysis_group.str.startswith("Calibration")).sum())
CAL_NOTE = (f"Calibration on {_cal_done} pilot runs (NB18)." if _cal_done else "Calibration against the ResFinder resource was NOT performed (0 of 6 calibration runs processed; NB18, table GI01).")
extra = pd.DataFrame([
 ("Dhaka reprocessing depth", "Dhaka and comparator runs were analysed at 4,000,000 read pairs (head-of-file subsample), 4-100% of each library.", "Rare genes may be missed; equal depth makes runs comparable.", "Full-depth reprocessing on a larger compute allocation."),
 ("Dhaka sample size and design", "16 Dhaka samples from 3 studies, 2017-2025, 3 instruments; sample type, study, year and site are partly confounded; DS84 is the only within-study type contrast (4 vs 3).", "Within-Dhaka contrasts are descriptive or underpowered (DS84 type contrast: exact p = 0.057 with 35 possible arrangements).", "A designed Dhaka survey with replicated sites per matrix."),
 ("Taxonomic resolution (reprocessed)", "Taxonomy from 16S rRNA fragments (Kraken 2 / SILVA 138): genus resolution; species not reported.", "Composition at phylum-genus level only.", "Genome-database profiling (e.g. Kraken 2 standard, MetaPhlAn) with adequate memory."),
 ("Cross-city comparability", "Dhaka vs Delhi compares different studies, matrices (sewer/STP inlet vs open drain) and periods with 4-5 independent units per city.", "Descriptive matched-matrix comparison, not a city effect.", "Concurrent, protocol-matched sampling in both cities."),
 ("Pipeline integration", "The reprocessed set uses newer KMA/ResFinder releases and a Kraken 2 16S denominator. " + (CAL_NOTE), "Kept as a separate analytical group; never pooled with the ResFinder resource or CARD/RGI values; Dhaka is not ranked against global cities.", "Process the 6 calibration runs (or reprocess the full resource with one pipeline)."),
 ("Processing environment", "Raw reads were processed in a temporary Linux cloud workspace (2 CPU cores, 7 GB RAM) with the archived scripts (02_Scripts/pipeline: process_sample.sh; process_sample_v2.sh differs only by a download-retry loop; collect_outputs.py). Per-run tool logs and outputs are archived.", "Reproducible from the archived scripts, pinned tool versions and database checksums (versions.json).", "Re-run on institutional hardware at full depth."),
 ("Karachi comparator provenance", "4 of the 7 Karachi DS01 runs come from ENA project PRJEB51229 (Global Sewage Surveillance 2016-2019 sewage metagenomics; monthly samples July-October 2017), the other 3 from PRJEB27621/PRJEB40815/PRJEB40816; all carry DTU Global Sewage Surveillance sample identifiers.", "Karachi contributes a monthly series from one sampling programme; samples are city-rounds, not independent sites.", "Confirm sampling site(s) from the Global Sewage Surveillance metadata."),
 ("DS01 Colombo site codes", "The three 2017 Colombo DS01 samples carry different site codes (LK_BO, LK_CO, LK_MA) on the same date; the harmonized master assigns all to Colombo.", "Colombo n = 3 may represent three sewage sites in one round rather than three rounds.", "Check Hendriksen et al. / Munk et al. site metadata."),
 ("DS84 collection year", "BioProject PRJNA1475805 states the samples were collected in 2025, whereas the BioSample records give collection dates in 2024 (e.g., 2024-09-17, 2024-11-01).", "Temporal comparisons that use DS84 dates carry a year-level uncertainty.", "Confirm with the submitting institution (icddr,b) or the publication."),
], columns=["domain", "limitation", "consequence_for_inference", "remedy_or_future_work"])
if N_NOT_DONE:
    extra = pd.concat([extra, pd.DataFrame([dict(domain="Unprocessed runs", limitation=f"{N_NOT_DONE} of {N_SEL} selected runs were not processed: {NOT_DONE_TXT}.", consequence_for_inference="Comparator groups are smaller than designed.", remedy_or_future_work="Process the remaining runs with the same pinned pipeline.")])], ignore_index=True)
L = pd.concat([lim, extra], ignore_index=True)
save_table(L, "R04u_updated_limitations", "UPDATED limitations register (global + Dhaka)", NB, sheet="Limitations")

# %%
ov = pd.read_csv(TAB / "D01_evidence_base_overview.csv").set_index("Quantity").iloc[:, 0]
dq = pd.read_csv(TAB / "DH06_dhaka_sequencing_qc.csv"); d16 = pd.read_csv(TAB / "DH16_dhaka_clinical_determinants.csv"); T29 = pd.read_csv(TAB / "DH29_dhaka_comparative_statistical_tests.csv")
summary = [f"{UPDATE_LABEL}. Inputs: harmonized public-data collection with sample-level corrections (AUD-04, AUD-11, AUD-18); producer-processed CARD/RGI (DS02) and ResFinder/KMA (Martiny et al.) collections analysed separately; {len(pd.read_csv(DER / 'reprocessed' / 'run_qc.csv'))} public runs uniformly reprocessed from raw reads.",
           f"Evidence base: {int(ov['Samples']):,} included urban metagenomes from {int(ov['Countries (reported)'])} countries. All 16 Dhaka samples now have taxonomic (16S-fragment, genus-level) and AMR (ResFinder/KMA) profiles from a uniform 4,000,000-read-pair subsample per run; {N_SEL - N_NOT_DONE} of {N_SEL} selected runs (16 Dhaka, comparators and calibration units) were processed."] + [f"[{r.id}] {r.finding}" for r in F.itertuples()]
save_table(pd.DataFrame({"line": range(1, len(summary) + 1), "summary": summary}), "R08u_updated_executive_summary", "UPDATED executive summary of verified findings", NB, sheet="Executive_Summary")
sm_ = load_csv("sample_master.csv"); inc_ = sm_[sm_.Inclusion_Status == "Included"]
qc_ = pd.read_csv(DER / "reprocessed" / "run_qc.csv"); rep_n = qc_[qc_.processing_status == "completed"].merge(sm_[["Sample_ID", "Dataset_ID"]], left_on="sample_accession", right_on="Sample_ID").groupby("Dataset_ID").size()
inv = inc_.groupby(["Dataset_ID", "Dataset_Name"]).agg(included_samples=("Sample_ID", "size"), countries=("Country", lambda s: "; ".join(sorted(set(s) - {"Not reported"}))[:200]),
        sample_types=("Sample_Type", lambda s: "; ".join(f"{k} {v}" for k, v in s.value_counts().items())), sample_type_corrections=("Sample_Type_Correction", lambda s: int(s.astype(str).str.startswith("AUD-04").sum())),
        site_ids_recovered=("Site_ID", lambda s: int(s.nunique()) if (s != "Not reported").any() else 0),
        CARD_RGI_profiles=("AMR_Reference_Database", lambda s: int(s.astype(str).str.contains("CARD").sum())), ResFinder_KMA_profiles=("AMR_Reference_Database", lambda s: int(s.astype(str).str.contains("ResFinder").sum()))).reset_index()
inv["uniformly_reprocessed_runs_this_update"] = inv.Dataset_ID.map(rep_n).fillna(0).astype(int)
save_table(inv.sort_values("included_samples", ascending=False), "R09u_updated_dataset_inventory", f"UPDATED dataset inventory: {len(inv)} datasets with included samples, verified sample types (AUD-04), recovered site identifiers (AUD-11), producer-processed AMR coverage and runs reprocessed in this update", NB, sheet="Dataset_Inventory")
tabs2 = load_catalog("tables")
save_table(tabs2[["file", "title", "sheet", "notebook", "rows", "cols", "description"]].sort_values("file"), "R07u_updated_table_catalog", "UPDATED catalogue of all analytical tables (CSV) and their workbook sheets", NB, sheet="Verification_Report")

# %% [markdown]
# ## Verification checks (computed from the outputs of this run)

# %%
chk = []
def c(name, ok, detail): chk.append(dict(check=name, result="PASS" if ok else "FAIL", detail=detail))
sm = load_csv("sample_master.csv"); inc = sm[sm.Inclusion_Status == "Included"]
c("AUD-04 sample types: hospital wastewater = 237 included", int((inc.Sample_Type == "Hospital wastewater").sum()) == 237, int((inc.Sample_Type == "Hospital wastewater").sum()))
c("AUD-11 DS02 site identifiers present for all 447 samples (19 sites)", (sm[sm.Dataset_ID == "DS02"].Site_ID != "Not reported").all() and sm[sm.Dataset_ID == "DS02"].Site_ID.nunique() == 19, sm[sm.Dataset_ID == "DS02"].Site_ID.nunique())
c("16 Dhaka samples processed with taxonomy and AMR profiles", len(dq) == 16 and (dq.prokaryotic_16S_fragments > 0).all() and (dq.resfinder_aligned_fragments > 0).all(), f"{len(dq)} runs")
c("All Dhaka runs analysed at 4,000,000 read pairs", (dq.acquired_pairs == 4_000_000).all(), dq.acquired_pairs.unique().tolist())
c("Figure catalogue complete: every number 1-max produced or explicitly reserved", figs.number.is_unique and all(n in RESERVED for n in missing), f"{len(figs)} produced; reserved {missing}")
c("All figures <= 190 mm wide (Elsevier double column)", bool((figs.width_mm <= 190).all()), f"max {figs.width_mm.max():.0f} mm")
c("No figure lettering below 5.5 pt", bool((figs.min_font_pt >= 5.5).all()), f"min {figs.min_font_pt.min():.1f} pt")
c("Findings register: unique IDs, every finding names a table or figure", F.id.is_unique and F.table_or_figure.str.len().min() > 3, len(F))
c("No pooled CARD/RGI and ResFinder/KMA values (reprocessed set separate)", True, "by design: separate tables and analytical groups; GI01-GI07 calibration only")
c("Comparative tests state the independent unit", T29.independent_unit.notna().all(), len(T29))
c("Dhaka clinical table covers 11 determinants x 1 overall group", d16[d16.group == "All Dhaka"].determinant.nunique() == 11, "")
CHK = pd.DataFrame(chk)
save_table(CHK, "V04_updated_analysis_checks", "Verification checks of the updated analysis (computed in notebook 20)", NB, sheet="Verification_Report")
assert (CHK.result == "PASS").all(), CHK[CHK.result != "PASS"]
CHK
