# %% [markdown]
# # 10 · Final Results Synthesis and Visualisation
# Compiles the outputs of notebooks 01–09 into objective-level results: the objective-to-method framework (with actual findings and the degree to which each objective was addressed), the consolidated register of statistical tests (with multiple-testing adjustment within pre-defined families), the findings register used for the Results section, limitations, the exclusion log, the figure and table catalogues, and a summary dashboard. No new statistics are estimated here except the family-wise Benjamini–Hochberg adjustment of the test register.

# %%
#BOOTSTRAP
NB = "NB10"
nb_header(NB)
import textwrap, json
figs = load_catalog("figures").sort_values("number"); tabs = load_catalog("tables"); tests = load_catalog("tests"); finds = load_catalog("findings"); excl = load_catalog("exclusions")
_g = [f"NB{i:02d}" for i in range(1, 11)]  # UPDATE: NB10 compiles the global notebooks 01-10 only; NB20 compiles the full updated set
figs = figs[figs.notebook.isin(_g)]; tabs = tabs[tabs.notebook.isin(_g)]; tests = tests[tests.notebook.isin(_g)]; finds = finds[finds.notebook.isin(_g)]; excl = excl[excl.notebook.isin(_g)]
print(len(figs), "figures;", len(tabs), "tables;", len(tests), "tests;", len(finds), "findings;", len(excl), "exclusion records")
(LOG / "environment.json").write_text(json.dumps(env_info(), indent=1), encoding="utf-8")

# %% [markdown]
# ## 1. Consolidated statistical test register

# %%
T = tests.copy()
T["p_value"] = pd.to_numeric(T.p_value, errors="coerce")
T["p_adj_BH_within_family"] = np.nan
for fam_, ix in T.groupby("family").groups.items():
    T.loc[ix, "p_adj_BH_within_family"] = bh(T.loc[ix, "p_value"].values)
T["significant_at_0.05_after_adjustment"] = np.where(T.p_adj_BH_within_family.isna(), "n/a", np.where(T.p_adj_BH_within_family < .05, "Yes", "No"))
cols = ["notebook", "objective", "family", "analysis", "test", "n", "statistic", "df", "p_value", "p_adj_BH_within_family", "significant_at_0.05_after_adjustment", "effect_size_name", "effect_size", "ci_low", "ci_high", "note"]
T = T[cols].sort_values(["notebook", "family", "analysis"])
save_table(T, "R02_statistical_tests_register", "Register of all statistical tests and model metrics (BH adjustment within test family)", NB, sheet="Statistical_Tests")
T.groupby(["notebook", "family"]).size().rename("tests").reset_index()

# %% [markdown]
# ## 2. Findings register (Results section), exclusion log and catalogues

# %%
F = finds.sort_values("id", key=lambda s: s.str.extract(r"(\d+)")[0].astype(int) + s.str.contains("b").astype(float) / 10)
F = F[["id", "objective", "evidence_type", "finding", "sample_size", "denominator", "statistical_evidence", "table_or_figure", "interpretation", "limitations", "notebook"]]
save_table(F, "R03_findings_register", "Verified findings by objective: result, denominator, statistical evidence, interpretation and limitations", NB, sheet="Results_Interpretation")
E = excl[["notebook", "analysis", "unit", "n_excluded", "reason", "example_ids"]]
save_table(E, "R05_exclusion_log", "Exclusion and data-handling log (all notebooks)", NB, sheet="Data_Quality_Audit")
F[["id", "objective", "finding"]]

# %% [markdown]
# ## 3. Objective-to-method framework with actual findings and status

# %%
def fids(*ids): return " | ".join(f"[{i}] " + F.set_index("id").loc[i, "finding"] for i in ids if i in set(F.id))
obj = pd.DataFrame([
 dict(objective="1. Urban microbial diversity", research_question="How do microbial composition and diversity vary across urban environments, cities and seasons?",
      dataset="DS02 Bracken phylum table (447 samples, 4 Indian cities); Martiny et al. SILVA rRNA diversity indices (551 units, 6 matrices)",
      variables="Prokaryotic phylum relative abundance; observed/rarefied richness, Shannon, Gini-Simpson, Pielou; genus Shannon/richness; city, season, year, matrix, study",
      methods="Relative abundance; alpha diversity; Kruskal-Wallis + Holm post hoc + epsilon-squared; OLS (HC3); linear mixed model (study random effect); Bray-Curtis/Aitchison PERMANOVA + PERMDISP; CLR-PCA",
      software="Python 3.12: pandas, numpy, scipy, statsmodels, scikit-bio, scikit-learn, matplotlib, seaborn",
      outputs="Tables T00-T12, C01, P01-P03c; Figures 5-8, 19", findings=fids("F04", "F05", "F06", "F10", "F13", "F13b"),
      limitations="Phylum resolution only for DS02; genus composition tables not downloaded; one study per Indian city set; matrix confounded with study; city inference rests on 19 DS02 sites",
      status="Partially addressed"),
 dict(objective="2. AMR patterns and microbiome associations", research_question="Which ARGs, classes and mechanisms occur, how abundant and diverse are they, and how do they relate to community characteristics?",
      dataset="CARD/RGI gene table (DS02); ResFinder/KMA templates (Martiny et al.) with rRNA-fragment denominators",
      variables="RPKM / FPKM per gene family, class group, mechanism; richness; Shannon; clinically important determinants; phylum CLR; genus diversity",
      methods="Length- and depth-normalised abundance; Wilson CIs; chi-square; Kruskal-Wallis; OLS; mixed models; PERMANOVA; Mantel; Procrustes; Spearman (BH); co-occurrence network; random-effects heterogeneity",
      software="As above + networkx", outputs="Tables A01-A15, C01-C15, P04-P08; Figures 9-18, 20-21", findings=fids("F07", "F08", "F09", "F11", "F12", "F14"),
      limitations="Two non-equivalent pipelines analysed separately; gene detection is not phenotype, expression, viability or host; CARD release unstated; detection depends on depth and thresholds",
      status="Addressed within the two processed collections"),
 dict(objective="3. Environmental and public-health implications", research_question="Which environmental factors are associated with AMR patterns, and how far does the evidence support public-health concern?",
      dataset="DS02 (season, year, city); ResFinder units (matrix, continent); MetaSUB city covariates (no AMR data)",
      variables="Season, year, city, matrix, continent; detection of clinically important determinants",
      methods="OLS with season/year; PERMANOVA season | city; monthly series; chi-square for determinant detection by city",
      software="As above", outputs="Tables A05-A11, C01; Figures 11-12", findings=fids("F08", "F10"),
      limitations="No pollutant, antibiotic-residue, sanitation, population or clinical covariates linked to AMR samples; public-health relevance limited to determinant presence",
      status="Partially addressed"),
 dict(objective="4. Predictive and spatial assessment", research_question="Do the data support predictive modelling, spatial pattern identification and forecasting?",
      dataset="ResFinder units (matrix labels); DS02 taxonomy + normalized ARG abundance; 2016 sewage pilot coordinates; all included-sample coordinates",
      variables="CLR class composition; CLR phylum composition; log10 normalized ARG abundance; coordinates; continent",
      methods="Leave-one-study-out and temporal/leave-one-city-out validation; RF, logistic, ridge, gradient boosting vs baselines; Moran's I; Mantel distance decay; observed-value maps",
      software="scikit-learn, geopandas (Natural Earth 1:50m)", outputs="Tables M01-M07, S01-S07; Figures 22-28", findings=fids("F15", "F16", "F17", "F18"),
      limitations="No outcome data for disease/outbreak/food-business prediction; spatial inference limited to one 2016 round; no hotspot or interpolation analysis justified",
      status="Partially addressed (exploratory models and observed-value spatial analysis; forecasting not feasible)"),
 dict(objective="5. Food systems and policy implications", research_question="What does the available evidence say about environmental AMR and food safety/food systems?",
      dataset="Food-related contextual units (155; 4 studies) in the ResFinder resource; FDA retail lamb (DS05; metadata only)",
      variables="Total ARG FPKM, class composition, clinically important determinants",
      methods="Descriptive comparison; mixed model food-related vs municipal wastewater (contextual); ordination",
      software="As above", outputs="Tables A12-A15, C10, SA03; Figures 13-14, 21", findings=fids("F09"),
      limitations="Food-related samples are a separate contextual dataset with different designs; no food-chain transmission, contamination, recall or business data; no causal or economic claims possible",
      status="Partially addressed (contextual comparison only); policy implications are interpretive"),
])
save_table(obj, "R01_objective_method_framework", "Objective-to-method analytical framework with actual findings and status", NB, sheet="Method_Selection")
ro = obj[["objective", "research_question", "status"]].copy()
save_table(ro, "R00_research_objectives", "Research objectives and the degree to which the available data addressed them", NB, sheet="Research_Objectives")
obj[["objective", "status"]]

# %% [markdown]
# ## 4. Limitations register

# %%
lim = pd.DataFrame([
 ("Data acquisition", "Global analyses (notebooks 01-10) use producer-processed tables (two pipelines); 51 public runs (16 Dhaka, 29 South Asian comparators, 6 calibration runs) were uniformly reprocessed from raw reads at 4 M read pairs (notebooks 12-18).", "Cross-study comparison of abundance is valid only within one pipeline; the reprocessed set is a third, separate analytical group.", "Uniform reprocessing of a larger, approved raw-read subset."),
 ("Pipeline heterogeneity", "CARD/RGI (DS02) and ResFinder/KMA (Martiny et al.) differ in database, alignment and count definition.", "Results are reported side by side, never pooled.", "Re-annotate both collections with one database/version."),
 ("Coverage", f"Producer-processed AMR tables cover 1,072 of 8,053 included samples; taxonomy covers 447. Dhaka (16) and Karachi (7) now have uniformly reprocessed profiles (4 M read pairs per run).", "MetaSUB transit surfaces and most other cities still contribute metadata only.", "Reprocess further cities and matrices with the same pipeline."),
 ("Confounding", "Matrix is confounded with study, KMA version and sequencing depth in the ResFinder set (ICC ~0.65).", "Matrix contrasts are study-limited.", "Multi-matrix studies with common protocols."),
 ("Pseudo-replication (resolved, AUD-11)", "DS02 site identities (19 sites) were recovered from ENA/SRA experiment titles; city inference now uses site-level permutation tests and site random-intercept models.", "City contrasts rest on 19 sites (Kolkata 2); sample-level p-values are descriptive only.", "Confirm site labels against Singh et al. Supplementary Data S7 when accessible."),
 ("Unexplained structure", "DS02 taxonomic profiles form discrete clusters not explained by available metadata.", "Possible site or batch effect; city contrasts re-estimated with cluster adjustment.", "Obtain site and processing-batch metadata."),
 ("Taxonomy", "DS02 Kraken 2 core_nt profile contains ~31% eukaryotic reads; only phylum level available.", "Analyses restricted to prokaryotic phyla; low resolution.", "Genus/species profiling of reprocessed reads."),
 ("Detection thresholds (AUD-12)", "Detection frequency of clinically important genes depends strongly on the threshold; identity and MAPQ were not deposited for DS02.", "High-confidence detection is the primary result; any-hit detection is an unconfirmed upper bound.", "Assembly-based or variant-aware detection; deposit identity/MAPQ."),
 ("Sequencing depth (AUD-16)", "231 included samples have < 10 Mbp and 642 < 100 Mbp; they enter descriptive denominators only.", "Depth-threshold sensitivity analysis (notebook 19) shows which results change.", "Set a minimum depth for future reprocessing."),
 ("Interpretation", "ARG detection does not establish phenotype, expression, viability, host, transmission or clinical infection.", "No public-health burden, transmission or food-chain claims are made.", "Linked clinical/food surveillance and culture-based confirmation."),
 ("Outcome data", "No matched clinical, epidemiological, food-safety or economic outcomes.", "Disease, outbreak and food-business prediction not feasible.", "Link to GLASS/NARMS-type indicators with matched place and time."),
 ("Spatial", "Coordinates for 76% of included samples; spatial inference limited to the 2016 pilot (one round).", "No hotspots, interpolation or risk maps.", "Denser, repeated sampling with coordinates."),
 ("Duplicates", "PRJEB13831 likely duplicates the 2016 pilot physical samples; PRJEB27054 accessions are technical replicates.", "Handled by exclusion/aggregation; sensitivity analysis provided.", "Confirm sample identity with the data producers."),
 ("Pending eligibility", "25 studies await author eligibility decisions and are not analysed.", "Coverage may change after decisions.", "Author review of curation_decisions.py."),
], columns=["domain", "limitation", "consequence_for_inference", "remedy_or_future_work"])
save_table(lim, "R04_limitations", "Limitations register", NB, sheet="Limitations")

# %% [markdown]
# ## 5. Summary dashboard

# %%
D02 = pd.read_csv(TAB / "D02_by_sample_type.csv"); a4 = pd.read_csv(TAB / "A04_ds02_arg_load_by_city.csv"); a12 = pd.read_csv(TAB / "A12_rf_arg_load_by_sample_type.csv")
ov = pd.read_csv(TAB / "D01_evidence_base_overview.csv").set_index("Quantity").iloc[:, 0]
m2s = pd.read_csv(TAB / "M02_modelA_performance.csv"); auc = m2s[(m2s.model == "Random forest") & m2s.validation.str.startswith("Stratified group")].roc_auc_mean.iloc[0]
s6 = pd.read_csv(TAB / "S06_spatial_autocorrelation.csv"); mi = s6.iloc[0]
fig = plt.figure(figsize=(6.9, 6.4)); gs = fig.add_gridspec(3, 4, height_ratios=[.55, 1, 1], hspace=.55, wspace=.35)
kpis = [(f"{int(ov['Samples']):,}", "included urban\nmetagenomes"), (f"{int(ov['Countries (reported)'])}", "countries"), (f"{int(ov['Source-labelled cities (City|Country)'])}", "source-labelled\ncities"),
        (f"{int(ov['Samples with any AMR gene table']):,}", "samples with a\nprocessed AMR table")]
for i, (v, l) in enumerate(kpis):
    ax = fig.add_subplot(gs[0, i]); ax.axis("off"); ax.text(.5, .62, v, ha="center", va="center", fontsize=15, weight="bold", color="#1F4E79"); ax.text(.5, .12, l, ha="center", va="center", fontsize=7.5)
ax = fig.add_subplot(gs[1, :2]); d = D02.sort_values("samples")
ax.barh(d.Sample_Type, d.samples, color=[SAMPLE_TYPE_COLORS.get(t) for t in d.Sample_Type]); ax.barh(d.Sample_Type, d.with_AMR, color="k", alpha=.35, height=.4, label="with AMR table")
ax.set_xlabel("Included samples"); ax.set_title("Samples by matrix (dark = with processed AMR table)", fontsize=7.5); ax.legend(frameon=False, fontsize=6.5)
ax = fig.add_subplot(gs[1, 2:]); ax.bar(a4.City, a4.median_total_RPKM, color=[CITY_COLORS[c] for c in a4.City])
ax.errorbar(a4.City, a4.median_total_RPKM, yerr=[a4.median_total_RPKM - a4.Q1, a4.Q3 - a4.median_total_RPKM], fmt="none", color="k", lw=.8, capsize=3)
ax.set_ylabel("Median total ARG RPKM (IQR)"); ax.set_title("Indian open-drain wastewater: normalized ARG abundance by city (CARD/RGI, n = 447)", fontsize=7.5)
ax = fig.add_subplot(gs[2, :2]); d = a12.sort_values("median_total_FPKM")
# AUDIT FIX AUD-08: bars on a log axis exaggerate differences; show medians as points with IQR on the log axis
ax.errorbar(d.median_total_FPKM, d.Sample_Type, xerr=[d.median_total_FPKM - d.Q1, d.Q3 - d.median_total_FPKM], fmt="none", ecolor="k", lw=.8, capsize=3)
ax.scatter(d.median_total_FPKM, d.Sample_Type, s=45, zorder=3, color=[SAMPLE_TYPE_COLORS.get(t, "#882255") for t in d.Sample_Type], edgecolors="k", linewidths=.4); ax.set_xscale("log")
ax.set_xlabel("Median total ARG FPKM (IQR), log scale"); ax.set_title("normalized ARG abundance by matrix (ResFinder/KMA; food-related = contextual)", fontsize=7.5)
ax = fig.add_subplot(gs[2, 2:]); ax.axis("off")
txt = [f"• Wastewater-vs-other-matrix resistome signature: mean fold ROC-AUC {auc:.2f} with whole studies held out (grouped CV).",
       f"• 2016 sewage pilot: Moran's I = {mi.morans_I:.2f} (p = {mi.p_value_two_sided:.3f}) for normalized ARG abundance across {int(mi.n)} pilot units (mostly one per city).",
       "• Two AMR pipelines analysed separately; detection ≠ phenotype or clinical risk.",
       "• No outcome data: disease, outbreak and food-business forecasting not feasible.",
       "• See Results_Interpretation and Limitations sheets for denominators and caveats."]
ax.text(0, 1, "\n".join(textwrap.fill(t, 58) for t in txt), va="top", fontsize=6.0, linespacing=1.35); ax.set_title("Key verified messages", fontsize=7.5, loc="left")
fig.suptitle("Urban microbiome and AMR: verified analysis dashboard (updated)", fontsize=8, weight="bold")
save_fig(fig, 30, "summary_dashboard", "Summary dashboard of verified findings",
         "Key counts and results from notebooks 02, 04, 07 and 08. All values are traceable to the tables named on each panel's source sheet; see the findings register for denominators and limitations.",
         NB, "All objectives", "tables/D01, D02, A04, A12, M02, S06")

# %%
# AUD-20 fix: the figure catalogue is written after Figure 30 so that it is complete
figs = load_catalog("figures").sort_values("number"); figs = figs[figs.number <= 30]  # global figure set (Dhaka figures 31+ are catalogued in NB20)
fcols = [c for c in ["number", "file", "tiff_file", "title", "caption", "objective", "notebook", "source_data", "width_mm", "height_mm", "print_layout", "min_font_pt", "text_overlaps_flagged"] if c in figs.columns]
save_table(figs[fcols], "R06_figure_catalog", f"Figure catalogue ({len(figs)} figures; final print size and legibility audit)", NB, sheet="Figure_Gallery")
assert len(figs) == figs.number.nunique() == 30, len(figs)
tabs2 = load_catalog("tables"); tabs2 = tabs2[tabs2.notebook.isin(_g)]
save_table(tabs2[["file", "title", "sheet", "notebook", "rows", "cols", "description"]].sort_values("file"), "R07_table_catalog", "Catalogue of analytical tables (CSV) and their workbook sheets", NB, sheet="Verification_Report")
summary_lines = [
 f"Analysis date: {ANALYSIS_DATE}. Inputs: harmonized public-data collection (read only). Software: Python {env_info()['python']} with pandas, scipy, statsmodels, scikit-learn, scikit-bio, geopandas.",
 f"Evidence base: {int(ov['Samples']):,} included urban shotgun metagenomes ({int(ov['Datasets'])} datasets, {int(ov['Countries (reported)'])} countries, {int(ov['Source-labelled cities (City|Country)'])} source-labelled cities) plus 840 food-related contextual samples; processed AMR tables exist for {int(ov['Samples with any AMR gene table']):,} included samples (447 CARD/RGI, 625 ResFinder/KMA) and taxonomy for 447.",
 "Two processed collections support inference: DS02 Indian open-drain wastewater (Delhi, Mumbai, Chennai, Kolkata; 2022-2024) and the Martiny et al. global resistome resource (ResFinder/KMA; six matrices).",
] + [f"[{r.id}] {r.finding}" for _, r in F.iterrows()]
save_table(pd.DataFrame({"line": range(1, len(summary_lines) + 1), "summary": summary_lines}), "R08_executive_summary", "Executive summary of verified findings", NB, sheet="Executive_Summary")
print("\n".join(summary_lines[:6]))
