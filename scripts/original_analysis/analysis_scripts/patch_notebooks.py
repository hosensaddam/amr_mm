"""Apply the 31 Oct 2026 update patches to notebook_sources of the corrected package.
Each content patch must match exactly once (asserted). A patch log is written as JSON."""
import sys, re, json
from pathlib import Path
SRC = Path(sys.argv[1]) / "notebook_sources"
LOG = []
files = {p.name: p.read_text(encoding="utf-8") for p in sorted(SRC.glob("*.py"))}

def rep(fn, old, new, issue, count=1):
    s = files[fn]; n = s.count(old)
    assert n == count, (fn, n, old[:160])
    files[fn] = s.replace(old, new); LOG.append(dict(issue=issue, file=fn, old=old[:200], new=new[:200]))

# ============================================================== AUD-11: DS02 site identifiers into ds02_meta
rep("01_Data_Inspection_and_QC.py",
    'd2[["Sample_ID", "Run_ID", "Source_Sample_Code", "City", "Collection_Date", "year", "month", "Season", "Sample_Type", "reads_total", "Base_Count", "bracken_phylum_reads"]].to_csv(DER / "ds02_meta.csv", index=False)',
    '# UPDATE (AUD-11): site identifiers recovered from ENA/SRA experiment titles (Site_ID in the updated sample master)\n'
    'd2["Site"] = d2.Site_ID if "Site_ID" in d2.columns else "Not reported"\n'
    'print("DS02 sites:", d2.Site.nunique(), d2.groupby("City").Site.nunique().to_dict())\n'
    'd2[["Sample_ID", "Run_ID", "Source_Sample_Code", "City", "Site", "Collection_Date", "year", "month", "Season", "Sample_Type", "reads_total", "Base_Count", "bracken_phylum_reads"]].to_csv(DER / "ds02_meta.csv", index=False)',
    "AUD-11")

# ============================================================== NB03 (AUD-11): site-aware alpha-diversity inference
f3 = "03_Taxonomic_and_Diversity_Analysis.py"
rep(f3, "the 19 DS02 sampling sites are not identified in the deposited metadata, so repeated monthly samples from the same site cannot be modelled as clusters (city-month aggregation is used as a sensitivity check).",
    "UPDATE (AUD-11): the 19 DS02 open-drain sites were identified from the ENA/SRA experiment titles of PRJNA1415085; the site is the independent unit for city contrasts, and site-aware tests (site-level Kruskal-Wallis with whole-site permutation; mixed models with a site random intercept) are the primary inferential evidence. Sample-level tests are retained as descriptive only.", "AUD-11")
rep(f3, 'family="DS02 alpha diversity", objective="Objective 1")',
    'family="DS02 alpha diversity [sample-level, descriptive only: pseudo-replicated]", note="Superseded for inference by site-aware tests (T06s); AUD-11", objective="Objective 1")', "AUD-11")
rep(f3, '"Kruskal-Wallis tests of DS02 alpha diversity by city and season (BH-adjusted across the family)"',
    '"DESCRIPTIVE (sample-level; 447 repeated samples of 19 sites, pseudo-replicated - see T06s for inference): Kruskal-Wallis of DS02 alpha diversity by city and season"', "AUD-11")
rep(f3, 'family="DS02 alpha models", objective="Objective 1")',
    'family="DS02 alpha models [sample-level OLS, descriptive only: pseudo-replicated]", objective="Objective 1")', "AUD-11")
rep(f3, "print(an.round(4)); print(an_cm.round(4))\n",
    '''print(an.round(4)); print(an_cm.round(4))

# %% [markdown]
# ### A2b. Site-aware inference (AUD-11 correction; primary inferential evidence)
# The independent experimental unit for a city contrast is the **open-drain site** (Delhi 4, Mumbai 8, Chennai 5, Kolkata 2 sites), not the monthly sample. Two complementary analyses are reported: (i) Kruskal-Wallis on site medians with a permutation p-value in which city labels are permuted among whole sites (9,999 permutations; exact-style reference distribution for 19 units); (ii) linear mixed models with a site random intercept and city, season and year fixed effects, each term tested by a maximum-likelihood ratio test (season and year vary within sites and are estimable within site). With only 19 sites (2 in Kolkata) the LRT chi-square reference is approximate, so the permutation p-value is preferred for the city term.

# %%
alpha["Site"] = m2.Site.reindex(alpha.index).astype(str).values
a_s = alpha.assign(City=alpha.City.astype(str), Season=alpha.Season.astype(str))
srows = []
for idx in ["rarefied_richness", "shannon", "gini_simpson", "pielou_evenness"]:
    k = site_level_kw(a_s, idx, "Site", "City")
    srows.append(dict(index=idx, factor="City", analysis="Site-level Kruskal-Wallis (site medians; whole-site permutation)", n_units=k["n_sites"], samples=k["samples"], H=k["H"], df=k["df"], p=k["p"], p_perm=k["p_perm"], epsilon2=k["epsilon2"], LR=np.nan, ICC_site=np.nan, sites_per_city=k["sites_per_group"]))
    record_test(NB, f"DS02 {idx} by City [site level]", "Kruskal-Wallis on site medians; permutation p (city labels permuted among sites)", k["n_sites"], k["H"], k["df"], k["p_perm"], "epsilon-squared (site level)", k["epsilon2"], family="DS02 alpha diversity (site-aware, primary)", note=f"{k['n_sites']} sites ({k['sites_per_group']}); asymptotic p = {k['p']:.3g}", objective="Objective 1")
    mm = mixed_site_lrt(a_s, f"{idx} ~ C(City) + C(Season) + C(year)", ["C(City)", "C(Season)", "C(year)"], "Site")
    for term, rr in mm["lrt"].items():
        srows.append(dict(index=idx, factor=term, analysis="Mixed model, site random intercept; ML likelihood-ratio test", n_units=mm["n_groups"], samples=mm["n"], H=np.nan, df=rr["df"], p=rr["p"], p_perm=np.nan, epsilon2=np.nan, LR=rr["LR"], ICC_site=mm["icc"], sites_per_city=""))
        record_test(NB, f"DS02 {idx} ~ City + Season + Year + (1|Site)", f"Mixed-model LRT (ML): {term}", mm["n"], rr["LR"], rr["df"], rr["p"], "ICC (site)", mm["icc"], family="DS02 alpha diversity (site-aware, primary)", note=f"{mm['n_groups']} sites", objective="Objective 1")
site_tests = pd.DataFrame(srows)
site_tests["p_primary"] = site_tests.p_perm.where(site_tests.p_perm.notna(), site_tests.p)
site_tests["p_primary_BH"] = bh(site_tests.p_primary)
save_table(site_tests.pipe(rnd, 6), "T06s_ds02_alpha_site_aware_tests", "PRIMARY INFERENCE (AUD-11 correction): DS02 alpha diversity - site-level Kruskal-Wallis on 19 site medians (whole-site permutation p) and site random-intercept mixed-model LRTs for city, season and year", NB, sheet="Statistical_Tests")
mm_sh = mixed_site_lrt(a_s, "shannon ~ C(City, Treatment('Delhi')) + C(Season, Treatment('Winter (Dec-Feb)')) + C(year)", [], "Site")
save_table(mm_sh["coef"].pipe(rnd, 5), "T08s_ds02_shannon_mixed_site", f"PRIMARY (AUD-11 correction): mixed model of DS02 phylum Shannon on city, season and year with a site random intercept (REML; {mm_sh['n']} samples, {mm_sh['n_groups']} sites; ICC = {mm_sh['icc']:.2f}; reference Delhi, winter, 2022)", NB, sheet="Comparative_Analysis")
site_tests
''', "AUD-11")
rep(f3, '''    k = tests[(tests["index"] == col) & (tests.factor == "City")].iloc[0]
    ax.set_title(f"{'A' if col=='shannon' else 'B'}  {lab.split(' (')[0]} by city\\nKruskal-Wallis H={k.H:.1f}, p={k.p:.1e}, ε²={k.epsilon2:.2f}", fontsize=8.5)''',
    '''    ks = site_tests[(site_tests["index"] == col) & (site_tests.factor == "City")].iloc[0]
    ax.set_title(f"{'A' if col=='shannon' else 'B'}  {lab.split(' (')[0]} by city\\nsite-level KW ({int(ks.n_units)} sites): H={ks.H:.1f}, perm. p={ks.p_perm:.3f}", fontsize=8.5)''', "AUD-11")
rep(f3, '''k = tests[(tests["index"] == "shannon") & (tests.factor == "Season")].iloc[0]
axes[2].set_xticks(range(4))''', '''kss = site_tests[(site_tests["index"] == "shannon") & (site_tests.factor == "C(Season)")].iloc[0]
axes[2].set_xticks(range(4))''', "AUD-11")
rep(f3, '''axes[2].set_title(f"C  Shannon by season\\nKruskal-Wallis H={k.H:.1f}, p={k.p:.1e}, ε²={k.epsilon2:.2f}", fontsize=8.5)''',
    '''axes[2].set_title(f"C  Shannon by season\\nmixed model (site RE) LRT p={kss.p:.1e}", fontsize=8.5)''', "AUD-11")
rep(f3, "Kruskal-Wallis statistics with epsilon-squared effect sizes are shown; pairwise Holm-adjusted tests are",
    "Panel statistics are site-aware (AUD-11): city contrasts use Kruskal-Wallis on the 19 site medians with a whole-site permutation p-value, and the season contrast a mixed-model likelihood-ratio test with a site random intercept (Table T06s); sample-level pairwise tests (descriptive only) are", "AUD-11")
rep(f3, "Samples from the same (unidentified) site are repeated monthly and are not independent.",
    "Points are repeated monthly samples of 19 identified open-drain sites and are not independent observations.", "AUD-11")
rep(f3, '''record_finding(NB, "F05", "Objective 1", "Phylum-level Shannon diversity differed between cities (median " + ", ".join(f"{c} {v:.2f}" for c, v in med.items()) + f"; Kruskal-Wallis H = {city_p.H:.1f}, df = {int(city_p.df)}, p = {city_p.p:.1e}, epsilon-squared = {city_p.epsilon2:.2f}) and seasons (H = {seas_p.H:.1f}, p = {seas_p.p:.1e}, epsilon-squared = {seas_p.epsilon2:.2f}).",
               len(alpha), "DS02 samples", f"Kruskal-Wallis; BH-adjusted p (city) = {city_p.p_BH:.1e}; adjusted OLS and city-month sensitivity in T08-T10", "Tables T04-T10; Figure 6",
               "City-level differences in diversity are detectable within one standardized study; effect sizes indicate how much variance city explains.", "Sites are unidentified, so repeated sampling of the same site inflates the effective sample size; city-month sensitivity addresses this only partly.")''',
    '''ksite = site_tests[(site_tests["index"] == "shannon") & (site_tests.factor == "City")].iloc[0]
kmc = site_tests[(site_tests["index"] == "shannon") & (site_tests.factor == "C(City)")].iloc[0]
kms = site_tests[(site_tests["index"] == "shannon") & (site_tests.factor == "C(Season)")].iloc[0]
record_finding(NB, "F05", "Objective 1", "Median phylum-level Shannon diversity was " + ", ".join(f"{c} {v:.2f}" for c, v in med.items()) + f". With the open-drain site as the independent unit (19 sites), the city difference was {'statistically supported' if ksite.p_perm < .05 else 'not statistically supported'} (Kruskal-Wallis on site medians H = {ksite.H:.1f}, df = {int(ksite.df)}, whole-site permutation p = {ksite.p_perm:.3f}; site-level epsilon-squared = {ksite.epsilon2:.2f}; mixed-model LRT p = {kmc.p:.2g}); within sites, season was {'associated' if kms.p < .05 else 'not associated'} with Shannon diversity (mixed-model LRT p = {kms.p:.2g}).",
               len(alpha), f"DS02: {len(alpha)} samples from {int(ksite.n_units)} sites", "Site-aware (AUD-11): site-level Kruskal-Wallis with whole-site permutation; site random-intercept mixed models (T06s, T08s). Sample-level p-values (T06-T09) are descriptive only.", "Tables T04, T06s, T08s; Figure 6",
               f"Effect directions agree with the sample-level description (sample-level epsilon-squared = {city_p.epsilon2:.2f}, descriptive), but the evidence for a city effect rests on 19 sites, not 447 samples; the site ICC was {mm_sh['icc']:.2f}.", "Few sites per city (Kolkata 2); city is confounded with site selection; phylum resolution only.")''', "AUD-11")

# ============================================================== NB04 (AUD-11 + AUD-12): site-aware ARG inference, detection tiers
f4 = "04_AMR_Gene_Analysis.py"
rep(f4, 'family="DS02 ARG load/diversity", objective="Objective 2")',
    'family="DS02 ARG abundance/diversity [sample-level, descriptive only: pseudo-replicated]", note="Superseded for inference by site-aware tests (A05s); AUD-11", objective="Objective 2")', "AUD-11")
rep(f4, '"Kruskal-Wallis tests: DS02 ARG load and diversity by city and season (BH-adjusted)"',
    '"DESCRIPTIVE (sample-level, pseudo-replicated - see A05s for inference): Kruskal-Wallis of DS02 normalized ARG abundance and diversity by city and season"', "AUD-11")
rep(f4, 'family="DS02 ARG models", objective="Objective 2/3")',
    'family="DS02 ARG models [sample-level OLS, descriptive only: pseudo-replicated]", objective="Objective 2/3")', "AUD-11")
rep(f4, "print(kt.round(4)); co.round(3)\n", '''print(kt.round(4)); co.round(3)

# %% [markdown]
# ### A3b. Site-aware inference for ARG abundance and diversity (AUD-11 correction; primary)
# Same design as notebook 03: Kruskal-Wallis on the 19 site medians with whole-site permutation of city labels, and mixed models with a site random intercept (city, season and year fixed; ML likelihood-ratio tests). Fold changes for the city contrasts come from the REML mixed model on log10 RPKM.

# %%
s2c["Site"] = m2.Site.reindex(s2c.index).astype(str).values
arows = []
for v in ["log10_total_RPKM", "family_richness_cov50", "family_shannon"]:
    k = site_level_kw(s2c, v, "Site", "City")
    arows.append(dict(variable=v, factor="City", analysis="Site-level Kruskal-Wallis (site medians; whole-site permutation)", n_units=k["n_sites"], samples=k["samples"], H=k["H"], df=k["df"], p=k["p"], p_perm=k["p_perm"], epsilon2=k["epsilon2"], LR=np.nan, ICC_site=np.nan))
    record_test(NB, f"DS02 {v} by City [site level]", "Kruskal-Wallis on site medians; permutation p (city labels permuted among sites)", k["n_sites"], k["H"], k["df"], k["p_perm"], "epsilon-squared (site level)", k["epsilon2"], family="DS02 ARG abundance/diversity (site-aware, primary)", note=f"{k['n_sites']} sites; asymptotic p = {k['p']:.3g}", objective="Objective 2")
    mm = mixed_site_lrt(s2c, f"{v} ~ C(City) + C(Season) + C(year)", ["C(City)", "C(Season)", "C(year)"], "Site")
    for term, rr in mm["lrt"].items():
        arows.append(dict(variable=v, factor=term, analysis="Mixed model, site random intercept; ML likelihood-ratio test", n_units=mm["n_groups"], samples=mm["n"], H=np.nan, df=rr["df"], p=rr["p"], p_perm=np.nan, epsilon2=np.nan, LR=rr["LR"], ICC_site=mm["icc"]))
        record_test(NB, f"DS02 {v} ~ City + Season + Year + (1|Site)", f"Mixed-model LRT (ML): {term}", mm["n"], rr["LR"], rr["df"], rr["p"], "ICC (site)", mm["icc"], family="DS02 ARG abundance/diversity (site-aware, primary)", note=f"{mm['n_groups']} sites", objective="Objective 2/3")
kt_site = pd.DataFrame(arows)
kt_site["p_primary"] = kt_site.p_perm.where(kt_site.p_perm.notna(), kt_site.p); kt_site["p_primary_BH"] = bh(kt_site.p_primary)
save_table(kt_site.pipe(rnd, 6), "A05s_ds02_arg_site_aware_tests", "PRIMARY INFERENCE (AUD-11 correction): DS02 normalized ARG abundance (log10 RPKM) and ARG diversity - site-level Kruskal-Wallis (19 sites, whole-site permutation) and site random-intercept mixed-model LRTs", NB, sheet="Statistical_Tests")
mm_arg = mixed_site_lrt(s2c, "log10_total_RPKM ~ C(City, Treatment('Delhi')) + C(Season, Treatment('Winter (Dec-Feb)')) + C(year)", [], "Site")
co_s = mm_arg["coef"].copy()
fe_ = co_s.term.str.startswith("C(") | (co_s.term == "Intercept")
co_s["fold_change"] = np.where(fe_ & (co_s.term != "Intercept"), 10 ** co_s.estimate, np.nan)
co_s["fc_ci_low"] = np.where(fe_ & (co_s.term != "Intercept"), 10 ** co_s.ci_low, np.nan); co_s["fc_ci_high"] = np.where(fe_ & (co_s.term != "Intercept"), 10 ** co_s.ci_high, np.nan)
save_table(co_s.pipe(rnd, 5), "A07s_ds02_arg_mixed_site", f"PRIMARY (AUD-11 correction): mixed model of log10 total normalized ARG abundance (RPKM) on city, season and year with a site random intercept (REML; {mm_arg['n']} samples, {mm_arg['n_groups']} sites; ICC = {mm_arg['icc']:.2f}); fold changes vs Delhi", NB, sheet="Comparative_Analysis")
kt_site
''', "AUD-11")
rep(f4, '''k = kt[(kt.variable == "log10_total_RPKM") & (kt.factor == "City")].iloc[0]
axes[1].set_ylabel("Total ARG abundance (RPKM, log scale)"); axes[1].set_xlabel("")''',
    '''k = kt_site[(kt_site.variable == "log10_total_RPKM") & (kt_site.factor == "City")].iloc[0]
axes[1].set_ylabel("Total normalized ARG abundance (RPKM, log scale)"); axes[1].set_xlabel("")''', "AUD-11/19")
rep(f4, '''axes[1].set_title(f"B  Total ARG load by city\\nKruskal-Wallis H={k.H:.1f}, p={k.p:.1e}, ε²={k.epsilon2:.2f}", fontsize=8.5)''',
    '''axes[1].set_title(f"B  Normalized ARG abundance by city\\nsite-level KW ({int(k.n_units)} sites): H={k.H:.1f}, perm. p={k.p_perm:.3f}", fontsize=8.5)''', "AUD-11/19")
rep(f4, "log scale, with Kruskal-Wallis statistics. n = 447.",
    "log scale; statistic = Kruskal-Wallis on 19 site medians with whole-site permutation p (Table A05s). n = 447 samples from 19 sites. RPKM is a depth-normalized relative abundance, not an absolute quantity per volume or per cell.", "AUD-11/19")
# AUD-12: detection tiers - strict (high-confidence, primary) and lenient (unconfirmed, sensitivity)
rep(f4, "# Presence is defined at allele coverage ≥ 50 % (primary) and at any deposited hit (lenient).",
    "# **UPDATE (AUD-12).** The deposited RGI bwt table provides, per CARD allele, all mapped reads, completely mapped reads and average percent coverage; it does not provide percent identity, MAPQ or the CARD release, and RGI bwt documents no detection threshold. Detection is therefore reported in three tiers: **(1) high-confidence (primary): at least one allele with ≥ 50 % reference coverage and ≥ 2 completely mapped reads;** (2) a stricter sensitivity tier (≥ 80 % coverage); (3) a lenient tier (any deposited read hit) reported only as an **unconfirmed upper bound** (cross-mapping from homologous genes cannot be excluded). Any-hit detection is never described as confirmed resistance. Presence is defined at allele coverage ≥ 50 % (primary) and at any deposited hit (lenient).", "AUD-12")
rep(f4, '''card["cov50"] = card.Average_Percent_Coverage >= 50''',
    '''card["cov50_only"] = card.Average_Percent_Coverage >= 50
card["completely_mapped"] = pd.to_numeric(card.Completely_Mapped_Reads, errors="coerce").fillna(0)
card["cov50"] = card.cov50_only & (card.completely_mapped >= 2)   # AUD-12 high-confidence tier (primary)
card["cov80"] = (card.Average_Percent_Coverage >= 80) & (card.completely_mapped >= 2)
print("alleles cov>=50%:", int(card.cov50_only.sum()), "| high-confidence (cov>=50% & >=2 complete reads):", int(card.cov50.sum()), "| cov>=80%:", int(card.cov80.sum()))''', "AUD-12")
rep(f4, '''clin2 = panel_detect(card, "Sample_ID", m2.City.astype(str), pipeline="CARD/RGI (DS02)")''',
    '''clin2 = panel_detect(card, "Sample_ID", m2.City.astype(str), pipeline="CARD/RGI (DS02)")
clin2_80 = panel_detect(card, "Sample_ID", m2.City.astype(str), strict_col="cov80", pipeline="CARD/RGI (DS02)")
clin2["detected_cov80"] = clin2_80.detected_cov50.values; clin2["pct_cov80"] = clin2_80.pct_cov50.values
clin2 = clin2.rename(columns={"detected_cov50": "detected_high_confidence", "pct_cov50": "pct_high_confidence", "ci_low": "hc_ci_low", "ci_high": "hc_ci_high",
                              "detected_any_hit": "detected_any_hit_unconfirmed", "pct_any_hit": "pct_any_hit_unconfirmed"})
clin2["uncertainty_range_pp"] = clin2.pct_any_hit_unconfirmed - clin2.pct_high_confidence
# keep the original column names used by the plotting helper
clin2_plot = clin2.rename(columns={"detected_high_confidence": "detected_cov50", "pct_high_confidence": "pct_cov50", "hc_ci_low": "ci_low", "hc_ci_high": "ci_high", "pct_any_hit_unconfirmed": "pct_any_hit"})''', "AUD-12")
rep(f4, '''save_table(clin2.pipe(rnd, 2), "A09_ds02_clinical_determinants", "DS02: detection of clinically important ARG determinants by city (Wilson 95% CI; coverage >= 50% and any hit)", NB, sheet="AMR_Gene_Summary")''',
    '''save_table(clin2.pipe(rnd, 2), "A09_ds02_clinical_determinants", "DS02 clinically important ARG determinants (AUD-12): high-confidence detection (primary; allele coverage >= 50% and >= 2 completely mapped reads; Wilson 95% CI), stricter >= 80% coverage tier, and any-hit detection (unconfirmed upper bound, sensitivity only)", NB, sheet="AMR_Gene_Summary")''', "AUD-12")
rep(f4, '''    t = clin2[(clin2.determinant == lab) & (clin2.group != "All")]''', '''    t = clin2_plot[(clin2_plot.determinant == lab) & (clin2_plot.group != "All")]''', "AUD-12")
rep(f4, '''    record_test(NB, f"DS02 detection of {r.determinant} (cov>=50%) by city", "Chi-square test of independence", len(m2), r.chi2, int(r.df), r.p, "Cramer's V", r.cramers_v, family="DS02 clinical determinants", note=r.note, objective="Objective 2/3")
save_table(chi.pipe(rnd, 6), "A10_ds02_clinical_by_city_tests", "Chi-square tests: city differences in detection of clinically important determinants (DS02; BH-adjusted)", NB, sheet="Statistical_Tests")''',
    '''    record_test(NB, f"DS02 detection of {r.determinant} (high-confidence) by city [sample level]", "Chi-square test of independence", len(m2), r.chi2, int(r.df), r.p, "Cramer's V", r.cramers_v, family="DS02 clinical determinants [sample-level, descriptive only: pseudo-replicated]", note=r.note + "; superseded by A10s (AUD-11)", objective="Objective 2/3")
save_table(chi.pipe(rnd, 6), "A10_ds02_clinical_by_city_tests", "DESCRIPTIVE (sample-level chi-square, pseudo-replicated - see A10s): city differences in high-confidence detection of clinically important determinants (DS02)", NB, sheet="Statistical_Tests")
# AUD-11: cluster (whole-site) permutation chi-square - primary
hc_units = {lab: set(card.loc[(card.Clinical == lab) & card.cov50, "Sample_ID"]) for lab in CLINICAL_PANEL}
chs = []
for lab in CLINICAL_PANEL:
    det = m2.index.isin(hc_units[lab])
    r_ = site_perm_chisq(det, m2.Site.astype(str).values, m2.City.astype(str).values)
    chs.append(dict(determinant=lab, detection_rule="high-confidence (cov >= 50% and >= 2 completely mapped reads)", samples=r_["n"], sites=r_["n_sites"], detected_samples=int(det.sum()), chi2=r_["chi2"], p_perm_whole_sites=r_["p_perm"], site_prevalence_by_city=r_.get("site_prevalence_by_group", "constant (all or none detected)")))
    if not np.isnan(r_["p_perm"]):
        record_test(NB, f"DS02 detection of {lab} (high-confidence) by city [site-aware]", "Chi-square statistic; permutation of city labels among whole sites", r_["n"], r_["chi2"], 3, r_["p_perm"], family="DS02 clinical determinants (site-aware, primary)", note=f"{r_['n_sites']} sites", objective="Objective 2/3")
chs = pd.DataFrame(chs); chs["p_BH"] = bh(chs.p_perm_whole_sites)
save_table(chs.pipe(rnd, 6), "A10s_ds02_clinical_by_city_site_aware", "PRIMARY (AUD-11/AUD-12): city differences in high-confidence detection of clinically important determinants, tested by permuting city labels among the 19 whole sites; site-level prevalence summaries", NB, sheet="Statistical_Tests")''', "AUD-11/12")
rep(f4, '''clin2[clin2.group == "All"].round(1)''', '''clin2[clin2.group == "All"].round(1)''', "AUD-12")
rep(f4, '''panel_plot(ax, clin2, cities, CITY_COLORS, "Clinically important determinants, CARD/RGI (DS02)")''',
    '''panel_plot(ax, clin2_plot, cities, CITY_COLORS, "Clinically important determinants, CARD/RGI (DS02)")''', "AUD-12")
rep(f4, '''ax.set_xlabel("Units with determinant detected (%): filled = coverage ≥ 50% (lower bound, Wilson 95% CI); open = any hit (upper bound)", fontsize=7.5)''',
    '''ax.set_xlabel("Units with determinant detected (%)\\nfilled = high-confidence (primary; Wilson 95% CI); open = any hit (unconfirmed upper bound)", fontsize=7.5)''', "AUD-12")
rep(f4, "Filled points: at least one CARD allele with ≥ 50% reference coverage (Wilson 95% CI); open points: any deposited RGI hit. The strict rule underestimates multi-allele families because RGI bwt splits reads among near-identical alleles; the lenient rule overestimates because of cross-mapping from homologous genes, so the two values bracket the detection frequency.",
    "Filled points (primary, high-confidence): at least one CARD allele with ≥ 50% reference coverage and ≥ 2 completely mapped reads (Wilson 95% CI). Open points (sensitivity only): any deposited RGI bwt read hit, an unconfirmed upper bound that is vulnerable to cross-mapping from homologous genes. The high-confidence rule may underestimate multi-allele families because RGI bwt splits reads among near-identical alleles. Percent identity, MAPQ and the CARD release were not deposited.", "AUD-12")
rep(f4, '''kc = kt[(kt.variable == "log10_total_RPKM") & (kt.factor == "City")].iloc[0]
record_finding(NB, "F07", "Objective 2", f"In Indian open-drain wastewater (CARD/RGI), total ARG load differed between cities (Kruskal-Wallis epsilon-squared = {kc.epsilon2:.2f}, p = {kc.p:.1e}); the leading class groups by share of abundance were {', '.join(top_cls)}.",
               len(s2), "DS02 samples", "Kruskal-Wallis; adjusted OLS (A07-A08)", "Tables A02, A04-A08; Figure 9",
               "City-level resistome differences exist within one harmonised study; adjusted city contrasts are given as fold changes.",
               "RPKM uses total reads (constant factor for paired data); site identity unavailable; CARD release unstated; efflux families may include housekeeping genes.")''',
    '''kc = kt[(kt.variable == "log10_total_RPKM") & (kt.factor == "City")].iloc[0]
kcs = kt_site[(kt_site.variable == "log10_total_RPKM") & (kt_site.factor == "City")].iloc[0]
kcm = kt_site[(kt_site.variable == "log10_total_RPKM") & (kt_site.factor == "C(City)")].iloc[0]
ksm = kt_site[(kt_site.variable == "log10_total_RPKM") & (kt_site.factor == "C(Season)")].iloc[0]
fcs = co_s[co_s.term.str.contains("City")]
record_finding(NB, "F07", "Objective 2", f"In Indian open-drain wastewater (CARD/RGI), median total normalized ARG abundance was " + ", ".join(f"{c} {summ_city.loc[c, 'median_total_RPKM']:.0f}" for c in summ_city.index) + f" RPKM. With the site as the unit (19 sites), the city difference was {'statistically supported' if kcs.p_perm < .05 else 'not statistically supported'} (site-level Kruskal-Wallis H = {kcs.H:.1f}, whole-site permutation p = {kcs.p_perm:.3f}; mixed-model LRT p = {kcm.p:.2g}); mixed-model fold changes vs Delhi: " + "; ".join(f"{term_label(t)} {f:.2f} (95% CI {l:.2f}-{h:.2f})" for t, f, l, h in zip(fcs.term, fcs.fold_change, fcs.fc_ci_low, fcs.fc_ci_high)) + f"; season within sites: LRT p = {ksm.p:.2g}. Leading class groups by share of abundance: {', '.join(top_cls)}.",
               len(s2), f"DS02: {len(s2)} samples from {int(kcs.n_units)} sites", "Site-aware (AUD-11): site-level Kruskal-Wallis with whole-site permutation; site random-intercept mixed model (A05s, A07s). Sample-level p-values (A05-A08) are descriptive only.", "Tables A02, A04, A05s, A07s; Figure 9",
               f"City differences in normalized ARG abundance are estimated from 19 sites; the site ICC was {mm_arg['icc']:.2f}, so sample-level tests overstated certainty (descriptive sample-level epsilon-squared {kc.epsilon2:.2f}).",
               "RPKM is a depth-normalized relative abundance (not an absolute quantity); few sites per city; CARD release unstated; efflux families may include housekeeping genes.")''', "AUD-11/19")
rep(f4, '''record_finding(NB, "F08", "Objective 2/3", "Detection frequencies of clinically important determinants in DS02 wastewater (strict coverage >= 50% to lenient any-hit bounds): " + "; ".join(f"{d} {allp.loc[d,'pct_cov50']:.0f}-{allp.loc[d,'pct_any_hit']:.0f}%" for d in CLINICAL_PANEL) + ".",
               len(m2), "DS02 samples", "Wilson 95% CIs for strict frequencies; chi-square by city (A10)", "Table A09; Figure 11",
               f"Detected in >= 80% of samples even under the strict rule: {', '.join(allp.index[allp.pct_cov50 >= 80]) or 'none'}; determinants whose frequency changes by > 20 percentage points between rules ({', '.join(allp.index[(allp.pct_any_hit - allp.pct_cov50) > 20])}) should be reported as a range.",
               "Gene-level read mapping without host, plasmid context or viability; allele-level coverage is diluted across near-identical alleles.")''',
    '''sig_s = chs[chs.p_BH < .05].determinant.tolist()
record_finding(NB, "F08", "Objective 2/3", "High-confidence detection (primary; allele coverage >= 50% and >= 2 completely mapped reads) of clinically important determinants in DS02 wastewater samples: " + "; ".join(f"{d} {allp.loc[d,'pct_high_confidence']:.0f}% (95% CI {allp.loc[d,'hc_ci_low']:.0f}-{allp.loc[d,'hc_ci_high']:.0f})" for d in CLINICAL_PANEL) + ". Any-hit detection (unconfirmed upper bound; sensitivity only): " + "; ".join(f"{d} {allp.loc[d,'pct_any_hit_unconfirmed']:.0f}%" for d in CLINICAL_PANEL) + ".",
               len(m2), f"DS02: {len(m2)} samples from 19 sites", "Wilson 95% CIs (high-confidence tier); city differences by whole-site permutation chi-square (A10s; BH)", "Tables A09, A10s; Figure 11",
               f"High-confidence detection in >= 80% of samples: {', '.join(allp.index[allp.pct_high_confidence >= 80]) or 'none'}. Determinants with > 20 percentage-point uncertainty between tiers ({', '.join(allp.index[allp.uncertainty_range_pp > 20]) or 'none'}) must be cited at the high-confidence value, with the lenient value described as unconfirmed. City differences supported at the site level (BH q < 0.05): {', '.join(sig_s) or 'none'}.",
               "Read mapping only (no assembly, host, plasmid context, expression or viability); identity/MAPQ not deposited; detection is not phenotypic or clinical resistance.")''', "AUD-12")
rep(f4, '''record_finding(NB, "F09", "Objective 2/5", "In the ResFinder resource, median total ARG load was "''',
    '''record_finding(NB, "F09", "Objective 2/5", "In the ResFinder resource, median total normalized ARG abundance (FPKM, a relative measure per million rRNA-aligned bacterial fragments) was "''', "AUD-19")

# ============================================================== NB05 (AUD-11, AUD-13)
f5 = "05_Comparative_and_Multivariate_Analysis.py"
rep(f5, 'meta2 = m2[["City", "Season", "year"]].reset_index()', 'meta2 = m2[["City", "Season", "year", "Site"]].reset_index()', "AUD-11")
rep(f5, '''    for term, cov, strata in [("City", (), None), ("Season", ("City",), "City"), ("year", ("City", "Season"), "City")]:
        r = permanova_dbrda(D, meta2, term, cov, strata=strata, n_perm=N_PERM)''',
    '''    # UPDATE (AUD-11): city is constant within site, so city labels are permuted among whole sites (block = Site);
    # season and year vary within site and are permuted within site. A free-permutation city row is kept only as a descriptive/cross-check value.
    rf_ = permanova_dbrda(D, meta2, "City", (), n_perm=N_PERM); rf_.update(features=feat, distance=dist, covariates="(none; free permutation, descriptive)"); res.append(rf_)
    for term, cov, strata, block in [("City", (), None, "Site"), ("Season", ("City",), "Site", None), ("year", ("City", "Season"), "Site", None)]:
        r = permanova_dbrda(D, meta2, term, cov, strata=strata, block=block, n_perm=N_PERM)''', "AUD-11")
rep(f5, 'family="DS02 PERMANOVA", note=r["permutation_scheme"], objective="Objective 1/2/3")',
    'family="DS02 PERMANOVA (site-aware, primary)", note=r["permutation_scheme"], objective="Objective 1/2/3")', "AUD-11")
rep(f5, 'record_test(NB, f"DS02 {feat} ({dist}) dispersion by {fac}", "PERMDISP", d["n"], d["F"], f"{d[\'groups\']-1}, {d[\'n\']-d[\'groups\']}", d["p_value"], family="DS02 PERMDISP", objective="Objective 1/2")',
    'record_test(NB, f"DS02 {feat} ({dist}) dispersion by {fac}", "PERMDISP", d["n"], d["F"], f"{d[\'groups\']-1}, {d[\'n\']-d[\'groups\']}", d["p_value"], family="DS02 PERMDISP [sample-level, descriptive only]", objective="Objective 1/2")', "AUD-11")
rep(f5, '"DS02 PERMANOVA (sequential, restricted permutations) and PERMDISP for taxonomic and resistome distances"',
    '"DS02 PERMANOVA (sequential): city tested by permuting whole sites (19), season/year permuted within site (AUD-11, primary); free-permutation city rows and PERMDISP are descriptive"', "AUD-11")
rep(f5, '''own = perm[(perm.distance == "Bray-Curtis") & (perm.term == "City") & (perm.covariates == "(none)")].iloc[0]''',
    '''own = perm[(perm.distance == "Bray-Curtis") & (perm.term == "City") & (perm.covariates == "(none; free permutation, descriptive)")].iloc[0]''', "AUD-11")
rep(f5, '''"PERMANOVA (999 permutations; season permuted within city) with PERMDISP", "Table C01",''',
    '''"PERMANOVA (999 permutations): city labels permuted among the 19 whole sites; season permuted within site (AUD-11); PERMDISP descriptive", "Table C01",''', "AUD-11")
rep(f5, '''"Repeated sampling of unidentified sites; season and year are partly collinear within a 25-month series.")''',
    '''"Inference on city rests on 19 sites; season and year are partly collinear within a 25-month series.")''', "AUD-11")
rep(f5, '''save_table(t1.pipe(rnd, 5), "C08_rf_mixed_arg_load",''', '''# UPDATE (AUD-13): identifiability of matrix effects given the study structure
sm_struct = urb.groupby("Sample_Type").agg(units=("Project_Accession", "size"), studies=("Project_Accession", "nunique"),
                                           study_list=("Project_Accession", lambda s: "; ".join(f"{k} ({v})" for k, v in s.value_counts().items()))).reset_index()
sm_struct["largest_study_share_pct"] = [100 * urb[urb.Sample_Type == t].Project_Accession.value_counts().iloc[0] / n for t, n in zip(sm_struct.Sample_Type, sm_struct.units)]
sm_struct["matrix_effect_identifiable_from_study_effect"] = np.where(sm_struct.studies >= 3, "Partly (>= 3 studies; still study-limited)", np.where(sm_struct.studies == 2, "Weakly (2 studies)", "No - single-study matrix: matrix and study effects are not separable"))
t1["single_study_matrix"] = [("Yes" if sm_struct.set_index("Sample_Type").studies.get(term_label(x), 99) == 1 else "No") if "Sample_Type" in str(x) else "" for x in t1.term]
t1["studies_for_matrix"] = [sm_struct.set_index("Sample_Type").studies.get(term_label(x), np.nan) if "Sample_Type" in str(x) else np.nan for x in t1.term]
save_table(sm_struct.pipe(rnd, 3), "C08b_rf_matrix_study_structure", "AUD-13: study structure of each urban matrix in the ResFinder resource and identifiability of matrix effects", NB, sheet="Comparative_Analysis")
# Within-study matrix contrasts: only studies containing >= 2 matrices with >= 3 units each
ws_rows = []
for prj, d in urb.groupby("Project_Accession"):
    cnt = d.Sample_Type.value_counts(); ok = cnt[cnt >= 3].index.tolist()
    if len(ok) < 2: continue
    ref = REF if REF in ok else ok[0]
    for t_ in [x for x in ok if x != ref]:
        a_ = d.loc[d.Sample_Type == t_, "log10_total_FPKM"]; b_ = d.loc[d.Sample_Type == ref, "log10_total_FPKM"]
        u_, p_ = stats.mannwhitneyu(a_, b_, alternative="two-sided")
        ws_rows.append(dict(study=prj, matrix=t_, reference=ref, n_matrix=len(a_), n_reference=len(b_), median_log10_FPKM_matrix=a_.median(), median_log10_FPKM_reference=b_.median(),
                            fold_change_of_medians=10 ** (a_.median() - b_.median()), mann_whitney_U=u_, p_value=p_, cliffs_delta=cliffs_delta_from_u(u_, len(a_), len(b_))))
ws = pd.DataFrame(ws_rows, columns=["study", "matrix", "reference", "n_matrix", "n_reference", "median_log10_FPKM_matrix", "median_log10_FPKM_reference", "fold_change_of_medians", "mann_whitney_U", "p_value", "cliffs_delta"])
if len(ws): ws["p_holm"] = holm(ws.p_value)
save_table(ws.pipe(rnd, 5), "C08c_rf_within_study_matrix_contrasts", "AUD-13: within-study matrix contrasts of log10 total normalized ARG abundance (FPKM), studies with >= 2 matrices of >= 3 units (no between-study confounding)", NB, sheet="Comparative_Analysis",
           description="Empty if no study sampled two or more matrices with at least three units each.")
print(sm_struct); print(ws)
save_table(t1.pipe(rnd, 5), "C08_rf_mixed_arg_load",''', "AUD-13")
rep(f5, '''record_finding(NB, "F12", "Objective 2", "After accounting for study as a random effect, urban matrices differed from municipal wastewater in total ARG load by fold changes of " + "; ".join(f"{term_label(x)} {10**e:.2f} (95% CI {10**l:.2f}-{10**h:.2f})" for x, e, l, h in zip(fc.term, fc.estimate, fc.ci_low, fc.ci_high)) + f"; the study ICC was {t1.set_index('term').loc['ICC (study)','estimate']:.2f}.",''',
    '''ss_ = sm_struct.set_index("Sample_Type").studies
record_finding(NB, "F12", "Objective 2", "Study-adjusted (random-intercept) fold changes in total normalized ARG abundance (FPKM) relative to municipal wastewater were " + "; ".join(f"{term_label(x)} {10**e:.2f} (95% CI {10**l:.2f}-{10**h:.2f}; {int(ss_.get(term_label(x), 0))} {'study - single-study matrix, not separable from study effect' if ss_.get(term_label(x), 0) == 1 else 'studies'})" for x, e, l, h in zip(fc.term, fc.estimate, fc.ci_low, fc.ci_high)) + f"; the study ICC was {t1.set_index('term').loc['ICC (study)','estimate']:.2f}. Within-study contrasts were possible for {ws.study.nunique() if len(ws) else 0} stud{'y' if (ws.study.nunique() if len(ws) else 0) == 1 else 'ies'} (Table C08c).",''', "AUD-13/19")
rep(f5, '''f"Matrices with credibly lower load than municipal wastewater (95% CI excludes 1): {', '.join(lower) or 'none'}; credibly higher: {', '.join(higher) or 'none'}. Most matrices are represented by one or two studies, so these are study-limited contrasts.",
               "Matrix contrasts rest on few studies; KMA versions and depth differ by study.")''',
    '''f"Single-study matrices ({', '.join(ss_[ss_ == 1].index) or 'none'}) cannot be interpreted as independent environmental-matrix effects; their fold changes describe one study and depend on the random-effect assumption. Estimates for the remaining matrices are study-limited and not generalizable (AUD-13).",
               "Matrix is confounded with study, KMA version and sequencing depth; few studies per matrix.")''', "AUD-13")

# ============================================================== NB08 (AUD-17): fold-internal zero replacement
f8 = "08_Predictive_Modeling.py"
rep(f8, "Pz, _ = mult_replace(Xc.values); XA = pd.DataFrame(clr(Pz), index=Xc.index, columns=Xc.columns)",
    '''from sklearn.base import BaseEstimator, TransformerMixin
class FoldwiseZeroReplaceCLR(BaseEstimator, TransformerMixin):
    """AUD-17 fix: multiplicative zero replacement whose constant (frac x smallest non-zero share) is
    estimated on the training rows of each fold only, followed by the centred log-ratio transform."""
    def __init__(self, frac=0.65): self.frac = frac
    def fit(self, X, y=None):
        P = closure(np.asarray(X, float)); self.delta_ = float(self.frac * P[P > 0].min()); return self
    def transform(self, X):
        Pz, _ = mult_replace(np.asarray(X, float), delta=self.delta_); return clr(Pz)
XA = pd.DataFrame(closure(Xc.values), index=Xc.index, columns=Xc.columns)   # closed class shares; CLR is fitted inside each training fold
delta_all = 0.65 * XA.values[XA.values > 0].min()
print("zero-replacement constant on all units (previous, transductive):", delta_all)''', "AUD-17")
rep(f8, '''dict(item="Features", value=f"CLR of {XA.shape[1]} class-group shares of total FPKM"),''',
    '''dict(item="Features", value=f"CLR of {XA.shape[1]} class-group shares of total FPKM"),
    dict(item="Pre-processing (AUD-17 fix)", value="Multiplicative zero-replacement constant (0.65 x smallest non-zero share) and CLR estimated on the training fold only, inside each model pipeline; test folds never inform pre-processing"),''', "AUD-17")
rep(f8, '''models = {"Baseline (prior)": DummyClassifier(strategy="prior"),
          "Logistic regression (L2)": make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=5000)),
          "Random forest": RandomForestClassifier(n_estimators=500, min_samples_leaf=2, class_weight="balanced", random_state=SEED, n_jobs=-1),
          "Gradient boosting (HistGB)": HistGradientBoostingClassifier(max_depth=3, learning_rate=.05, max_iter=300, random_state=SEED)}''',
    '''models = {"Baseline (prior)": make_pipeline(FoldwiseZeroReplaceCLR(), DummyClassifier(strategy="prior")),
          "Logistic regression (L2)": make_pipeline(FoldwiseZeroReplaceCLR(), StandardScaler(), LogisticRegression(C=1.0, max_iter=5000)),
          "Random forest": make_pipeline(FoldwiseZeroReplaceCLR(), RandomForestClassifier(n_estimators=500, min_samples_leaf=2, class_weight="balanced", random_state=SEED, n_jobs=-1)),
          "Gradient boosting (HistGB)": make_pipeline(FoldwiseZeroReplaceCLR(), HistGradientBoostingClassifier(max_depth=3, learning_rate=.05, max_iter=300, random_state=SEED))}''', "AUD-17")
rep(f8, '''rf.fit(XA.values, yA); gini = pd.Series(rf.feature_importances_, index=XA.columns, name="impurity_importance_full_fit")''',
    '''rf.fit(XA.values, yA); gini = pd.Series(rf[-1].feature_importances_, index=XA.columns, name="impurity_importance_full_fit")''', "AUD-17")
rep(f8, '''save_table(feasP, "M07_predictive_feasibility",''', '''# AUD-17 documentation: Model B uses a fixed pseudo-count (0.5) for CLR, which is not estimated from data; its prevalence filter uses no outcome information.
save_table(feasP, "M07_predictive_feasibility",''', "AUD-17")

# ============================================================== NB10 (AUD-20, AUD-11/16 limitation text)
f10 = "10_Final_Results_and_Visualization.py"
rep(f10, '''save_table(figs[["number", "file", "title", "caption", "objective", "notebook", "source_data"]], "R06_figure_catalog", "Figure catalogue", NB, sheet="Figure_Gallery")\n''', "", "AUD-20")
rep(f10, '''tabs2 = load_catalog("tables")''', '''# AUD-20 fix: the figure catalogue is written after Figure 30 so that it is complete
figs = load_catalog("figures").sort_values("number")
fcols = [c for c in ["number", "file", "tiff_file", "title", "caption", "objective", "notebook", "source_data", "width_mm", "height_mm", "print_layout", "min_font_pt", "text_overlaps_flagged"] if c in figs.columns]
save_table(figs[fcols], "R06_figure_catalog", f"Figure catalogue ({len(figs)} figures; final print size and legibility audit)", NB, sheet="Figure_Gallery")
assert len(figs) == figs.number.nunique() == 30, len(figs)
tabs2 = load_catalog("tables")''', "AUD-20")
rep(f10, ''' ("Data acquisition", "No raw reads were reprocessed; all AMR and taxonomic results come from producer-processed tables (two pipelines).", "Cross-study comparison of absolute abundance is not possible; results are within-pipeline.", "Uniform reprocessing of an approved raw-read subset (est. 24 TB total)."),''',
    ''' ("Data acquisition", "Global analyses (notebooks 01-10) use producer-processed tables (two pipelines); only 45 runs (16 Dhaka + 29 South Asian comparators) were uniformly reprocessed from raw reads (notebooks 12-16).", "Cross-study comparison of abundance is valid only within one pipeline; the reprocessed set is a third, separate analytical group.", "Uniform reprocessing of a larger, approved raw-read subset."),''', "UPDATE")
rep(f10, ''' ("Coverage", f"Processed AMR tables cover 1,072 of 8,053 included samples; taxonomy covers 447.", "Most cities and matrices (MetaSUB transit surfaces, Dhaka, Karachi) contribute metadata only.", "Reprocess MetaSUB, Dhaka and Karachi raw reads."),''',
    ''' ("Coverage", f"Producer-processed AMR tables cover 1,072 of 8,053 included samples; taxonomy covers 447. Dhaka (16) and Karachi (7) now have uniformly reprocessed profiles (4 M read pairs per run).", "MetaSUB transit surfaces and most other cities still contribute metadata only.", "Reprocess further cities and matrices with the same pipeline."),''', "UPDATE")
rep(f10, ''' ("Pseudo-replication", "DS02 sites (19) are not identified; monthly repeats are treated as independent.", "Effective sample size is overstated; p-values are optimistic (city-month sensitivity provided).", "Obtain site IDs (Singh et al. Supplementary Data S7)."),''',
    ''' ("Pseudo-replication (resolved, AUD-11)", "DS02 site identities (19 sites) were recovered from ENA/SRA experiment titles; city inference now uses site-level permutation tests and site random-intercept models.", "City contrasts rest on 19 sites (Kolkata 2); sample-level p-values are descriptive only.", "Confirm site labels against Singh et al. Supplementary Data S7 when accessible."),''', "AUD-11")
rep(f10, ''' ("Detection thresholds", "Detection frequency of clinically important genes depends strongly on the coverage threshold.", "Frequencies reported as strict-lenient ranges.", "Assembly-based or variant-aware detection."),''',
    ''' ("Detection thresholds (AUD-12)", "Detection frequency of clinically important genes depends strongly on the threshold; identity and MAPQ were not deposited for DS02.", "High-confidence detection is the primary result; any-hit detection is an unconfirmed upper bound.", "Assembly-based or variant-aware detection; deposit identity/MAPQ."),
 ("Sequencing depth (AUD-16)", "231 included samples have < 10 Mbp and 642 < 100 Mbp; they enter descriptive denominators only.", "Depth-threshold sensitivity analysis (notebook 19) shows which results change.", "Set a minimum depth for future reprocessing."),''', "AUD-16")
files[f10] = files[f10].replace("matrix confounded with study; sites unidentified", "matrix confounded with study; city inference rests on 19 DS02 sites")

# ============================================================== AUD-15: final-size figure layout
FIGSIZE = {
 "01_Data_Inspection_and_QC.py": [("figsize=(9, 5.2)", "figsize=(6.9, 4.6)")],
 "02_Descriptive_Statistical_Analysis.py": [('figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.25, 1]}', 'figsize=(6.9, 3.4), gridspec_kw={"width_ratios": [1.25, 1]}'), ("figsize=(9, 3.8)", "figsize=(6.9, 3.0)"), ("figsize=(11, 4.2))", "figsize=(6.9, 3.2))")],
 "03_Taxonomic_and_Diversity_Analysis.py": [('figsize=(11, 4.4), gridspec_kw={"width_ratios": [1, 1.1]}', 'figsize=(6.9, 3.5), gridspec_kw={"width_ratios": [1, 1.1]}'), ("figsize=(12, 3.9)", "figsize=(6.9, 3.0)"), ("figsize=(10, 6)", "figsize=(6.9, 4.8)"), ("figsize=(11, 4.2)", "figsize=(6.9, 3.2)")],
 "04_AMR_Gene_Analysis.py": [('figsize=(11.5, 4.4), gridspec_kw={"width_ratios": [1.1, 1]}', 'figsize=(6.9, 3.6), gridspec_kw={"width_ratios": [1.1, 1]}'), ("figsize=(10.5, 8)", "figsize=(6.9, 8.4)"), ("figsize=(7, 5)", "figsize=(6.9, 5.0)"), ("figsize=(10, 3.8)", "figsize=(6.9, 2.9)"), ('figsize=(12, 4.6), gridspec_kw={"width_ratios": [1, 1.25]}', 'figsize=(6.9, 3.6), gridspec_kw={"width_ratios": [1, 1.25]}'), ("figsize=(7.5, 5.2)", "figsize=(6.9, 5.0)")],
 "05_Comparative_and_Multivariate_Analysis.py": [("figsize=(9, 6.5)", "figsize=(6.9, 5.4)"), ("figsize=(9, 7.5)", "figsize=(6.9, 6.2)"), ("figsize=(8.5, 7)", "figsize=(6.9, 7.2)"), ("figsize=(11, 4.3)", "figsize=(6.9, 3.2)")],
 "06_PCA_and_Ordination.py": [("figsize=(12, 4.4)", "figsize=(6.9, 3.0)"), ('figsize=(10, 4.3), gridspec_kw={"width_ratios": [1.4, .8]}', 'figsize=(6.9, 3.3), gridspec_kw={"width_ratios": [1.4, .8]}'), ("figsize=(15, 4.6)", "figsize=(6.9, 2.9)")],
 "07_Spatial_Analysis.py": [("figsize=(11, 10)", "figsize=(6.9, 6.8)"), ('figsize=(12, 4.2), gridspec_kw={"width_ratios": [.9, 1.3]}', 'figsize=(6.9, 3.0), gridspec_kw={"width_ratios": [.9, 1.3]}'), ("figsize=(11, 4))", "figsize=(6.9, 3.0))")],
 "08_Predictive_Modeling.py": [("figsize=(14, 4.2)", "figsize=(6.9, 2.9)")],
 "09_Sensitivity_and_Validation.py": [("figsize=(14, 6.2)", "figsize=(6.9, 5.6)")],
 "10_Final_Results_and_Visualization.py": [("figsize=(14, 8.5)", "figsize=(6.9, 6.4)")],
}
for fn, pairs in FIGSIZE.items():
    for old, new in pairs:
        rep(fn, old, new, "AUD-15", count=max(1, files[fn].count(old)))
files["07_Spatial_Analysis.py"] = files["07_Spatial_Analysis.py"].replace("figsize=(12, 6)); basemap(ax)", "figsize=(6.9, 3.6)); basemap(ax)")
# Fig 27: 1 x 4 panels -> 2 x 2 at double-column width
rep(f8, 'fig, axes = plt.subplots(1, 4, figsize=(16, 3.9))', 'fig, axes = plt.subplots(2, 2, figsize=(6.9, 5.8)); axes = axes.ravel()', "AUD-15")
# Dashboard typography
rep(f10, 'fontsize=24, weight="bold"', 'fontsize=15, weight="bold"', "AUD-15")
rep(f10, 'fig.suptitle("Urban microbiome and AMR: verified analysis dashboard", fontsize=13, weight="bold")', 'fig.suptitle("Urban microbiome and AMR: verified analysis dashboard (updated)", fontsize=9.5, weight="bold")', "AUD-15")
rep(f10, 'ax.text(0, 1, "\\n".join(textwrap.fill(t, 70) for t in txt), va="top", fontsize=8.5, linespacing=1.5)', 'ax.text(0, 1, "\\n".join(textwrap.fill(t, 58) for t in txt), va="top", fontsize=6.3, linespacing=1.35)', "AUD-15")
# Fig 10 long family labels: wrap instead of truncating
rep(f4, 'ax.set_yticklabels([f"{f[:55]} [{grp_of.get(f)}]" for f in top30], fontsize=6.5)', 'ax.set_yticklabels([textwrap.fill(f"{f} [{grp_of.get(f)}]", 52) for f in top30], fontsize=5.8)', "AUD-15")
rep(f4, "import seaborn as sns\namr = load_csv", "import seaborn as sns, textwrap\namr = load_csv", "AUD-15")
# Scale explicit font sizes to final print size (base 7 pt; nothing below 6 pt)
FS = {"24": "15", "13": "9.5", "10": "8", "9.5": "8", "9": "7.5", "8.5": "7.2", "8": "7", "7.5": "6.8", "7": "6.5", "6.5": "6.2", "6.3": "6.0"}
for fn in files:
    files[fn] = re.sub(r"fontsize=(\d+(?:\.\d+)?)", lambda m: "fontsize=" + FS.get(m.group(1), m.group(1)), files[fn])

# ============================================================== AUD-19: terminology
WORD = [(r"\btotal ARG load\b", "total normalized ARG abundance"), (r"\bTotal ARG load\b", "Total normalized ARG abundance"),
        (r"\bARG load\b", "normalized ARG abundance"), (r"\bARG-load\b", "normalized-ARG-abundance"), (r"\bresistome load\b", "normalized resistome abundance")]
n_word = 0
for fn in files:
    for pat, new in WORD:
        files[fn], k = re.subn(pat, new, files[fn]); n_word += k
LOG.append(dict(issue="AUD-19", file="all notebook sources", old="ARG load (all variants)", new=f"{n_word} replacements"))
for fn, s in files.items():
    (SRC / fn).write_text(s, encoding="utf-8")
(Path(sys.argv[1]) / "logs").mkdir(exist_ok=True)
(Path(sys.argv[1]) / "logs" / "update_patch_log.json").write_text(json.dumps(LOG, indent=1), encoding="utf-8")
print(len(LOG), "patches applied;", n_word, "terminology replacements")
