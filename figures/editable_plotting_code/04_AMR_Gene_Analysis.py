# %% [markdown]
# # 04 · AMR Gene Analysis
# **Research objective:** Objective 2 — occurrence, abundance, diversity and patterns of antimicrobial resistance genes (ARGs) in urban metagenomes; Objective 5 (food-related contextual comparison).
#
# **Two pipelines, analysed separately (never pooled):**
#
# | | Pipeline A | Pipeline B |
# |---|---|---|
# | Data | DS02, Singh et al. 2026 (Indian open drains), 447 samples | Martiny et al. 2022 resource (Zenodo 6919377), 551 analytic units, 6 matrices |
# | Method | RGI bwt (Bowtie2) against CARD (release not stated) | KMA read mapping against ResFinder |
# | Raw unit | reads mapped per CARD allele | fragments aligned per ResFinder template |
# | Normalisation used here | RPKM: reads / (allele length kb) / (total sequenced reads, millions) | FPKM: fragments / (template length kb) / (rRNA-aligned bacterial fragments, millions), as defined by the resource's `bacterial_fragment` column |
#
# **Primary detection rule:** all deposited hits are retained for abundance (as in both source studies); a stricter rule (allele/template coverage ≥ 50 %) is used for presence/absence of clinically important genes and in the sensitivity notebook. Detection of a gene is not evidence of phenotypic resistance, expression, viability, host organism, transmission or clinical infection. Absence of a hit in a sample with no ARG table is *missing*, not zero.

# %%
#BOOTSTRAP
NB = "NB04"
nb_header(NB)
import seaborn as sns, textwrap
amr = load_csv("amr_genes_long.csv")
m2 = pd.read_csv(DER / "ds02_meta.csv").set_index("Sample_ID")

# %% [markdown]
# ## A. CARD/RGI — DS02 Indian open-drain wastewater
# ### A1. Allele-level normalisation and aggregation

# %%
card = amr[amr.AMR_Reference_Database == "CARD"].copy()
card["reads"] = card.Detection_Value.astype(float)
card["len_kb"] = card.Reference_Length_bp / 1000
card["reads_total_M"] = card.Sample_ID.map(m2.reads_total) / 1e6
card["RPKM"] = card.reads / card.len_kb / card.reads_total_M
card["Class_Group"] = card.Resistance_Class.map(lambda c: class_group(c, "CARD"))
card["Clinical"] = [clinical_label(g, f, "CARD") for g, f in zip(card.AMR_Gene, card.AMR_Gene_Family)]
card["cov50_only"] = card.Average_Percent_Coverage >= 50
card["completely_mapped"] = pd.to_numeric(card.Completely_Mapped_Reads, errors="coerce").fillna(0)
card["cov50"] = card.cov50_only & (card.completely_mapped >= 2)   # AUD-12 high-confidence tier (primary)
card["cov80"] = (card.Average_Percent_Coverage >= 80) & (card.completely_mapped >= 2)
print("alleles cov>=50%:", int(card.cov50_only.sum()), "| high-confidence (cov>=50% & >=2 complete reads):", int(card.cov50.sum()), "| cov>=80%:", int(card.cov80.sum()))
print(card.RPKM.describe()); print(card.Class_Group.value_counts())
fam = card.pivot_table(index="Sample_ID", columns="AMR_Gene_Family", values="RPKM", aggfunc="sum", fill_value=0).reindex(m2.index, fill_value=0)
fam50 = card[card.cov50].pivot_table(index="Sample_ID", columns="AMR_Gene_Family", values="RPKM", aggfunc="sum", fill_value=0).reindex(m2.index, fill_value=0)
cls = card.pivot_table(index="Sample_ID", columns="Class_Group", values="RPKM", aggfunc="sum", fill_value=0).reindex(m2.index, fill_value=0)
mech = card.assign(Mechanism=card.Resistance_Mechanism.where(~card.Resistance_Mechanism.str.contains(";"), "multiple mechanisms")).pivot_table(index="Sample_ID", columns="Mechanism", values="RPKM", aggfunc="sum", fill_value=0).reindex(m2.index, fill_value=0)
fam.to_csv(DER / "ds02_card_family_rpkm.csv"); fam50.to_csv(DER / "ds02_card_family_rpkm_cov50.csv"); cls.to_csv(DER / "ds02_card_class_rpkm.csv"); mech.to_csv(DER / "ds02_card_mechanism_rpkm.csv")
card[["Sample_ID", "AMR_Gene", "AMR_Gene_Family", "Class_Group", "Resistance_Mechanism", "Clinical", "reads", "Reference_Length_bp", "Average_Percent_Coverage", "RPKM"]].to_csv(DER / "ds02_card_allele_long.csv", index=False)
s2 = pd.DataFrame({"total_RPKM": fam.sum(axis=1), "total_RPKM_cov50": fam50.sum(axis=1),
                   "family_richness_all": (fam > 0).sum(axis=1), "family_richness_cov50": (fam50 > 0).sum(axis=1),
                   "allele_richness_cov50": card[card.cov50].groupby("Sample_ID").AMR_Gene.nunique().reindex(m2.index, fill_value=0),
                   "family_shannon": [shannon(r) for r in fam.values], "total_reads_mapped": card.groupby("Sample_ID").reads.sum().reindex(m2.index)}, index=m2.index)
s2["bact_norm_total"] = card.groupby("Sample_ID").reads.sum() / (m2.bracken_phylum_reads / 1e6)  # sensitivity denominator
s2 = s2.join(m2[["City", "Season", "year", "month", "Collection_Date", "reads_total"]])
s2["log10_total_RPKM"] = np.log10(s2.total_RPKM)
s2.to_csv(DER / "ds02_arg_summary.csv")
s2.describe().T.round(3)

# %% [markdown]
# ### A2. Dominant gene families, classes and mechanisms

# %%
prev_all = (fam > 0).mean(); prev50 = (fam50 > 0).mean()
fam_tab = pd.DataFrame({"mean_RPKM": fam.mean(), "median_RPKM": fam.median(), "sd_RPKM": fam.std(),
                        "share_of_total_pct": 100 * fam.sum() / fam.values.sum(),
                        "detected_any_hit_n": (fam > 0).sum(), "detected_cov50_n": (fam50 > 0).sum()}).sort_values("mean_RPKM", ascending=False)
ci = [wilson(k, len(fam)) for k in fam_tab.detected_cov50_n]
fam_tab["prevalence_cov50_pct"] = 100 * fam_tab.detected_cov50_n / len(fam)
fam_tab["prev_cov50_ci_low"] = [100 * a for a, b in ci]; fam_tab["prev_cov50_ci_high"] = [100 * b for a, b in ci]
grp_of = card.drop_duplicates("AMR_Gene_Family").set_index("AMR_Gene_Family").Class_Group
fam_tab.insert(0, "Class_Group", fam_tab.index.map(grp_of)); fam_tab.index.name = "AMR_Gene_Family"
save_table(fam_tab.reset_index().pipe(rnd, 4), "A01_ds02_card_gene_families", "DS02 CARD gene families: abundance (RPKM) and prevalence (coverage >= 50%, Wilson 95% CI)", NB, sheet="AMR_Gene_Summary")
cl_tab = pd.DataFrame({"share_of_total_RPKM_pct": 100 * cls.sum() / cls.values.sum(), "median_RPKM": cls.median(), "mean_RPKM": cls.mean()}).sort_values("share_of_total_RPKM_pct", ascending=False)
cl_city = (100 * cls.div(cls.sum(axis=1), axis=0)).groupby(m2.City).mean().T
cl_tab = cl_tab.join(cl_city.add_prefix("mean_share_pct_"))
save_table(cl_tab.reset_index().rename(columns={"index": "Class_Group"}).pipe(rnd, 3), "A02_ds02_card_class_groups", "DS02 CARD resistance-class groups: share of total RPKM, overall and by city", NB, sheet="AMR_Resistance_Classes")
me_tab = pd.DataFrame({"share_of_total_RPKM_pct": 100 * mech.sum() / mech.values.sum()}).join((100 * mech.div(mech.sum(axis=1), axis=0)).groupby(m2.City).mean().T.add_prefix("mean_share_pct_")).sort_values("share_of_total_RPKM_pct", ascending=False)
save_table(me_tab.reset_index().rename(columns={"index": "Resistance_Mechanism"}).pipe(rnd, 3), "A03_ds02_card_mechanisms", "DS02 CARD resistance mechanisms: share of total RPKM, overall and by city", NB, sheet="AMR_Resistance_Classes")
fam_tab.head(15).round(2)

# %% [markdown]
# ### A3. City and season comparisons of total normalized ARG abundance and diversity

# %%
import statsmodels.formula.api as smf, statsmodels.api as sm_api
s2c = s2.assign(City=s2.City.astype(str), Season=s2.Season.astype(str))
summ_city = s2c.groupby("City").agg(n=("total_RPKM", "size"), median_total_RPKM=("total_RPKM", "median"),
        Q1=("total_RPKM", lambda s: s.quantile(.25)), Q3=("total_RPKM", lambda s: s.quantile(.75)),
        median_family_richness_cov50=("family_richness_cov50", "median"), median_family_shannon=("family_shannon", "median"))
summ_city["median_CI95"] = [str(tuple(float(round(v, 1)) for v in boot_ci(s2c.loc[s2c.City == c, "total_RPKM"]))) for c in summ_city.index]
save_table(summ_city.reset_index().pipe(rnd, 3), "A04_ds02_arg_load_by_city", "DS02 total normalized ARG abundance (RPKM) and family richness by city (bootstrap 95% CI of median)", NB, sheet="AMR_Gene_Summary")
kt = []
for v in ["log10_total_RPKM", "family_richness_cov50", "family_shannon"]:
    for f in ["City", "Season"]:
        k = kruskal_eps2(s2c, v, f); kt.append(dict(variable=v, factor=f, **k))
        record_test(NB, f"DS02 {v} by {f}", "Kruskal-Wallis", k["n"], k["H"], k["df"], k["p"], "epsilon-squared", k["epsilon2"], family="DS02 ARG abundance/diversity [sample-level, descriptive only: pseudo-replicated]", note="Superseded for inference by site-aware tests (A05s); AUD-11", objective="Objective 2")
kt = pd.DataFrame(kt); kt["p_BH"] = bh(kt.p)
save_table(kt.pipe(rnd, 6), "A05_ds02_arg_kruskal", "DESCRIPTIVE (sample-level, pseudo-replicated - see A05s for inference): Kruskal-Wallis of DS02 normalized ARG abundance and diversity by city and season", NB, sheet="Statistical_Tests")
pw = pairwise_mwu(s2c, "total_RPKM", "City")
save_table(pw.pipe(rnd, 6), "A06_ds02_arg_pairwise_city", "Pairwise Mann-Whitney tests of total ARG RPKM between cities (Holm)", NB, sheet="Statistical_Tests")
ols = smf.ols("log10_total_RPKM ~ C(City, Treatment('Delhi')) + C(Season, Treatment('Winter (Dec-Feb)')) + C(year)", s2c).fit(cov_type="HC3")
co = pd.DataFrame({"estimate_log10": ols.params, "SE_HC3": ols.bse, "ci_low": ols.conf_int()[0], "ci_high": ols.conf_int()[1], "p_value": ols.pvalues})
co["fold_change"] = 10 ** co.estimate_log10; co["fc_ci_low"] = 10 ** co.ci_low; co["fc_ci_high"] = 10 ** co.ci_high
save_table(co.reset_index().rename(columns={"index": "term"}).pipe(rnd, 5), "A07_ds02_arg_load_ols", "OLS of log10 total ARG RPKM on city, season and year (HC3; reference Delhi, winter, 2022); fold changes", NB, sheet="Comparative_Analysis")
an = sm_api.stats.anova_lm(smf.ols("log10_total_RPKM ~ C(City) + C(Season) + C(year)", s2c).fit(), typ=2)
an["eta2_partial"] = an.sum_sq / (an.sum_sq + an.loc["Residual", "sum_sq"])
for t in ["C(City)", "C(Season)", "C(year)"]:
    record_test(NB, "DS02 log10 total ARG RPKM ~ City + Season + Year (type II)", f"OLS F-test: {t}", len(s2c), an.loc[t, "F"], f"{int(an.loc[t,'df'])}, {int(an.loc['Residual','df'])}", an.loc[t, "PR(>F)"], "partial eta-squared", an.loc[t, "eta2_partial"], family="DS02 ARG models [sample-level OLS, descriptive only: pseudo-replicated]", objective="Objective 2/3")
save_table(an.reset_index().rename(columns={"index": "term"}).pipe(rnd, 6), "A08_ds02_arg_load_anova", "Type II ANOVA: log10 total ARG RPKM ~ city + season + year (DS02)", NB, sheet="Statistical_Tests")
print(kt.round(4)); co.round(3)

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

# %%
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.6), gridspec_kw={"width_ratios": [1.1, 1]})
cities = list(CITY_COLORS)
sh = (100 * cls.div(cls.sum(axis=1), axis=0)).groupby(m2.City).mean().reindex(cities)
sh = sh[[c for c in CLASS_ORDER if c in sh.columns]]
bottom = np.zeros(len(sh))
for c in sh.columns:
    axes[0].bar(range(len(sh)), sh[c], bottom=bottom, color=CLASS_COLORS[c], label=c, edgecolor="white", linewidth=.3); bottom += sh[c].values
axes[0].set_xticks(range(len(sh))); axes[0].set_xticklabels([f"{c}\n(n={int((m2.City==c).sum())})" for c in cities]); axes[0].set_ylabel("Mean share of total ARG RPKM (%)"); axes[0].set_ylim(0, 100)
axes[0].legend(frameon=False, fontsize=6.2, bbox_to_anchor=(1.0, 1), loc="upper left"); axes[0].set_title("A  Resistance-class composition (CARD/RGI)")
sns.boxplot(data=s2c, x="City", y="total_RPKM", order=cities, hue="City", palette=CITY_COLORS, legend=False, showfliers=False, linewidth=.7, ax=axes[1], log_scale=True)
sns.stripplot(data=s2c, x="City", y="total_RPKM", order=cities, color="k", size=1.6, alpha=.35, jitter=.25, ax=axes[1])
k = kt_site[(kt_site.variable == "log10_total_RPKM") & (kt_site.factor == "City")].iloc[0]
axes[1].set_ylabel("Total normalized ARG abundance (RPKM, log scale)"); axes[1].set_xlabel("")
axes[1].set_title(f"B  Normalized ARG abundance by city\nsite-level KW ({int(k.n_units)} sites): H={k.H:.1f}, perm. p={k.p_perm:.3f}", fontsize=7.2)
fig.tight_layout()
save_fig(fig, 9, "ds02_resistome_class_and_load", "Resistome composition and total normalized ARG abundance in Indian open-drain wastewater (DS02, CARD/RGI)",
         "(A) Mean share of total ARG abundance contributed by each harmonised resistance-class group (genes spanning several groups = 'Multi-class'). (B) Total ARG abundance per sample (RPKM: reads per kb of CARD reference per million sequenced reads; all deposited RGI bwt hits), log scale; statistic = Kruskal-Wallis on 19 site medians with whole-site permutation p (Table A05s). n = 447 samples from 19 sites. RPKM is a depth-normalized relative abundance, not an absolute quantity per volume or per cell. ARG detection does not establish phenotypic resistance or human exposure.",
         NB, "Objective 2", "tables/A02_ds02_card_class_groups.csv; derived/ds02_arg_summary.csv")

# %%
top30 = fam_tab.index[:30]
hm = np.log10(fam[top30].groupby([m2.City, m2.Season]).mean() + 1e-3)
hm = hm.reindex(pd.MultiIndex.from_product([cities, SEASON_ORDER])).dropna(how="all")
fig, ax = plt.subplots(figsize=(6.9, 8.4))
sns.heatmap(hm.T, cmap="viridis", cbar_kws={"label": "log10 mean RPKM"}, linewidths=.25, ax=ax)
ax.set_xticklabels([f"{c} · {s.split(' (')[0]}" for c, s in hm.index], rotation=60, ha="right", fontsize=6.5)
ax.set_yticklabels([textwrap.fill(f"{f} [{grp_of.get(f)}]", 52) for f in top30], fontsize=5.8); ax.set_xlabel(""); ax.set_ylabel("")
ax.set_title("Thirty most abundant CARD gene families by city and season (DS02)")
save_fig(fig, 10, "ds02_gene_family_heatmap", "Abundance of the 30 most abundant CARD gene families by city and season (DS02)",
         "log10 mean RPKM of the 30 CARD AMR gene families with the highest mean abundance, for each city-season cell (n = 447 samples). Bracketed labels give the harmonised class group. Efflux-pump families can include genes with housekeeping functions; abundance is not resistance phenotype.",
         NB, "Objective 2", "derived/ds02_card_family_rpkm.csv; tables/A01_ds02_card_gene_families.csv")

# %% [markdown]
# ### A4. Clinically important resistance determinants (detection frequency)
# **UPDATE (AUD-12).** The deposited RGI bwt table provides, per CARD allele, all mapped reads, completely mapped reads and average percent coverage; it does not provide percent identity, MAPQ or the CARD release, and RGI bwt documents no detection threshold. Detection is therefore reported in three tiers: **(1) high-confidence (primary): at least one allele with ≥ 50 % reference coverage and ≥ 2 completely mapped reads;** (2) a stricter sensitivity tier (≥ 80 % coverage); (3) a lenient tier (any deposited read hit) reported only as an **unconfirmed upper bound** (cross-mapping from homologous genes cannot be excluded). Any-hit detection is never described as confirmed resistance. Presence is defined at allele coverage ≥ 50 % (primary) and at any deposited hit (lenient). Wilson 95 % confidence intervals are reported. Panel: NDM, KPC, OXA-48-like, VIM, IMP and CTX-M β-lactamases; MCR; tet(X); qnr; vanA; sul.

# %%
def panel_detect(df, unit_col, group_series, label_col="Clinical", strict_col="cov50", pipeline=""):
    rows = []
    units = group_series.index
    for lab in CLINICAL_PANEL:
        d = df[df[label_col] == lab]
        anyhit = set(d[unit_col]); strict = set(d.loc[d[strict_col], unit_col])
        for g in ["All"] + sorted(group_series.dropna().unique()):
            u = units if g == "All" else group_series.index[group_series == g]
            n = len(u); k1 = len(set(u) & strict); k2 = len(set(u) & anyhit)
            lo, hi = wilson(k1, n)
            rows.append(dict(pipeline=pipeline, determinant=lab, group=g, n=n, detected_cov50=k1, pct_cov50=100 * k1 / n if n else np.nan,
                             ci_low=100 * lo, ci_high=100 * hi, detected_any_hit=k2, pct_any_hit=100 * k2 / n if n else np.nan))
    return pd.DataFrame(rows)
clin2 = panel_detect(card, "Sample_ID", m2.City.astype(str), pipeline="CARD/RGI (DS02)")
clin2_80 = panel_detect(card, "Sample_ID", m2.City.astype(str), strict_col="cov80", pipeline="CARD/RGI (DS02)")
clin2["detected_cov80"] = clin2_80.detected_cov50.values; clin2["pct_cov80"] = clin2_80.pct_cov50.values
clin2 = clin2.rename(columns={"detected_cov50": "detected_high_confidence", "pct_cov50": "pct_high_confidence", "ci_low": "hc_ci_low", "ci_high": "hc_ci_high",
                              "detected_any_hit": "detected_any_hit_unconfirmed", "pct_any_hit": "pct_any_hit_unconfirmed"})
clin2["uncertainty_range_pp"] = clin2.pct_any_hit_unconfirmed - clin2.pct_high_confidence
# keep the original column names used by the plotting helper
clin2_plot = clin2.rename(columns={"detected_high_confidence": "detected_cov50", "pct_high_confidence": "pct_cov50", "hc_ci_low": "ci_low", "hc_ci_high": "ci_high", "pct_any_hit_unconfirmed": "pct_any_hit"})
save_table(clin2.pipe(rnd, 2), "A09_ds02_clinical_determinants", "DS02 clinically important ARG determinants (AUD-12): high-confidence detection (primary; allele coverage >= 50% and >= 2 completely mapped reads; Wilson 95% CI), stricter >= 80% coverage tier, and any-hit detection (unconfirmed upper bound, sensitivity only)", NB, sheet="AMR_Gene_Summary")
chi = []
from scipy.stats import chi2_contingency, fisher_exact
for lab in CLINICAL_PANEL:
    t = clin2_plot[(clin2_plot.determinant == lab) & (clin2_plot.group != "All")]
    tab = np.column_stack([t.detected_cov50, t.n - t.detected_cov50])
    if tab[:, 0].sum() == 0 or tab[:, 1].sum() == 0: chi.append(dict(determinant=lab, chi2=np.nan, df=np.nan, p=np.nan, cramers_v=np.nan, note="constant (all or none detected)")); continue
    c2, p, dfree, exp = chi2_contingency(tab)
    v = np.sqrt(c2 / (tab.sum() * (min(tab.shape) - 1)))
    chi.append(dict(determinant=lab, chi2=c2, df=dfree, p=p, cramers_v=v, min_expected=exp.min(), note="chi-square (expected >=5)" if exp.min() >= 5 else "chi-square; some expected counts <5, interpret with caution"))
chi = pd.DataFrame(chi); chi["p_BH"] = bh(chi.p)
for _, r in chi.dropna(subset=["p"]).iterrows():
    record_test(NB, f"DS02 detection of {r.determinant} (high-confidence) by city [sample level]", "Chi-square test of independence", len(m2), r.chi2, int(r.df), r.p, "Cramer's V", r.cramers_v, family="DS02 clinical determinants [sample-level, descriptive only: pseudo-replicated]", note=r.note + "; superseded by A10s (AUD-11)", objective="Objective 2/3")
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
save_table(chs.pipe(rnd, 6), "A10s_ds02_clinical_by_city_site_aware", "PRIMARY (AUD-11/AUD-12): city differences in high-confidence detection of clinically important determinants, tested by permuting city labels among the 19 whole sites; site-level prevalence summaries", NB, sheet="Statistical_Tests")
clin2[clin2.group == "All"].round(1)

# %%
def panel_plot(ax, tab, groups, colors, title):
    labs = CLINICAL_PANEL; y = np.arange(len(labs)); off = np.linspace(-.3, .3, len(groups))
    for o, g in zip(off, groups):
        t = tab[tab.group == g].set_index("determinant").reindex(labs)
        ax.hlines(y + o, t.pct_cov50, t.pct_any_hit, color=colors.get(g, "k"), lw=.9, alpha=.6)
        ax.errorbar(t.pct_cov50, y + o, xerr=[(t.pct_cov50 - t.ci_low).clip(lower=0), (t.ci_high - t.pct_cov50).clip(lower=0)], fmt="o", ms=3.5, lw=.8, capsize=1.5, color=colors.get(g, "k"), label=f"{g} (n={int(t.n.iloc[0])})")
        ax.plot(t.pct_any_hit, y + o, "o", ms=3.5, mfc="white", mec=colors.get(g, "k"), mew=.8)
    ax.set_yticks(y); ax.set_yticklabels(labs); ax.invert_yaxis(); ax.set_xlim(-2, 102)
    ax.set_xlabel("Units with determinant detected (%)\nfilled = high-confidence (primary; Wilson 95% CI); open = any hit (unconfirmed upper bound)", fontsize=6.8)
    ax.set_title(title, fontsize=7.5); ax.legend(frameon=False, fontsize=6.2, loc="upper center", bbox_to_anchor=(.5, -.13), ncol=3)
fig, ax = plt.subplots(figsize=(6.9, 5.0))
panel_plot(ax, clin2_plot, cities, CITY_COLORS, "Clinically important determinants, CARD/RGI (DS02)")
save_fig(fig, 11, "ds02_clinical_determinants", "Detection frequency of clinically important ARG determinants by city (DS02, CARD/RGI)",
         "Percentage of samples in which each determinant was detected. Filled points (primary, high-confidence): at least one CARD allele with ≥ 50% reference coverage and ≥ 2 completely mapped reads (Wilson 95% CI). Open points (sensitivity only): any deposited RGI bwt read hit, an unconfirmed upper bound that is vulnerable to cross-mapping from homologous genes. The high-confidence rule may underestimate multi-allele families because RGI bwt splits reads among near-identical alleles. Percent identity, MAPQ and the CARD release were not deposited. Detection reflects DNA presence in wastewater, not carriage, infection or phenotypic resistance.",
         NB, "Objective 2/3", "tables/A09_ds02_clinical_determinants.csv")

# %%
mon = s2c.groupby(["City", "Collection_Date"]).total_RPKM.agg(["median", lambda s: s.quantile(.25), lambda s: s.quantile(.75), "size"]).reset_index()
mon.columns = ["City", "Month", "median_RPKM", "Q1", "Q3", "n"]
save_table(mon.pipe(rnd, 3), "A11_ds02_monthly_arg_load", "DS02 monthly total normalized ARG abundance by city (median, IQR)", NB, sheet="Environmental_Analysis")
fig, ax = plt.subplots(figsize=(6.9, 2.9))
for c in cities:
    t = mon[mon.City == c]; x = pd.to_datetime(t.Month + "-15")
    ax.plot(x, t.median_RPKM, "-o", ms=2.5, lw=1, color=CITY_COLORS[c], label=c); ax.fill_between(x, t.Q1, t.Q3, color=CITY_COLORS[c], alpha=.12)
ax.set_yscale("log"); ax.set_ylabel("Total ARG RPKM (median, IQR)"); ax.legend(frameon=False, ncol=4, fontsize=6.8)
ax.set_title("Monthly total normalized ARG abundance in open-drain wastewater, March 2022 - March 2024 (DS02)")
save_fig(fig, 12, "ds02_monthly_arg_load", "Monthly total normalized ARG abundance by city (DS02 longitudinal surveillance)",
         "Median (line) and interquartile range (band) of total ARG RPKM across the 2-8 sites sampled per city each month. The series is observed, not modelled or forecast; site identities are not available, so between-month changes may partly reflect site composition.",
         NB, "Objective 3/4 (monitoring)", "tables/A11_ds02_monthly_arg_load.csv")

# %% [markdown]
# ## B. ResFinder/KMA — Martiny et al. global resistome resource
# ### B1. Aggregation of runs to analytic units and FPKM normalisation

# %%
rmap = pd.read_csv(DER / "rf_run_to_unit.csv")
arg = pd.read_csv(HARM / "subset_6919377_ARG_core_runs.tsv", sep="\t").merge(rmap[["Run_ID", "unit_id", "primary"]], left_on="run_accession", right_on="Run_ID")
cls_rf = amr[amr.AMR_Reference_Database == "ResFinder"].drop_duplicates("Gene_Identifier").set_index("Gene_Identifier").Resistance_Class
arg["gene"] = arg.refSequence.map(rf_gene_name); arg["allele"] = arg.refSequence.str.split("_").str[0]
arg["Class_Group"] = arg.refSequence.map(cls_rf).map(lambda c: class_group(c, "ResFinder"))
arg["Clinical"] = [clinical_label(a, "", "ResFinder") for a in arg.allele]
arg["run_cov"] = 100 * arg.refCoveredPositions / arg.refSequence_length
runden = arg.drop_duplicates("run_accession")[["run_accession", "unit_id", "primary", "bacterial_fragment", "trimmed_fragments"]]
def build_units(a, rd):
    den = rd.groupby("unit_id")[["bacterial_fragment", "trimmed_fragments"]].sum()
    t = a.groupby(["unit_id", "refSequence", "gene", "allele", "Class_Group", "Clinical"], dropna=False, as_index=False).agg(
        fragments=("fragmentCountAln", "sum"), length=("refSequence_length", "first"), max_run_cov=("run_cov", "max"))
    t = t.join(den, on="unit_id")
    t["FPKM"] = t.fragments / (t.length / 1000) / (t.bacterial_fragment / 1e6)
    t["FPKM_trimmed"] = t.fragments / (t.length / 1000) / (t.trimmed_fragments / 1e6)
    return t
rfl = build_units(arg, runden)
rfl.to_csv(DER / "rf_template_long_all_units.csv", index=False)
umeta_cols = ["Dataset_ID", "Project_Accession", "Sample_Type", "Inclusion_Status"]
ru = rmap.groupby("unit_id").agg(n_runs=("Run_ID", "size"), n_accessions=("Sample_ID", "nunique"), Sample_ID_first=("Sample_ID", "first"), primary=("primary", "all"),
        has_ARG=("has_ARG_rows", "any"), kma_versions=("kma_version", lambda s: ";".join(sorted(s.dropna().unique()))), **{c: (c, "first") for c in umeta_cols})
meta = pd.read_csv(DER / "metadata_set.csv", low_memory=False).set_index("Sample_ID")
geo = rmap.merge(meta[["City", "Country", "ISO2", "Continent", "Subregion", "lat", "lon", "year"]], left_on="Sample_ID", right_index=True, how="left").groupby("unit_id").agg(
        City=("City", "first"), Country=("Country", "first"), ISO2=("ISO2", "first"), Continent=("Continent", "first"), Subregion=("Subregion", "first"),
        lat=("lat", "mean"), lon=("lon", "mean"), lat_sd=("lat", "std"), year=("year", "min"))
ru = ru.join(geo).join(runden.groupby("unit_id")[["bacterial_fragment", "trimmed_fragments"]].sum())
ru["Stratum"] = np.where(ru.Inclusion_Status == "Included", "Urban environmental", "Food-related (contextual)")
ru["kma_major"] = ru.kma_versions.str.extract(r"(\d+\.\d+)")[0]
ru.to_csv(DER / "rf_units_meta.csv")
P = rfl[rfl.unit_id.isin(ru.index[ru.primary & ru.has_ARG])]
gfp = P.pivot_table(index="unit_id", columns="gene", values="FPKM", aggfunc="sum", fill_value=0)
gfr = P.pivot_table(index="unit_id", columns="gene", values="fragments", aggfunc="sum", fill_value=0)
cfp = P.pivot_table(index="unit_id", columns="Class_Group", values="FPKM", aggfunc="sum", fill_value=0)
gfp.to_csv(DER / "rf_gene_fpkm.csv"); gfr.to_csv(DER / "rf_gene_fragments.csv"); cfp.to_csv(DER / "rf_class_fpkm.csv")
us = pd.DataFrame({"total_FPKM": gfp.sum(axis=1), "gene_richness": (gfp > 0).sum(axis=1),
                   "template_richness": P.groupby("unit_id").refSequence.nunique(), "gene_shannon_fragments": [shannon(r) for r in gfr.loc[gfp.index].values],
                   "total_FPKM_trimmed": P.groupby("unit_id").FPKM_trimmed.sum()}).join(ru)
us["log10_total_FPKM"] = np.log10(us.total_FPKM); us["log10_bacterial_fragments"] = np.log10(us.bacterial_fragment)
rr = pd.read_csv(DER / "rf_unit_rrna_diversity.csv").set_index("unit_id")[["genus_richness", "genus_shannon", "genus_simpson", "phylum_richness", "phylum_shannon", "total_fragments"]]
us = us.join(rr)
us.to_csv(DER / "rf_unit_summary.csv")
print("primary units with ARG data:", len(us)); print(us.groupby(["Stratum", "Sample_Type"]).size())
record_exclusion(NB, "ResFinder resistome analyses", "unit", int((ru.primary & ~ru.has_ARG).sum()), "Primary unit without ARG rows (missing, not zero)")

# %% [markdown]
# ### B2. normalized ARG abundance, richness and composition by matrix

# %%
rf_type = us.groupby(["Stratum", "Sample_Type"]).agg(units=("total_FPKM", "size"), studies=("Project_Accession", "nunique"),
        countries=("Country", lambda s: s[~is_missing(s)].nunique()), median_total_FPKM=("total_FPKM", "median"),
        Q1=("total_FPKM", lambda s: s.quantile(.25)), Q3=("total_FPKM", lambda s: s.quantile(.75)), mean_log10_FPKM=("log10_total_FPKM", "mean"),
        sd_log10_FPKM=("log10_total_FPKM", "std"), median_gene_richness=("gene_richness", "median"), median_bacterial_fragments=("bacterial_fragment", "median")).reset_index()
rf_type["median_CI95"] = [str(tuple(float(round(v, 0)) for v in boot_ci(us.loc[(us.Stratum == a) & (us.Sample_Type == b), "total_FPKM"]))) for a, b in zip(rf_type.Stratum, rf_type.Sample_Type)]
save_table(rf_type.pipe(rnd, 3), "A12_rf_arg_load_by_sample_type", "ResFinder/KMA: total ARG FPKM and gene richness by sample type (bootstrap 95% CI of median)", NB, sheet="AMR_Gene_Summary")
gtab = pd.DataFrame({"Class_Group": P.drop_duplicates("gene").set_index("gene").Class_Group.reindex(gfp.columns), "mean_FPKM_all_units": gfp.mean(),
                     "prevalence_any_hit_pct": 100 * (gfp > 0).mean()})
for st, ix in us.groupby("Sample_Type").groups.items():
    gtab[f"prev_pct_{st}"] = 100 * (gfp.loc[ix] > 0).mean(); gtab[f"mean_FPKM_{st}"] = gfp.loc[ix].mean()
gtab = gtab.sort_values("mean_FPKM_all_units", ascending=False); gtab.index.name = "gene"
save_table(gtab.reset_index().pipe(rnd, 3), "A13_rf_gene_summary", "ResFinder genes (allele-collapsed for bla/mcr): mean FPKM and prevalence by sample type", NB, sheet="AMR_Gene_Summary")
csh = (100 * cfp.div(cfp.sum(axis=1), axis=0)).join(us.Sample_Type).groupby("Sample_Type").mean().T
save_table(csh.pipe(rnd, 3).reset_index().rename(columns={"index": "Class_Group"}), "A14_rf_class_share_by_sample_type", "ResFinder: mean share (%) of total ARG FPKM by class group and sample type", NB, sheet="AMR_Resistance_Classes")
clin_rf = panel_detect(P.assign(cov50=P.max_run_cov >= 50), "unit_id", us.Sample_Type, pipeline="ResFinder/KMA (Martiny et al.)")
save_table(clin_rf.pipe(rnd, 2), "A15_rf_clinical_determinants", "ResFinder: detection of clinically important determinants by sample type (Wilson 95% CI)", NB, sheet="AMR_Gene_Summary")
rf_type

# %%
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.6), gridspec_kw={"width_ratios": [1, 1.25]})
order = ["Municipal wastewater", "Surface water", "Air", "Built environment", "Transit surface", "Food-related environment"]
order = [o for o in order if o in us.Sample_Type.unique()]
sns.boxplot(data=us, y="Sample_Type", x="total_FPKM", order=order, hue="Sample_Type", palette=SAMPLE_TYPE_COLORS, legend=False, showfliers=False, linewidth=.7, ax=axes[0], log_scale=True)
sns.stripplot(data=us, y="Sample_Type", x="total_FPKM", order=order, color="k", size=1.7, alpha=.35, jitter=.25, ax=axes[0])
nn = us.Sample_Type.value_counts(); sn = us.groupby("Sample_Type").Project_Accession.nunique()
axes[0].set_yticks(range(len(order))); axes[0].set_yticklabels([f"{t}{' (contextual)' if t.startswith('Food') else ''}\nunits={nn[t]}, studies={sn[t]}" for t in order], fontsize=6.8)
axes[0].axhline(len(order) - 1.5, color="grey", lw=.6, ls="--"); axes[0].set_ylabel(""); axes[0].set_xlabel("Total ARG abundance (FPKM, log scale)"); axes[0].set_title("A  Total normalized ARG abundance (ResFinder/KMA)")
cs = csh[order]; cs = cs.loc[[c for c in CLASS_ORDER if c in cs.index]]
left = np.zeros(len(order))
for c in cs.index:
    axes[1].barh(range(len(order)), cs.loc[c].values, left=left, color=CLASS_COLORS[c], label=c, edgecolor="white", linewidth=.3); left += cs.loc[c].values
axes[1].set_yticks(range(len(order))); axes[1].set_yticklabels([]); axes[1].invert_yaxis(); axes[1].set_xlim(0, 100); axes[1].set_xlabel("Mean share of total ARG FPKM (%)")
axes[1].legend(frameon=False, fontsize=6.2, bbox_to_anchor=(1.0, 1), loc="upper left"); axes[1].set_title("B  Resistance-class composition")
fig.tight_layout()
save_fig(fig, 13, "rf_arg_load_and_classes_by_matrix", "Total normalized ARG abundance and resistance-class composition across urban matrices (ResFinder/KMA)",
         f"(A) Total ARG abundance per analytic unit (FPKM: fragments per kb of ResFinder template per million rRNA-aligned bacterial fragments) by matrix; food-related environments are a separate contextual dataset (below dashed line). (B) Mean share of total FPKM by harmonised class group. {len(us)} units from {us.Project_Accession.nunique()} studies. Units within a matrix often come from one or two studies, so matrix and study effects are confounded (see notebook 05). ARG abundance does not establish resistance phenotype, viability or human exposure.",
         NB, "Objective 2 / 5", "tables/A12_rf_arg_load_by_sample_type.csv; tables/A14_rf_class_share_by_sample_type.csv")

# %%
fig, ax = plt.subplots(figsize=(6.9, 5.0))
panel_plot(ax, clin_rf, order, SAMPLE_TYPE_COLORS, "Clinically important determinants, ResFinder/KMA")
save_fig(fig, 14, "rf_clinical_determinants", "Detection frequency of clinically important ARG determinants by matrix (ResFinder/KMA)",
         "Filled points: at least one ResFinder template detected with ≥ 50% coverage in at least one run (Wilson 95% CI); open points: any aligned fragment. Detection depends strongly on sequencing depth, which differs by matrix and study; frequencies are not adjusted and are not comparable with the CARD/RGI panel (Figure 11).",
         NB, "Objective 2 / 5", "tables/A15_rf_clinical_determinants.csv")

# %%
top_cls = cl_tab.index[:3].tolist()
kc = kt[(kt.variable == "log10_total_RPKM") & (kt.factor == "City")].iloc[0]
kcs = kt_site[(kt_site.variable == "log10_total_RPKM") & (kt_site.factor == "City")].iloc[0]
kcm = kt_site[(kt_site.variable == "log10_total_RPKM") & (kt_site.factor == "C(City)")].iloc[0]
ksm = kt_site[(kt_site.variable == "log10_total_RPKM") & (kt_site.factor == "C(Season)")].iloc[0]
fcs = co_s[co_s.term.str.contains("City")]
record_finding(NB, "F07", "Objective 2", f"In Indian open-drain wastewater (CARD/RGI), median total normalized ARG abundance was " + ", ".join(f"{c} {summ_city.loc[c, 'median_total_RPKM']:.0f}" for c in summ_city.index) + f" RPKM. With the site as the unit (19 sites), the city difference was {'statistically supported' if kcs.p_perm < .05 else 'not statistically supported'} (site-level Kruskal-Wallis H = {kcs.H:.1f}, whole-site permutation p = {kcs.p_perm:.3f}; mixed-model LRT p = {kcm.p:.2g}); mixed-model fold changes vs Delhi: " + "; ".join(f"{term_label(t)} {f:.2f} (95% CI {l:.2f}-{h:.2f})" for t, f, l, h in zip(fcs.term, fcs.fold_change, fcs.fc_ci_low, fcs.fc_ci_high)) + f"; season within sites: LRT p = {ksm.p:.2g}. Leading class groups by share of abundance: {', '.join(top_cls)}.",
               len(s2), f"DS02: {len(s2)} samples from {int(kcs.n_units)} sites", "Site-aware (AUD-11): site-level Kruskal-Wallis with whole-site permutation; site random-intercept mixed model (A05s, A07s). Sample-level p-values (A05-A08) are descriptive only.", "Tables A02, A04, A05s, A07s; Figure 9",
               f"City differences in normalized ARG abundance are estimated from 19 sites; the site ICC was {mm_arg['icc']:.2f}, so sample-level tests overstated certainty (descriptive sample-level epsilon-squared {kc.epsilon2:.2f}).",
               "RPKM is a depth-normalized relative abundance (not an absolute quantity); few sites per city; CARD release unstated; efflux families may include housekeeping genes.")
allp = clin2[(clin2.group == "All")].set_index("determinant")
sig_s = chs[chs.p_BH < .05].determinant.tolist()
record_finding(NB, "F08", "Objective 2/3", "High-confidence detection (primary; allele coverage >= 50% and >= 2 completely mapped reads) of clinically important determinants in DS02 wastewater samples: " + "; ".join(f"{d} {allp.loc[d,'pct_high_confidence']:.0f}% (95% CI {allp.loc[d,'hc_ci_low']:.0f}-{allp.loc[d,'hc_ci_high']:.0f})" for d in CLINICAL_PANEL) + ". Any-hit detection (unconfirmed upper bound; sensitivity only): " + "; ".join(f"{d} {allp.loc[d,'pct_any_hit_unconfirmed']:.0f}%" for d in CLINICAL_PANEL) + ".",
               len(m2), f"DS02: {len(m2)} samples from 19 sites", "Wilson 95% CIs (high-confidence tier); city differences by whole-site permutation chi-square (A10s; BH)", "Tables A09, A10s; Figure 11",
               f"High-confidence detection in >= 80% of samples: {', '.join(allp.index[allp.pct_high_confidence >= 80]) or 'none'}. Determinants with > 20 percentage-point uncertainty between tiers ({', '.join(allp.index[allp.uncertainty_range_pp > 20]) or 'none'}) must be cited at the high-confidence value, with the lenient value described as unconfirmed. City differences supported at the site level (BH q < 0.05): {', '.join(sig_s) or 'none'}.",
               "Read mapping only (no assembly, host, plasmid context, expression or viability); identity/MAPQ not deposited; detection is not phenotypic or clinical resistance.")
r = rf_type.set_index("Sample_Type")
rank = r.median_total_FPKM.sort_values(ascending=False).index.tolist()
record_finding(NB, "F09", "Objective 2/5", "In the ResFinder resource, median total normalized ARG abundance (FPKM, a relative measure per million rRNA-aligned bacterial fragments) was " + "; ".join(f"{t} {r.loc[t,'median_total_FPKM']:,.0f} FPKM (units={int(r.loc[t,'units'])}, studies={int(r.loc[t,'studies'])})" for t in order) + ".",
               len(us), "Primary ResFinder analytic units", "Descriptive (median, IQR, bootstrap CI); model-based contrasts in notebook 05", "Table A12; Figure 13",
               f"Unadjusted medians rank {rank[0]} and {rank[1]} highest and {rank[-1]} lowest; between-study variation is large, so the ranking is descriptive only.",
               "Matrix is confounded with study, KMA version and depth; FPKM is per rRNA-aligned fragment, which differs from CARD/RGI RPKM.")
