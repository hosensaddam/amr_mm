# %% [markdown]
# # 03 · Taxonomic Relative Abundance and Alpha Diversity
# **Research objective:** Objective 1 — variation in microbial community composition and diversity across urban environments and locations.
#
# **Analytic sets and rationale**
# * **Part A — DS02 (Indian metropolitan open-drain wastewater; Singh et al. 2026; n = 447; Delhi, Mumbai, Chennai, Kolkata; Mar 2022–Mar 2024).** The only included collection with a sample-level taxonomic abundance table (Kraken 2 + Bracken, phylum level). It allows within-study, within-matrix comparison of cities and seasons, which avoids cross-pipeline confounding.
# * **Part B — Martiny et al. global resistome resource (Zenodo 6919377).** Provides genus- and phylum-level richness, Shannon and Simpson indices from 16S/18S rRNA read mapping (KMA, SILVA) for 551 analytic units across six matrices. Composition tables (rRNA.tsv, 6.4 GB) were not downloaded, so only diversity indices can be compared.
#
# **Key constraints:** relative abundance is compositional and is not microbial load; phylum is the only rank available for DS02; UPDATE (AUD-11): the 19 DS02 open-drain sites were identified from the ENA/SRA experiment titles of PRJNA1415085; the site is the independent unit for city contrasts, and site-aware tests (site-level Kruskal-Wallis with whole-site permutation; mixed models with a site random intercept) are the primary inferential evidence. Sample-level tests are retained as descriptive only.

# %%
#BOOTSTRAP
NB = "NB03"
nb_header(NB)
import seaborn as sns
import statsmodels.formula.api as smf
import statsmodels.api as sm_api
m2 = pd.read_csv(DER / "ds02_meta.csv").set_index("Sample_ID")
ph_all = pd.read_csv(DER / "ds02_phylum_counts.csv", index_col=0).loc[m2.index]
domain = pd.Series({c: phylum_domain(c) for c in ph_all.columns}, name="Domain")
ph = ph_all.loc[:, domain.isin(["Bacteria", "Archaea"])]
ph.to_csv(DER / "ds02_prokaryote_phylum_counts.csv"); domain.rename_axis("Phylum").reset_index().to_csv(DER / "ds02_phylum_domain_map.csv", index=False)
m2["City"] = pd.Categorical(m2.City, ["Delhi", "Mumbai", "Chennai", "Kolkata"])
m2["Season"] = pd.Categorical(m2.Season, SEASON_ORDER)
print(ph_all.shape, "->", ph.shape, domain.value_counts().to_dict(), m2.City.value_counts().to_dict())

# %% [markdown]
# ## A0. Domain composition of the phylum table (DS02)
# The deposited Bracken table was produced against the Kraken 2 core_nt index, which contains eukaryotic and viral genomes. Phyla were assigned to domains (NCBI taxonomy names; rule list in `amr_common.phylum_domain`). Large eukaryotic fractions (e.g. Arthropoda, Chordata) in wastewater are plausibly driven by host DNA and database composition and cannot be validated here, so **all community analyses below use prokaryotic phyla only (Bacteria + Archaea), re-closed to 100%.**

# %%
ra_all = ph_all.div(ph_all.sum(axis=1), axis=0)
rep = load_csv("taxonomy_long.csv").pivot_table(index="Sample_ID", columns="Taxon_Name", values="Relative_Abundance", aggfunc="sum", fill_value=0).loc[ra_all.index, ra_all.columns]
print("max |recomputed RA - reported RA| =", float((ra_all - rep).abs().max().max()))
dom = (100 * ra_all.T.groupby(domain).sum().T)
dsum = dom.groupby(m2.City.astype(str)).agg(["mean", "std"]).round(2)
dsum.columns = [f"{a}_{b}" for a, b in dsum.columns]
dsum.loc["All samples"] = [f(dom[c]) for c in dom.columns for f in (lambda x: round(x.mean(), 2), lambda x: round(x.std(), 2))]
dcount = domain.value_counts().rename("phyla")
save_table(dsum.reset_index().rename(columns={"City": "Group"}), "T00_ds02_domain_composition", "DS02: share (%) of phylum-classified reads by domain (mean, SD)", NB,
           f"Phyla per domain: {dcount.to_dict()}", sheet="Taxonomic_Abundance")
record_exclusion(NB, "DS02 community analyses", "phylum", int((~domain.isin(['Bacteria','Archaea'])).sum()), "Eukaryotic and viral phyla removed; prokaryotic phyla re-closed")
dsum

# %% [markdown]
# ## A1. Relative abundance and dominant prokaryotic phyla (DS02)

# %%
tot = ph.sum(axis=1)
ra = ph.div(tot, axis=0)
prev = (ph > 0).mean()
summ = pd.DataFrame({"mean_RA_pct": 100 * ra.mean(), "sd_RA_pct": 100 * ra.std(), "median_RA_pct": 100 * ra.median(),
                     "Q1_pct": 100 * ra.quantile(.25), "Q3_pct": 100 * ra.quantile(.75), "prevalence_pct": 100 * prev,
                     "samples_detected": (ph > 0).sum()}).sort_values("mean_RA_pct", ascending=False)
summ.index.name = "Phylum"
save_table(summ.assign(Domain=summ.index.map(domain)).reset_index().pipe(rnd, 4), "T01_ds02_phylum_relative_abundance", "DS02 prokaryotic phylum relative abundance (Bracken; % of prokaryote-classified reads), n = 447", NB, sheet="Taxonomic_Abundance")
top10 = summ.index[:10].tolist()
bycity = (100 * ra[top10].groupby(m2.City, observed=True).mean()).T
bycity.loc["Other phyla"] = 100 - bycity.sum()
byseason = (100 * ra[top10].groupby(m2.Season, observed=True).mean()).T
byseason.loc["Other phyla"] = 100 - byseason.sum()
save_table(bycity.pipe(rnd, 3).reset_index().rename(columns={"index": "Phylum"}), "T02_ds02_top_phyla_by_city", "Mean relative abundance (%) of the ten dominant phyla by city (DS02)", NB, sheet="Taxonomic_Abundance")
save_table(byseason.pipe(rnd, 3).reset_index().rename(columns={"index": "Phylum"}), "T03_ds02_top_phyla_by_season", "Mean relative abundance (%) of the ten dominant phyla by season (DS02)", NB, sheet="Taxonomic_Abundance")
summ.head(12).round(3)

# %%
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.5), gridspec_kw={"width_ratios": [1, 1.1]})
pal = dict(zip(top10 + ["Other phyla"], sns.color_palette("tab10", 10) + [(.8, .8, .8)]))
for ax, tab, ttl in [(axes[0], bycity, "A  By city"), (axes[1], byseason, "B  By season")]:
    bottom = np.zeros(tab.shape[1])
    for p in tab.index:
        ax.bar(range(tab.shape[1]), tab.loc[p].values, bottom=bottom, color=pal[p], label=p, edgecolor="white", linewidth=.3); bottom += tab.loc[p].values
    ax.set_xticks(range(tab.shape[1])); ax.set_xticklabels([str(c).replace(" (", "\n(") for c in tab.columns], fontsize=6.8)
    ax.set_ylabel("Mean relative abundance (%)"); ax.set_ylim(0, 100); ax.set_title(ttl)
ns = m2.City.value_counts(); axes[0].set_xticklabels([f"{c}\n(n={ns[c]})" for c in bycity.columns], fontsize=6.8)
axes[1].legend(frameon=False, fontsize=6.5, bbox_to_anchor=(1.01, 1), loc="upper left")
fig.tight_layout()
save_fig(fig, 5, "ds02_phylum_composition", "Prokaryotic phylum-level composition of Indian open-drain wastewater by city and season (DS02)",
         "Mean Bracken relative abundance (percentage of prokaryote-classified reads; eukaryotic and viral phyla removed) of the ten most abundant prokaryotic phyla (remaining phyla pooled as 'Other') in 447 shotgun metagenomes from 19 open-drain sites in Delhi, Mumbai, Chennai and Kolkata (March 2022-March 2024). Relative abundances are proportions of phylum-classified reads, not absolute loads. Seasons follow India Meteorological Department definitions.",
         NB, "Objective 1", "tables/T02_ds02_top_phyla_by_city.csv; tables/T03_ds02_top_phyla_by_season.csv")

# %% [markdown]
# ## A2. Alpha diversity (DS02, phylum level)
# Observed richness, Shannon (natural log), Gini–Simpson (1 − Σp²) and Pielou evenness are computed from Bracken counts. Because observed richness depends on sequencing depth, a rarefied richness (random subsampling without replacement to the smallest library, seed fixed) is also reported.

# %%
rng = np.random.default_rng(SEED)
depth_min = int(tot.min())
rare = np.array([rng.multivariate_hypergeometric(ph.loc[s].values.astype(np.int64), depth_min) for s in ph.index])
alpha = pd.DataFrame({"observed_richness": (ph > 0).sum(axis=1), "rarefied_richness": (rare > 0).sum(axis=1),
                      "shannon": [shannon(r) for r in ph.values], "gini_simpson": [gini_simpson(r) for r in ph.values],
                      "bracken_phylum_reads": tot}, index=ph.index)
alpha["pielou_evenness"] = alpha.shannon / np.log(alpha.observed_richness)
alpha = alpha.join(m2[["City", "Season", "year", "month", "Collection_Date", "reads_total"]])
alpha.to_csv(DER / "ds02_alpha_diversity.csv")
print("rarefaction depth:", depth_min)
r_depth = spearman_ci(np.log10(alpha.bracken_phylum_reads), alpha.observed_richness)
r_depth_r = spearman_ci(np.log10(alpha.bracken_phylum_reads), alpha.rarefied_richness)
record_test(NB, "DS02 observed phylum richness vs log10 classified reads", "Spearman", r_depth["n"], r_depth["rho"], None, r_depth["p"], "rho", r_depth["rho"], r_depth["ci_low"], r_depth["ci_high"], "depth check", objective="Objective 1")
record_test(NB, "DS02 rarefied phylum richness vs log10 classified reads", "Spearman", r_depth_r["n"], r_depth_r["rho"], None, r_depth_r["p"], "rho", r_depth_r["rho"], r_depth_r["ci_low"], r_depth_r["ci_high"], "depth check", objective="Objective 1")
print(r_depth, r_depth_r)
desc = alpha.groupby("City", observed=True)[["observed_richness", "rarefied_richness", "shannon", "gini_simpson", "pielou_evenness"]].agg(["median", lambda s: s.quantile(.25), lambda s: s.quantile(.75), "mean", "std"])
desc.columns = [f"{a}_{b if isinstance(b, str) else ('Q1' if i % 5 == 1 else 'Q3')}" for i, (a, b) in enumerate(desc.columns)]
desc = desc.rename(columns=lambda c: c.replace("<lambda_0>", "Q1").replace("<lambda_1>", "Q3"))
save_table(desc.pipe(rnd, 4).reset_index(), "T04_ds02_alpha_by_city", "DS02 phylum-level alpha diversity by city (median, IQR, mean, SD)", NB, sheet="Alpha_Diversity")
descs = alpha.groupby("Season", observed=True)[["rarefied_richness", "shannon", "gini_simpson"]].agg(["median", "mean", "std", "count"])
descs.columns = ["_".join(c) for c in descs.columns]
save_table(descs.pipe(rnd, 4).reset_index(), "T05_ds02_alpha_by_season", "DS02 phylum-level alpha diversity by season", NB, sheet="Alpha_Diversity")
desc.T.round(3)

# %%
tests = []
for idx in ["rarefied_richness", "shannon", "gini_simpson", "pielou_evenness"]:
    for fac in ["City", "Season"]:
        k = kruskal_eps2(alpha.assign(**{fac: alpha[fac].astype(str)}), idx, fac)
        tests.append(dict(index=idx, factor=fac, **k))
        record_test(NB, f"DS02 {idx} by {fac}", "Kruskal-Wallis", k["n"], k["H"], k["df"], k["p"], "epsilon-squared", k["epsilon2"], family="DS02 alpha diversity [sample-level, descriptive only: pseudo-replicated]", note="Superseded for inference by site-aware tests (T06s); AUD-11", objective="Objective 1")
tests = pd.DataFrame(tests); tests["p_BH"] = bh(tests.p)
save_table(tests.pipe(rnd, 6), "T06_ds02_alpha_kruskal", "DESCRIPTIVE (sample-level; 447 repeated samples of 19 sites, pseudo-replicated - see T06s for inference): Kruskal-Wallis of DS02 alpha diversity by city and season", NB, sheet="Statistical_Tests")
pw = pd.concat([pairwise_mwu(alpha.assign(City=alpha.City.astype(str)), i, "City").assign(index=i) for i in ["rarefied_richness", "shannon"]])
save_table(pw.pipe(rnd, 6), "T07_ds02_alpha_pairwise_city", "Pairwise Mann-Whitney tests between cities (Holm-adjusted), DS02 alpha diversity", NB, sheet="Statistical_Tests")
tests

# %% [markdown]
# **Multivariable model.** Shannon diversity is modelled jointly on city and season (OLS with heteroscedasticity-consistent HC3 standard errors; type II Wald tests), so that the city contrast is adjusted for the uneven seasonal mix and vice versa. Sensitivity: the same model on city–month means, which removes within-month replication of sites.

# %%
ols = smf.ols("shannon ~ C(City, Treatment('Delhi')) + C(Season, Treatment('Winter (Dec-Feb)')) + C(year)", data=alpha.assign(City=alpha.City.astype(str), Season=alpha.Season.astype(str))).fit(cov_type="HC3")
an = sm_api.stats.anova_lm(smf.ols("shannon ~ C(City) + C(Season) + C(year)", data=alpha.assign(City=alpha.City.astype(str), Season=alpha.Season.astype(str))).fit(), typ=2)
coefs = pd.DataFrame({"estimate": ols.params, "SE_HC3": ols.bse, "ci_low": ols.conf_int()[0], "ci_high": ols.conf_int()[1], "p_value": ols.pvalues}).reset_index().rename(columns={"index": "term"})
save_table(coefs.pipe(rnd, 5), "T08_ds02_shannon_ols", "OLS of DS02 phylum Shannon diversity on city, season and year (HC3 SEs; reference Delhi, winter, 2022)", NB, sheet="Comparative_Analysis")
an["eta2_partial"] = an.sum_sq / (an.sum_sq + an.loc["Residual", "sum_sq"])
save_table(an.reset_index().rename(columns={"index": "term"}).pipe(rnd, 6), "T09_ds02_shannon_anova", "Type II ANOVA for DS02 Shannon model (classical SEs; partial eta-squared)", NB, sheet="Statistical_Tests")
for t in ["C(City)", "C(Season)"]:
    record_test(NB, "DS02 Shannon ~ City + Season + Year (type II)", f"OLS F-test: {t}", len(alpha), an.loc[t, "F"], f"{int(an.loc[t,'df'])}, {int(an.loc['Residual','df'])}", an.loc[t, "PR(>F)"], "partial eta-squared", an.loc[t, "eta2_partial"], family="DS02 alpha models [sample-level OLS, descriptive only: pseudo-replicated]", objective="Objective 1")
cm = alpha.assign(City=alpha.City.astype(str), Season=alpha.Season.astype(str), ym=alpha.Collection_Date).groupby(["City", "ym", "Season", "year"], as_index=False).shannon.mean()
an_cm = sm_api.stats.anova_lm(smf.ols("shannon ~ C(City) + C(Season) + C(year)", data=cm).fit(), typ=2)
an_cm["eta2_partial"] = an_cm.sum_sq / (an_cm.sum_sq + an_cm.loc["Residual", "sum_sq"])
save_table(an_cm.reset_index().rename(columns={"index": "term"}).pipe(rnd, 6), "T10_ds02_shannon_anova_citymonth", f"Sensitivity: type II ANOVA on city-month mean Shannon (n = {len(cm)} city-months)", NB, sheet="Sensitivity_Analysis")
for t in ["C(City)", "C(Season)"]:
    record_test(NB, "DS02 city-month mean Shannon ~ City + Season + Year (sensitivity)", f"OLS F-test: {t}", len(cm), an_cm.loc[t, "F"], f"{int(an_cm.loc[t,'df'])}, {int(an_cm.loc['Residual','df'])}", an_cm.loc[t, "PR(>F)"], "partial eta-squared", an_cm.loc[t, "eta2_partial"], family="DS02 alpha models (sensitivity)", objective="Objective 1")
print(an.round(4)); print(an_cm.round(4))

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

# %%
fig, axes = plt.subplots(1, 3, figsize=(6.9, 3.0))
a2 = alpha.assign(City=alpha.City.astype(str), Season=alpha.Season.astype(str))
for ax, (col, lab) in zip(axes[:2], [("shannon", "Shannon diversity (phylum)"), ("rarefied_richness", f"Rarefied phylum richness ({depth_min:,} reads)")]):
    sns.boxplot(data=a2, x="City", y=col, order=list(CITY_COLORS), hue="City", palette=CITY_COLORS, legend=False, showfliers=False, linewidth=.7, ax=ax)
    sns.stripplot(data=a2, x="City", y=col, order=list(CITY_COLORS), color="k", size=1.6, alpha=.35, jitter=.25, ax=ax)
    ks = site_tests[(site_tests["index"] == col) & (site_tests.factor == "City")].iloc[0]
    ax.set_title(f"{'A' if col=='shannon' else 'B'}  {lab.split(' (')[0]} by city\nsite-level KW ({int(ks.n_units)} sites): H={ks.H:.1f}, perm. p={ks.p_perm:.3f}", fontsize=7.2)
    ax.set_xlabel(""); ax.set_ylabel(lab)
sns.boxplot(data=a2, x="Season", y="shannon", order=SEASON_ORDER, hue="Season", palette=SEASON_COLORS, legend=False, showfliers=False, linewidth=.7, ax=axes[2])
sns.stripplot(data=a2, x="Season", y="shannon", order=SEASON_ORDER, color="k", size=1.6, alpha=.35, jitter=.25, ax=axes[2])
kss = site_tests[(site_tests["index"] == "shannon") & (site_tests.factor == "C(Season)")].iloc[0]
axes[2].set_xticks(range(4)); axes[2].set_xticklabels(["Winter\nDec-Feb", "Pre-mons.\nMar-May", "Monsoon\nJun-Sep", "Post-mons.\nOct-Nov"], fontsize=6.2); axes[2].set_xlabel(""); axes[2].set_ylabel("Shannon diversity (phylum)")
axes[2].set_title(f"C  Shannon by season\nmixed model (site RE) LRT p={kss.p:.2g}", fontsize=7.2)
fig.tight_layout()
save_fig(fig, 6, "ds02_alpha_diversity", "Phylum-level alpha diversity of Indian open-drain wastewater by city and season (DS02)",
         f"Shannon diversity (natural log) and rarefied observed richness (subsampled to {depth_min:,} phylum-classified reads) for 447 samples. Boxes: median and IQR; points: individual samples. Panel statistics are site-aware (AUD-11): city contrasts use Kruskal-Wallis on the 19 site medians with a whole-site permutation p-value, and the season contrast a mixed-model likelihood-ratio test with a site random intercept (Table T06s); sample-level pairwise tests (descriptive only) are in Table T07. Points are repeated monthly samples of 19 identified open-drain sites and are not independent observations.",
         NB, "Objective 1", "derived/ds02_alpha_diversity.csv; tables/T06_ds02_alpha_kruskal.csv")

# %%
hm = 100 * ra[summ.index[:20]].groupby([m2.City.astype(str), m2.Season.astype(str)]).mean()
hm = hm.reindex(pd.MultiIndex.from_product([list(CITY_COLORS), SEASON_ORDER])).dropna(how="all")
z = np.log10(hm.T + 1e-4); z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(6.9, 4.8))
sns.heatmap(z, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5, cbar_kws={"label": "Row z-score of log10 mean RA"}, linewidths=.3, ax=ax)
ax.set_xticklabels([f"{c} · {s.split(' (')[0]}" for c, s in hm.index], rotation=60, ha="right", fontsize=6.5)
ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title("Top 20 prokaryotic phyla: city-by-season pattern (DS02)")
save_fig(fig, 7, "ds02_phylum_heatmap", "Relative abundance of the 20 most abundant prokaryotic phyla by city and season (DS02)",
         "Heatmap of row-standardised log10 mean relative abundance (Bracken) of the 20 most abundant phyla for each city-season combination (n = 447 samples). Colours compare a phylum across columns only; they do not compare phyla with each other.",
         NB, "Objective 1", "tables/T01_ds02_phylum_relative_abundance.csv")

# %% [markdown]
# ## B. Diversity across urban matrices (Martiny et al. resource; SILVA rRNA mapping)
# Unit = analytic sample (PRJEB27054 accessions merged to biological samples; PRJEB13831 held out). For units with several runs, genus-level indices come from the deepest run (indices are not additive across runs). Food-related records are a separate contextual dataset and are shown apart from urban matrices.

# %%
rmap = pd.read_csv(DER / "rf_run_to_unit.csv")
dv = pd.read_csv(HARM / "subset_6919377_diversity_core_runs.tsv", sep="\t")
gen = dv[dv.category == "Genera"].merge(rmap.drop(columns="kma_version"), left_on="run_accession", right_on="Run_ID")
gen = gen[gen.primary].sort_values("total_fragments", ascending=False).drop_duplicates("unit_id")
phy = dv[dv.category == "Phyla"].set_index("run_accession")
gen["phylum_richness"] = gen.run_accession.map(phy.n); gen["phylum_shannon"] = gen.run_accession.map(phy.Shannon)
gen = gen.rename(columns={"n": "genus_richness", "Shannon": "genus_shannon", "Simpson": "genus_simpson"})
gen["log10_total_fragments"] = np.log10(gen.total_fragments)
gen["Stratum"] = np.where(gen.Inclusion_Status == "Included", "Urban environmental", "Food-related (contextual)")
gen[["unit_id", "Run_ID", "Dataset_ID", "Project_Accession", "Sample_Type", "Stratum", "genus_richness", "genus_shannon", "genus_simpson", "phylum_richness", "phylum_shannon", "total_fragments", "log10_total_fragments", "kma_version"]].to_csv(DER / "rf_unit_rrna_diversity.csv", index=False)
gs = gen.groupby(["Stratum", "Sample_Type"]).agg(units=("unit_id", "size"), studies=("Project_Accession", "nunique"),
        genus_shannon_median=("genus_shannon", "median"), genus_shannon_Q1=("genus_shannon", lambda s: s.quantile(.25)), genus_shannon_Q3=("genus_shannon", lambda s: s.quantile(.75)),
        genus_richness_median=("genus_richness", "median"), genus_richness_Q1=("genus_richness", lambda s: s.quantile(.25)), genus_richness_Q3=("genus_richness", lambda s: s.quantile(.75)),
        total_fragments_median=("total_fragments", "median")).reset_index()
save_table(gs.pipe(rnd, 3), "T11_rrna_diversity_by_sample_type", "Genus-level (SILVA rRNA) diversity by sample type, Martiny et al. resource units", NB, sheet="Alpha_Diversity")
gs

# %%
urb = gen[gen.Stratum == "Urban environmental"].copy()
k = kruskal_eps2(urb, "genus_shannon", "Sample_Type")
record_test(NB, "Genus Shannon by urban sample type (unit level)", "Kruskal-Wallis", k["n"], k["H"], k["df"], k["p"], "epsilon-squared", k["epsilon2"], family="Cross-matrix diversity", objective="Objective 1")
st_med = urb.groupby(["Project_Accession", "Sample_Type"], as_index=False).genus_shannon.median()
k2 = kruskal_eps2(st_med, "genus_shannon", "Sample_Type")
record_test(NB, "Genus Shannon by urban sample type (study-level medians; sensitivity)", "Kruskal-Wallis", k2["n"], k2["H"], k2["df"], k2["p"], "epsilon-squared", k2["epsilon2"], family="Cross-matrix diversity", note="unit = study x sample type", objective="Objective 1")
md = smf.mixedlm("genus_shannon ~ C(Sample_Type, Treatment('Municipal wastewater')) + log10_total_fragments", urb, groups=urb.Project_Accession).fit(reml=True)
mdt = pd.DataFrame({"estimate": md.fe_params, "SE": md.bse_fe, "ci_low": md.conf_int().loc[md.fe_params.index, 0], "ci_high": md.conf_int().loc[md.fe_params.index, 1], "p_value": md.pvalues[md.fe_params.index]}).reset_index().rename(columns={"index": "term"})
mdt.loc[len(mdt)] = dict(term="Study random-effect variance", estimate=float(md.cov_re.iloc[0, 0]), SE=np.nan, ci_low=np.nan, ci_high=np.nan, p_value=np.nan)
mdt.loc[len(mdt)] = dict(term="Residual variance", estimate=float(md.scale), SE=np.nan, ci_low=np.nan, ci_high=np.nan, p_value=np.nan)
icc = float(md.cov_re.iloc[0, 0]) / (float(md.cov_re.iloc[0, 0]) + md.scale)
mdt.loc[len(mdt)] = dict(term="ICC (study)", estimate=icc, SE=np.nan, ci_low=np.nan, ci_high=np.nan, p_value=np.nan)
save_table(mdt.pipe(rnd, 5), "T12_genus_shannon_mixed_model", f"Linear mixed model: genus Shannon ~ sample type + log10 mapped fragments + (1 | study); urban units n = {len(urb)}, studies = {urb.Project_Accession.nunique()}", NB, sheet="Comparative_Analysis")
print(k, k2, "ICC", round(icc, 3)); mdt

# %%
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.2))
order = gen.groupby("Sample_Type").genus_shannon.median().sort_values(ascending=False).index.tolist()
for ax, (col, lab) in zip(axes, [("genus_shannon", "Genus Shannon diversity"), ("genus_richness", "Observed genus richness")]):
    sns.boxplot(data=gen, x=col, y="Sample_Type", order=order, hue="Sample_Type", palette=SAMPLE_TYPE_COLORS, legend=False, showfliers=False, linewidth=.7, ax=ax)
    sns.stripplot(data=gen, x=col, y="Sample_Type", order=order, color="k", size=1.8, alpha=.35, jitter=.25, ax=ax)
    ax.set_ylabel(""); ax.set_xlabel(lab)
nn = gen.Sample_Type.value_counts(); sn = gen.groupby("Sample_Type").Project_Accession.nunique()
axes[0].set_yticklabels([f"{t}\n(units={nn[t]}, studies={sn[t]})" for t in order], fontsize=6.8); axes[1].set_yticklabels([])
axes[0].set_title(f"A  Genus Shannon (urban: KW ε²={k['epsilon2']:.2f}; study-level ε²={k2['epsilon2']:.2f})", fontsize=7.2); axes[1].set_title("B  Genus richness", fontsize=7.2)
fig.tight_layout()
save_fig(fig, 8, "rrna_diversity_by_matrix", "Bacterial genus diversity across urban matrices and food-related environments (Martiny et al. resource)",
         f"Genus-level Shannon diversity and observed richness from SILVA rRNA read mapping (KMA) for {len(gen)} analytic units ({gen.Project_Accession.nunique()} studies). Food-related environments are a separate contextual dataset. Richness depends on sequencing depth; the mixed model (Table T12) adjusts for log10 mapped fragments and treats study as a random effect.",
         NB, "Objective 1", "derived/rf_unit_rrna_diversity.csv; tables/T11, T12")

# %%
city_p = tests[(tests["index"] == "shannon") & (tests.factor == "City")].iloc[0]; seas_p = tests[(tests["index"] == "shannon") & (tests.factor == "Season")].iloc[0]
med = alpha.groupby("City", observed=True).shannon.median()
record_finding(NB, "F04", "Objective 1", f"Prokaryotic communities of open-drain wastewater in the four Indian cities were dominated by {', '.join(top10[:3])} (mean {summ.mean_RA_pct.iloc[0]:.1f}%, {summ.mean_RA_pct.iloc[1]:.1f}% and {summ.mean_RA_pct.iloc[2]:.1f}% of prokaryote-classified reads); eukaryotic and viral phyla accounted for a mean {dom['Eukaryota'].mean():.1f}% and {dom['Viruses'].mean():.1f}% of all phylum-classified reads and were excluded.",
               len(ra), "DS02 samples with Bracken phylum profiles", "Descriptive (mean, SD, prevalence)", "Table T01; Figure 5",
               "Phylum structure is broadly shared across cities; differences are of degree, not of dominant taxa.", "Phylum resolution only; relative not absolute abundance; single pipeline (Kraken 2/Bracken, core_nt); large eukaryotic fraction unexplained.")
ksite = site_tests[(site_tests["index"] == "shannon") & (site_tests.factor == "City")].iloc[0]
kmc = site_tests[(site_tests["index"] == "shannon") & (site_tests.factor == "C(City)")].iloc[0]
kms = site_tests[(site_tests["index"] == "shannon") & (site_tests.factor == "C(Season)")].iloc[0]
record_finding(NB, "F05", "Objective 1", "Median phylum-level Shannon diversity was " + ", ".join(f"{c} {v:.2f}" for c, v in med.items()) + f". With the open-drain site as the independent unit (19 sites), the city difference was {'statistically supported' if ksite.p_perm < .05 else 'not statistically supported'} (Kruskal-Wallis on site medians H = {ksite.H:.1f}, df = {int(ksite.df)}, whole-site permutation p = {ksite.p_perm:.3f}; site-level epsilon-squared = {ksite.epsilon2:.2f}; mixed-model LRT p = {kmc.p:.2g}); within sites, season was {'associated' if kms.p < .05 else 'not associated'} with Shannon diversity (mixed-model LRT p = {kms.p:.2g}).",
               len(alpha), f"DS02: {len(alpha)} samples from {int(ksite.n_units)} sites", "Site-aware (AUD-11): site-level Kruskal-Wallis with whole-site permutation; site random-intercept mixed models (T06s, T08s). Sample-level p-values (T06-T09) are descriptive only.", "Tables T04, T06s, T08s; Figure 6",
               f"Effect directions agree with the sample-level description (sample-level epsilon-squared = {city_p.epsilon2:.2f}, descriptive), but the evidence for a city effect rests on 19 sites, not 447 samples; the site ICC was {mm_sh['icc']:.2f}.", "Few sites per city (Kolkata 2); city is confounded with site selection; phylum resolution only.")
record_finding(NB, "F06", "Objective 1", f"Across urban matrices in the Martiny et al. resource, genus Shannon diversity varied by sample type at unit level (epsilon-squared = {k['epsilon2']:.2f}) and at study level (epsilon-squared = {k2['epsilon2']:.2f}); the study random effect accounted for an ICC of {icc:.2f} in the adjusted mixed model.",
               len(urb), f"Urban analytic units ({urb.Project_Accession.nunique()} studies)", f"Kruskal-Wallis p = {k['p']:.1e} (unit), {k2['p']:.2f} (study); mixed model T12", "Tables T11-T12; Figure 8",
               "Matrix differences in diversity are confounded with study; most variance is between studies, so matrix contrasts must be read cautiously.", "Few studies per matrix (air, built environment, transit each from 1-3 studies); genus composition tables not downloaded.")
