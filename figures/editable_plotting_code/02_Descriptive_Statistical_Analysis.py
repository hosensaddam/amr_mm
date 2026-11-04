# %% [markdown]
# # 02 · Descriptive Statistical Analysis
# **Research objectives:** context for Objectives 1–4 (coverage of cities, matrices, platforms and annotations) and Objective 5 (food-related contextual records).
#
# **Rationale:** before any comparison, the analysable evidence base must be described honestly: how many samples, studies, countries and cities are represented, in which matrices, with which sequencing platforms, over which years, and how many carry taxonomic or AMR annotations. These denominators constrain every later claim.
#
# **Input:** `derived/metadata_set.csv` (notebook 01): 8,053 included samples and 840 food-related samples held as a separate contextual dataset.
#
# **Methods:** frequencies, percentages, medians, interquartile ranges (IQR), means and standard deviations. Sequencing depth is summarised on the log10 scale because it spans five orders of magnitude. No inferential tests are performed in this notebook.

# %%
#BOOTSTRAP
NB = "NB02"
nb_header(NB)
import seaborn as sns
meta = pd.read_csv(DER / "metadata_set.csv", low_memory=False)
amr = load_csv("amr_genes_long.csv", usecols=["Sample_ID", "AMR_Reference_Database"])
tax_ids = set(load_csv("taxonomy_long.csv", usecols=["Sample_ID"]).Sample_ID)
div_ids = set(load_csv("diversity_long.csv", usecols=["Sample_ID"]).Sample_ID)
card_ids = set(amr.loc[amr.AMR_Reference_Database == "CARD", "Sample_ID"]); rf_ids = set(amr.loc[amr.AMR_Reference_Database == "ResFinder", "Sample_ID"])
meta["has_coords"] = meta.lat.notna() & meta.lon.notna()
meta["has_CARD_RGI"] = meta.Sample_ID.isin(card_ids); meta["has_ResFinder_KMA"] = meta.Sample_ID.isin(rf_ids)
meta["has_any_AMR"] = meta.has_CARD_RGI | meta.has_ResFinder_KMA
meta["has_taxonomy"] = meta.Sample_ID.isin(tax_ids); meta["has_diversity"] = meta.Sample_ID.isin(div_ids)
meta["city_label"] = np.where(is_missing(meta.City) | is_missing(meta.Country), np.nan, meta.City + " | " + meta.Country)
inc = meta[meta.Inclusion_Status == "Included"].copy(); ctx = meta[meta.Inclusion_Status != "Included"].copy()
print(len(inc), len(ctx))

# %% [markdown]
# ## 1. Overall evidence base

# %%
def overview(d, label):
    yrs = d.year.dropna()
    return {"Set": label, "Samples": len(d), "Datasets": d.Dataset_ID.nunique(), "Projects (accessions)": d.Project_Accession.nunique(),
            "Countries (reported)": d.loc[~is_missing(d.Country), "Country"].nunique(),
            "Source-labelled cities (City|Country)": d.city_label.nunique(),
            "Continents": d.loc[d.Continent != "Not reported", "Continent"].nunique(),
            "Sample types": d.Sample_Type.nunique(), "Samples with coordinates": int(d.has_coords.sum()),
            "Samples with any AMR gene table": int(d.has_any_AMR.sum()), "  CARD/RGI": int(d.has_CARD_RGI.sum()),
            "  ResFinder/KMA": int(d.has_ResFinder_KMA.sum()), "Samples with taxonomy table": int(d.has_taxonomy.sum()),
            "Samples with diversity indices": int(d.has_diversity.sum()),
            "Collection years (min-max)": f"{int(yrs.min())}-{int(yrs.max())}" if len(yrs) else "Not reported",
            "Samples with collection year": int(d.year.notna().sum())}
ov = pd.DataFrame([overview(inc, "Included urban environmental samples"), overview(ctx, "Food-related contextual samples")]).T
ov.columns = ov.iloc[0]; ov = ov.iloc[1:].reset_index().rename(columns={"index": "Quantity"})
save_table(ov, "D01_evidence_base_overview", "Overview of the analysable evidence base", NB, sheet="Descriptive_Statistics")
ov

# %% [markdown]
# ## 2. Distribution by sample type, dataset, country and continent

# %%
def depth_stats(s):
    g = s.dropna() / 1e9; lg = np.log10(s.dropna()[s.dropna() > 0])
    return pd.Series({"depth_n": len(g), "median_Gbp": g.median(), "Q1_Gbp": g.quantile(.25), "Q3_Gbp": g.quantile(.75),
                      "mean_log10_bases": lg.mean(), "sd_log10_bases": lg.std()})
bytype = inc.groupby("Sample_Type").agg(samples=("Sample_ID", "size"), datasets=("Dataset_ID", "nunique"),
        countries=("Country", lambda s: s[~is_missing(s)].nunique()), cities=("city_label", "nunique"),
        with_coordinates=("has_coords", "sum"), with_AMR=("has_any_AMR", "sum"), with_taxonomy=("has_taxonomy", "sum"),
        year_min=("year", "min"), year_max=("year", "max")).join(inc.groupby("Sample_Type").bases.apply(depth_stats).unstack())
bytype.insert(1, "pct_of_included", (100 * bytype.samples / len(inc)).round(1))
bytype = bytype.sort_values("samples", ascending=False).reset_index()
save_table(bytype.pipe(rnd, 3), "D02_by_sample_type", "Included samples by sample type: coverage and sequencing depth", NB, sheet="Sample_Distribution")
bytype

# %%
names = meta.drop_duplicates("Dataset_ID").set_index("Dataset_ID").Dataset_Name
byds = meta.groupby(["Dataset_ID", "Inclusion_Status"]).agg(samples=("Sample_ID", "size"),
        sample_types=("Sample_Type", lambda s: "; ".join(sorted(s.unique()))), countries=("Country", lambda s: s[~is_missing(s)].nunique()),
        cities=("city_label", "nunique"), platform=("Platform_Family", lambda s: "; ".join(sorted(s.unique()))),
        year_min=("year", "min"), year_max=("year", "max"), CARD_RGI=("has_CARD_RGI", "sum"), ResFinder_KMA=("has_ResFinder_KMA", "sum"),
        taxonomy=("has_taxonomy", "sum")).reset_index()
byds.insert(1, "Dataset_Name", byds.Dataset_ID.map(names))
byds = byds.sort_values("samples", ascending=False)
save_table(byds, "D03_by_dataset", "Samples by dataset (included and contextual)", NB, sheet="Sample_Distribution")
byds.head(15)

# %%
bycountry = inc[~is_missing(inc.Country)].groupby(["Continent", "Country"]).agg(samples=("Sample_ID", "size"),
        sample_types=("Sample_Type", "nunique"), cities=("city_label", "nunique"), datasets=("Dataset_ID", "nunique"),
        with_AMR=("has_any_AMR", "sum")).reset_index().sort_values("samples", ascending=False)
save_table(bycountry, "D04_by_country", "Included samples by country", NB, sheet="Geographic_Coverage")
ct = pd.crosstab(inc.Continent, inc.Sample_Type, margins=True, margins_name="Total")
save_table(ct.reset_index(), "D05_continent_by_sample_type", "Included samples: continent by sample type", NB, sheet="Geographic_Coverage")
top_c = bycountry.head(10); print(f"Top 10 countries hold {top_c.samples.sum()} of {len(inc)} included samples ({100*top_c.samples.sum()/len(inc):.1f}%)")
ct

# %%
bycity = inc[inc.city_label.notna()].groupby("city_label").agg(samples=("Sample_ID", "size"), sample_types=("Sample_Type", lambda s: "; ".join(sorted(s.unique()))),
        datasets=("Dataset_ID", lambda s: "; ".join(sorted(s.unique()))), with_AMR=("has_any_AMR", "sum"), with_taxonomy=("has_taxonomy", "sum")).reset_index().sort_values("samples", ascending=False)
save_table(bycity, "D06_by_city", "Included samples by source-labelled city", NB, sheet="Geographic_Coverage")
city_n = bycity.samples
city_sum = pd.DataFrame([{"Quantity": "Source-labelled cities", "Value": len(bycity)},
    {"Quantity": "Cities with >=10 samples", "Value": int((city_n >= 10).sum())},
    {"Quantity": "Cities with >=30 samples", "Value": int((city_n >= 30).sum())},
    {"Quantity": "Cities with 1 sample", "Value": int((city_n == 1).sum())},
    {"Quantity": "Median samples per city (IQR)", "Value": f"{city_n.median():.0f} ({city_n.quantile(.25):.0f}-{city_n.quantile(.75):.0f})"},
    {"Quantity": "Cities represented in >1 sample type", "Value": int((bycity.sample_types.str.count(";") >= 1).sum())},
    {"Quantity": "Cities with >=1 AMR-annotated sample", "Value": int((bycity.with_AMR > 0).sum())},
    {"Quantity": "Cities with >=1 taxonomy-annotated sample", "Value": int((bycity.with_taxonomy > 0).sum())},
    {"Quantity": "Included samples without a City|Country pair (city missing, or country missing/unverified)", "Value": int(inc.city_label.isna().sum())}])
save_table(city_sum, "D07_city_coverage_summary", "City coverage summary (relevance to the proposal's 35-city framing)", NB, sheet="Geographic_Coverage")
city_sum

# %% [markdown]
# ## 3. Proposal target cities (Dhaka, Delhi, Karachi)

# %%
loc = load_csv("sample_master.csv", usecols=["Sample_ID", "Locality_As_Reported"]).set_index("Sample_ID").Locality_As_Reported
meta["locality"] = meta.Sample_ID.map(loc)
tg = []
for c in ["Dhaka", "Delhi", "Karachi"]:
    m = meta[(meta.City.astype(str).str.contains(c, case=False)) | (meta.locality.astype(str).str.contains(c, case=False))]
    for (ds, st, inc_s), g in m.groupby(["Dataset_ID", "Sample_Type", "Inclusion_Status"]):
        tg.append(dict(target_city=c, Dataset_ID=ds, Dataset_Name=names.get(ds), Sample_Type=st, Inclusion_Status=inc_s, samples=len(g),
                       source_labelled_city=int(g.City.astype(str).str.contains(c, case=False).sum()),
                       years=f"{g.year.min():.0f}-{g.year.max():.0f}" if g.year.notna().any() else "Not reported",
                       CARD_RGI=int(g.has_CARD_RGI.sum()), ResFinder_KMA=int(g.has_ResFinder_KMA.sum()), taxonomy=int(g.has_taxonomy.sum())))
tg = pd.DataFrame(tg)
save_table(tg, "D08_target_cities", "Records for the proposal's named cities (Dhaka, Delhi, Karachi)", NB, sheet="Geographic_Coverage")
tg

# %% [markdown]
# ## 4. Sequencing platforms, read types and temporal coverage

# %%
plat = pd.crosstab(inc.Platform_Family, inc.Sample_Type, margins=True, margins_name="Total")
save_table(plat.reset_index(), "D09_platform_by_sample_type", "Platform family by sample type (included samples)", NB, sheet="Descriptive_Statistics")
rt = inc.Read_Type.fillna("Not reported").value_counts().rename_axis("Read_Type").reset_index(name="samples")
save_table(rt, "D10_read_type", "Read layout (included samples)", NB, sheet="Descriptive_Statistics")
yr = pd.crosstab(inc.year.fillna(-1).astype(int).replace(-1, "Not reported").astype(str), inc.Sample_Type)
save_table(yr.reset_index().rename(columns={"year": "Collection_Year"}), "D11_year_by_sample_type", "Collection year by sample type (included samples)", NB, sheet="Descriptive_Statistics")
plat

# %%
order = bytype.Sample_Type.tolist()
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.4), gridspec_kw={"width_ratios": [1.25, 1]})
cont_order = inc.Continent.value_counts().index.tolist()
cpal = dict(zip(cont_order, sns.color_palette("Set2", len(cont_order))))
tab = pd.crosstab(inc.Sample_Type, inc.Continent).reindex(order)[cont_order]
left = np.zeros(len(tab))
for c in cont_order:
    axes[0].barh(tab.index, tab[c], left=left, color=cpal[c], label=c, edgecolor="white", linewidth=.3); left += tab[c].values
for i, v in enumerate(left): axes[0].text(v + 15, i, f"{int(v):,}", va="center", fontsize=6.8)
axes[0].invert_yaxis(); axes[0].set_xlabel("Included samples (n)"); axes[0].set_title("A  Sample type by continent")
axes[0].legend(frameon=False, fontsize=6.5, loc="lower right")
ptab = pd.crosstab(inc.Sample_Type, inc.Platform_Family).reindex(order)
ppct = ptab.div(ptab.sum(axis=1), axis=0) * 100
pcols = ppct.columns.tolist(); ppal = dict(zip(pcols, sns.color_palette("Dark2", len(pcols))))
left = np.zeros(len(ppct))
for c in pcols:
    axes[1].barh(ppct.index, ppct[c], left=left, color=ppal[c], label=c, edgecolor="white", linewidth=.3); left += ppct[c].values
axes[1].invert_yaxis(); axes[1].set_yticklabels([]); axes[1].set_xlabel("% of samples"); axes[1].set_title("B  Sequencing platform family")
axes[1].legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(.5, -.15), ncol=3)
fig.tight_layout()
save_fig(fig, 2, "sample_type_continent_platform", "Included samples by sample type, continent and sequencing platform",
         f"(A) Number of included samples (n = {len(inc):,}) per sample type, stacked by continent (Natural Earth continent of the reported country; 'Not reported' where country is missing). (B) Percentage of each sample type sequenced on each platform family.",
         NB, "Objective 1 (context)", "tables/D02_by_sample_type.csv; tables/D09_platform_by_sample_type.csv")

# %%
fig, ax = plt.subplots(figsize=(6.9, 3.0))
y = inc[inc.year.notna()]; years = np.arange(int(y.year.min()), int(y.year.max()) + 1)
bottom = np.zeros(len(years))
for st in order:
    v = y[y.Sample_Type == st].year.value_counts().reindex(years, fill_value=0).values
    ax.bar(years, v, bottom=bottom, color=SAMPLE_TYPE_COLORS.get(st), label=st, width=.8, edgecolor="white", linewidth=.3); bottom += v
ax.set_xlabel("Collection year"); ax.set_ylabel("Included samples (n)"); ax.set_xticks(years[::2]); ax.tick_params(axis="x", rotation=45)
ax.legend(frameon=False, ncol=2, fontsize=6.5); ax.set_title("Temporal coverage of included samples")
save_fig(fig, 3, "temporal_coverage", "Collection year of included samples by sample type",
         f"Number of included samples per reported collection year, stacked by sample type (n = {len(y):,} of {len(inc):,} with a collection year; placeholder dates excluded).",
         NB, "Objective 1 / 3 (context)", "tables/D11_year_by_sample_type.csv")

# %% [markdown]
# ## 5. Annotation coverage and sequencing depth

# %%
cov = inc.groupby("Sample_Type")[["has_coords", "has_any_AMR", "has_CARD_RGI", "has_ResFinder_KMA", "has_taxonomy", "has_diversity"]].sum()
cov.insert(0, "samples", inc.groupby("Sample_Type").size()); cov = cov.reindex(order)
cov.loc["Total"] = cov.sum()
covp = cov.copy()
for c in cov.columns[1:]: covp[c + "_pct"] = (100 * cov[c] / cov.samples).round(1)
save_table(covp.reset_index(), "D12_annotation_coverage_by_type", "Coordinates, AMR, taxonomy and diversity coverage by sample type (included samples)", NB, sheet="Descriptive_Statistics")
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.2))
c2 = cov.drop(index="Total")
w = .16; x = np.arange(len(c2))
for i, (col, lab) in enumerate([("samples", "All included"), ("has_coords", "Coordinates"), ("has_any_AMR", "AMR gene table"), ("has_taxonomy", "Taxonomy table"), ("has_diversity", "Diversity indices")]):
    axes[0].bar(x + (i - 2) * w, c2[col].replace(0, np.nan), w, label=lab, color=sns.color_palette("Blues_r", 6)[i] if i else "#333333")
axes[0].set_yscale("log"); axes[0].set_xticks(x); axes[0].set_xticklabels(c2.index, rotation=35, ha="right"); axes[0].set_ylabel("Samples (log scale)")
axes[0].legend(frameon=False, fontsize=6.5); axes[0].set_title("A  Annotation coverage")
dd = inc[inc.bases > 0].assign(Gbp=lambda d: d.bases / 1e9)
sns.boxplot(data=dd, y="Sample_Type", x="Gbp", order=order, hue="Sample_Type", palette=SAMPLE_TYPE_COLORS, legend=False, fliersize=1.2, linewidth=.7, ax=axes[1], log_scale=True)
axes[1].set_xlabel("Sequenced bases per sample (Gbp, log scale)"); axes[1].set_ylabel(""); axes[1].set_title("B  Sequencing depth")
fig.tight_layout()
save_fig(fig, 4, "annotation_coverage_depth", "Annotation coverage and sequencing depth by sample type",
         "(A) Included samples per sample type with source-reported coordinates, a gene-level AMR table (CARD/RGI or ResFinder/KMA), a taxonomic table, or diversity indices (log scale; absent bars = zero). (B) Total sequenced bases per sample (Gbp, log scale); boxes show median and IQR. Most samples have metadata and raw reads only; processed AMR or taxonomic tables exist for a minority.",
         NB, "Objectives 1-2 (data availability)", "tables/D12_annotation_coverage_by_type.csv; tables/D02_by_sample_type.csv")
covp

# %%
n_amr = int(inc.has_any_AMR.sum()); n_tax = int(inc.has_taxonomy.sum())
record_finding(NB, "F01", "Objective 1-4 (evidence base)",
    f"The analysable evidence base comprises {len(inc):,} included urban environmental shotgun metagenomes from {inc.Dataset_ID.nunique()} datasets, {inc.loc[~is_missing(inc.Country),'Country'].nunique()} countries and {inc.city_label.nunique()} source-labelled cities, but only {n_amr:,} ({100*n_amr/len(inc):.1f}%) carry a processed AMR gene table and {n_tax} ({100*n_tax/len(inc):.1f}%) a taxonomic table.",
    len(inc), "Included samples", "Descriptive counts", "Table D01, D12; Figure 4",
    "Most records support only metadata-level description; taxonomic and AMR inference rests on two processed sub-collections.",
    "Raw reads (approx. 24 TB) were not downloaded or reprocessed; processed tables come from two non-equivalent pipelines.", "Descriptive")
record_finding(NB, "F02", "Objective 1 (geographic coverage)",
    f"Coverage is geographically uneven: the ten best-sampled countries contribute {top_c.samples.sum():,} of {len(inc):,} included samples ({100*top_c.samples.sum()/len(inc):.1f}%); {int((city_n == 1).sum())} of {len(bycity)} source-labelled cities contribute a single sample.",
    len(inc), "Included samples with a reported country / city", "Descriptive counts", "Tables D04, D06, D07",
    "City-level comparisons and rankings are not supported for most cities.", "City labels are source-reported; 625 locality strings were not standardized.", "Descriptive")
tgs = tg.groupby("target_city")[["samples", "CARD_RGI", "ResFinder_KMA", "taxonomy"]].sum()
desc = "; ".join(f"{c}: {int(r.samples)} samples ({int(r.CARD_RGI)} CARD/RGI, {int(r.ResFinder_KMA)} ResFinder/KMA, {int(r.taxonomy)} taxonomy)" for c, r in tgs.iterrows())
record_finding(NB, "F03", "Objective 1 (proposal target cities)",
    f"All three named cities are present, but processed annotation is sparse: {desc}.",
    int(tg.samples.sum()), "Included and contextual samples labelled with a target city (City or locality text)", "Descriptive counts", "Table D08",
    "Target-city AMR and taxonomy statements must be restricted to the cities with processed tables; the rest require raw-read reprocessing.",
    "Small and heterogeneous samples; locality strings are not standardized.", "Descriptive")
print(tgs)
