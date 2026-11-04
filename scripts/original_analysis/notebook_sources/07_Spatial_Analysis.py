# %% [markdown]
# # 07 · Geographic and Spatial Analysis
# **Research objective:** Objective 4 — whether the available data support geographical analysis and spatial pattern identification; Objective 1–2 geographic variation.
#
# **Rules applied.** Only producer-reported coordinates are mapped (no geocoding of city names, no imputation). Where coordinates are absent (e.g., DS02, which reports only the city), the country or city is the spatial unit. Maps show observed values with their sample sizes; no interpolated surfaces or "risk" maps are produced. Spatial autocorrelation is tested only for one internally comparable set — the 2016 Global Sewage Surveillance pilot (PRJEB27054; one protocol, one sampling round, ResFinder/KMA) — using global Moran's I with k-nearest-neighbour weights and a Mantel test of distance decay. "Risk" is not inferred: the mapped quantity is ARG abundance in sewage, not human-health risk.

# %%
#BOOTSTRAP
NB = "NB07"
nb_header(NB)
import geopandas as gpd, seaborn as sns
from skbio import DistanceMatrix
from skbio.stats.distance import mantel
world = load_world()
meta = pd.read_csv(DER / "metadata_set.csv", low_memory=False)
inc = meta[meta.Inclusion_Status == "Included"].copy()
us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0)

def basemap(ax):
    world.plot(ax=ax, color="#EDEDED", edgecolor="white", linewidth=.3)
    ax.set_xlim(-180, 180); ax.set_ylim(-60, 85); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel("")
    for s in ax.spines.values(): s.set_visible(False)

# %% [markdown]
# ## 1. Global distribution of included samples (source-reported coordinates)

# %%
pts = inc.dropna(subset=["lat", "lon"])
cov = pd.DataFrame([dict(item="Included samples", n=len(inc)), dict(item="with source-reported coordinates (mapped)", n=len(pts)),
                    dict(item="without coordinates (not mapped)", n=len(inc) - len(pts)),
                    dict(item="distinct coordinate pairs", n=pts[["lat", "lon"]].drop_duplicates().shape[0]),
                    dict(item="flagged coordinate-country mismatch >25 km (mapped as reported)", n=int(pts.coord_country_flag.sum()))])
by_prec = pts.Coordinate_Precision.value_counts().rename_axis("Coordinate_Precision").reset_index(name="samples")
save_table(cov, "S01_mapping_coverage", "Coverage of the global sample map", NB, sheet="Spatial_Analysis")
save_table(by_prec, "S02_coordinate_precision", "Precision of source-reported coordinates (mapped included samples)", NB, sheet="Spatial_Analysis")
fig, ax = plt.subplots(figsize=(6.9, 3.6)); basemap(ax)
for st in pts.Sample_Type.value_counts().index:
    d = pts[pts.Sample_Type == st].drop_duplicates(["lat", "lon"])
    ax.scatter(d.lon, d.lat, s=9, color=SAMPLE_TYPE_COLORS.get(st), alpha=.75, edgecolors="k", linewidths=.15,
               label=f"{st} ({(pts.Sample_Type == st).sum():,} samples, {len(d):,} sites)")
ax.legend(frameon=False, fontsize=6.5, loc="lower left", markerscale=1.6)
ax.set_title(f"Included urban metagenomes with source-reported coordinates (n = {len(pts):,} of {len(inc):,})", fontsize=8)
save_fig(fig, 22, "global_sample_map", "Global distribution of included urban metagenomes with source-reported coordinates",
         f"Each point is a distinct reported coordinate pair (overlapping samples drawn once); colour = sample type. {len(inc)-len(pts):,} included samples without coordinates (including all DS02 Indian samples, which report city only) are not shown. Coordinates were not geocoded or corrected; {int(pts.coord_country_flag.sum())} flagged mismatches are plotted as reported. Basemap: Natural Earth 1:50m.",
         NB, "Objective 1 / 4", "tables/S01, S02; derived/metadata_set.csv")

# %%
cc = inc[inc.ISO2.notna()].groupby("ISO2").agg(samples=("Sample_ID", "size"), sample_types=("Sample_Type", "nunique"), with_AMR=("Dataset_ID", lambda s: 0)).reset_index()
amr_ids = set(load_csv("amr_genes_long.csv", usecols=["Sample_ID"]).Sample_ID)
cc["with_AMR_table"] = cc.ISO2.map(inc[inc.Sample_ID.isin(amr_ids)].groupby("ISO2").size()).fillna(0).astype(int); cc = cc.drop(columns="with_AMR")
w = world.merge(cc, on="ISO2", how="left")
fig, axes = plt.subplots(2, 1, figsize=(6.9, 6.8))
for ax, col, lab in [(axes[0], "samples", "Included samples per country (log scale)"), (axes[1], "with_AMR_table", "Included samples with a processed AMR gene table")]:
    basemap(ax)
    v = w[w[col] > 0]
    v.plot(ax=ax, column=np.log10(v[col]), cmap="viridis", edgecolor="white", linewidth=.2, legend=False)
    sm_ = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(np.log10(v[col]).min(), np.log10(v[col]).max()))
    cb = fig.colorbar(sm_, ax=ax, fraction=.025, pad=.01); ticks = cb.get_ticks(); cb.set_ticks(ticks); cb.set_ticklabels([f"{10**t:,.0f}" for t in ticks]); cb.set_label(lab, fontsize=7)
    ax.set_title(f"{'A' if col=='samples' else 'B'}  {lab}: {int(w.loc[w[col]>0, 'ISO2'].nunique())} countries", fontsize=8)  # AUDIT FIX AUD-03: count countries (ISO2), not map polygons
save_fig(fig, 23, "country_coverage_choropleth", "Country-level coverage of included samples and of processed AMR tables",
         f"(A) Number of included samples per reported country (log colour scale; grey = none). (B) Included samples with a processed AMR gene table (CARD/RGI or ResFinder/KMA). Countries: {cc.shape[0]} with samples; {int((cc.with_AMR_table>0).sum())} with AMR tables. Small countries may not be visible at this scale; counts are in Table D04.",
         NB, "Objective 4", "tables/D04_by_country.csv")
save_table(cc.sort_values("samples", ascending=False), "S03_country_coverage", "Country coverage: included samples, sample types and samples with AMR tables", NB, sheet="Geographic_Coverage")

# %% [markdown]
# ## 2. Observed normalized ARG abundance at sampling locations (ResFinder/KMA; municipal wastewater)

# %%
ww = us[(us.Sample_Type == "Municipal wastewater") & us.lat.notna()].copy()
pil = ww[ww.Project_Accession == "PRJEB27054"].copy()
print("wastewater units with coordinates:", len(ww), "| 2016 pilot units:", len(pil), "| pilot units whose accessions disagree on coordinates:", int((pil.lat_sd > 1e-6).sum()))
fig, ax = plt.subplots(figsize=(6.9, 3.6)); basemap(ax)
sc = ax.scatter(ww.lon, ww.lat, c=ww.log10_total_FPKM, cmap="magma_r", s=22, edgecolors="k", linewidths=.25, vmin=ww.log10_total_FPKM.quantile(.02), vmax=ww.log10_total_FPKM.quantile(.98))
cb = fig.colorbar(sc, ax=ax, fraction=.025, pad=.01); cb.set_label("log10 total ARG FPKM (observed)")
ax.set_title(f"Observed normalized ARG abundance in municipal wastewater at reported locations (ResFinder/KMA; {len(ww)} units, {ww.Project_Accession.nunique()} studies)", fontsize=8)
save_fig(fig, 24, "wastewater_arg_load_map", "Observed ARG abundance in municipal wastewater at reported sampling locations (ResFinder/KMA)",
         f"Each point is one analytic unit (biological sample) with source-reported coordinates; colour = log10 total ARG FPKM (colour range clipped to the 2nd-98th percentiles). Values are observed, not interpolated; no surface between points is implied. {ww.Project_Accession.value_counts().get('PRJEB27054', 0)} points come from the 2016 Global Sewage Surveillance pilot (mostly one sample per city); other studies contribute clustered points.",
         NB, "Objective 4", "derived/rf_unit_summary.csv")

# %%
pil["Continent"] = pil.Continent.fillna("Not reported")
cont = pil.groupby("Continent").agg(units=("total_FPKM", "size"), median_FPKM=("total_FPKM", "median"), Q1=("total_FPKM", lambda s: s.quantile(.25)),
        Q3=("total_FPKM", lambda s: s.quantile(.75)), median_gene_richness=("gene_richness", "median")).reset_index().sort_values("median_FPKM", ascending=False)
save_table(cont.pipe(rnd, 2), "S04_pilot_arg_load_by_continent", "2016 sewage pilot (PRJEB27054): total ARG FPKM by continent", NB, sheet="Spatial_Analysis")
pk = pil[pil.Continent.isin(cont.Continent[cont.units >= 5])]
k = kruskal_eps2(pk, "log10_total_FPKM", "Continent")
record_test(NB, "2016 sewage pilot: log10 total ARG FPKM by continent (continents with >=5 units)", "Kruskal-Wallis", k["n"], k["H"], k["df"], k["p"], "epsilon-squared", k["epsilon2"], family="Spatial", objective="Objective 4")
cfp = pd.read_csv(DER / "rf_class_fpkm.csv", index_col=0)
csh = (100 * cfp.loc[pil.index].div(cfp.loc[pil.index].sum(axis=1), axis=0)).join(pil.Continent).groupby("Continent").mean()
save_table(csh.pipe(rnd, 2).reset_index(), "S05_pilot_class_share_by_continent", "2016 sewage pilot: mean resistance-class share (%) by continent", NB, sheet="Spatial_Analysis")
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.0), gridspec_kw={"width_ratios": [.9, 1.3]})
order = cont.Continent.tolist()
sns.boxplot(data=pil, y="Continent", x="total_FPKM", order=order, color="#9ecae1", showfliers=False, linewidth=.7, ax=axes[0], log_scale=True)
sns.stripplot(data=pil, y="Continent", x="total_FPKM", order=order, color="k", size=2.5, alpha=.6, ax=axes[0])
axes[0].set_yticks(range(len(order))); axes[0].set_yticklabels([f"{c} (n={int(cont.set_index('Continent').loc[c,'units'])})" for c in order]); axes[0].set_ylabel("")
axes[0].set_xlabel("Total ARG FPKM (log scale)"); axes[0].set_title(f"A  normalized ARG abundance by continent\nKruskal-Wallis H={k['H']:.1f}, p={k['p']:.2g}, ε²={k['epsilon2']:.2f} (continents n≥5)", fontsize=7.2)
cs = csh.loc[order, [c for c in CLASS_ORDER if c in csh.columns]]
left = np.zeros(len(cs))
for c in cs.columns:
    axes[1].barh(range(len(cs)), cs[c], left=left, color=CLASS_COLORS[c], label=c, edgecolor="white", linewidth=.3); left += cs[c].values
axes[1].set_yticks(range(len(cs))); axes[1].set_yticklabels([]); axes[1].invert_yaxis(); axes[1].set_xlim(0, 100); axes[1].set_xlabel("Mean share of total ARG FPKM (%)")
axes[1].legend(frameon=False, fontsize=6.2, bbox_to_anchor=(1, 1), loc="upper left"); axes[1].set_title("B  Resistance-class composition", fontsize=7.5)
fig.tight_layout()
save_fig(fig, 25, "pilot_continent_resistome", "Regional variation in the sewage resistome, 2016 Global Sewage Surveillance pilot (ResFinder/KMA)",
         f"(A) Total ARG FPKM per city sample by continent (PRJEB27054; technical accessions merged to biological samples, mostly one per city; {int((~pil.index.str.startswith('BIO:')).sum())} accessions without a biological-sample label kept as separate units; n = {len(pil)}). (B) Mean class composition. Within-study comparison: one protocol and pipeline. Continents are coarse aggregates of heterogeneous countries.",
         NB, "Objective 2 / 4", "tables/S04, S05")
cont

# %% [markdown]
# ## 3. Spatial autocorrelation and distance decay (2016 pilot)

# %%
def haversine(lat, lon):
    la, lo = np.radians(lat)[:, None], np.radians(lon)[:, None]
    d = np.sin((la - la.T) / 2) ** 2 + np.cos(la) * np.cos(la.T) * np.sin((lo - lo.T) / 2) ** 2
    return 2 * 6371 * np.arcsin(np.sqrt(np.clip(d, 0, 1)))
def morans_i(x, W):
    z = x - x.mean(); n = len(x); S0 = W.sum()
    return (n / S0) * (z @ W @ z) / (z @ z)
pp = pil.dropna(subset=["lat", "lon"]).copy()
Dg = haversine(pp.lat.values, pp.lon.values)
rng = np.random.default_rng(SEED); mres = []
cls_p = cfp.loc[pp.index]; share_bl = (cls_p.get("Beta-lactam", 0) / cls_p.sum(axis=1)).values
for var, x in [("log10 total ARG FPKM", pp.log10_total_FPKM.values), ("ARG gene richness", pp.gene_richness.values.astype(float)), ("Beta-lactam share of FPKM", share_bl)]:
    for kk in [5, 8]:
        W = np.zeros_like(Dg)
        for i in range(len(Dg)):
            d_ = Dg[i].copy(); d_[i] = np.inf  # AUDIT FIX AUD-02: exclude self explicitly (co-located units tie at distance 0)
            nb = np.argsort(d_, kind="stable")[:kk]; W[i, nb] = 1
        W = W / W.sum(axis=1, keepdims=True)
        I0 = morans_i(x, W); null = np.array([morans_i(rng.permutation(x), W) for _ in range(9999)])
        p = (np.sum(np.abs(null - null.mean()) >= abs(I0 - null.mean())) + 1) / (len(null) + 1)
        mres.append(dict(variable=var, weights=f"k-nearest neighbours (k={kk}), row-standardised", n=len(x), morans_I=I0, expected_I=-1 / (len(x) - 1),
                         null_mean=null.mean(), null_sd=null.std(), z=(I0 - null.mean()) / null.std(), p_value_two_sided=p, permutations=9999))
        record_test(NB, f"2016 pilot Moran's I: {var} (kNN k={kk})", "Global Moran's I (permutation)", len(x), I0, None, p, "Moran's I", I0, family="Spatial autocorrelation", objective="Objective 4")
mo = pd.DataFrame(mres)
gfp = pd.read_csv(DER / "rf_gene_fpkm.csv", index_col=0).loc[pp.index]
from scipy.spatial.distance import pdist, squareform
Dr = squareform(pdist(closure(gfp.values), "braycurtis"))
ids = [str(i) for i in range(len(pp))]
mt = mantel(DistanceMatrix(Dr, ids), DistanceMatrix(Dg, ids), method="spearman", permutations=9999, seed=SEED)
mo.loc[len(mo)] = dict(variable="Resistome Bray-Curtis vs great-circle distance", weights="Mantel (Spearman)", n=int(mt[2]), morans_I=np.nan, expected_I=np.nan, null_mean=np.nan, null_sd=np.nan, z=np.nan, p_value_two_sided=mt[1], permutations=9999)
mo["mantel_r"] = [np.nan] * (len(mo) - 1) + [mt[0]]
record_test(NB, "2016 pilot distance decay: resistome Bray-Curtis vs geographic distance", "Mantel (Spearman)", int(mt[2]), mt[0], None, mt[1], "Mantel r", mt[0], family="Spatial autocorrelation", objective="Objective 4")
save_table(mo.pipe(rnd, 5), "S06_spatial_autocorrelation", f"Spatial autocorrelation (global Moran's I) and distance decay (Mantel) in the 2016 sewage pilot (n = {len(pp)} cities)", NB, sheet="Spatial_Analysis")
iu = np.triu_indices(len(pp), 1)
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.0))
axes[0].scatter(Dg[iu] / 1000, Dr[iu], s=3, alpha=.25, color="#4C72B0", edgecolors="none")
bins = np.arange(0, 20.5, 1.0); mids = []; med = []
for a, b in zip(bins[:-1], bins[1:]):
    m_ = (Dg[iu] / 1000 >= a) & (Dg[iu] / 1000 < b)
    if m_.sum() > 10: mids.append((a + b) / 2); med.append(np.median(Dr[iu][m_]))
axes[0].plot(mids, med, "o-", color="#C44E52", ms=3, label="median per 1,000 km bin"); axes[0].legend(frameon=False, fontsize=6.5)
axes[0].set_xlabel("Great-circle distance between cities (1,000 km)"); axes[0].set_ylabel("Resistome Bray-Curtis dissimilarity")
axes[0].set_title(f"A  Distance decay (Mantel r = {mt[0]:.2f}, p = {mt[1]:.4f})", fontsize=7.5)
r0 = mo[(mo.variable == "log10 total ARG FPKM") & mo.weights.str.contains("k=5")].iloc[0]
W = np.zeros_like(Dg)
for i in range(len(Dg)):
    d_ = Dg[i].copy(); d_[i] = np.inf; W[i, np.argsort(d_, kind="stable")[:5]] = 1  # AUDIT FIX AUD-02
W = W / W.sum(axis=1, keepdims=True)
x = pp.log10_total_FPKM.values; z = (x - x.mean()) / x.std()
axes[1].scatter(z, W @ z, s=12, color="#55A868", edgecolors="k", linewidths=.2)
b = np.polyfit(z, W @ z, 1); xx = np.linspace(z.min(), z.max(), 10); axes[1].plot(xx, np.polyval(b, xx), color="k", lw=.8)
axes[1].axhline(0, color="grey", lw=.4); axes[1].axvline(0, color="grey", lw=.4)
axes[1].set_xlabel("Standardised log10 total ARG FPKM"); axes[1].set_ylabel("Spatial lag (mean of 5 nearest units)")
axes[1].set_title(f"B  Moran scatterplot (I = {r0.morans_I:.2f}, permutation p = {r0.p_value_two_sided:.4f})", fontsize=7.5)
fig.tight_layout()
save_fig(fig, 26, "pilot_spatial_autocorrelation", "Spatial structure of the sewage resistome among cities in the 2016 pilot",
         f"(A) Pairwise resistome Bray-Curtis dissimilarity versus great-circle distance for {len(pp)} cities, with bin medians; Mantel test with 9,999 permutations. (B) Moran scatterplot of standardised log10 total ARG FPKM against its spatial lag (row-standardised 5-nearest-neighbour weights). Mostly one sample per city; results describe between-city spatial structure in 2016 and are not a predictive surface.",
         NB, "Objective 4", "tables/S06_spatial_autocorrelation.csv")
mo.round(4)

# %% [markdown]
# ## 4. Feasibility assessment of spatial methods

# %%
feas = pd.DataFrame([
    dict(method="Point maps of observed values", data_requirement="Producer-reported coordinates", available=f"{len(pts):,} included samples; {len(ww)} wastewater ARG units", decision="Performed (Figures 22, 24)"),
    dict(method="Country-level choropleth", data_requirement="Country per sample", available=f"{cc.shape[0]} countries", decision="Performed (Figure 23)"),
    dict(method="Global Moran's I / Mantel distance decay", data_requirement="Comparable measurements at >=30 locations from one protocol", available=f"{len(pp)} cities in the 2016 pilot (one study)", decision="Performed for the pilot only (Figure 26)"),
    dict(method="Local hotspot statistics (LISA, Getis-Ord Gi*)", data_requirement="Dense, comparable sampling within regions; multiple-testing control", available="One sample per city; sparse in Africa, South America, Oceania", decision="Not performed: density insufficient, hotspots would be artefacts of sampling gaps"),
    dict(method="Interpolated surfaces (kriging, IDW) / risk maps", data_requirement="Dense sampling, validated variogram, outcome data defining 'risk'", available="No such data available", decision="Not performed: would imply measurements where none exist"),
    dict(method="City-level mapping of DS02 Indian cities", data_requirement="Site coordinates", available="City labels only (19 sites unidentified)", decision="City used as unit (bar/box plots); site coordinates require Singh et al. Supplementary Data S7"),
])
save_table(feas, "S07_spatial_feasibility", "Feasibility of spatial analyses with the available data", NB, sheet="Spatial_Analysis")
record_finding(NB, "F15", "Objective 4", f"Among {len(pp)} units in the 2016 sewage pilot (mostly one per city), total normalized ARG abundance showed {'positive' if r0.morans_I>0 else 'negative'} spatial autocorrelation (Moran's I = {r0.morans_I:.2f}, permutation p = {r0.p_value_two_sided:.4f}, k = 5) and resistome composition showed distance decay (Mantel r = {mt[0]:.2f}, p = {mt[1]:.4f}); median normalized ARG abundance differed by continent (epsilon-squared = {k['epsilon2']:.2f}, p = {k['p']:.2g}).",
               len(pp), "2016 pilot units with coordinates (mostly one biological sample per city)", "Moran's I and Mantel with 9,999 permutations; Kruskal-Wallis", "Tables S04-S06; Figures 24-26",
               ("Geographically closer cities tend to have more similar sewage resistomes, consistent with regional structuring; regional summaries are more defensible than city rankings." if (mt[0] > 0 and mt[1] < .05) else "No clear distance-decay of resistome composition was detected at this sampling density.")
               + (" normalized ARG abundance is spatially clustered among neighbouring cities." if (r0.morans_I > 0 and r0.p_value_two_sided < .05) else " normalized ARG abundance shows no significant spatial clustering at k = 5."),
               "Single sampling round (2016), mostly one sample per city, uneven continental coverage; no outcome data; hotspot detection and interpolation not justified.")
record_finding(NB, "F16", "Objective 4", f"Only {len(pts):,} of {len(inc):,} included samples ({100*len(pts)/len(inc):.0f}%) have producer-reported coordinates, and processed AMR tables exist for samples in {int((cc.with_AMR_table>0).sum())} of {cc.shape[0]} countries; DS02 has no site coordinates.",
               len(inc), "Included samples", "Descriptive", "Tables S01, S03, S07; Figures 22-23",
               "Global AMR 'risk mapping' at city resolution is not supported by the processed data; observed-value maps and regional summaries are the defensible outputs.", "Coordinates are source-reported and may be approximate (precision in Table S02).", "Descriptive")
