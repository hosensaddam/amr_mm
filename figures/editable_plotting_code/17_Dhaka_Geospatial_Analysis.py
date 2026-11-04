# %% [markdown]
# # 17 · Dhaka geographic verification and sampling coverage (AUD-18)
# **Research objective:** verify where the Dhaka samples were taken, distinguish Dhaka-labelled, geographically verified and boundary-confirmed samples, and map Dhaka-centred sampling coverage in South Asia.
#
# **Inputs.** Updated sample master (both coordinate sources retained for SAMEA4777503: ENA sample-level record and the publication city-level value — never averaged); official administrative boundaries of Bangladesh (BBS/OCHA common operational dataset, COD-AB, admin level 3 v03, valid from 21 May 2023; Dhaka North and Dhaka South City Corporation polygons) obtained from the Humanitarian Data Exchange (point-in-polygon results computed in the processing workspace and stored in `derived/reprocessed/dhaka_geographic_verification.csv`); Natural Earth 1:50m countries.
#
# **Definitions.** *Dhaka-labelled*: source geo_loc_name/locality or city = Dhaka. *Geographically verified*: sample-level coordinates reported in the source record with >= 3 decimals. *Confirmed within Dhaka City*: verified coordinates inside the DNCC or DSCC polygon and farther from the boundary than the coordinate precision. *Uncertain*: city-level or 2-decimal coordinates, or a point closer to the boundary than its precision.
#
# **Limitations.** No sample location is inferred; the coordinate precision of 2-decimal values is about ±0.55 km; spatial statistics are not attempted (16 samples, at most a handful of distinct sites per dataset). No hotspot or risk surface is produced.

# %%
#BOOTSTRAP
NB = "NB17"
nb_header(NB)
from dhaka_common import *
import geopandas as gpd
from shapely.geometry import Point
meta_all = load_run_meta(); meta = meta_all[meta_all.processing_status == "completed"]
geo = pd.read_csv(REP / "dhaka_geographic_verification.csv")
dh = meta_all[meta_all.analysis_group == "Dhaka"].reset_index().set_index("Sample_ID")
print(geo.groupby("coordinate_source").size())

# %%
prim = geo[geo.coordinate_source == "primary (harmonized)"].set_index("Sample_ID")
def prec_km(p):
    return 0.55 if "<=2" in str(p) else 0.06
cls = []
for s, r in prim.iterrows():
    exact = "Exact" in str(r.precision); within = bool(r.within_Dhaka_city_corporations_DNCC_DSCC); pk = prec_km(r.precision)
    city_level = s in ("SAMEA7428327", "SAMEA7428511", "SAMEA7426329")
    if city_level: status = "Dhaka-labelled; location uncertain (publication city-level coordinate; sampling site not reported)"
    elif within and exact and r.distance_to_DNCC_DSCC_boundary_km > pk: status = "Confirmed within Dhaka City (DNCC/DSCC), sample-level coordinate"
    elif within and r.distance_to_DNCC_DSCC_boundary_km > pk: status = "Within Dhaka City (DNCC/DSCC); approximate (2-decimal) coordinate"
    elif within: status = "Within DNCC/DSCC but closer to the boundary than the coordinate precision - uncertain"
    else: status = "Outside DNCC/DSCC"
    alt = geo[(geo.Sample_ID == s) & (geo.coordinate_source == "alternative source")]
    cls.append(dict(Sample_ID=s, Dataset_ID=r.Dataset_ID, Sample_Type=dh.loc[s, "Sample_Type"], Dhaka_labelled=True, latitude=r.latitude, longitude=r.longitude, coordinate_precision=r.precision,
                    coordinate_source=dh.loc[s, "Coordinate_Provenance_Note"] if str(dh.loc[s, "Coordinate_Provenance_Note"]) not in ("Not applicable", "nan") else "NCBI BioSample lat_lon (sample-level)",
                    admin2_district=r.admin2_district, admin3_unit=r.admin3_upazila_or_city_corporation, within_DNCC_DSCC=within, distance_to_city_boundary_km=r.distance_to_DNCC_DSCC_boundary_km,
                    alternative_coordinate=(f"{alt.latitude.iloc[0]}, {alt.longitude.iloc[0]} -> {alt.admin3_upazila_or_city_corporation.iloc[0]}" if len(alt) else "none"),
                    geographic_status=status))
G = pd.DataFrame(cls).sort_values(["Dataset_ID", "Sample_ID"])
save_table(G.pipe(rnd, 5), "DH32_dhaka_geographic_verification", "AUD-18: Dhaka geographic verification - coordinate provenance and precision, point-in-polygon against official DNCC/DSCC boundaries (BBS/OCHA COD-AB adm3 v03), distance to the city boundary and classification", NB, sheet="Dhaka_Spatial_Analysis")
gs = G.geographic_status.value_counts().rename_axis("geographic_status").reset_index(name="samples")
save_table(gs, "DH33_dhaka_geographic_status_summary", "AUD-18: number of Dhaka samples by geographic status (16 Dhaka-labelled samples)", NB, sheet="Dhaka_Spatial_Analysis")
gs

# %%
# Pairwise distances between Dhaka sampling points (verified coordinates only; co-location <= 0.5 km)
pts = G[~G.geographic_status.str.startswith("Dhaka-labelled; location uncertain")].reset_index(drop=True)
def hav(a, b, c, d):
    R = 6371.0; p1, p2 = np.radians(a), np.radians(c); dp = p2 - p1; dl = np.radians(d - b)
    return 2 * R * np.arcsin(np.sqrt(np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2))
pr = []
for i in range(len(pts)):
    for j in range(i + 1, len(pts)):
        dkm = hav(pts.latitude[i], pts.longitude[i], pts.latitude[j], pts.longitude[j])
        if dkm <= 0.5: pr.append(dict(sample_a=pts.Sample_ID[i], dataset_a=pts.Dataset_ID[i], type_a=pts.Sample_Type[i], sample_b=pts.Sample_ID[j], dataset_b=pts.Dataset_ID[j], type_b=pts.Sample_Type[j], distance_km=dkm))
pr = pd.DataFrame(pr)
save_table(pr.pipe(rnd, 4), "DH34_dhaka_colocated_sampling_points", "Dhaka sampling points within 0.5 km of each other (verified coordinates): related site-level measurements across datasets", NB, sheet="Dhaka_Spatial_Analysis")
pr

# %%
adm = gpd.read_file(REP / "dhaka_area_admin3_subset.geojson")
city = adm[adm.adm3_name.str.contains("Dhaka") & adm.adm3_name.str.contains("City Corporation")]
fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.9), gridspec_kw={"width_ratios": [1, 1.05]})
ax = axes[0]
adm.plot(ax=ax, color="#F2F2F2", edgecolor="#BBBBBB", linewidth=.3)
city.plot(ax=ax, color="#FCE9DC", edgecolor="#D55E00", linewidth=.8)
for t_, col in TYPE_COL.items():
    g = G[G.Sample_Type == t_]
    for d_, mk in [("DS01", "o"), ("DS03", "s"), ("DS84", "^")]:
        gg = g[g.Dataset_ID == d_]
        if len(gg): ax.scatter(gg.longitude, gg.latitude, marker=mk, s=26, color=col, edgecolors="k", linewidths=.4, zorder=3)
alt = geo[geo.coordinate_source == "alternative source"]
ax.scatter(alt.longitude, alt.latitude, marker="x", s=22, color="k", lw=.8, zorder=4)
ax.annotate("DS01 city-level coordinate\n(3 samples; site not reported)", (90.38, 23.76), xytext=(90.275, 23.86), fontsize=6, arrowprops=dict(arrowstyle="-", lw=.5))
ax.annotate("Pagla sewage-treatment area\n(DS01 ENA + DS84 UWW1)", (90.4515, 23.683), xytext=(90.47, 23.615), fontsize=6, arrowprops=dict(arrowstyle="-", lw=.5))
ax.set_xlim(90.25, 90.56); ax.set_ylim(23.6, 23.92); ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)"); ax.set_title("A  Dhaka sampling points (DNCC/DSCC outlined)")
from matplotlib.lines import Line2D
hd = [Line2D([], [], marker=m, ls="", mfc="w", mec="k", ms=4, label=l) for m, l in [("o", "DS01"), ("s", "DS03"), ("^", "DS84"), ("x", "alternative source value")]] + [Line2D([], [], marker="o", ls="", mfc=c, mec="k", ms=4, label=t) for t, c in TYPE_COL.items()]
ax.legend(handles=hd, frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(.5, -.16), ncol=2)
ax = axes[1]
world = load_world(); sa = world[world.ISO2.isin(["BD", "IN", "PK", "NP", "LK", "BT", "MM", "CN", "AF"])]
sa.plot(ax=ax, color="#F2F2F2", edgecolor="#AAAAAA", linewidth=.3)
cm = meta[~meta.analysis_group.str.startswith("Calibration")].copy()
cm["lat"] = pd.to_numeric(cm.Latitude, errors="coerce"); cm["lon"] = pd.to_numeric(cm.Longitude, errors="coerce")
locs = {"Delhi": (28.61, 77.21)}
for c, g in cm.groupby("City_group"):
    la, lo = (g.lat.median(), g.lon.median()) if g.lat.notna().any() else locs.get(c, (np.nan, np.nan))
    ax.scatter(lo, la, s=12 + 6 * len(g), color=CITY_COL.get(c, "k"), edgecolors="k", linewidths=.4, zorder=3, alpha=.9)
    ax.annotate(f"{c} ({len(g)})", (lo, la), xytext=(4, 4), textcoords="offset points", fontsize=6)
ax.set_xlim(60, 97); ax.set_ylim(5, 36); ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)"); ax.set_title("B  Uniformly reprocessed runs by city")
fig.tight_layout()
save_fig(fig, 43, "dhaka_sampling_map", "Dhaka sampling locations against official city boundaries and South Asian comparator coverage",
         "(A) Coordinates of the 16 Dhaka samples as reported in the source records (marker = dataset; colour = verified sample type), over BBS/OCHA COD-AB admin-3 units with Dhaka North and Dhaka South City Corporations outlined; x = the alternative coordinate retained for SAMEA4777503 (publication city-level value; both values kept, not averaged). Three DS01 samples have only the city-level coordinate. (B) Cities contributing uniformly reprocessed runs (symbol size and label = number of runs; Delhi plotted at the city centre because DS02 records give site names but no coordinates). Observed sampling locations only; no interpolation, hotspot or risk surface.",
         NB, "Dhaka spatial (AUD-18)", "tables/DH32-DH34")

# %%
n_conf = int(G.geographic_status.str.startswith("Confirmed").sum()); n_unc = int(G.geographic_status.str.contains("uncertain").sum())
record_finding(NB, "DF09", "Dhaka spatial (AUD-18)", f"All 16 Dhaka-labelled samples have reported coordinates inside the official Dhaka North/South City Corporation boundaries (COD-AB adm3 v03); {n_conf} are confirmed with sample-level coordinates, {int(G.geographic_status.str.contains('approximate').sum())} lie inside with 2-decimal coordinates, and {n_unc} remain uncertain because only the publication city-level coordinate exists. SAMN60675067 (DS84 UWW1) lies {G.set_index('Sample_ID').loc['SAMN60675067','distance_to_city_boundary_km']:.2f} km inside DSCC near the Pagla sewage-treatment area, about 0.13 km from the ENA coordinate of DS01 SAMEA4777503.",
               16, "Dhaka-labelled samples", "Point-in-polygon against official boundaries; haversine distances", "Tables DH32-DH34; Figure 43",
               "Dhaka sample locations are consistent with Dhaka City; DS03 and DS84 re-sampled two hospital sites and DS01/DS84 the Pagla sewage-treatment area, so these are related site-level measurements.",
               "Boundary polygons are administrative; 2-decimal coordinates carry ~0.55 km uncertainty; three DS01 sites unknown.", kind="Empirical (verification)")
