# %% [markdown]
# # 01 · Data Inspection and Quality Control
# **Chapter:** Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping (Academic Press, Elsevier)
#
# **Research objective:** prerequisite for Objectives 1–5. Establishes which records are analysable, how analytic units are defined, and which data-quality problems constrain every later inference.
#
# **Analytical rationale:** the harmonized tables were assembled from several repositories and two independent AMR pipelines. Before any statistic is computed, the audit checks structure, missingness, duplicates, identifier linkage, geographic validity, sample-type classification, sequencing characteristics, annotation completeness and agreement between the CSV tables, the Excel workbooks and the counts reported in the data-preparation README.
#
# **Inputs:** `AMR_Chapter_Reproducibility_Bundle/AMR_Chapter_Public_Data_Collection/harmonized/*.csv` and the five project workbooks (read only).
#
# **Integrity rules applied:** original files are never modified; source missing-value codes (*Not reported*, *Not available*, *Not yet verified*, *Not applicable*, *Restricted access*) are treated as missing and never converted to zero; no coordinate, date, city or AMR value is imputed.

# %%
#BOOTSTRAP
NB = "NB01"
nb_header(NB)
import json, zipfile, urllib.request, hashlib
import seaborn as sns

# %% [markdown]
# ## 1. Directory inventory and file integrity
# Every file under the research directory is listed with size and MD5 checksum. Byte-identical copies are identified so that duplicated workbooks (files with the `_1` suffix) are not counted twice.

# %%
rows = []
for p in sorted(ROOT.rglob("*")):
    if p.is_file() and "Analysis_Results" not in p.parts:
        rows.append(dict(relative_path=str(p.relative_to(ROOT)), extension=p.suffix.lower() or "(none)",
                         size_bytes=p.stat().st_size, md5=md5(p)))
inv = pd.DataFrame(rows)
dup_groups = inv[inv.duplicated("md5", keep=False)].groupby("md5").relative_path.apply(lambda s: " | ".join(sorted(s)))
inv["identical_copy_group"] = inv.md5.map(dup_groups).fillna("")
top = inv.relative_path.str.split(r"[\\/]").str[0]
inv["location"] = np.where(inv.relative_path.str.contains(r"[\\/]"), top, "(root)")
save_table(inv, "QC01_file_inventory", "Inventory of all files in the research directory", NB,
           "Path, type, size and MD5 of every file (Analysis_Results excluded).", sheet="Dataset_Inventory")
ftype = inv.groupby("extension").agg(files=("relative_path", "size"), total_MB=("size_bytes", lambda s: round(s.sum() / 1e6, 2))).reset_index()
save_table(ftype, "QC02_file_types", "File counts by type", NB, sheet="Dataset_Inventory")
n_dupsets = inv.loc[inv.identical_copy_group != "", "md5"].nunique()
print(f"{len(inv)} files inspected; {n_dupsets} sets of byte-identical copies")
print(inv.loc[inv.identical_copy_group != "", ["relative_path", "size_bytes"]].to_string())

# %% [markdown]
# ## 2. Structure of the harmonized tables

# %%
tables = {n: load_csv(n) for n in ["sample_master.csv", "run_master.csv", "taxonomy_long.csv", "amr_genes_long.csv",
                                   "diversity_long.csv", "env_covariates_long.csv", "biosample_attributes.csv"]}
tables["subset_6919377_ARG_core_runs.tsv"] = pd.read_csv(HARM / "subset_6919377_ARG_core_runs.tsv", sep="\t")
tables["subset_6919377_diversity_core_runs.tsv"] = pd.read_csv(HARM / "subset_6919377_diversity_core_runs.tsv", sep="\t")
struct = pd.DataFrame([dict(table=k, rows=len(v), columns=v.shape[1], memory_MB=round(v.memory_usage(deep=True).sum() / 1e6, 1))
                       for k, v in tables.items()])
save_table(struct, "QC03_table_structure", "Dimensions of harmonized analysis tables", NB, sheet="Dataset_Inventory")
struct

# %%
sm, rm = tables["sample_master.csv"], tables["run_master.csv"]
dd = []
for tname in ["sample_master.csv", "run_master.csv", "taxonomy_long.csv", "amr_genes_long.csv", "diversity_long.csv", "env_covariates_long.csv"]:
    t = tables[tname]
    for c in t.columns:
        s = t[c]; miss = is_missing(s)
        num = pd.to_numeric(s.where(~miss), errors="coerce")
        inferred = "numeric" if (num.notna().sum() == (~miss).sum() and (~miss).sum() > 0) else "text/categorical"
        codes = {k: int((s.astype(str) == k).sum()) for k in MISSING_CODES if (s.astype(str) == k).any()}
        ex = s[~miss].astype(str).unique()[:3]
        dd.append(dict(table=tname, column=c, pandas_dtype=str(s.dtype), inferred_type=inferred, n_rows=len(s),
                       n_observed=int((~miss).sum()), n_missing=int(miss.sum()), pct_missing=round(100 * miss.mean(), 2),
                       missing_codes=json.dumps(codes), n_unique=int(s[~miss].nunique()), example_values=" | ".join(ex)[:200]))
dd = pd.DataFrame(dd)
save_table(dd, "QC04_data_dictionary", "Data dictionary with observed/missing counts per column", NB,
           "Missing = blank or one of the five source missing-value codes.", sheet="Data_Dictionary")

# %% [markdown]
# ## 3. Inclusion status, pending studies and exclusion log

# %%
incl = sm.Inclusion_Status.value_counts().rename_axis("Inclusion_Status").reset_index(name="samples")
reasons = sm.groupby(["Inclusion_Status", "Inclusion_Reason"]).size().rename("samples").reset_index().sort_values("samples", ascending=False)
save_table(reasons, "QC05_inclusion_reasons", "Sample-level inclusion status and reasons (sample master)", NB, sheet="Data_Quality_Audit")
src = (BUNDLE / "scripts" / "curation_decisions.py").read_text(encoding="utf-8")
blk = src[src.index("PENDING = {"):]; blk = blk[:blk.index("\n}")]
pending = re.findall(r'"((?:PRJ[A-Z]{2}\d+)|(?:mgp\d+))":\s*"([^"]+)"', blk)
mg = re.search(r'"(mgp87891)": \(P, [A-Z]+, "([^"]+)"\)', src)
if mg: pending.append((mg.group(1), mg.group(2)))
pend = pd.DataFrame(pending, columns=["study_accession", "reason_pending"])
ids = json.loads((HARM / "dataset_ids.json").read_text())
pend["present_in_sample_master"] = pend.study_accession.isin(ids.keys())
pend["treatment_in_this_analysis"] = "Not analysed (eligibility pending author decision; no samples harmonized)"
save_table(pend, "QC06_pending_studies", "Studies with pending eligibility decisions and their treatment", NB, sheet="Data_Quality_Audit")
for _, r in reasons[reasons.Inclusion_Status != "Included"].iterrows():
    record_exclusion(NB, "All analyses (sample master)", "sample", r.samples, f"{r.Inclusion_Status}: {r.Inclusion_Reason}")
record_exclusion(NB, "All analyses", "study", len(pend), "Eligibility pending author decision (not harmonized)", ";".join(pend.study_accession[:5]))
print(incl.to_string(index=False)); print(f"Pending studies: {len(pend)}; present in sample master: {pend.present_in_sample_master.sum()}")

# %% [markdown]
# ## 4. Duplicate samples, runs and probable duplicate physical samples

# %%
dup = []
dup.append(dict(check="Duplicate Sample_ID rows in sample master", n=int(sm.Sample_ID.duplicated().sum())))
dup.append(dict(check="Duplicate Run_ID rows in run master", n=int(rm.Run_ID.duplicated().sum())))
dup.append(dict(check="Runs linked to >1 sample", n=int((rm.groupby("Run_ID").Sample_ID.nunique() > 1).sum())))
dup.append(dict(check="Samples assigned to >1 dataset", n=int((sm.groupby("Sample_ID").Dataset_ID.nunique() > 1).sum())))
dup.append(dict(check="Duplicate Secondary_Sample_ID (observed values)", n=int(sm.loc[~is_missing(sm.Secondary_Sample_ID), "Secondary_Sample_ID"].duplicated().sum())))
d1 = sm[sm.Dataset_ID == "DS01"]
bio = d1[~is_missing(d1.Biological_Sample_ID)]
dup.append(dict(check="DS01 accessions sharing a Biological_Sample_ID (technical/reproducibility accessions)", n=int(bio.Biological_Sample_ID.duplicated(keep=False).sum())))
dup.append(dict(check="DS01 distinct Biological_Sample_ID", n=int(bio.Biological_Sample_ID.nunique())))
# PRJEB13831 (Hendriksen et al. 2019 pilot) vs PRJEB27054 (Munk et al. 2022 2016 pilot): match by site code from sample alias
alias = rm.drop_duplicates("Sample_ID").set_index("Sample_ID").Sample_Alias
def site_code(a):
    m = re.match(r"([A-Z]{3})[.\-_ ]?(\d+)", str(a)); return f"{m.group(1)}-{int(m.group(2))}" if m else None
pp = sm[sm.Project_Accession.isin(["PRJEB13831", "PRJEB27054"])][["Sample_ID", "Project_Accession", "Country", "Collection_Date", "Latitude", "Longitude", "Biological_Sample_ID"]].copy()
pp["alias"] = pp.Sample_ID.map(alias); pp["site_code"] = pp.alias.map(site_code)
a13 = pp[pp.Project_Accession == "PRJEB13831"]; a27 = pp[pp.Project_Accession == "PRJEB27054"]
mm = a13.merge(a27, on="site_code", suffixes=("_PRJEB13831", "_PRJEB27054"))
mm["same_country"] = mm.Country_PRJEB13831 == mm.Country_PRJEB27054
def dm_swap(a, b):
    try: return a[5:7] == b[8:10] and a[8:10] == b[5:7] and a[:4] == b[:4]
    except Exception: return False
mm["date_identical"] = mm.Collection_Date_PRJEB13831 == mm.Collection_Date_PRJEB27054
mm["date_day_month_transposed"] = [dm_swap(a, b) for a, b in zip(mm.Collection_Date_PRJEB13831, mm.Collection_Date_PRJEB27054)]
save_table(mm, "QC07_probable_duplicate_physical_samples", "PRJEB13831 samples matching PRJEB27054 2016-pilot site codes", NB,
           "Site code parsed from sample alias (e.g., AUS.18a / AUS-18_1 -> AUS-18).", sheet="Data_Quality_Audit")
n13 = a13.Sample_ID.nunique(); n13m = mm.Sample_ID_PRJEB13831.nunique()
dup.append(dict(check="PRJEB13831 samples whose site code matches a PRJEB27054 pilot sample (probable same physical sewage sample)", n=n13m))
dup.append(dict(check="  ...of which same country", n=int(mm.drop_duplicates("Sample_ID_PRJEB13831").same_country.sum())))
dup.append(dict(check="  ...of which identical collection date", n=int(mm.drop_duplicates("Sample_ID_PRJEB13831").date_identical.sum())))
dup.append(dict(check="  ...of which day and month transposed between projects", n=int(mm.drop_duplicates("Sample_ID_PRJEB13831").date_day_month_transposed.sum())))
dup = pd.DataFrame(dup)
save_table(dup, "QC08_duplicate_checks", "Duplicate and identifier-uniqueness checks", NB, sheet="Data_Quality_Audit")
print(dup.to_string(index=False)); print(f"PRJEB13831 samples: {n13}")

# %% [markdown]
# **Decision.** PRJEB13831 and PRJEB27054 share site codes for most 2016 pilot sites, so PRJEB13831 very probably re-uses the same physical sewage samples. To avoid double counting, PRJEB27054 (which carries the Munk et al. city labels and biological-sample grouping) is used in primary analyses and PRJEB13831 is held out for a sensitivity analysis. Accessions sharing a `Biological_Sample_ID` are aggregated to one analytic unit by summing read counts.

# %% [markdown]
# ## 5. Identifier linkage across tables

# %%
sm_ids, rm_ids = set(sm.Sample_ID), set(rm.Run_ID)
amr, tax, div = tables["amr_genes_long.csv"], tables["taxonomy_long.csv"], tables["diversity_long.csv"]
sub_arg, sub_div = tables["subset_6919377_ARG_core_runs.tsv"], tables["subset_6919377_diversity_core_runs.tsv"]
sm_runs = set(";".join(sm.Run_ID.astype(str)).split(";")) - set(MISSING_CODES)
link = []
def L(check, n_bad, n_total): link.append(dict(check=check, n_unlinked=int(n_bad), n_total=int(n_total), status="Pass" if n_bad == 0 else "Fail"))
L("run_master Sample_ID present in sample_master", (~rm.Sample_ID.isin(sm_ids)).sum(), len(rm))
L("sample_master Run_ID entries present in run_master (MG-RAST included)", len(sm_runs - rm_ids), len(sm_runs))
for nm, t in [("amr_genes_long", amr), ("taxonomy_long", tax), ("diversity_long", div)]:
    L(f"{nm} Sample_ID present in sample_master", (~t.Sample_ID.isin(sm_ids)).sum(), len(t))
    L(f"{nm} Run_ID present in run_master", (~t.Run_ID.isin(rm_ids)).sum(), len(t))
    pair = t[["Sample_ID", "Run_ID"]].drop_duplicates().merge(rm[["Sample_ID", "Run_ID"]], how="left", indicator=True)
    L(f"{nm} Sample_ID-Run_ID pair matches run_master", (pair._merge != "both").sum(), len(pair))
    dsm = t[["Sample_ID", "Dataset_ID"]].drop_duplicates().merge(sm[["Sample_ID", "Dataset_ID"]], on="Sample_ID", suffixes=("", "_sm"))
    L(f"{nm} Dataset_ID agrees with sample_master", (dsm.Dataset_ID != dsm.Dataset_ID_sm).sum(), len(dsm))
k_runs = set(amr.loc[amr.AMR_Reference_Database == "ResFinder", "Run_ID"])
L("ResFinder runs in amr_genes_long also in Zenodo ARG subset", len(k_runs - set(sub_arg.run_accession)), len(k_runs))
L("Zenodo ARG subset runs also in amr_genes_long", len(set(sub_arg.run_accession) - k_runs), sub_arg.run_accession.nunique())
core = set((HARM / "global_resistome_core_runs.txt").read_text().split())
L("Core resistome run list runs with ARG rows", len(core - set(sub_arg.run_accession)), len(core))
L("Core resistome run list runs with diversity rows", len(core - set(sub_div.run_accession)), len(core))
dsamp = set(tax.Sample_ID); asamp = set(amr.loc[amr.Dataset_ID == "DS02", "Sample_ID"])
L("DS02 taxonomy samples with CARD/RGI rows", len(dsamp - asamp), len(dsamp))
L("DS02 CARD/RGI samples with taxonomy rows", len(asamp - dsamp), len(asamp))
link = pd.DataFrame(link)
save_table(link, "QC09_identifier_linkage", "Identifier linkage checks across tables", NB, sheet="Data_Quality_Audit")
missing_arg_run = sorted(core - set(sub_arg.run_accession))
print(link.to_string(index=False)); print("Run(s) in core list without ARG rows:", missing_arg_run)

# %% [markdown]
# ## 6. Reference geography and geographic validation
# Country polygons: Natural Earth 1:50m Admin-0 countries (public domain), retrieved once and checksummed. They are used only to (i) check that source-reported coordinates fall in the stated country and (ii) assign continents for aggregation. No coordinates are created.

# %%
ne_zip = EXT / "ne_50m_admin_0_countries.zip"
if not ne_zip.exists():
    url = "https://naciscdn.org/naturalearth/50m/cultural/ne_50m_admin_0_countries.zip"
    urllib.request.urlretrieve(url, ne_zip)
    (EXT / "PROVENANCE.txt").write_text(f"Natural Earth 1:50m Admin 0 Countries\nURL: {url}\nRetrieved: {dt.datetime.now().isoformat(timespec='seconds')}\n"
                                        f"SHA256: {hashlib.sha256(ne_zip.read_bytes()).hexdigest()}\nLicence: public domain (naturalearthdata.com/about/terms-of-use)\n", encoding="utf-8")
print((EXT / "PROVENANCE.txt").read_text())
import geopandas as gpd
world = load_world()
iso_cont = world.set_index("ISO2")[["CONTINENT", "SUBREGION", "NAME"]]
iso_cont = iso_cont[~iso_cont.index.duplicated()]

# %%
g = sm.copy()
g["lat"] = to_num(g.Latitude); g["lon"] = to_num(g.Longitude)
has = g.lat.notna() & g.lon.notna()
geo = [dict(check="Samples with numeric latitude and longitude", n=int(has.sum()), of=len(g))]
geo.append(dict(check="Latitude outside [-90, 90]", n=int((g.lat.abs() > 90).sum()), of=int(has.sum())))
geo.append(dict(check="Longitude outside [-180, 180]", n=int((g.lon.abs() > 180).sum()), of=int(has.sum())))
geo.append(dict(check="Coordinates exactly (0, 0)", n=int(((g.lat == 0) & (g.lon == 0)).sum()), of=int(has.sum())))
geo.append(dict(check="Only one of latitude/longitude reported", n=int((g.lat.notna() ^ g.lon.notna()).sum()), of=len(g)))
pts = gpd.GeoDataFrame(g[has], geometry=gpd.points_from_xy(g.loc[has, "lon"], g.loc[has, "lat"]), crs="EPSG:4326")
j = gpd.sjoin(pts[["Sample_ID", "Dataset_ID", "Country", "ISO_Country_Code", "Inclusion_Status", "geometry"]],
              world[["ISO2", "NAME", "geometry"]], how="left", predicate="within")
j = j[~j.index.duplicated()]
stated = j.ISO_Country_Code.where(~is_missing(j.ISO_Country_Code))
j["status"] = np.select([stated.isna(), j.ISO2.isna(), stated == j.ISO2], ["Country not stated", "Outside all land polygons", "Inside stated country"], "Inside another country")
# distance from the stated country polygon for non-matching points (equal-area projection, km)
wproj = world.to_crs("EPSG:6933").set_index("ISO2"); jp = j.to_crs("EPSG:6933")
chk = j.status.isin(["Outside all land polygons", "Inside another country"]) & stated.notna()
dist = []
for ix in j.index[chk]:
    iso = stated.loc[ix]
    dist.append(wproj.loc[[iso]].geometry.distance(jp.loc[ix, "geometry"]).min() / 1000 if iso in wproj.index else np.nan)
j.loc[chk, "km_from_stated_country"] = dist
j["flag"] = np.where(chk & (j.km_from_stated_country > 25), "Coordinate-country mismatch >25 km", np.where(chk, "Within 25 km of stated country (coastline/simplification tolerance)", ""))
geo_status = j.groupby(["Inclusion_Status", "status"]).size().rename("samples").reset_index()
mis = j[j.flag == "Coordinate-country mismatch >25 km"].drop(columns="geometry")
save_table(geo_status, "QC10_coordinate_country_consistency", "Source-reported coordinates versus stated country (Natural Earth 1:50m)", NB, sheet="Data_Quality_Audit")
save_table(mis, "QC11_coordinate_mismatches", "Samples whose coordinates fall >25 km from the stated country", NB, sheet="Data_Quality_Audit")
geo.append(dict(check="Coordinate-country mismatch >25 km (flagged, not corrected)", n=len(mis), of=int(has.sum())))
kyiv = sm[sm.City.astype(str).str.contains("Kyiv|Kiev", case=False)]
geo.append(dict(check="Kyiv records with country 'Not yet verified' (source listed Poland or blank)", n=int((kyiv.Country == "Not yet verified").sum()), of=len(kyiv)))
geo.append(dict(check="Collection date 1970-01-01 placeholder", n=int(sm.Collection_Date.astype(str).str.startswith("1970-01-01").sum()), of=len(sm)))
geo = pd.DataFrame(geo)
save_table(geo, "QC12_geographic_validation", "Geographic and date validation checks", NB, sheet="Data_Quality_Audit")
print(geo.to_string(index=False)); print(mis.groupby(["Dataset_ID", "Country", "NAME"]).size().head(20))

# %% [markdown]
# ## 7. Sample-type classification audit
# Sample types were assigned at study level and refined by keyword rules. Source matrix text is screened for terms that would contradict the assigned type; hits are flagged for author review, not reclassified.

# %%
st = sm[sm.Inclusion_Status != "Excluded"]
rule = st.groupby(["Sample_Type", "Sample_Type_Rule"]).size().rename("samples").reset_index()
kw = {"Hospital wastewater": r"hospital", "Surface water": r"river|lake|stream|creek|surface water|estuar",
      "Air": r"\bair\b|aerosol|dust|pm2\.?5", "Urban soil": r"\bsoil\b", "Transit surface": r"subway|metro|station|railing|kiosk|bench",
      "Municipal wastewater": r"influent|municipal|wwtp|treatment plant"}
aud = []
for t_assigned, grp in st.groupby("Sample_Type"):
    srcv = grp.Sample_Type_Source_Value.astype(str).str.lower()
    for t_kw, pat in kw.items():
        if t_kw == t_assigned: continue
        hit = srcv.str.contains(pat, regex=True)
        if hit.any():
            aud.append(dict(assigned_type=t_assigned, keyword_suggests=t_kw, n_samples=int(hit.sum()),
                            example_source_values=" | ".join(srcv[hit].value_counts().index[:4])))
aud = pd.DataFrame(aud)
save_table(rule, "QC13_sample_type_rules", "Sample types by classification rule (included and contextual samples)", NB, sheet="Data_Quality_Audit")
save_table(aud, "QC14_sample_type_keyword_audit", "Source matrix text that may conflict with assigned sample type (for author review)", NB, sheet="Data_Quality_Audit")
aud.head(20)

# %% [markdown]
# ## 8. Sequencing platform, data type and depth

# %%
def platform_family(s):
    s = str(s).lower()
    for k, v in [("nanopore", "Oxford Nanopore"), ("promethion", "Oxford Nanopore"), ("minion", "Oxford Nanopore"), ("gridion", "Oxford Nanopore"),
                 ("pacbio", "PacBio"), ("sequel", "PacBio"), ("dnbseq", "DNBSEQ/MGI/BGISEQ"), ("bgiseq", "DNBSEQ/MGI/BGISEQ"), ("mgiseq", "DNBSEQ/MGI/BGISEQ"),
                 ("ion", "Ion Torrent"), ("illumina", "Illumina"), ("hiseq", "Illumina"), ("novaseq", "Illumina"), ("nextseq", "Illumina"), ("miseq", "Illumina")]:
        if k in s: return v
    return "Not reported" if s in [m.lower() for m in MISSING_CODES] or s in ("", "nan") else "Other"
sm["Platform_Family"] = sm.Sequencing_Platform.map(platform_family)
sm["bases"] = to_num(sm.Base_Count_Total)
inc = sm[sm.Inclusion_Status == "Included"].copy()
inc["log10_bases"] = np.log10(inc.bases.where(inc.bases > 0))
q = inc.groupby("Dataset_ID").log10_bases.quantile([.25, .75]).unstack()
inc = inc.join(q, on="Dataset_ID")
iqr = inc[0.75] - inc[0.25]
inc["depth_outlier_within_dataset"] = (inc.log10_bases < inc[0.25] - 1.5 * iqr) | (inc.log10_bases > inc[0.75] + 1.5 * iqr)
inc["shallow_lt_100Mbp"] = inc.bases < 1e8
depth = inc.groupby("Platform_Family").agg(samples=("Sample_ID", "size"), bases_reported=("bases", "count"),
        median_Gbp=("bases", lambda s: s.median() / 1e9), q1_Gbp=("bases", lambda s: s.quantile(.25) / 1e9),
        q3_Gbp=("bases", lambda s: s.quantile(.75) / 1e9), shallow_lt_100Mbp=("shallow_lt_100Mbp", "sum"),
        within_dataset_outliers=("depth_outlier_within_dataset", "sum")).reset_index()
save_table(depth, "QC15_sequencing_depth_by_platform", "Sequencing depth (total bases) by platform family, included samples", NB, sheet="Data_Quality_Audit")
outl = inc.loc[inc.depth_outlier_within_dataset | inc.shallow_lt_100Mbp, ["Sample_ID", "Dataset_ID", "Sample_Type", "Platform_Family", "bases", "depth_outlier_within_dataset", "shallow_lt_100Mbp"]]
save_table(outl, "QC16_depth_outliers", "Included samples flagged for extreme sequencing depth (flagged, not removed)", NB, sheet="Data_Quality_Audit")
print(depth.to_string(index=False))

# %% [markdown]
# ## 9. Annotation completeness and internal consistency of processed tables

# %%
card = amr[amr.AMR_Reference_Database == "CARD"]; rf = amr[amr.AMR_Reference_Database == "ResFinder"]
ann = []
def A(check, value, note=""): ann.append(dict(check=check, value=value, note=note))
A("CARD/RGI rows (DS02)", len(card)); A("CARD/RGI samples", card.Sample_ID.nunique())
A("CARD/RGI reference model types", "; ".join(card.Reference_Model_Type.unique()))
A("CARD/RGI rows with allele coverage < 50%", int((card.Average_Percent_Coverage < 50).sum()), f"{100 * (card.Average_Percent_Coverage < 50).mean():.1f}% of rows")
A("CARD database release stated", card.AMR_Database_Version.iloc[0])
A("ResFinder/KMA rows", len(rf)); A("ResFinder/KMA runs", rf.Run_ID.nunique()); A("ResFinder/KMA samples", rf.Sample_ID.nunique())
A("ResFinder/KMA KMA versions", "; ".join(sorted(rf.AMR_Detection_Type.unique())))
A("ResFinder rows with class 'Missing'", int((rf.Resistance_Class == "Missing").sum()))
A("Taxonomy (Bracken phylum) samples", tax.Sample_ID.nunique()); A("Taxonomy distinct phyla", tax.Taxon_Name.nunique())
rs = tax.groupby("Sample_ID").Relative_Abundance.sum()
A("Taxonomy: samples whose reported relative abundances sum to 1 +/- 0.01", int(((rs - 1).abs() <= 0.01).sum()), f"range {rs.min():.4f}-{rs.max():.4f}")
argd = sub_div[sub_div.category == "ARG"]
A("Zenodo ARG diversity rows with Simpson = 0 but Shannon > 0 (internally inconsistent)", int(((argd.Simpson == 0) & (argd.Shannon > 0)).sum()), f"of {len(argd)} ARG rows; ARG diversity is recomputed from counts")
gend = sub_div[sub_div.category == "Genera"]
A("Zenodo genus diversity rows with Simpson = 0 but Shannon > 0", int(((gend.Simpson == 0) & (gend.Shannon > 0)).sum()), f"of {len(gend)}")
A("Runs in core list without ARG rows (treated as missing, not zero)", len(missing_arg_run), ";".join(missing_arg_run))
ann = pd.DataFrame(ann)
save_table(ann, "QC17_annotation_completeness", "AMR and taxonomic annotation completeness and consistency", NB, sheet="Data_Quality_Audit")
cov = sm[sm.Inclusion_Status != "Excluded"].assign(
    has_card=lambda d: d.Sample_ID.isin(card.Sample_ID), has_resfinder=lambda d: d.Sample_ID.isin(rf.Sample_ID),
    has_taxonomy=lambda d: d.Sample_ID.isin(tax.Sample_ID), has_diversity=lambda d: d.Sample_ID.isin(div.Sample_ID),
    has_coords=lambda d: to_num(d.Latitude).notna() & to_num(d.Longitude).notna()
).groupby(["Inclusion_Status", "Dataset_ID", "Dataset_Name", "Sample_Type"]).agg(
    samples=("Sample_ID", "size"), CARD_RGI=("has_card", "sum"), ResFinder_KMA=("has_resfinder", "sum"),
    taxonomy=("has_taxonomy", "sum"), diversity=("has_diversity", "sum"), coordinates=("has_coords", "sum")).reset_index().sort_values("samples", ascending=False)
save_table(cov, "QC18_annotation_coverage_by_dataset", "Annotation and coordinate coverage by dataset and sample type", NB, sheet="Data_Quality_Audit")
ann

# %% [markdown]
# ## 10. Cross-file consistency: CSV tables, Excel workbooks and README counts

# %%
import openpyxl
wbchk = []
for wbn in ["00_Master_Data_Index.xlsx", "01_Public_Dataset_Inventory.xlsx", "02_Sample_Metadata_Master.xlsx",
            "03_Analysis_Ready_Metagenomic_Data.xlsx", "04_AMR_Reference_Database_Register.xlsx"]:
    wb = openpyxl.load_workbook(ROOT / wbn, read_only=True)
    for ws in wb.worksheets:
        wbchk.append(dict(workbook=wbn, sheet=ws.title, max_row=ws.max_row, max_col=ws.max_column))
    wb.close()
wbchk = pd.DataFrame(wbchk)
save_table(wbchk, "QC19_workbook_sheets", "Sheets and dimensions of the project workbooks", NB, sheet="Dataset_Inventory")
print(wbchk.to_string(index=False))

# %%
def find_sheet_rows(wbn, pattern):
    r = wbchk[(wbchk.workbook == wbn) & wbchk.sheet.str.contains(pattern, case=False, regex=True)]
    return (r.sheet.iloc[0], int(r.max_row.iloc[0])) if len(r) else (None, None)
cons = []
for wbn, pat, csvn in [("02_Sample_Metadata_Master.xlsx", "^Sample_Master$|Sample_Metadata|^Samples", "sample_master.csv"),
                       ("02_Sample_Metadata_Master.xlsx", "Run", "run_master.csv"),
                       ("03_Analysis_Ready_Metagenomic_Data.xlsx", "Taxonom", "taxonomy_long.csv"),
                       ("03_Analysis_Ready_Metagenomic_Data.xlsx", "AMR_Genes|AMR_Gene", "amr_genes_long.csv"),
                       ("03_Analysis_Ready_Metagenomic_Data.xlsx", "Divers", "diversity_long.csv"),
                       ("03_Analysis_Ready_Metagenomic_Data.xlsx", "Covariat|Environ", "env_covariates_long.csv")]:
    sh, mr = find_sheet_rows(wbn, pat)
    n_csv = len(tables[csvn])
    cons.append(dict(workbook=wbn, sheet=sh, csv=csvn, csv_rows=n_csv, sheet_max_row=mr,
                     note="sheet max_row includes header/title rows; split sheets possible" if sh else "no matching sheet found"))
cons = pd.DataFrame(cons)
idx = json.loads((HARM / "index_numbers.json").read_text())["index"]
incl_s = sm[sm.Inclusion_Status == "Included"]
recomp = {
    "Unique samples in sample master": len(sm), "Included unique samples": len(incl_s),
    "Countries represented (included samples)": incl_s.loc[~is_missing(incl_s.Country), "Country"].nunique(),
    "Included samples with source-reported coordinates": int((to_num(incl_s.Latitude).notna() & to_num(incl_s.Longitude).notna()).sum()),
    "Included samples with gene-level AMR table": int(incl_s.Sample_ID.isin(amr.Sample_ID).sum()),
    "Included samples with taxonomy table (phylum)": int(incl_s.Sample_ID.isin(tax.Sample_ID).sum()),
    "Included samples with diversity indices": int(incl_s.Sample_ID.isin(div.Sample_ID).sum()),
    "Verified unique samples (Verification_Status = Verified)": int((sm.Verification_Status == "Verified").sum()),
    "Datasets in sample master": sm.Dataset_ID.nunique(),
    "Cities represented (distinct City|Country pairs, source-labelled city, country known)":
        incl_s.loc[~is_missing(incl_s.City) & ~is_missing(incl_s.Country), ["City", "Country"]].drop_duplicates().shape[0],
}
rc = pd.DataFrame([dict(quantity=k, reported_in_index=idx.get(k), recomputed=v) for k, v in recomp.items()])
rc["agrees"] = [str(a).replace(",", "") == str(b) for a, b in zip(rc.reported_in_index, rc.recomputed)]
tc = idx["Proposal target cities present (Dhaka / Delhi / Karachi)"]
for city in ["Dhaka", "Delhi", "Karachi"]:
    lab = incl_s.City.astype(str).str.contains(city, case=False)
    loc = incl_s.Locality_As_Reported.astype(str).str.contains(city, case=False)
    n_any = int((lab | loc).sum())
    rc.loc[len(rc)] = dict(quantity=f"Included samples with '{city}' in City or Locality_As_Reported (index definition)", reported_in_index=tc,
                           recomputed=n_any, agrees=(f"{city}: {n_any}" in tc) or (f"{city}: {n_any} " in tc))
    rc.loc[len(rc)] = dict(quantity=f"  of which source-labelled City = {city}", reported_in_index="(not reported separately)", recomputed=int(lab.sum()), agrees="n/a")
save_table(cons, "QC20_csv_vs_workbook_rows", "Row counts: harmonized CSV tables versus workbook sheets", NB, sheet="Verification_Report")
save_table(rc, "QC21_reported_vs_recomputed_counts", "Counts in the data-preparation index recomputed from the CSV tables", NB, sheet="Verification_Report")
print(cons.to_string(index=False)); rc

# %% [markdown]
# ## 11. Missingness of core metadata (included samples)

# %%
core_fields = ["City", "Country", "Latitude", "Longitude", "Collection_Date", "Sample_Type", "Environmental_Matrix",
               "Sequencing_Platform", "Read_Type", "Base_Count_Total", "Sampling_Design", "Sampling_Location", "Urban_Setting", "Publication_ID"]
incl_s = sm[sm.Inclusion_Status == "Included"]
top_ds = incl_s.Dataset_ID.value_counts().index[:15]
grp = np.where(incl_s.Dataset_ID.isin(top_ds), incl_s.Dataset_ID, "Other (" + str(incl_s.Dataset_ID.nunique() - 15) + " datasets)")
M = pd.DataFrame({f: is_missing(incl_s[f]).values for f in core_fields}); M["group"] = grp
miss_ds = M.groupby("group").mean().mul(100)
miss_ds["n"] = M.groupby("group").size()
miss_ds = miss_ds.sort_values("n", ascending=False)
overall = pd.DataFrame({"field": core_fields, "n_missing": [int(M[f].sum()) for f in core_fields], "pct_missing": [round(100 * M[f].mean(), 1) for f in core_fields]})
save_table(overall, "QC22_metadata_missingness_overall", "Missingness of core metadata fields, included samples (n = %d)" % len(incl_s), NB, sheet="Metadata_Missingness")
save_table(miss_ds.reset_index().pipe(rnd, 1), "QC23_metadata_missingness_by_dataset", "Percentage missing by field and dataset, included samples", NB, sheet="Metadata_Missingness")
fig, ax = plt.subplots(figsize=(6.9, 4.6))
sns.heatmap(miss_ds[core_fields], cmap="Reds", vmin=0, vmax=100, linewidths=.4, linecolor="white", cbar_kws={"label": "% missing"}, ax=ax)
ax.set_yticklabels([f"{i} (n={n})" for i, n in zip(miss_ds.index, miss_ds.n)], rotation=0)
ax.set_xticklabels([c.replace("_", " ") for c in core_fields], rotation=40, ha="right"); ax.set_xlabel(""); ax.set_ylabel("")
ax.set_title("Missingness of core metadata among included samples")
save_fig(fig, 1, "metadata_missingness_heatmap", "Missingness of core metadata fields by dataset (included samples)",
         f"Percentage of included samples (n = {len(incl_s):,}) with a missing value for each core field, by dataset (15 largest datasets shown individually). Missing = blank or a source missing-value code; no values were imputed.",
         NB, "Prerequisite (all objectives)", "tables/QC23_metadata_missingness_by_dataset.csv")

# %% [markdown]
# ## 12. Analytic datasets derived for downstream notebooks
# Three analytic sets are defined here so that all notebooks use the same units:
# 1. **Metadata set** – included samples plus food-related contextual samples (flagged), with numeric coordinates, years, platform family and continent.
# 2. **DS02 set** – Indian metropolitan open-drain wastewater (Singh et al. 2026): Bracken phylum counts and CARD/RGI gene counts for 447 samples (one run each).
# 3. **ResFinder set** – runs covered by the Martiny et al. global resistome resource (Zenodo 6919377, ResFinder/KMA), aggregated to analytic units.

# %%
meta = sm[sm.Inclusion_Status != "Excluded"].copy()
meta["lat"] = to_num(meta.Latitude); meta["lon"] = to_num(meta.Longitude)
bad = ~(meta.lat.between(-90, 90) & meta.lon.between(-180, 180)) | ((meta.lat == 0) & (meta.lon == 0))
meta.loc[bad, ["lat", "lon"]] = np.nan
meta["year"] = to_num(meta.Collection_Year)
meta.loc[meta.Collection_Date.astype(str).str.startswith("1970-01-01"), "year"] = np.nan
meta["ISO2"] = meta.ISO_Country_Code.where(~is_missing(meta.ISO_Country_Code))
meta["Continent"] = meta.ISO2.map(iso_cont.CONTINENT).fillna("Not reported")
meta["Subregion"] = meta.ISO2.map(iso_cont.SUBREGION).fillna("Not reported")
meta["Stratum"] = np.where(meta.Inclusion_Status == "Included", "Urban environmental (included)", "Food-related (separate contextual dataset)")
meta["coord_country_flag"] = meta.Sample_ID.isin(mis.Sample_ID)
keep = ["Dataset_ID", "Dataset_Name", "Project_Accession", "Sample_ID", "Biological_Sample_ID", "Run_ID", "N_Runs", "Repository", "City", "Country", "ISO2",
        "Continent", "Subregion", "lat", "lon", "Coordinate_Precision", "Sample_Type", "Sample_Subtype", "Collection_Date", "year", "Collection_Month",
        "Sequencing_Platform", "Platform_Family", "Sequencing_Type", "Read_Type", "bases", "Verification_Status", "Inclusion_Status", "Stratum", "QC_Status", "QC_Flags", "coord_country_flag"]
meta[keep].to_csv(DER / "metadata_set.csv", index=False)
print("metadata set:", meta.shape, meta.Stratum.value_counts().to_dict())
print("unmatched ISO2 codes:", sorted(set(meta.ISO2.dropna()) - set(iso_cont.index)))

# %%
# DS02 set ---------------------------------------------------------------
d2 = sm[sm.Dataset_ID == "DS02"].copy()
d2["month"] = d2.Collection_Date.str[5:7].astype(int); d2["year"] = d2.Collection_Date.str[:4].astype(int)
d2["Season"] = d2.month.map(season_india)
d2 = d2.merge(rm[["Run_ID", "Read_Count", "Base_Count", "Read_Type"]].rename(columns={"Read_Type": "RT"}), on="Run_ID", how="left")
d2["reads_total"] = to_num(d2.Read_Count) * np.where(d2.RT == "PAIRED", 2, 1)
code = tax.drop_duplicates("Sample_ID").set_index("Sample_ID").Source_Sample_Code
kassigned = tax.drop_duplicates("Sample_ID").set_index("Sample_ID").Kraken_Assigned_Reads
d2["Source_Sample_Code"] = d2.Sample_ID.map(code)
d2["bracken_phylum_reads"] = d2.Sample_ID.map(tax.groupby("Sample_ID").Abundance_Value.sum())
# UPDATE (AUD-11): site identifiers recovered from ENA/SRA experiment titles (Site_ID in the updated sample master)
d2["Site"] = d2.Site_ID if "Site_ID" in d2.columns else "Not reported"
print("DS02 sites:", d2.Site.nunique(), d2.groupby("City").Site.nunique().to_dict())
d2[["Sample_ID", "Run_ID", "Source_Sample_Code", "City", "Site", "Collection_Date", "year", "month", "Season", "Sample_Type", "reads_total", "Base_Count", "bracken_phylum_reads"]].to_csv(DER / "ds02_meta.csv", index=False)
ph = tax.pivot_table(index="Sample_ID", columns="Taxon_Name", values="Abundance_Value", aggfunc="sum", fill_value=0)
ph.to_csv(DER / "ds02_phylum_counts.csv")
print("DS02:", d2.shape, "phylum matrix", ph.shape, d2.City.value_counts().to_dict(), d2.RT.value_counts().to_dict())

# %%
# ResFinder set: map runs to analytic units -----------------------------------
runs = pd.DataFrame({"Run_ID": sorted(set(sub_arg.run_accession) | set(sub_div.run_accession))})
runs = runs.merge(rm[["Run_ID", "Sample_ID"]], on="Run_ID", how="left")
runs = runs.merge(sm[["Sample_ID", "Dataset_ID", "Project_Accession", "Biological_Sample_ID", "Inclusion_Status", "Sample_Type"]], on="Sample_ID", how="left")
ptr = rm.set_index("Run_ID")
runs["unit_id"] = np.where((runs.Project_Accession == "PRJEB27054") & ~is_missing(runs.Biological_Sample_ID),
                           "BIO:" + runs.Biological_Sample_ID.astype(str), runs.Sample_ID)
runs["has_ARG_rows"] = runs.Run_ID.isin(sub_arg.run_accession)
runs["primary"] = runs.Project_Accession != "PRJEB13831"
runs["exclusion_reason"] = np.where(~runs.primary, "PRJEB13831: probable duplicate physical samples of PRJEB27054 pilot (sensitivity only)",
                                    np.where(~runs.has_ARG_rows, "No ARG rows in resource (missing, not zero)", ""))
kv = sub_arg.drop_duplicates("run_accession").set_index("run_accession").kma_version
runs["kma_version"] = runs.Run_ID.map(kv).fillna(runs.Run_ID.map(sub_div.drop_duplicates("run_accession").set_index("run_accession").kma_version))
runs.to_csv(DER / "rf_run_to_unit.csv", index=False)
record_exclusion(NB, "ResFinder resistome analyses", "run", (~runs.primary).sum(), "PRJEB13831 probable duplicate physical samples of PRJEB27054 (held out for sensitivity analysis)")
record_exclusion(NB, "ResFinder resistome analyses", "run", (~runs.has_ARG_rows).sum(), "Run listed in resource without ARG rows (treated as missing, not zero)", ";".join(runs.Run_ID[~runs.has_ARG_rows]))
nunits = runs.groupby("primary").unit_id.nunique()
record_exclusion(NB, "ResFinder resistome analyses", "accession", int(runs[runs.Project_Accession == "PRJEB27054"].Sample_ID.nunique() - runs[runs.Project_Accession == "PRJEB27054"].unit_id.nunique()),
                 "PRJEB27054 accessions merged into their biological sample (counts summed)")
print(runs.groupby(["primary", "Inclusion_Status", "Sample_Type"]).agg(runs=("Run_ID", "size"), samples=("Sample_ID", "nunique"), units=("unit_id", "nunique")))
print("units by primary flag:", nunits.to_dict())

# %% [markdown]
# ## 13. Data-quality summary report

# %%
lines = ["# Data Quality and Exclusion Report", f"Generated by notebook 01 on {dt.date.today().isoformat()} from the harmonized tables (read only).", "",
         f"- Files inspected in the research directory: **{len(inv)}** ({n_dupsets} sets of byte-identical copies, i.e. the `_1` duplicates).",
         f"- Sample master: **{len(sm):,}** samples; included **{(sm.Inclusion_Status=='Included').sum():,}**, contextual **{(sm.Inclusion_Status=='Separate contextual dataset').sum():,}**, excluded **{(sm.Inclusion_Status=='Excluded').sum():,}**.",
         f"- Studies with pending eligibility: **{len(pend)}**; none were harmonized, so none enter any analysis.",
         f"- Identifier linkage: {int((link.status=='Pass').sum())} of {len(link)} checks passed.", "", "## Linkage checks", df_to_md(link), "",
         "## Duplicate checks", df_to_md(dup), "", "## Geographic checks", df_to_md(geo), "",
         "## Annotation completeness", df_to_md(ann), "", "## Reported vs recomputed counts", df_to_md(rc)]
(QCD / "Data_Quality_Report.md").write_text("\n".join(lines), encoding="utf-8")
print("written qc/Data_Quality_Report.md")
