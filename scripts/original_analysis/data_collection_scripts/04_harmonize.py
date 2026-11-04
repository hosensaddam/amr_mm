"""
04_harmonize.py
Build harmonized, analysis-ready CSV tables from preserved raw downloads.
No values are imputed. Every harmonized field keeps its source value alongside.
Missing-value codes: "Not reported", "Not available", "Not yet verified",
"Not applicable", "Restricted access".

Outputs (harmonized/):
  run_master.csv, sample_master.csv, study_inventory.csv, accession_register.csv,
  taxonomy_long.csv, diversity_long.csv, amr_genes_long.csv, env_covariates_long.csv,
  transformation_log.csv
"""
import csv, glob, json, os, re, importlib.util
import numpy as np, pandas as pd, pycountry

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H = os.path.join(B, "harmonized"); R = os.path.join(B, "raw_downloads")
spec = importlib.util.spec_from_file_location("cd", os.path.join(B, "scripts", "curation_decisions.py"))
cd = importlib.util.module_from_spec(spec); spec.loader.exec_module(cd)
NR, NAV, NYV, NAP, RA = "Not reported", "Not available", "Not yet verified", "Not applicable", "Restricted access"
RETRIEVED = "2026-10-31"
TLOG = []
def tlog(step, desc, n=None): TLOG.append({"Step": step, "Description": desc, "Records_Affected": n})

SEWAGE = ["PRJEB13831", "PRJEB27054", "PRJEB27621", "PRJEB40798", "PRJEB40815", "PRJEB40816", "PRJEB51229"]
CORE_DS = {"DS01": SEWAGE, "DS02": ["PRJNA1415085"], "DS03": ["PRJNA1399457"], "DS04": ["PRJNA732392"], "DS05": ["PRJNA1522440"]}
DS_NAME = {"DS01": "Global Sewage Surveillance Project (Munk et al. 2022; incl. 2016 pilot)",
           "DS02": "Indian metropolitan open-drain wastewater (Singh et al. 2026)",
           "DS03": "Dhaka hospital and community wastewater (BCSIR)",
           "DS04": "MetaSUB City Sampling Day 2016-2017 (SRA subset)",
           "DS05": "FDA retail lamb food metagenome (contextual)"}
DS_PUB = {"DS01": "10.1038/s41467-022-34312-7", "DS02": "10.1038/s41467-026-70702-x", "DS04": "10.1016/j.cell.2021.05.002"}
DS_PUBID = {"DS01": "PMID:36456547; PMC9715550", "DS02": "PMID:41844673; PMC13144492", "DS04": "PMC8238498"}

# ----------------------------------------------------------------- helpers
def clean(v):
    if v is None or (isinstance(v, float) and np.isnan(v)): return ""
    return str(v).strip()

MISSING_TOKENS = {"", "missing", "not collected", "not applicable", "na", "n/a", "nan", "none", "unknown",
                  "not provided", "restricted access", "-", "missing: data agreement established pre-2023"}
def is_missing(v): return clean(v).lower() in MISSING_TOKENS

COUNTRY_FIX = {"USA": "US", "Russia": "RU", "South Korea": "KR", "Viet Nam": "VN", "Vietnam": "VN", "Iran": "IR",
               "Czech Republic": "CZ", "Hong Kong": "HK", "Taiwan": "TW", "Bolivia": "BO", "Tanzania": "TZ",
               "Turkey": "TR", "Cote d'Ivoire": "CI", "Côte d'Ivoire": "CI", "Ivory Coast": "CI", "Laos": "LA",
               "Moldova": "MD", "Venezuela": "VE", "Syria": "SY", "Palestine": "PS", "Kosovo": "XK",
               "Democratic Republic of the Congo": "CD", "Republic of the Congo": "CG", "Brunei": "BN",
               "Macedonia": "MK", "North Macedonia": "MK", "United Kingdom": "GB", "UK": "GB"}
def iso2(country):
    c = clean(country)
    if not c or c.startswith("Not "): return NAP if c.startswith("Not ") else NR
    if c in COUNTRY_FIX: return COUNTRY_FIX[c]
    try: return pycountry.countries.lookup(c).alpha_2
    except LookupError: return NYV

def std_country(c):
    c = clean(c).replace("_", " ")
    if is_missing(c): return NR
    special = {"USA": "United States", "UK": "United Kingdom", "Viet Nam": "Viet Nam"}
    return special.get(c, c[:1].upper() + c[1:])

def parse_date(v):
    """Return (iso, year, month, note). Ranges keep the start date and are flagged."""
    s = clean(v)
    if is_missing(s): return NR, NR, NR, ""
    note = ""
    if "/" in s:
        s, note = s.split("/")[0], "Date range reported; start date retained"
    s = s.split("T")[0]
    m = re.match(r"^(\d{4})(?:-(\d{1,2}))?(?:-(\d{1,2}))?$", s)
    if not m:
        m3 = re.match(r"^(\d{4})-([A-Za-z]{3})-(\d{1,2})$", s)
        if m3:
            mon = ["jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"].index(m3.group(2).lower()) + 1
            return f"{m3.group(1)}-{mon:02d}-{int(m3.group(3)):02d}", int(m3.group(1)), mon, "Reformatted from YYYY-Mon-DD"
        m2 = re.match(r"^(\d{1,2})-([A-Za-z]{3})-(\d{4})$", s)
        if m2:
            mon = ["jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"].index(m2.group(2).lower()) + 1
            return f"{m2.group(3)}-{mon:02d}-{int(m2.group(1)):02d}", int(m2.group(3)), mon, "Reformatted from DD-Mon-YYYY"
        return NYV, NYV, NYV, f"Unparsed date string: {clean(v)}"
    y = int(m.group(1)); mo = int(m.group(2)) if m.group(2) else None; d = int(m.group(3)) if m.group(3) else None
    if not 1950 <= y <= 2026: return NYV, NYV, NYV, f"Implausible year in source date '{clean(v)}'"
    if s == "1970-01-01": return NYV, NYV, NYV, "Epoch placeholder date 1970-01-01 in source"
    if mo and not 1 <= mo <= 12: return NYV, y, NYV, f"Invalid month in {clean(v)}"
    iso = f"{y}" + (f"-{mo:02d}" if mo else "") + (f"-{d:02d}" if d else "")
    return iso, y, (mo if mo else NR), note

def decimals(x):
    s = clean(x)
    return len(s.split(".")[1]) if "." in s else 0

def parse_coord(lat, lon):
    la, lo = clean(lat), clean(lon)
    if is_missing(la) or is_missing(lo): return None, None, ""
    try: fa, fo = float(la), float(lo)
    except ValueError: return None, None, f"Non-numeric coordinates ({la}, {lo})"
    if not (-90 <= fa <= 90 and -180 <= fo <= 180): return None, None, f"Coordinates out of range ({la}, {lo})"
    return fa, fo, ""

def parse_latlon_string(s):
    m = re.match(r"^\s*([\d.]+)\s*([NS])\s+([\d.]+)\s*([EW])\s*$", clean(s))
    if not m: return "", ""
    la = float(m.group(1)) * (1 if m.group(2) == "N" else -1); lo = float(m.group(3)) * (1 if m.group(4) == "E" else -1)
    return m.group(1) if m.group(2) == "N" else "-" + m.group(1), m.group(3) if m.group(4) == "E" else "-" + m.group(3)

def data_type(strategy, source):
    s, so = clean(strategy).upper(), clean(source).upper()
    if s == "WGS" and so == "METAGENOMIC": return "Shotgun metagenomic"
    if s == "AMPLICON": return "Amplicon"
    if so == "METATRANSCRIPTOMIC" or s == "RNA-SEQ": return "Metatranscriptomic"
    if s == "WGA": return "Single-cell/WGA genomic"
    if s == "WGS" and so == "GENOMIC": return "Isolate/genomic WGS"
    if so == "SYNTHETIC": return "Synthetic/mock"
    return f"Other ({s}/{so})"

def repo_of(run):
    return {"SRR": "NCBI SRA (INSDC)", "ERR": "ENA (INSDC)", "DRR": "DDBJ DRA (INSDC)"}.get(run[:3], "INSDC")

# ----------------------------------------------------------------- load ENA runs
sel = json.load(open(os.path.join(H, "selected_studies.json")))
frames = []
for acc in sel:
    p = os.path.join(R, "ena_sample_metadata", f"{acc}.tsv")
    d = pd.read_csv(p, sep="\t", dtype=str).fillna("")
    d["Curated_Study"] = acc
    frames.append(d)
ena = pd.concat(frames, ignore_index=True)
dups = ena.duplicated("run_accession").sum()
ena = ena.drop_duplicates("run_accession")
tlog("Load", f"Loaded ENA run reports for {len(sel)} selected studies; removed {dups} duplicate run rows across study files", dups)

# dataset ids
ds_of = {}
for ds, accs in CORE_DS.items():
    for a in accs: ds_of[a] = ds
others = sorted(a for a in sel if a not in ds_of)
for i, a in enumerate(others, start=6): ds_of[a] = f"DS{i:02d}"
MGR_INCL = ["mgp16346", "mgp89046", "mgp90475", "mgp90740", "mgp92067", "mgp92140", "mgp93181"]
nxt = 6 + len(others)
ds_of["mgp16346"] = f"DS{nxt:02d}"; ds_of["mgp89046"] = f"DS{nxt+1:02d}"
for p in ["mgp90475", "mgp90740", "mgp92067", "mgp92140", "mgp93181"]: ds_of[p] = f"DS{nxt+2:02d}"
DS_NAME[ds_of["mgp16346"]] = "GLUSEEN global urban soils (MG-RAST)"
DS_NAME[ds_of["mgp89046"]] = "Swiss wastewater and river bacteriomes, NRP72 (MG-RAST)"
DS_NAME[ds_of["mgp90475"]] = "Mobile antibiotic resistome in Polish WWTPs, Tychy/Olsztyn (MG-RAST series)"
cand = json.load(open(os.path.join(H, "study_candidates.json")))
title_of = {v["accession"]: (v.get("ena") or {}).get("study_title") or (v["titles"][0] if v["titles"] else "") for v in cand.values()}
mgnify_of = {v["accession"]: ";".join(v.get("mgnify") or []) for v in cand.values()}
for a in others: DS_NAME[ds_of[a]] = title_of.get(a, a)

# ----------------------------------------------------------------- auxiliary sources
munk = pd.read_excel(os.path.join(R, "munk2022_supplementary", "extracted", "41467_2022_34312_MOESM3_ESM.xlsx"),
                     sheet_name="Supplementary Data 1", dtype=str).fillna("")
munk_by_run = munk.set_index("ena_run_acc")
metasub = pd.read_csv(os.path.join(R, "metasub_metadata", "complete_metadata.csv"), dtype=str).fillna("")
ms_by_uuid = metasub.drop_duplicates("uuid").set_index("uuid")
india_meta = pd.read_csv(os.path.join(R, "zenodo_17590646_india", "extracted", "data", "meta.tsv"), sep="\t", dtype=str).fillna("")
india_by_code = india_meta.set_index("sample_code")
bios = pd.read_csv(os.path.join(H, "biosample_attributes.csv"), dtype=str).fillna("").drop_duplicates("sample_accession").set_index("sample_accession")
g214 = pd.read_csv(os.path.join(R, "zenodo_6919377_global_resistome", "metadata.tsv"), sep="\t", dtype=str, usecols=["run_accession"])
in214 = set(g214.run_accession)

# MetaSUB: country by city from ENA records of the same dataset (used only when unambiguous)
ms_ena = ena[ena.Curated_Study == "PRJNA732392"].copy()
ms_ena["ms_city"] = ms_ena.sample_alias.map(lambda u: ms_by_uuid.city.get(u, ""))
ms_city_country = {}
for city, g in ms_ena.groupby("ms_city"):
    cs = {clean(x).split(":")[0] for x in g.country if not is_missing(x)}
    ms_city_country[city] = cs

# ----------------------------------------------------------------- run-level harmonization
rows = []
for _, r in ena.iterrows():
    study = r.Curated_Study; ds = ds_of[study]; run = r.run_accession; smp = r.sample_accession
    dt = data_type(r.library_strategy, r.library_source)
    out = {"Dataset_ID": ds, "Study_ID": study, "Project_Accession": r.study_accession,
           "Secondary_Study_Accession": r.secondary_study_accession or NR, "Sample_ID": smp,
           "Secondary_Sample_ID": r.secondary_sample_accession or NR, "Run_ID": run,
           "Experiment_ID": r.experiment_accession, "Repository": repo_of(run),
           "Sample_Alias": r.sample_alias or NR, "Library_Name": r.library_name or NR,
           "Sequencing_Method": f"{r.library_strategy}/{r.library_source}/{r.library_selection}",
           "Sequencing_Platform": f"{r.instrument_platform}: {r.instrument_model}" if r.instrument_model else NR,
           "Sequencing_Type": dt, "Read_Type": r.library_layout or NR,
           "Base_Count": int(r.base_count) if r.base_count.isdigit() else NR,
           "Read_Count": int(r.read_count) if r.read_count.isdigit() else NR,
           "FASTQ_Bytes": sum(int(b) for b in r.fastq_bytes.split(";") if b.isdigit()) if r.fastq_bytes else NR,
           "FASTQ_MD5": r.fastq_md5 or NR,
           "Raw_Read_Availability": "Public FASTQ (ENA)" if r.fastq_ftp else NAV,
           "Raw_Data_URL": ";".join("ftp://" + u for u in r.fastq_ftp.split(";")) if r.fastq_ftp else NAV,
           "First_Public": r.first_public or NR, "Center_Name": r.center_name or NR,
           "Source_URL": f"https://www.ebi.ac.uk/ena/browser/view/{run}",
           # raw geography/date/matrix (source values)
           "SRC_country": r.country, "SRC_lat": r.lat, "SRC_lon": r.lon, "SRC_collection_date": r.collection_date,
           "SRC_isolation_source": r.isolation_source, "SRC_env_medium": r.environmental_medium or r.environment_material,
           "SRC_env_local": r.local_environmental_context or r.environment_feature,
           "SRC_env_broad": r.broad_scale_environmental_context or r.environment_biome,
           "SRC_sample_title": r.sample_title, "SRC_scientific_name": r.scientific_name}
    rows.append(out)
runs = pd.DataFrame(rows)

# MG-RAST metagenomes as run-equivalents
mg_rows = []
for p in MGR_INCL:
    d = json.load(open(os.path.join(R, "mgrast_metagenomes", f"{p}.json")))
    for m in d["data"]:
        st = clean(m.get("sequence_type"))
        dt = {"WGS": "Shotgun metagenomic", "MT": "Metatranscriptomic", "Amplicon": "Amplicon"}.get(st, f"Other ({st})")
        mg = m["metagenome_id"]
        mg_rows.append({"Dataset_ID": ds_of[p], "Study_ID": p, "Project_Accession": p, "Secondary_Study_Accession": NR,
                        "Sample_ID": m.get("sample_id") or mg, "Secondary_Sample_ID": mg, "Run_ID": mg, "Experiment_ID": NAP,
                        "Repository": "MG-RAST", "Sample_Alias": clean(m.get("library_name")) or NR, "Library_Name": clean(m.get("library_name")) or NR,
                        "Sequencing_Method": f"{st}/{clean(m.get('seq_meth'))}", "Sequencing_Platform": clean(m.get("seq_meth")) or NR,
                        "Sequencing_Type": dt, "Read_Type": NR, "Base_Count": m.get("bp_count_raw") or NR,
                        "Read_Count": m.get("sequence_count_raw") or NR, "FASTQ_Bytes": NR, "FASTQ_MD5": NR,
                        "Raw_Read_Availability": "Public (MG-RAST download)" if m.get("public") else RA,
                        "Raw_Data_URL": f"https://api.mg-rast.org/download/{mg}", "First_Public": NR, "Center_Name": clean(m.get("PI_organization")) or NR,
                        "Source_URL": f"https://www.mg-rast.org/mgmain.html?mgpage=overview&metagenome={mg}",
                        "SRC_country": clean(m.get("country")), "SRC_lat": clean(m.get("latitude")), "SRC_lon": clean(m.get("longitude")),
                        "SRC_collection_date": clean(m.get("collection_date")), "SRC_isolation_source": clean(m.get("material")),
                        "SRC_env_medium": clean(m.get("material")), "SRC_env_local": clean(m.get("feature")), "SRC_env_broad": clean(m.get("biome")),
                        "SRC_sample_title": clean(m.get("location")), "SRC_scientific_name": clean(m.get("env_package_type"))})
runs = pd.concat([runs, pd.DataFrame(mg_rows)], ignore_index=True)
tlog("Load", f"Added {len(mg_rows)} MG-RAST metagenome records (7 projects) as run-equivalent rows", len(mg_rows))

# ----------------------------------------------------------------- sample-level geography, date, matrix
def wastewater_like(cat): return cat in (cd.MW, cd.HW, cd.CW)

def classify(study_cat, text, ds, extra):
    """Return (Sample_Type, rule) using documented keyword rules within the study category."""
    t = text.lower()
    if ds == "DS04":
        if extra.get("control"): return cd.OTH, "MetaSUB control sample"
        if extra.get("ms_sample_type") == "air_microbiome": return cd.AIR, "MetaSUB sample_type=air_microbiome"
        return cd.TS, "MetaSUB transit/urban surface (study design)"
    if ds_of.get("mgp89046") == ds:
        if re.search(r"river|sediment|lake", t): return cd.SW, "MG-RAST material/feature indicates river/sediment"
        return cd.MW, "MG-RAST wastewater"
    if wastewater_like(study_cat):
        if re.search(r"hospital", t): return cd.HW, "Sample text mentions hospital"
        if re.search(r"\briver|\blake|stream|estuar", t) and not re.search(r"sewage|wastewater|influent|effluent", t):
            return cd.SW, "Sample text indicates receiving surface water"
        if re.search(r"tap water|drinking", t): return cd.OTH, "Sample text indicates tap/drinking water"
        if re.search(r"community|sewer|drain", t) and study_cat == cd.HW: return cd.CW, "Sample text indicates community sewer/drain"
        return study_cat, "Study-level category"
    if study_cat == cd.TS and re.search(r"\bair\b|aerosol", t): return cd.AIR, "Sample text indicates air"
    return study_cat or "Unclassified", "Study-level category"

samples = []
for (repo_sample), g in runs.groupby("Sample_ID", sort=False):
    r0 = g.iloc[0]; study = r0.Study_ID; ds = r0.Dataset_ID
    dec = cd.decide(study) if not study.startswith("mgp") else cd.MGRAST[study]
    study_decision, study_cat = dec[0], dec[1]
    dts = sorted(set(g.Sequencing_Type))
    shotgun = "Shotgun metagenomic" in dts
    flags, remarks = [], []
    other = g[g.Sequencing_Type != "Shotgun metagenomic"]
    other_runs = "; ".join(f"{a} ({b})" for a, b in zip(other.Run_ID, other.Sequencing_Type)) or NAP
    if shotgun and len(other):
        g = g[g.Sequencing_Type == "Shotgun metagenomic"]
        flags.append(f"Sample also has {len(other)} non-shotgun run(s), excluded from Run_ID and totals (see Other_Run_IDs)")
    src = {k: next((clean(x) for x in g[k] if not is_missing(x)), "") for k in
           ["SRC_country", "SRC_lat", "SRC_lon", "SRC_collection_date", "SRC_isolation_source", "SRC_env_medium",
            "SRC_env_local", "SRC_env_broad", "SRC_sample_title", "SRC_scientific_name"]}
    bs = bios.loc[repo_sample] if repo_sample in bios.index else None
    geo_src, date_src = "", ""
    city, country, locality, lat, lon, date_raw = NR, NR, NR, "", "", ""
    country_src_val, city_src_val = src["SRC_country"], ""
    geounit, bio_id, extra = "", NAP, {}
    run_ids = list(g.Run_ID)
    # --- source priority
    if ds == "DS01":
        m = munk_by_run.loc[[x for x in run_ids if x in munk_by_run.index]] if any(x in munk_by_run.index for x in run_ids) else None
        if m is not None and len(m):
            m0 = m.iloc[0]
            city, country, lat, lon = std_country(m0.city) if m0.city else NR, std_country(m0.country), m0.latitude, m0.longtitude
            city_src_val = m0.city; country_src_val = m0.country
            date_raw = clean(m0.collection_date)[:10]; bio_id = m0.unique_sample_id
            geo_src = date_src = "Munk et al. 2022, Supplementary Data 1 (author-reported)"; geounit = "City (source-labelled)"
            extra["batch"] = m0.batch_name; extra["year"] = m0.year
        else:
            flags.append("Run not listed in Munk et al. 2022 Supplementary Data 1; repository geography used")
            extra["generic"] = True
    if ds == "DS04":
        u = clean(g.Sample_Alias.iloc[0])
        if u in ms_by_uuid.index:
            ms = ms_by_uuid.loc[u]
            city_src_val = ms.city
            ctrl = bool(clean(ms.control_type)) or ms.city in ("control", "neg_control", "pos_control", "other_control")
            extra.update({"control": ctrl, "ms_sample_type": ms.sample_type, "ms": ms})
            if not ctrl and ms.city:
                city = ms.city.replace("_", " ").title(); geounit = "City (source-labelled)"
            lat, lon = ms.latitude, ms.longitude
            geo_src = "MetaSUB-metadata complete_metadata.csv (GitHub commit 1d73874)"
            if not is_missing(src["SRC_country"]):
                country = std_country(src["SRC_country"].split(":")[0]); geo_src += "; country from SRA geo_loc_name"
            else:
                cs = ms_city_country.get(ms.city, set())
                if len(cs) == 1 and not ctrl:
                    country = std_country(next(iter(cs)))
                    flags.append("Country propagated from same-city records in this dataset (SRA field missing)")
                    geo_src += "; country propagated within dataset"
            if city.lower() == "kyiv" and country == "Poland":
                if not is_missing(src["SRC_country"]):
                    flags.append("City-country inconsistency in source (SRA lists Kyiv under Poland); country set to Not yet verified")
                else:
                    flags.append("City-country inconsistency (country propagated from same-city SRA records lists Poland for Kyiv); country set to Not yet verified")
                country = NYV
        else:
            flags.append("SRA sample alias not found in MetaSUB metadata")
        date_raw = src["SRC_collection_date"]; date_src = "SRA collection_date"
    elif ds == "DS02":
        code = clean(g.Sample_Alias.iloc[0])
        if code in india_by_code.index:
            city = india_by_code.loc[code].city_state; city_src_val = city; geounit = "City (source-labelled)"
            extra["india_code"] = code; extra["india_month"] = india_by_code.loc[code].month; extra["india_year"] = india_by_code.loc[code].year
        country = "India"; geo_src = "Zenodo 17590646 meta.tsv (city) + SRA geo_loc_name (country)"
        date_raw = src["SRC_collection_date"]; date_src = "SRA collection_date"
        ena_city = clean(src["SRC_country"]).split(":")[-1].strip()
        if ena_city and city != NR and ena_city.lower() != city.lower():
            flags.append(f"City mismatch: SRA '{ena_city}' vs Zenodo '{city}'")
    if ds not in ("DS01", "DS02", "DS04") or extra.get("generic"):
        loc = src["SRC_country"]
        if (is_missing(loc) or ":" not in loc) and bs is not None and not is_missing(bs.geo_loc_name):
            loc = bs.geo_loc_name; geo_src = "NCBI BioSample geo_loc_name"
        elif not is_missing(loc):
            geo_src = "Repository geo_loc_name (ENA/MG-RAST)"
        if not is_missing(loc):
            country = std_country(loc.split(":")[0]); country_src_val = loc
            if ":" in loc:
                locality = clean(loc.split(":", 1)[1]); geounit = "Locality as reported (not standardized to city)"
            else:
                geounit = "Country only"
        if study.startswith("mgp"):
            locality = src["SRC_sample_title"] or NR
            geounit = "Locality as reported (not standardized to city)" if locality != NR else "Country only"
        if ds == "DS03" and bs is not None and "Dhaka" in bs.geo_loc_name:
            city, geounit = "Dhaka", "City (source-labelled)"; locality = clean(bs.geo_loc_name.split(":", 1)[1])
        lat, lon = src["SRC_lat"], src["SRC_lon"]
        if (is_missing(lat) or is_missing(lon)) and bs is not None and bs.lat_lon:
            lat, lon = parse_latlon_string(bs.lat_lon)
            if lat: geo_src += "; coordinates from BioSample lat_lon"
        date_raw = src["SRC_collection_date"]; date_src = "Repository collection_date"
        if is_missing(date_raw) and bs is not None and not is_missing(bs.collection_date):
            date_raw = bs.collection_date; date_src = "NCBI BioSample collection_date"
    iso_date, yr, mo, dnote = parse_date(date_raw)
    if dnote: flags.append(dnote)
    fa, fo, cnote = parse_coord(lat, lon)
    if cnote: flags.append(cnote)
    if fa is not None:
        dec = min(decimals(lat), decimals(lon))
        precision = "Exact (source-reported, >=3 decimals)" if dec >= 3 else "Approximate (source-reported, <=2 decimals)"
        if fa == 0 and fo == 0: precision = "Invalid (0,0)"; flags.append("Coordinates 0,0"); fa = fo = None
    else:
        precision = "City-level location only" if (city != NR or locality not in (NR, "")) else "Missing geographic information"
    # matrix text and classification
    text = " ".join([src["SRC_isolation_source"], src["SRC_env_medium"], src["SRC_env_local"], src["SRC_sample_title"],
                     bs.isolation_source if bs is not None else "", bs.env_medium if bs is not None else ""])
    stype, rule = classify(study_cat, text, ds, extra)
    # inclusion at sample level
    if study_decision == cd.S:
        incl, why = cd.S, "Study classified as separate contextual dataset"
    elif study_decision != cd.I:
        incl, why = study_decision, "Study-level decision"
    elif not shotgun:
        incl, why = cd.E, f"No shotgun metagenomic run (data type: {', '.join(dts)})"
    elif extra.get("control"):
        incl, why = cd.E, "Negative/positive control sample"
    elif stype == cd.OTH and "tap" in rule:
        incl, why = cd.E, "Tap/drinking water outside urban environmental scope"
    else:
        incl, why = cd.I, "Meets eligibility: public shotgun metagenome from included urban study"
    # processed data availability
    ind_runs = [x for x in run_ids if x in in214]
    proc, purl, tax, amr, rank, func = [], [], NAV, "No processed AMR data collected (raw reads permit reanalysis)", NAP, NAV
    if ds == "DS02":
        proc.append("Bracken taxonomy + RGI/CARD ARG table (Zenodo 17590646, downloaded)"); purl.append("https://doi.org/10.5281/zenodo.17590646")
        tax = "Yes (Kraken2/Bracken read estimates; downloaded)"; rank = "Phylum in Taxonomy_Long; genus and species in source files"
        amr = "Yes (RGI read mapping against CARD; downloaded)"
    if ind_runs:
        proc.append("KMA/ResFinder ARG counts + genus/phylum/ARG diversity (Zenodo 6919377 resource; subset extracted)"); purl.append("https://doi.org/10.5281/zenodo.6919377")
        tax = "Partial (genus/phylum diversity indices downloaded; rRNA count table not downloaded, 6.4 GB)"; rank = "Genus and phylum (diversity indices only)"
        amr = "Yes (KMA read mapping against ResFinder; downloaded)"
    if study.startswith("mgp"):
        proc.append("MG-RAST taxonomic/functional annotations (online; not downloaded)"); purl.append(f"https://www.mg-rast.org/mgmain.html?mgpage=project&project={study}")
        tax = "Available in MG-RAST (not downloaded)"; func = "Available in MG-RAST (not downloaded)"
    if mgnify_of.get(study):
        proc.append(f"MGnify study {mgnify_of[study]} (analysis availability not checked per sample)")
        purl.append(f"https://www.ebi.ac.uk/metagenomics/studies/{mgnify_of[study].split(';')[0]}")
    if study == "PRJEB27054":
        proc.append("MGnify MGYS00005846 (assembly analyses; not downloaded)"); purl.append("https://www.ebi.ac.uk/metagenomics/studies/MGYS00005846")
    # QC
    miss = [n for n, v in [("City/locality", city if city != NR else locality), ("Country", country), ("Collection_Date", iso_date),
                            ("Coordinates", "x" if fa is not None else NR)] if v in (NR, "", NYV)]
    if incl == cd.I and stype in ("Unclassified", ""): flags.append("Sample type unclassified")
    qc = "Fail" if any(("out of range" in f or "Invalid" in f or "inconsistency" in f) for f in flags) else ("Warning" if flags or miss else "Pass")
    if ds in ("DS01", "DS02", "DS03", "DS04"):
        ver = "Verified" if (country not in (NR, NYV) and iso_date not in (NR, NYV)) else "Partially verified"
    else:
        ver = "Partially verified"
    samples.append({
        "Dataset_ID": ds, "Dataset_Name": DS_NAME.get(ds, ""), "Study_ID": study,
        "Project_Accession": r0.Project_Accession, "Sample_ID": repo_sample, "Secondary_Sample_ID": r0.Secondary_Sample_ID,
        "Biological_Sample_ID": bio_id, "Run_ID": ";".join(run_ids), "N_Runs": len(run_ids), "Other_Run_IDs": other_runs,
        "Analysis_ID": ";".join(x for x in [mgnify_of.get(study, ""), "MGYS00005846" if study == "PRJEB27054" else ""] if x) or NAV,
        "Repository": r0.Repository, "Publication_ID": DS_PUBID.get(ds, NYV), "DOI": DS_PUB.get(ds, NYV),
        "City": city, "Locality_As_Reported": locality if locality else NR, "Country": country, "ISO_Country_Code": iso2(country) if country not in (NR, NYV) else country,
        "Latitude": round(fa, 6) if fa is not None else NR, "Longitude": round(fo, 6) if fo is not None else NR,
        "Geographic_Unit": geounit or "Missing geographic information", "Coordinate_Precision": precision,
        "Urban_Setting": "Urban (study design, publication-confirmed)" if ds in ("DS01", "DS02", "DS03", "DS04") else ("Not applicable (food chain)" if ds == "DS05" else "Urban (study title; to confirm)"),
        "Geographic_Source": geo_src or NR,
        "Sample_Type": stype, "Sample_Type_Rule": rule,
        "Sample_Subtype": src["SRC_isolation_source"] or (bs.isolation_source if bs is not None and bs.isolation_source else NR),
        "Environmental_Matrix": src["SRC_env_medium"] or (bs.env_medium if bs is not None and bs.env_medium else NR),
        "Collection_Date": iso_date, "Collection_Year": yr, "Collection_Month": mo, "Collection_Date_Source": date_src or NR,
        "Sample_Size": NAP, "Sampling_Design": extra.get("batch", NR) if ds == "DS01" else ("Monthly grab samples, 19 open-drain sites (publication)" if ds == "DS02" else NR),
        "Sampling_Location": locality if ds == "DS03" else ((clean(extra["ms"].station) or clean(extra["ms"].surface) or NR) if "ms" in extra else NR),
        "Sequencing_Method": ";".join(sorted(set(g.Sequencing_Method))), "Sequencing_Platform": ";".join(sorted(set(g.Sequencing_Platform))),
        "Sequencing_Type": ";".join(sorted(set(g.Sequencing_Type))), "Read_Type": ";".join(sorted(set(map(str, g.Read_Type)))),
        "Base_Count_Total": int(sum(int(x) for x in g.Base_Count if str(x).isdigit())) if any(str(x).isdigit() for x in g.Base_Count) else NR,
        "FASTQ_Bytes_Total": int(sum(int(x) for x in g.FASTQ_Bytes if str(x).isdigit())) if any(str(x).isdigit() for x in g.FASTQ_Bytes) else NR,
        "Raw_Read_Availability": ";".join(sorted(set(g.Raw_Read_Availability))), "Raw_Data_URL": ";".join(g.Raw_Data_URL)[:2000],
        "Processed_Data_Availability": "; ".join(proc) or NAV, "Processed_Data_URL": "; ".join(purl) or NAV,
        "Taxonomic_Profile_Available": tax, "Taxonomic_Rank": rank,
        "Microbial_Taxa": "See 03_Analysis_Ready Taxonomy_Long" if ds == "DS02" else ("See Diversity_Long" if ind_runs else NAV),
        "Taxonomic_Abundance": "See 03_Analysis_Ready Taxonomy_Long" if ds == "DS02" else NAV,
        "Abundance_Unit": "Bracken estimated reads; fraction of classified reads" if ds == "DS02" else NAP,
        "Functional_Profile_Available": func,
        "AMR_Data_Available": amr,
        "AMR_Gene": "See AMR_Genes_Long" if (ds == "DS02" or ind_runs) else NAV,
        "AMR_Gene_Family": "See AMR_Genes_Long" if ds == "DS02" else (NAV if not ind_runs else "Not provided by source (ResFinder gene IDs only)"),
        "Resistance_Class": "See AMR_Genes_Long" if (ds == "DS02" or ind_runs) else NAV,
        "Resistance_Mechanism": "See AMR_Genes_Long" if ds == "DS02" else (NAV if not ind_runs else "Not provided by source"),
        "AMR_Detection_Type": "Read mapping, RGI bwt (Bowtie2)" if ds == "DS02" else ("Read mapping, KMA 1.3.0" if ind_runs else NAP),
        "AMR_Abundance": "See AMR_Genes_Long" if (ds == "DS02" or ind_runs) else NAV,
        "AMR_Abundance_Unit": "All mapped reads (count)" if ds == "DS02" else ("Aligned read fragments (count)" if ind_runs else NAP),
        "AMR_Reference_Database": "CARD" if ds == "DS02" else ("ResFinder" if ind_runs else NAP),
        "AMR_Database_Version": "Not reported in deposited files" if ds == "DS02" else ("Not stated in table; see Martiny et al. resource documentation" if ind_runs else NAP),
        "Data_Access_Status": "Public" if "Restricted" not in ";".join(g.Raw_Read_Availability) else RA,
        "Licence": ("Raw reads: INSDC unrestricted-access policy; processed data CC BY 4.0 (Zenodo)" if (ds == "DS02" or ind_runs) else
                    ("MG-RAST public data; no licence stated" if study.startswith("mgp") else "INSDC unrestricted-access policy; no study-specific licence stated")),
        "Date_Collected": RETRIEVED, "Verification_Status": ver, "Inclusion_Status": incl, "Inclusion_Reason": why,
        "Missing_Data_Status": "Complete (core fields)" if not miss else "Missing: " + "; ".join(miss),
        "QC_Status": qc, "QC_Flags": "; ".join(flags) or "None",
        "Limitations": "", "Source_URL": (f"https://www.ebi.ac.uk/ena/browser/view/{repo_sample}" if not study.startswith("mgp") else r0.Source_URL),
        "Remarks": "",
        # source values kept alongside harmonized values
        "City_Source_Value": city_src_val or NR, "Country_Source_Value": country_src_val or NR,
        "Collection_Date_Source_Value": date_raw or NR, "Sample_Type_Source_Value": (src["SRC_isolation_source"] or src["SRC_env_medium"] or NR),
        "Latitude_Source_Value": lat or NR, "Longitude_Source_Value": lon or NR,
        "_extra": extra})
sm = pd.DataFrame(samples)
LIMIT = {"DS01": "Sewage resistome reflects pooled community inputs; read-mapping counts not comparable with other pipelines; 2016 pilot round sits in PRJEB27054.",
         "DS02": "Open drains, not WWTP influent; site coordinates not in repository (publication Suppl. S7 not retrieved); CARD release not stated.",
         "DS03": "Four samples on one date; no linked publication found; not representative of Dhaka.",
         "DS04": "SRA holds 1,898 of 4,728 published samples; run platform (NovaSeq 6000) differs from publication (HiSeq X); surface swabs have low biomass.",
         "DS05": "Food product; metagenomic and enrichment libraries not separable from run metadata; not urban-located."}
sm["Limitations"] = sm.Dataset_ID.map(LIMIT).fillna("Study-level classification from title review; sample metadata not checked against publication.")
tlog("Harmonize", "Aggregated runs to one row per repository sample (BioSample/ENA sample or MG-RAST sample)", len(sm))
tlog("Harmonize", "DS01 geography/date taken from Munk et al. 2022 Suppl. Data 1 by run accession", int((sm.Dataset_ID == "DS01").sum()))
tlog("Harmonize", "DS04 city/coordinates joined from MetaSUB metadata on SRA sample_alias = uuid", int((sm.Dataset_ID == "DS04").sum()))
tlog("Harmonize", "DS04 country propagated from same-city SRA records where SRA country missing and unambiguous",
     int(sm.QC_Flags.str.contains("propagated").sum()))
tlog("Harmonize", "DS02 city joined from Zenodo meta.tsv on SRA sample_alias = sample_code", int((sm.Dataset_ID == "DS02").sum()))
tlog("Harmonize", "Location/date filled from NCBI BioSample attributes where ENA record lacked them", int(sm.Geographic_Source.str.contains("BioSample").sum()))
tlog("Harmonize", "Coordinates validated for numeric type and range; no coordinates assigned where absent", int((sm.Latitude != NR).sum()))
tlog("Harmonize", "Dates converted to ISO 8601 (partial dates kept at reported precision; ranges keep start date)", int((sm.Collection_Date != NR).sum()))

# ----------------------------------------------------------------- long tables
# Taxonomy (DS02 Bracken phylum)
alias_to_sample = sm[sm.Dataset_ID == "DS02"].set_index(sm[sm.Dataset_ID == "DS02"]._extra.map(lambda e: e.get("india_code")))
idx = {code: (row.Sample_ID, row.Run_ID) for code, row in alias_to_sample.iterrows()}
bp = pd.read_csv(os.path.join(R, "zenodo_17590646_india", "extracted", "data", "bracken_phylum.tsv"), sep="\t", dtype={"taxonomy_id": str})
bp = bp[bp.fileName.isin(idx)]
tax = pd.DataFrame({"Dataset_ID": "DS02", "Sample_ID": bp.fileName.map(lambda c: idx[c][0]), "Run_ID": bp.fileName.map(lambda c: idx[c][1]),
                    "Source_Sample_Code": bp.fileName, "Taxon_Name": bp["name"], "Taxon_ID_NCBI": bp.taxonomy_id,
                    "Taxonomic_Rank": bp.taxonomy_lvl.map({"P": "Phylum"}).fillna(bp.taxonomy_lvl),
                    "Abundance_Value": bp.new_est_reads, "Abundance_Unit": "Bracken estimated reads (count)",
                    "Relative_Abundance": bp.fraction_total_reads, "Relative_Abundance_Unit": "Fraction of Bracken-estimated reads at phylum level",
                    "Kraken_Assigned_Reads": bp.kraken_assigned_reads, "Annotation_Pipeline": "Kraken 2 + Bracken (Singh et al. 2026)",
                    "Reference_Database": "Kraken 2 core_nt index (k2_nt, downloaded 2023-05-02)", "Source_File": "Zenodo 17590646 data/bracken_phylum.tsv"})
tax.to_csv(os.path.join(H, "taxonomy_long.csv"), index=False)
tlog("Long tables", "Taxonomy_Long: DS02 Bracken phylum rows mapped to BioSample/run via sample_code", len(tax))

# Diversity (214K resource subset)
run_to_sample = {}
for _, s in sm.iterrows():
    for rr in s.Run_ID.split(";"): run_to_sample[rr] = (s.Sample_ID, s.Dataset_ID)
dv = pd.read_csv(os.path.join(H, "subset_6919377_diversity_core_runs.tsv"), sep="\t")
dv = dv[dv.run_accession.isin(run_to_sample)]
div = pd.DataFrame({"Dataset_ID": dv.run_accession.map(lambda x: run_to_sample[x][1]), "Sample_ID": dv.run_accession.map(lambda x: run_to_sample[x][0]),
                    "Run_ID": dv.run_accession, "Feature_Category": dv.category.map({"ARG": "ARG (ResFinder genes)", "Genera": "Bacterial genera (16S rRNA, Silva)", "Phyla": "Bacterial phyla (16S rRNA, Silva)"}),
                    "Richness_n": dv.n, "Shannon": dv.Shannon, "Simpson": dv.Simpson, "Total_Fragments": dv.total_fragments,
                    "Category_Fragments": dv.category_fragments, "Annotation_Pipeline": "KMA " + dv.kma_version.astype(str) + " read mapping (Martiny et al. resource)",
                    "Source_File": "Zenodo 6919377 diversity.tsv"})
div.to_csv(os.path.join(H, "diversity_long.csv"), index=False)
tlog("Long tables", "Diversity_Long: Zenodo 6919377 diversity rows for all included/contextual runs covered by the resource", len(div))

# AMR long: DS02 RGI
rgi = pd.read_csv(os.path.join(R, "zenodo_17590646_india", "extracted", "data", "rgi_bowtie.tsv"), sep="\t", dtype=str)
rgi = rgi[rgi.fileName.isin(idx)]
def num(s): return pd.to_numeric(s, errors="coerce")
a1 = pd.DataFrame({"Dataset_ID": "DS02", "Sample_ID": rgi.fileName.map(lambda c: idx[c][0]), "Run_ID": rgi.fileName.map(lambda c: idx[c][1]),
                   "Source_Sample_Code": rgi.fileName, "AMR_Gene": rgi["ARO Term"], "Gene_Identifier": "ARO:" + rgi["ARO Accession"],
                   "AMR_Gene_Family": rgi["AMR Gene Family"], "Resistance_Class": rgi["Drug Class"], "Resistance_Mechanism": rgi["Resistance Mechanism"],
                   "Reference_Model_Type": rgi["Reference Model Type"], "Detection_Value": num(rgi["All Mapped Reads"]),
                   "Detection_Unit": "All mapped reads (count)", "Completely_Mapped_Reads": num(rgi["Completely Mapped Reads"]),
                   "Average_Percent_Coverage": num(rgi["Average Percent Coverage"]), "Reference_Length_bp": num(rgi["Reference Length"]),
                   "AMR_Detection_Type": "Read mapping (RGI bwt, Bowtie2)", "AMR_Reference_Database": "CARD",
                   "AMR_Database_Version": "Not reported in deposited files", "Detection_Threshold": "Not reported; all RGI bwt output rows retained",
                   "Source_File": "Zenodo 17590646 data/rgi_bowtie.tsv"})
# 214K ResFinder
arg = pd.read_csv(os.path.join(H, "subset_6919377_ARG_core_runs.tsv"), sep="\t")
arg = arg[arg.run_accession.isin(run_to_sample)]
ann = pd.read_csv(os.path.join(R, "zenodo_6919377_global_resistome", "ResFinder_anno.tsv"), sep="\t")
cls = ann[ann.anno_type == "Class"].groupby("gene").anno_value.apply(lambda x: "; ".join(sorted(set(x)))).to_dict()
a2 = pd.DataFrame({"Dataset_ID": arg.run_accession.map(lambda x: run_to_sample[x][1]), "Sample_ID": arg.run_accession.map(lambda x: run_to_sample[x][0]),
                   "Run_ID": arg.run_accession, "Source_Sample_Code": NAP, "AMR_Gene": arg.refSequence.str.split("_").str[0],
                   "Gene_Identifier": arg.refSequence, "AMR_Gene_Family": "Not provided by source",
                   "Resistance_Class": arg.refSequence.map(cls).fillna(NR), "Resistance_Mechanism": "Not provided by source",
                   "Reference_Model_Type": NAP, "Detection_Value": arg.fragmentCountAln, "Detection_Unit": "Aligned read fragments (count)",
                   "Completely_Mapped_Reads": NAV, "Average_Percent_Coverage": (arg.refCoveredPositions / arg.refSequence_length * 100).round(2),
                   "Reference_Length_bp": arg.refSequence_length, "AMR_Detection_Type": "Read mapping (KMA " + arg.kma_version.astype(str) + ")",
                   "AMR_Reference_Database": "ResFinder", "AMR_Database_Version": "Not stated in table; see Martiny et al. resource documentation",
                   "Detection_Threshold": "As applied by resource (not restated in table)", "Source_File": "Zenodo 6919377 ARG.tsv (subset)"})
amrl = pd.concat([a1, a2], ignore_index=True)
amrl["Normalization_Note"] = np.where(amrl.Dataset_ID == "DS02", "Raw counts; normalize by gene length and sequencing depth before comparison",
                                      "Raw counts; trimmed_fragments and bacterial_fragment denominators available in source")
amrl.to_csv(os.path.join(H, "amr_genes_long.csv"), index=False)
tlog("Long tables", "AMR_Genes_Long: DS02 RGI rows + Zenodo 6919377 ResFinder rows; pipelines kept separate via AMR_Detection_Type", len(amrl))
tlog("Long tables", "Average_Percent_Coverage for ResFinder rows computed as refCoveredPositions/refSequence_length*100", len(a2))

# Environmental covariates (as reported)
cov = []
for _, s in sm.iterrows():
    e = s._extra
    if "ms" in e:
        ms = e["ms"]
        for k, unit, lvl in [("surface_material", "", "Sample"), ("surface_ontology_coarse", "", "Sample"), ("setting", "", "Sample"),
                             ("location_type", "", "Sample"), ("temperature", "As reported", "Sample"), ("elevation", "As reported", "Sample"),
                             ("coastal_city", "", "City"), ("city_total_population", "persons", "City"), ("city_population_density", "persons/km2", "City"),
                             ("city_ave_june_temp_c", "degC", "City"), ("city_koppen_climate", "", "City"), ("continent", "", "City")]:
            v = clean(ms.get(k, ""))
            if v: cov.append([s.Dataset_ID, s.Sample_ID, k, v, unit, lvl, "MetaSUB-metadata complete_metadata.csv"])
    if "batch" in e:
        cov.append([s.Dataset_ID, s.Sample_ID, "sampling_round", e["batch"], "", "Sample", "Munk et al. 2022 Suppl. Data 1"])
    if "india_code" in e:
        cov.append([s.Dataset_ID, s.Sample_ID, "sampling_month", e["india_month"], "month", "Sample", "Zenodo 17590646 meta.tsv"])
        cov.append([s.Dataset_ID, s.Sample_ID, "sampling_year_2digit", e["india_year"], "year (YY)", "Sample", "Zenodo 17590646 meta.tsv"])
for p in ["mgp16346", "mgp89046"]:
    d = json.load(open(os.path.join(R, "mgrast_metagenomes", f"{p}.json")))
    for m in d["data"]:
        sid = m.get("sample_id") or m["metagenome_id"]
        for k in ("biome", "feature", "material", "env_package_type", "depth", "temperature", "ph"):
            v = clean(m.get(k))
            if v: cov.append([ds_of[p], sid, k, v, "As reported", "Sample", "MG-RAST search API"])
covdf = pd.DataFrame(cov, columns=["Dataset_ID", "Sample_ID", "Covariate", "Value", "Unit", "Level", "Source"]).drop_duplicates()
covdf.to_csv(os.path.join(H, "env_covariates_long.csv"), index=False)
tlog("Long tables", "Environmental_Covariates: source-reported covariates only (no external covariates joined)", len(covdf))

# ----------------------------------------------------------------- save masters
sm = sm.drop(columns=["_extra"])
amr_n = amrl.groupby("Sample_ID").AMR_Gene.nunique()
sm["AMR_Gene"] = [f"See AMR_Genes_Long (n={amr_n[s]} distinct genes)" if s in amr_n.index else v for s, v in zip(sm.Sample_ID, sm.AMR_Gene)]
sm.to_csv(os.path.join(H, "sample_master.csv"), index=False)
runs.to_csv(os.path.join(H, "run_master.csv"), index=False)
pd.DataFrame(TLOG).to_csv(os.path.join(H, "transformation_log.csv"), index=False)
json.dump({k: v for k, v in ds_of.items()}, open(os.path.join(H, "dataset_ids.json"), "w"), indent=1)
json.dump(DS_NAME, open(os.path.join(H, "dataset_names.json"), "w"), indent=1)
print("samples", len(sm), "runs", len(runs), "tax", len(tax), "div", len(div), "amr", len(amrl), "cov", len(covdf))
print(sm.Inclusion_Status.value_counts().to_dict())
print(sm[sm.Inclusion_Status == "Included"].Sample_Type.value_counts().to_dict())
