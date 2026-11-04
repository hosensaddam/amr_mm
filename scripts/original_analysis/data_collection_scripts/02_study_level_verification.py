"""
02_study_level_verification.py
(1) Deduplicate potentially relevant records from all query logs into a
    study-level candidate set keyed by stable accession.
(2) For every INSDC study candidate (PRJNA/PRJEB/PRJDB/PRJDA), retrieve its
    public run records from the ENA Portal API (read_run) and summarise data
    type, platform, geography and dates. ENA mirrors all INSDC projects,
    so NCBI-, ENA- and MGnify-discovered projects are checked the same way.
(3) Add the brief-specified discovery leads explicitly, even if not retrieved
    by the queries.
Outputs: harmonized/study_candidates.json, raw_downloads/ena_study_runs/*.tsv
"""
import csv, io, json, glob, os, re, time, urllib.parse
import importlib.util
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("qs", os.path.join(B, "scripts", "01_query_search.py"))
qs = importlib.util.module_from_spec(spec); spec.loader.exec_module(qs)
OUT = os.path.join(B, "raw_downloads", "ena_study_runs"); os.makedirs(OUT, exist_ok=True)

LEADS = {  # brief-specified discovery leads (Section 3D of the brief)
    "PRJEB40798": "Munk et al. 2022 lead", "PRJEB40816": "Munk et al. 2022 lead", "PRJEB40815": "Munk et al. 2022 lead",
    "PRJEB27621": "Munk et al. 2022 lead", "PRJEB51229": "Munk et al. 2022 lead", "ERP015409": "Hendriksen/Munk lead (secondary study accession)",
    "PRJNA1415085": "Indian metro wastewater lead", "PRJNA1399457": "Dhaka hospital wastewater lead",
    "PRJNA1522440": "FDA food-associated lead", "PRJNA732392": "MetaSUB SRA project (identified 31 Oct 2026 in project files)",
    "PRJEB27054": "Underlies MGnify MGYS00005846 (identified 31 Oct 2026 in project files)",
}
FIELDS = ("run_accession,sample_accession,secondary_sample_accession,study_accession,secondary_study_accession,"
          "study_title,library_strategy,library_source,library_selection,library_layout,instrument_platform,"
          "instrument_model,country,location,lat,lon,collection_date,first_public,base_count,read_count,"
          "fastq_bytes,fastq_ftp,submitted_ftp,center_name,scientific_name,sample_title,broker_name")

def candidates():
    cand = {}
    for f in sorted(glob.glob(os.path.join(B, "logs", "query_log_*.jsonl"))):
        for l in open(f):
            r = json.loads(l)
            raw = json.load(open(os.path.join(B, r["Raw_Response_File"])))
            recs = {x["id"]: x for x in raw["records"]}
            for i in r["Relevant_IDs"]:
                x = recs.get(i, {}); db = r["Database"]
                if db.startswith("MG-RAST"): kind, key = "MGRAST_PROJECT", x.get("project")
                elif db.startswith("Europe PMC"): kind, key = "PUBLICATION", i
                elif db.startswith("Zenodo"): kind, key = "ZENODO_RECORD", i
                elif db.startswith("MGnify"):
                    ins = [a for a in (x.get("insdc") or "").split(";") if a.startswith("PRJ")]
                    kind, key = "INSDC_STUDY", (ins[0] if ins else i)
                else: kind, key = "INSDC_STUDY", i
                c = cand.setdefault(f"{kind}|{key}", {"kind": kind, "accession": key, "queries": set(),
                                                     "titles": set(), "dbs": set(), "mgnify": set(), "meta": {}})
                c["queries"].add(r["Query_ID"]); c["titles"].add(x.get("title") or ""); c["dbs"].add(db)
                if db.startswith("MGnify"): c["mgnify"].add(i)
                for k in ("year", "journal", "pmid", "pmcid", "licence", "date", "biome"):
                    if x.get(k): c["meta"][k] = x[k]
    for acc, why in LEADS.items():
        k = f"INSDC_STUDY|{acc}"
        c = cand.setdefault(k, {"kind": "INSDC_STUDY", "accession": acc, "queries": set(), "titles": set(),
                                "dbs": set(), "mgnify": set(), "meta": {}})
        c["dbs"].add("Brief discovery lead"); c["meta"]["lead"] = why
    return cand

def ena_runs(acc):
    p = os.path.join(OUT, f"{acc}.tsv")
    if os.path.exists(p) and os.path.getsize(p) > 0:
        return open(p).read(), None
    params = {"result": "read_run", "query": f'study_accession="{acc}" OR secondary_study_accession="{acc}"',
              "fields": FIELDS, "format": "tsv", "limit": 0}
    url = qs.ENA_PORTAL + urllib.parse.urlencode(params)
    txt, err = qs.get(url, as_json=False)
    if txt is not None:
        open(p, "w").write(txt)
    return txt, err

def summarise(txt):
    rows = list(csv.DictReader(io.StringIO(txt), delimiter="\t")) if txt else []
    s = {"n_runs": len(rows), "n_samples": len({r["sample_accession"] for r in rows}),
         "n_wgs_metagenomic_runs": sum(1 for r in rows if r["library_strategy"] == "WGS" and r["library_source"] == "METAGENOMIC"),
         "strategies": sorted({f'{r["library_strategy"]}/{r["library_source"]}' for r in rows}),
         "platforms": sorted({r["instrument_model"] for r in rows if r["instrument_model"]}),
         "countries": sorted({r["country"].split(":")[0].strip() for r in rows if r["country"]}),
         "locations": sorted({r["country"] for r in rows if r["country"]})[:50],
         "n_with_latlon": sum(1 for r in rows if r["lat"] and r["lon"]),
         "n_with_date": sum(1 for r in rows if r["collection_date"]),
         "dates": sorted({r["collection_date"][:7] for r in rows if r["collection_date"]}),
         "base_count_total": sum(int(r["base_count"]) for r in rows if r["base_count"].isdigit()),
         "fastq_bytes_total": sum(int(b) for r in rows for b in r["fastq_bytes"].split(";") if b.isdigit()),
         "study_title": rows[0]["study_title"] if rows else "", "center": rows[0]["center_name"] if rows else "",
         "study_accession": rows[0]["study_accession"] if rows else "",
         "secondary_study_accession": rows[0]["secondary_study_accession"] if rows else ""}
    return s

def main():
    cand = candidates()
    ins = [c for c in cand.values() if c["kind"] == "INSDC_STUDY" and re.match(r"^(PRJ|ERP|SRP|DRP)", c["accession"] or "")]
    print("candidates:", {k: sum(1 for c in cand.values() if c["kind"] == k) for k in {c["kind"] for c in cand.values()}})
    print("INSDC studies to check:", len(ins))
    for n, c in enumerate(ins, 1):
        txt, err = ena_runs(c["accession"])
        c["ena"] = summarise(txt); c["ena_error"] = err
        c["ena_url"] = f"https://www.ebi.ac.uk/ena/browser/view/{c['accession']}"
        if n % 50 == 0: print(n, flush=True)
        time.sleep(0.2)
    ser = {k: {**v, "queries": sorted(v["queries"]), "titles": sorted(t for t in v["titles"] if t),
               "dbs": sorted(v["dbs"]), "mgnify": sorted(v["mgnify"])} for k, v in cand.items()}
    json.dump(ser, open(os.path.join(B, "harmonized", "study_candidates.json"), "w"), indent=1)
    print("done")

if __name__ == "__main__":
    main()
