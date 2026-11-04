"""
01_query_search.py
Systematic, query-wise discovery of public metagenomic datasets.

For every (database x query) pair this script:
  * calls the official public API (no scraping),
  * saves the raw response under raw_downloads/search_responses/,
  * records the exact request URL, filters, returned count, and number of
    records retrieved and screened,
  * applies a documented programmatic title/metadata screen (SCREEN_RULE)
    to flag "potentially relevant" records.

The programmatic screen is a triage step only. "Verified" status is assigned
later, per accession, from official repository records (02_verify_accessions.py).

Run:  python3 -I scripts/01_query_search.py
"""
import json, os, re, sys, time, hashlib, urllib.parse, urllib.request
from datetime import datetime, timezone

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(BASE, "raw_downloads", "search_responses")
LOG = os.path.join(BASE, "logs")
os.makedirs(RAW, exist_ok=True)
SEARCH_DATE = datetime.now(timezone.utc).strftime("%Y-%m-%d")
UA = {"User-Agent": "AMR-chapter-public-data-audit/1.0 (academic research; metadata only)"}

QUERY_GROUPS = {
    "QG1_Urban_microbiome": ["urban microbiome", "urban metagenomics", "urban microbial diversity",
                             "urban environmental microbiome", "built environment metagenomics",
                             "urban microbial communities"],
    "QG2_Urban_AMR_resistome": ["urban antimicrobial resistance", "urban resistome",
                                "antimicrobial resistance genes metagenomics", "environmental resistome",
                                "urban AMR gene mapping"],
    "QG3_Wastewater_sewage": ["urban wastewater metagenomics", "sewage metagenomic analysis",
                              "wastewater resistome", "hospital wastewater AMR",
                              "municipal sewage antimicrobial resistance", "urban wastewater microbiome"],
    "QG4_Environmental": ["urban soil metagenomics", "urban surface water microbiome",
                          "airborne microbiome urban", "environmental shotgun metagenomics",
                          "urban environmental resistome"],
    "QG5_Food_environments": ["food environment metagenomics", "food-associated antimicrobial resistance",
                              "food contact surface microbiome", "food market microbiome",
                              "food processing environmental resistome"],
    "QG6_Geographic": ["Dhaka urban microbiome", "Dhaka wastewater metagenomics", "Delhi urban microbiome",
                       "Delhi wastewater AMR", "Karachi environmental microbiome",
                       "Karachi wastewater metagenomics", "South Asian urban resistome",
                       "global urban microbiome mapping"],
}

# ---- documented programmatic screen ---------------------------------------
URBAN_CTX = r"urban|city|cities|metropol|municipal|sewage|sewer|wastewater|waste water|influent|hospital|subway|metro\b|transit|mass.transit|built environment|public transport|street|household|indoor|market|dhaka|delhi|karachi|mumbai|kolkata|chennai|lahore|kathmandu|colombo"
FOOD_CTX = r"food|meat|retail|slaughter|abattoir|kitchen|dairy|vegetable|produce|processing plant|poultry|seafood"
OFF_TARGET = r"\bgut\b|faec|fecal microbiota|stool|rumen|mouse|mice|patient|clinical isolate|tumou?r|oral microbiome|vagin|skin of"
SCREEN_RULE = ("Potentially relevant if title/description matches URBAN_CTX (or FOOD_CTX for QG5) and does not "
               "match OFF_TARGET, case-insensitive. URBAN_CTX=" + URBAN_CTX + " | FOOD_CTX=" + FOOD_CTX +
               " | OFF_TARGET=" + OFF_TARGET)

def screen(text, group):
    t = (text or "").lower()
    ctx = FOOD_CTX if group == "QG5_Food_environments" else URBAN_CTX
    return bool(re.search(ctx, t)) and not re.search(OFF_TARGET, t)

def get(url, as_json=True, tries=3):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            return (json.loads(body) if as_json else body.decode("utf-8", "replace")), None
        except Exception as e:  # noqa
            last = repr(e)
            time.sleep(3 * (i + 1))
    return None, last

def save_raw(qid, payload):
    p = os.path.join(RAW, f"{qid}.json")
    with open(p, "w") as f:
        json.dump(payload, f)
    return os.path.relpath(p, BASE)

# ---- per-database query functions -----------------------------------------
# Each returns: url, filters, n_returned(total reported by API), records(list of dicts with id,title,text)

def q_mgnify(term):
    params = {"search": term, "page_size": 100}
    url = "https://www.ebi.ac.uk/metagenomics/api/v2/studies/?" + urllib.parse.urlencode(params)
    recs, total, nxt, page, err = [], None, url, 1, None
    while nxt:
        d, err = get(nxt)
        if d is None: break
        total = d.get("count", total)
        for it in d.get("items", []):
            md = it.get("metadata") or {}
            recs.append({"id": it.get("accession"),
                         "insdc": ";".join(it.get("ena_accessions") or []),
                         "title": it.get("title"),
                         "text": f"{it.get('title')} {md.get('study_description','')} {(it.get('biome') or {}).get('lineage','')}",
                         "biome": (it.get("biome") or {}).get("lineage")})
        page += 1
        if total and len(recs) < total and page <= 20:
            nxt = "https://www.ebi.ac.uk/metagenomics/api/v2/studies/?" + urllib.parse.urlencode({**params, "page": page})
        else:
            nxt = None
    return url, "studies endpoint; full-text 'search'; page_size=100; all pages (cap 20)", total, recs, err

def q_mgrast(term):
    # MG-RAST search API: 'all' = free-text across metagenome metadata; public data only.
    params = {"all": term, "limit": 1000}
    url = "https://api.mg-rast.org/search?" + urllib.parse.urlencode(params)
    d, err = get(url)
    if d is None: return url, "metagenome_index; all=<term>; limit=1000", None, [], err
    total = d.get("total_count")
    recs = []
    for it in d.get("data", []):
        recs.append({"id": it.get("metagenome_id"), "project": it.get("project_id"),
                     "title": it.get("project_name"),
                     "text": " ".join(str(it.get(k, "")) for k in ("project_name", "biome", "feature", "material",
                                                                    "env_package_type", "location", "country", "city")),
                     "sequence_type": it.get("sequence_type"), "public": it.get("public")})
    return url, "metagenome_index; all=<term>; limit=1000 (first page only)", total, recs, err

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/search?"
def q_ena(term):
    # Words AND-ed across study_title; restricted to shotgun metagenomic reads.
    words = [w for w in re.split(r"[\s\-]+", term) if len(w) > 2 and w.lower() not in ("and", "the")]
    w2 = []
    for w in words:
        wl = w.lower()
        stem = {"metagenomics": "metagenom", "metagenomic": "metagenom", "microbiome": "microbio",
                "microbial": "microb", "communities": "communit", "resistome": "resistom",
                "antimicrobial": "antimicrob", "airborne": "air"}.get(wl, wl)
        w2.append(stem)
    title_q = " AND ".join(f'study_title="*{w}*"' for w in w2)
    query = f'({title_q}) AND library_source="METAGENOMIC" AND library_strategy="WGS"'
    params = {"result": "read_run", "query": query, "fields": "study_accession,study_title,country",
              "format": "json", "limit": 0}
    url = ENA_PORTAL + urllib.parse.urlencode(params)
    d, err = get(url)
    if d is None: return url, query, None, [], err
    studies = {}
    for r in d:
        s = studies.setdefault(r.get("study_accession"), {"id": r.get("study_accession"), "title": r.get("study_title"),
                                                            "text": r.get("study_title"), "n_runs": 0, "countries": set()})
        s["n_runs"] += 1
        if r.get("country"): s["countries"].add(r["country"].split(":")[0])
    recs = []
    for s in studies.values():
        s["countries"] = ";".join(sorted(s["countries"]))
        recs.append(s)
    return url, "read_run; " + query + "; result aggregated to distinct study_accession", len(recs), recs, err

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def q_ncbi_bioproject(term):
    t = f'({term}) AND (metagenome[All Fields] OR metagenomic[All Fields]) AND "bioproject sra"[Filter]'
    params = {"db": "bioproject", "term": t, "retmax": 500, "retmode": "json"}
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params)
    d, err = get(url)
    if d is None: return url, t, None, [], err
    total = int(d["esearchresult"].get("count", 0)); ids = d["esearchresult"].get("idlist", [])
    recs = []
    for i in range(0, len(ids), 200):
        su = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode({"db": "bioproject", "id": ",".join(ids[i:i+200]), "retmode": "json"})
        s, e2 = get(su); time.sleep(0.4)
        if s is None: err = e2; continue
        for uid in s["result"].get("uids", []):
            r = s["result"][uid]
            recs.append({"id": r.get("project_acc"), "title": r.get("project_title"),
                         "text": f"{r.get('project_title')} {r.get('project_description','')}",
                         "organism": r.get("organism_name")})
    return url, t + "; retmax=500", total, recs, err

def q_europepmc(term):
    q = f'({term}) AND (ACCESSION_TYPE:bioproject OR ACCESSION_TYPE:ena) AND (metagenom* OR resistome)'
    params = {"query": q, "format": "json", "pageSize": 100, "resultType": "lite"}
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(params)
    d, err = get(url)
    if d is None: return url, q, None, [], err
    total = d.get("hitCount")
    recs = [{"id": r.get("doi") or r.get("pmid") or r.get("id"), "title": r.get("title"),
             "text": r.get("title"), "year": r.get("pubYear"), "journal": r.get("journalTitle"),
             "pmid": r.get("pmid"), "pmcid": r.get("pmcid")} for r in d.get("resultList", {}).get("result", [])]
    return url, q + "; first 100 records", total, recs, err

def q_zenodo(term):
    params = {"q": " AND ".join(re.split(r"[\s]+", term)), "resource_type.type": "dataset", "size": 25, "sort": "bestmatch"}
    url = "https://zenodo.org/api/records?" + urllib.parse.urlencode(params)
    d, err = get(url)
    if d is None: return url, "resource_type.type=dataset", None, [], err
    total = (d.get("hits") or {}).get("total")
    recs = []
    for r in (d.get("hits") or {}).get("hits", []):
        md = r.get("metadata", {})
        recs.append({"id": r.get("doi") or str(r.get("id")), "title": md.get("title"),
                     "text": f"{md.get('title')} {re.sub('<[^>]+>', ' ', md.get('description',''))[:1500]}",
                     "licence": (md.get("license") or {}).get("id"), "date": md.get("publication_date")})
    return url, "all words required (AND); resource_type.type=dataset; size=25 (anonymous API page cap); sort=bestmatch", total, recs, err

DBS = [("MGnify (API v2)", "MGN", q_mgnify), ("MG-RAST", "MGR", q_mgrast), ("ENA Portal", "ENA", q_ena),
       ("NCBI BioProject", "NCB", q_ncbi_bioproject), ("Europe PMC", "EPM", q_europepmc), ("Zenodo", "ZEN", q_zenodo)]

def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    out_path = os.path.join(LOG, f"query_log_{only or 'ALL'}.jsonl")
    done = set()
    if os.path.exists(out_path):
        for line in open(out_path):
            done.add(json.loads(line)["Query_ID"])
    out = open(out_path, "a")
    for dbname, code, fn in DBS:
        if only and code != only: continue
        n = 0
        for group, terms in QUERY_GROUPS.items():
            for term in terms:
                n += 1
                qid = f"{code}-{group[:3]}-{n:02d}"
                if qid in done: continue
                url, filters, total, recs, err = fn(term)
                flags = [screen(r.get("text"), group) for r in recs]
                raw_rel = save_raw(qid, {"query_id": qid, "url": url, "total": total, "records": recs,
                                         "error": err})
                row = {"Query_ID": qid, "Database": dbname, "Query_Group": group, "Search_String": term,
                       "Search_Date": SEARCH_DATE, "Search_URL_or_Endpoint": url, "Filters_Used": filters,
                       "N_Results_Returned": total, "N_Records_Inspected": len(recs),
                       "N_Potentially_Relevant": sum(flags), "Screen_Rule": "SCREEN_RULE v1 (see README)",
                       "Raw_Response_File": raw_rel, "API_Error": err,
                       "Relevant_IDs": [r["id"] for r, f in zip(recs, flags) if f]}
                out.write(json.dumps(row) + "\n"); out.flush()
                print(qid, dbname, repr(term), "total=", total, "inspected=", len(recs), "relevant=", sum(flags),
                      "ERR" if err else "", flush=True)
                time.sleep(0.5)
    out.close()
    with open(os.path.join(LOG, "screen_rule.txt"), "w") as f:
        f.write(SCREEN_RULE + "\n")

if __name__ == "__main__":
    main()
