"""
01b_mgnify_adapted_queries.py
Adapted MGnify v2 search strategy.

Reason: in the main run (01_query_search.py) every multi-word string returned
count=0 from /api/v2/studies/?search=..., whereas single words return results
(e.g. 'sewage' -> 22 studies). The v2 'search' parameter therefore behaves as a
contiguous phrase match. This script runs (a) single-keyword searches derived
from the six query groups and (b) biome-lineage browsing of the urban-relevant
MGnify biomes, logging each as its own query.
"""
import json, os, sys, time, urllib.parse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location("qs", os.path.join(os.path.dirname(os.path.abspath(__file__)), "01_query_search.py"))
qs = importlib.util.module_from_spec(spec); spec.loader.exec_module(qs)

KEYWORDS = [("QG1_Urban_microbiome", "urban"), ("QG1_Urban_microbiome", "city"),
            ("QG1_Urban_microbiome", "subway"), ("QG1_Urban_microbiome", "metasub"),
            ("QG2_Urban_AMR_resistome", "resistome"), ("QG2_Urban_AMR_resistome", "antibiotic resistance"),
            ("QG2_Urban_AMR_resistome", "antimicrobial resistance"),
            ("QG3_Wastewater_sewage", "sewage"), ("QG3_Wastewater_sewage", "wastewater"),
            ("QG3_Wastewater_sewage", "hospital"), ("QG3_Wastewater_sewage", "influent"),
            ("QG4_Environmental", "air"), ("QG4_Environmental", "river"),
            ("QG5_Food_environments", "food"), ("QG5_Food_environments", "meat"),
            ("QG6_Geographic", "Bangladesh"), ("QG6_Geographic", "India"), ("QG6_Geographic", "Pakistan"),
            ("QG6_Geographic", "Dhaka"), ("QG6_Geographic", "Delhi"), ("QG6_Geographic", "Karachi")]
BIOMES = [("QG1_Urban_microbiome", "root:Engineered:Built environment"),
          ("QG3_Wastewater_sewage", "root:Engineered:Wastewater:Water and sludge"),
          ("QG4_Environmental", "root:Environmental:Air:Outdoor Air"),
          ("QG4_Environmental", "root:Environmental:Air:Indoor Air"),
          ("QG5_Food_environments", "root:Engineered:Food production")]

def browse_biome(lineage):
    params = {"biome_lineage": lineage, "page_size": 100}
    url = "https://www.ebi.ac.uk/metagenomics/api/v2/studies/?" + urllib.parse.urlencode(params)
    recs, total, page, err = [], None, 1, None
    while True:
        u = "https://www.ebi.ac.uk/metagenomics/api/v2/studies/?" + urllib.parse.urlencode({**params, "page": page})
        d, err = qs.get(u)
        if d is None: break
        total = d.get("count")
        for it in d.get("items", []):
            md = it.get("metadata") or {}
            recs.append({"id": it.get("accession"), "insdc": ";".join(it.get("ena_accessions") or []),
                         "title": it.get("title"), "biome": (it.get("biome") or {}).get("lineage"),
                         "text": f"{it.get('title')} {md.get('study_description','')}"})
        page += 1
        if not total or len(recs) >= total or page > 30: break
    return url, f"studies endpoint; biome_lineage={lineage} (includes child biomes); all pages (cap 30)", total, recs, err

def main():
    out = open(os.path.join(qs.LOG, "query_log_MGNADJ.jsonl"), "w")
    n = 0
    jobs = [("kw", g, k) for g, k in KEYWORDS] + [("biome", g, b) for g, b in BIOMES]
    for kind, group, term in jobs:
        n += 1
        qid = f"MGN-ADJ-{n:02d}"
        url, filters, total, recs, err = qs.q_mgnify(term) if kind == "kw" else browse_biome(term)
        flags = [qs.screen(r.get("text"), group) for r in recs]
        raw = qs.save_raw(qid, {"query_id": qid, "url": url, "total": total, "records": recs, "error": err})
        row = {"Query_ID": qid, "Database": "MGnify (API v2) - adapted", "Query_Group": group,
               "Search_String": term if kind == "kw" else f"[biome browse] {term}", "Search_Date": qs.SEARCH_DATE,
               "Search_URL_or_Endpoint": url, "Filters_Used": filters, "N_Results_Returned": total,
               "N_Records_Inspected": len(recs), "N_Potentially_Relevant": sum(flags),
               "Screen_Rule": "SCREEN_RULE v1 (see README)", "Raw_Response_File": raw, "API_Error": err,
               "Relevant_IDs": [r["id"] for r, f in zip(recs, flags) if f]}
        out.write(json.dumps(row) + "\n"); out.flush()
        print(qid, repr(term), total, len(recs), sum(flags), err or "")
        time.sleep(0.5)

if __name__ == "__main__":
    main()
