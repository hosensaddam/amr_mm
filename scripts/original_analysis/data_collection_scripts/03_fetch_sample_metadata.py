"""
03_fetch_sample_metadata.py
Download run-level metadata (one row per sequencing run, extended fields) from
the ENA Portal API for every study listed in harmonized/selected_studies.json.
Raw responses are preserved unmodified in raw_downloads/ena_sample_metadata/.
No sequence data are downloaded; fastq_ftp/fastq_bytes are recorded only.
"""
import json, os, time, urllib.parse, importlib.util
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("qs", os.path.join(B, "scripts", "01_query_search.py"))
qs = importlib.util.module_from_spec(spec); spec.loader.exec_module(qs)
OUT = os.path.join(B, "raw_downloads", "ena_sample_metadata"); os.makedirs(OUT, exist_ok=True)
FIELDS = ("run_accession,experiment_accession,sample_accession,secondary_sample_accession,study_accession,"
          "secondary_study_accession,study_title,sample_alias,sample_title,sample_description,library_name,"
          "experiment_title,library_strategy,library_source,library_selection,library_layout,instrument_platform,"
          "instrument_model,scientific_name,tax_id,country,location,lat,lon,collection_date,isolation_source,"
          "environment_biome,environment_feature,environment_material,broad_scale_environmental_context,"
          "local_environmental_context,environmental_medium,host,first_public,last_updated,base_count,read_count,"
          "fastq_bytes,fastq_md5,fastq_ftp,center_name,status")

def main():
    sel = json.load(open(os.path.join(B, "harmonized", "selected_studies.json")))
    for acc in sel:
        p = os.path.join(OUT, f"{acc}.tsv")
        if os.path.exists(p) and os.path.getsize(p) > 0:
            continue
        params = {"result": "read_run", "query": f'study_accession="{acc}" OR secondary_study_accession="{acc}"',
                  "fields": FIELDS, "format": "tsv", "limit": 0}
        url = qs.ENA_PORTAL + urllib.parse.urlencode(params)
        txt, err = qs.get(url, as_json=False)
        if txt is None:
            print(acc, "ERROR", err); continue
        open(p, "w").write(txt)
        open(p + ".url", "w").write(url + "\n")
        print(acc, txt.count("\n") - 1, "runs", flush=True)
        time.sleep(0.3)

if __name__ == "__main__":
    main()
