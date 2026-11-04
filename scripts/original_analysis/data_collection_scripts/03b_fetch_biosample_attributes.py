"""
03b_fetch_biosample_attributes.py
For NCBI BioSamples (SAMN*) in the selected studies whose ENA run record lacks a
country/locality or collection date, retrieve the BioSample attributes from NCBI
E-utilities (official record). Raw XML is preserved; parsed attributes are written
to harmonized/biosample_attributes.csv. Nothing is inferred.
"""
import csv, glob, io, json, os, re, time, urllib.parse, urllib.request, html
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(B, "raw_downloads", "ncbi_biosample"); os.makedirs(OUT, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def need():
    acc = set()
    for f in glob.glob(os.path.join(B, "raw_downloads", "ena_sample_metadata", "*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            s = r["sample_accession"]
            miss_geo = (not r["country"]) or r["country"].lower().startswith("missing") or ":" not in r["country"]
            miss_date = (not r["collection_date"]) or r["collection_date"].lower().startswith("missing")
            if s.startswith("SAMN") and (miss_geo or miss_date):
                acc.add(s)
    return sorted(acc)

def fetch(batch, n):
    p = os.path.join(OUT, f"batch_{n:03d}.xml")
    if os.path.exists(p): return open(p).read()
    term = " OR ".join(f"{a}[accn]" for a in batch)
    d = json.load(urllib.request.urlopen(urllib.request.Request(EU + "esearch.fcgi",
        data=urllib.parse.urlencode({"db": "biosample", "term": term, "retmax": 500, "retmode": "json"}).encode()), timeout=120))
    uids = d["esearchresult"]["idlist"]; time.sleep(0.4)
    x = urllib.request.urlopen(urllib.request.Request(EU + "efetch.fcgi",
        data=urllib.parse.urlencode({"db": "biosample", "id": ",".join(uids)}).encode()), timeout=300).read().decode()
    open(p, "w").write(x); time.sleep(0.4)
    return x

def main():
    acc = need(); print("BioSamples needing attributes:", len(acc))
    rows = []
    for n, i in enumerate(range(0, len(acc), 150)):
        x = fetch(acc[i:i+150], n)
        for b in re.findall(r"<BioSample .*?</BioSample>", x, re.S):
            a = re.search(r'accession="([^"]+)"', b).group(1)
            attrs = {k: html.unescape(v) for k, v in re.findall(r'<Attribute attribute_name="([^"]+)"[^>]*>([^<]*)<', b)}
            rows.append({"sample_accession": a, "geo_loc_name": attrs.get("geo_loc_name", ""),
                         "collection_date": attrs.get("collection_date", ""), "lat_lon": attrs.get("lat_lon", ""),
                         "isolation_source": attrs.get("isolation_source", ""), "env_medium": attrs.get("env_medium", ""),
                         "env_local_scale": attrs.get("env_local_scale", ""), "env_broad_scale": attrs.get("env_broad_scale", "")})
        print("batch", n, "parsed", len(rows), flush=True)
    with open(os.path.join(B, "harmonized", "biosample_attributes.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

if __name__ == "__main__":
    main()
