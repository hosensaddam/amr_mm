"""02b: print a compact review table of INSDC candidate studies with public
WGS/METAGENOMIC runs, for documented manual eligibility review."""
import json, os
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
c = json.load(open(os.path.join(B, "harmonized", "study_candidates.json")))
rows = []
for k, v in c.items():
    if v["kind"] != "INSDC_STUDY" or "ena" not in v: continue
    e = v["ena"]
    if e["n_wgs_metagenomic_runs"] == 0: continue
    rows.append((v["accession"], e["n_wgs_metagenomic_runs"], e["n_samples"], ";".join(e["countries"])[:40],
                 e["n_with_latlon"], e["n_with_date"], (e["study_title"] or "")[:110]))
rows.sort(key=lambda r: r[0])
for r in rows: print("\t".join(str(x) for x in r))
print("N with WGS/METAGENOMIC runs:", len(rows))
