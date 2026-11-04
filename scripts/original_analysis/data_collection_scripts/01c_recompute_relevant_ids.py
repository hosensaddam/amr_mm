"""01c: recompute Relevant_IDs in every query log from the saved raw responses
using the same screen() as 01_query_search.py (corrects a list-filter bug in the
first run; counts N_Potentially_Relevant were already correct and are re-checked)."""
import json, glob, os, importlib.util
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("qs", os.path.join(B, "scripts", "01_query_search.py"))
qs = importlib.util.module_from_spec(spec); spec.loader.exec_module(qs)
mismatch = 0
for f in sorted(glob.glob(os.path.join(B, "logs", "query_log_*.jsonl"))):
    rows = [json.loads(l) for l in open(f)]
    for r in rows:
        raw = json.load(open(os.path.join(B, r["Raw_Response_File"])))
        ids = [x["id"] for x in raw["records"] if qs.screen(x.get("text"), r["Query_Group"])]
        if len(ids) != r["N_Potentially_Relevant"]: mismatch += 1
        r["Relevant_IDs"] = ids
    with open(f, "w") as o:
        for r in rows: o.write(json.dumps(r) + "\n")
print("count mismatches:", mismatch)
