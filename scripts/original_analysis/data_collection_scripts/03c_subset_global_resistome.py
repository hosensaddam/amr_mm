"""03c: subset Zenodo 6919377 ARG.tsv and diversity.tsv to every run in
harmonized/run_master.csv that the resource covers (written to
harmonized/global_resistome_core_runs.txt by this script). Originals untouched."""
import os, pandas as pd
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Z = os.path.join(B, "raw_downloads", "zenodo_6919377_global_resistome"); H = os.path.join(B, "harmonized")
runs = set(pd.read_csv(os.path.join(H, "run_master.csv"), dtype=str, usecols=["Run_ID"]).Run_ID)
meta = set(pd.read_csv(os.path.join(Z, "metadata.tsv"), sep="\t", usecols=["run_accession"], dtype=str).run_accession)
keep = runs & meta
open(os.path.join(H, "global_resistome_core_runs.txt"), "w").write("\n".join(sorted(keep)))
for src, dst in [("ARG.tsv", "subset_6919377_ARG_core_runs.tsv"), ("diversity.tsv", "subset_6919377_diversity_core_runs.tsv")]:
    n = k = 0
    with open(os.path.join(Z, src)) as f, open(os.path.join(H, dst), "w") as o:
        o.write(f.readline())
        for line in f:
            n += 1
            if line.split("\t", 1)[0] in keep: o.write(line); k += 1
    print(src, "scanned", n, "kept", k, "runs", len(keep))
