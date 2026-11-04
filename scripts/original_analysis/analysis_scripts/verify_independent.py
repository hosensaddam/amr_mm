"""Independent technical, statistical, data-integrity, Excel and scientific verification.
Re-derives key numbers from the original harmonized files with separate code paths (csv module /
plain numpy) and compares them with the reported outputs. Usage:
    python verify_independent.py full     -> verification/V01_verification_checks.csv (+ report)
    python verify_independent.py excel    -> verification/V02_final_workbook_checks.csv (final workbook only)"""
import csv, json, math, sys, re, hashlib, py_compile, random, datetime as dt
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent))
from amr_common import ROOT, HARM, OUT, TAB, DER, FIG, VER, LOG, CAT, md5, load_catalog
csv.field_size_limit(10 ** 9)
checks = []
def chk(cat, check, expected, observed, ok, note=""):
    checks.append(dict(category=cat, check=check, expected=str(expected)[:300], observed=str(observed)[:300], status="PASS" if ok else "FAIL", note=note))
def close(a, b, rtol=1e-6, atol=1e-9): return abs(float(a) - float(b)) <= atol + rtol * abs(float(b))

def code_checks():
    log = [json.loads(l) for l in open(LOG / "notebook_run_log.jsonl", encoding="utf-8")]
    last = {}
    for r in log: last[r["notebook"]] = r
    nbs = sorted(p.stem for p in (OUT / "notebook_sources").glob("*.py"))
    chk("A. Code", "Every notebook's most recent fresh-kernel run succeeded", "success x %d" % len(nbs), {k: v["status"] for k, v in last.items()}, all(last.get(n, {}).get("status") == "success" for n in nbs))
    nerr = 0
    for nb in (OUT / "notebooks").glob("*.ipynb"):
        j = json.load(open(nb, encoding="utf-8"))
        nerr += sum(o.get("output_type") == "error" for c in j["cells"] if c["cell_type"] == "code" for o in c.get("outputs", []))
    chk("A. Code", "No error outputs stored in executed notebooks", 0, nerr, nerr == 0)
    bad = []
    for f in list((OUT / "scripts").glob("*.py")) + list((OUT / "scripts" / "nb_exports").glob("*.py")):
        try: py_compile.compile(str(f), doraise=True)
        except Exception as e: bad.append(f.name)
    chk("A. Code", "All scripts and exported notebook scripts compile", "0 failures", bad or "0 failures", not bad)
    figs = load_catalog("figures"); missing = [f for f in figs.file if not (FIG / f).exists()]
    chk("A. Code", "Every catalogued figure file exists", len(figs), len(figs) - len(missing), not missing)
    nums = sorted(figs.number.astype(int)); chk("A. Code", "Figure numbers are unique and consecutive", f"1..{len(nums)}", nums[:3] + ["..."] + nums[-2:], nums == list(range(1, len(nums) + 1)))
    tabs = load_catalog("tables"); miss = [t for t in tabs.file if not (TAB / t).exists()]
    chk("A. Code", "Every catalogued table file exists", len(tabs), len(tabs) - len(miss), not miss)

def read_csv_rows(path, cols, filt=None, delim=","):
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter=delim)
        for row in r:
            if filt is None or filt(row): yield {c: row[c] for c in cols}

def stat_checks():
    MISS = {"Not reported", "Not available", "Not yet verified", "Not applicable", "Restricted access", ""}
    # 1. evidence base
    inc = list(read_csv_rows(HARM / "sample_master.csv", ["Sample_ID", "Country", "City", "Sample_Type", "Dataset_ID"], lambda r: r["Inclusion_Status"] == "Included"))
    d01 = pd.read_csv(TAB / "D01_evidence_base_overview.csv").set_index("Quantity").iloc[:, 0]
    chk("B. Statistics", "Included samples (csv module) = D01", d01["Samples"], len(inc), len(inc) == int(d01["Samples"]))
    nc = len({r["Country"] for r in inc if r["Country"] not in MISS})
    chk("B. Statistics", "Countries among included samples = D01", d01["Countries (reported)"], nc, nc == int(d01["Countries (reported)"]))
    ncity = len({(r["City"], r["Country"]) for r in inc if r["City"] not in MISS and r["Country"] not in MISS})
    chk("B. Statistics", "Source-labelled City|Country pairs = D01", d01["Source-labelled cities (City|Country)"], ncity, ncity == int(d01["Source-labelled cities (City|Country)"]))
    st = defaultdict(int)
    for r in inc: st[r["Sample_Type"]] += 1
    d02 = pd.read_csv(TAB / "D02_by_sample_type.csv").set_index("Sample_Type").samples
    chk("B. Statistics", "Sample-type counts = D02", dict(d02), dict(st), all(st[k] == v for k, v in d02.items()) and sum(st.values()) == d02.sum())
    # 2. DS02 RPKM recomputed independently
    reads = {r["Sample_ID"]: int(r["Read_Count"]) * (2 if r["Read_Type"] == "PAIRED" else 1) for r in read_csv_rows(HARM / "run_master.csv", ["Sample_ID", "Read_Count", "Read_Type"], lambda r: r["Dataset_ID"] == "DS02")}
    tot = defaultdict(float)
    for r in read_csv_rows(HARM / "amr_genes_long.csv", ["Sample_ID", "Detection_Value", "Reference_Length_bp"], lambda r: r["AMR_Reference_Database"] == "CARD"):
        tot[r["Sample_ID"]] += float(r["Detection_Value"]) / (float(r["Reference_Length_bp"]) / 1000) / (reads[r["Sample_ID"]] / 1e6)
    s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0)
    diff = max(abs(tot[k] - v) / v for k, v in s2.total_RPKM.items())
    chk("B. Statistics", "DS02 total ARG RPKM: independent recomputation vs analysis (max relative difference)", "< 1e-9", f"{diff:.2e} over {len(s2)} samples", diff < 1e-9 and len(tot) == len(s2))
    a04 = pd.read_csv(TAB / "A04_ds02_arg_load_by_city.csv").set_index("City")
    city = {r["Sample_ID"]: r["City"] for r in read_csv_rows(HARM / "sample_master.csv", ["Sample_ID", "City"], lambda r: r["Dataset_ID"] == "DS02")}
    med = {c: float(np.median([v for k, v in tot.items() if city[k] == c])) for c in a04.index}
    chk("B. Statistics", "DS02 median RPKM by city = A04", a04.median_total_RPKM.round(3).to_dict(), {k: round(v, 3) for k, v in med.items()}, all(close(med[c], a04.loc[c, "median_total_RPKM"], 1e-4) for c in a04.index), "relative tolerance 1e-4 = precision of saved table (5 significant digits)")
    # 3. Kruskal-Wallis H by own rank implementation
    x = np.array([math.log10(tot[k]) for k in s2.index]); g = np.array([city[k] for k in s2.index])
    rk = stats.rankdata(x); N = len(x); H = 12 / (N * (N + 1)) * sum(rk[g == c].sum() ** 2 / (g == c).sum() for c in set(g)) - 3 * (N + 1)
    _, cnt = np.unique(x, return_counts=True); H /= 1 - (cnt ** 3 - cnt).sum() / (N ** 3 - N)
    a05 = pd.read_csv(TAB / "A05_ds02_arg_kruskal.csv"); rep = a05[(a05.variable == "log10_total_RPKM") & (a05.factor == "City")].iloc[0]
    chk("B. Statistics", "Kruskal-Wallis H (DS02 load by city): own implementation vs A05", round(rep.H, 4), round(H, 4), close(H, rep.H, 1e-6))
    chk("B. Statistics", "Kruskal-Wallis p from chi-square(3) vs A05", f"{rep.p:.3e}", f"{stats.chi2.sf(H, 3):.3e}", close(stats.chi2.sf(H, 3), rep.p, 1e-4, 1e-300))
    # 4. ResFinder FPKM recomputed from the TSV with own aggregation
    rmap = {r["Run_ID"]: (r["unit_id"], r["primary"]) for r in read_csv_rows(DER / "rf_run_to_unit.csv", ["Run_ID", "unit_id", "primary"])}
    frag = defaultdict(float); bact = {}
    for r in read_csv_rows(HARM / "subset_6919377_ARG_core_runs.tsv", ["run_accession", "fragmentCountAln", "refSequence_length", "bacterial_fragment"], delim="\t"):
        u, p = rmap[r["run_accession"]]
        if p != "True": continue
        frag[u] += float(r["fragmentCountAln"]) / (float(r["refSequence_length"]) / 1000); bact[(u, r["run_accession"])] = float(r["bacterial_fragment"])
    bu = defaultdict(float)
    for (u, _), b in bact.items(): bu[u] += b
    fp = {u: frag[u] / (bu[u] / 1e6) for u in frag}
    us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0)
    d = max(abs(fp[u] - v) / v for u, v in us.total_FPKM.items())
    chk("B. Statistics", "ResFinder unit total FPKM: independent aggregation vs analysis (max relative difference)", "< 1e-9", f"{d:.2e} over {len(us)} units", d < 1e-9 and set(fp) == set(us.index))
    # 5. Wilson CI by closed-form formula
    a09 = pd.read_csv(TAB / "A09_ds02_clinical_determinants.csv"); r = a09[(a09.determinant == "OXA-48-like") & (a09.group == "All")].iloc[0]
    k, n, z = r.detected_cov50, r.n, 1.959963984540054; ph = k / n; den = 1 + z * z / n
    lo = (ph + z * z / (2 * n) - z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))) / den * 100
    chk("B. Statistics", "Wilson 95% CI lower bound (OXA-48-like) by closed form vs A09", round(r.ci_low, 2), round(lo, 2), abs(lo - r.ci_low) < 0.01)
    # 6. prokaryotic Shannon from taxonomy_long
    dm = pd.read_csv(DER / "ds02_phylum_domain_map.csv").set_index("Phylum").Domain.to_dict()
    cnts = defaultdict(lambda: defaultdict(float))
    for r in read_csv_rows(HARM / "taxonomy_long.csv", ["Sample_ID", "Taxon_Name", "Abundance_Value"]):
        if dm[r["Taxon_Name"]] in ("Bacteria", "Archaea"): cnts[r["Sample_ID"]][r["Taxon_Name"]] += float(r["Abundance_Value"])
    al = pd.read_csv(DER / "ds02_alpha_diversity.csv", index_col=0)
    def H_(v): v = np.array([x for x in v if x > 0]); p = v / v.sum(); return -(p * np.log(p)).sum()
    dd = max(abs(H_(cnts[s].values()) - v) for s, v in al.shannon.items())
    chk("B. Statistics", "DS02 prokaryotic Shannon: independent recomputation (max abs difference)", "< 1e-9", f"{dd:.2e}", dd < 1e-9)
    # 7. PERMANOVA pseudo-F by the classical sum-of-squares formula
    D = np.load(DER / "ds02_D_arg_bc.npy"); gg = pd.read_csv(DER / "ds02_meta.csv").City.values; N = len(gg); a = len(set(gg))
    SST = (D[np.triu_indices(N, 1)] ** 2).sum() / N
    SSW = sum((D[np.ix_(g == c, g == c)][np.triu_indices((g == c).sum(), 1)] ** 2).sum() / (g == c).sum() for c in set(gg) for g in [gg])
    F = ((SST - SSW) / (a - 1)) / (SSW / (N - a))
    c01 = pd.read_csv(TAB / "C01_ds02_permanova_permdisp.csv"); rep = c01[(c01.distance == "Bray-Curtis (relative RPKM)") & (c01.term == "City") & (c01.covariates == "(none)")].iloc[0]
    chk("B. Statistics", "PERMANOVA pseudo-F (DS02 resistome BC ~ City) by classical formula vs C01", round(rep.pseudo_F, 4), round(F, 4), close(F, rep.pseudo_F, 1e-6))
    # 8. mixed model refit with an alternative optimiser
    import statsmodels.formula.api as smf
    urb = us[us.Stratum == "Urban environmental"]
    m = smf.mixedlm("log10_total_FPKM ~ C(Sample_Type, Treatment('Municipal wastewater'))", urb, groups=urb.Project_Accession).fit(reml=True, method="powell")
    c08 = pd.read_csv(TAB / "C08_rf_mixed_arg_load.csv").set_index("term")
    t = [x for x in m.fe_params.index if "Surface water" in x][0]
    chk("B. Statistics", "Mixed model (surface water vs wastewater) refit with Powell optimiser vs C08", round(c08.loc[t, "estimate"], 4), round(m.fe_params[t], 4), abs(m.fe_params[t] - c08.loc[t, "estimate"]) < 1e-3)
    # 9. Moran's I recomputed with an explicit double loop
    pil = us[(us.Project_Accession == "PRJEB27054") & us.lat.notna()]
    lat, lon, xv = np.radians(pil.lat.values), np.radians(pil.lon.values), pil.log10_total_FPKM.values
    n = len(xv); Dg = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            Dg[i, j] = 2 * math.asin(math.sqrt(min(1, math.sin((lat[i] - lat[j]) / 2) ** 2 + math.cos(lat[i]) * math.cos(lat[j]) * math.sin((lon[i] - lon[j]) / 2) ** 2)))
    zx = xv - xv.mean(); num = 0.0
    for i in range(n):
        d_ = Dg[i].copy(); d_[i] = np.inf; nb = np.argsort(d_, kind="stable")[:5]  # AUDIT FIX AUD-02
        num += sum(zx[i] * zx[j] / 5 for j in nb)
    I = (n / n) * num / (zx ** 2).sum()
    s06 = pd.read_csv(TAB / "S06_spatial_autocorrelation.csv").iloc[0]
    chk("B. Statistics", "Moran's I (pilot normalized ARG abundance, k=5) by explicit loop vs S06", round(s06.morans_I, 4), round(I, 4), abs(I - s06.morans_I) < 1e-6)
    # 10. BH adjustment re-implemented
    reg = pd.read_csv(TAB / "R02_statistical_tests_register.csv"); fam = reg[reg.family == "DS02 phylum-class correlations (BH within family)"]
    p = fam.p_value.values; o = np.argsort(p); m_ = len(p); q = np.empty(m_); prev = 1.0
    for rank in range(m_ - 1, -1, -1):
        prev = min(prev, p[o[rank]] * m_ / (rank + 1)); q[o[rank]] = prev
    chk("B. Statistics", "Benjamini-Hochberg q-values re-implemented (phylum-class family) vs register", "max diff < 1e-12", f"{np.nanmax(np.abs(q - fam.p_adj_BH_within_family.values)):.1e}", np.nanmax(np.abs(q - fam.p_adj_BH_within_family.values)) < 1e-12)
    # 11. model metrics plausibility
    m2 = pd.read_csv(TAB / "M02_modelA_performance.csv")
    b = m2[m2.model == "Baseline (prior)"].roc_auc_mean
    chk("B. Statistics", "Model A prior baseline has fold-mean ROC-AUC = 0.5 (sanity)", 0.5, b.round(3).tolist(), bool(np.allclose(b, 0.5)))
    chk("B. Statistics", "All reported ROC-AUC values lie in [0, 1]", "[0,1]", f"{m2.roc_auc_mean.min():.3f}-{m2.roc_auc_mean.max():.3f}", bool(m2.roc_auc_mean.between(0, 1).all()))

def integrity_checks():
    inv = pd.read_csv(TAB / "QC01_file_inventory.csv").set_index("relative_path")
    changed = []
    for rel, r in inv.iterrows():
        p = ROOT / rel
        if not p.exists() or p.stat().st_size != r.size_bytes: changed.append(rel); continue
        if rel.startswith("AMR_Chapter_Reproducibility_Bundle") and "harmonized" in rel and md5(p) != r.md5: changed.append(rel)
    chk("C. Data integrity", "Original files unchanged since audit (size for all; MD5 for harmonized tables)", "0 changed", f"{len(changed)} changed of {len(inv)}", not changed)
    sm_ids = {r["Sample_ID"] for r in read_csv_rows(HARM / "sample_master.csv", ["Sample_ID"])}
    meta = pd.read_csv(DER / "metadata_set.csv", usecols=["Sample_ID", "Inclusion_Status"])
    chk("C. Data integrity", "Metadata set: unique sample IDs, all in sample master, no excluded samples", "8,893 unique", f"{meta.Sample_ID.nunique()} unique; {(~meta.Sample_ID.isin(sm_ids)).sum()} unknown; {(meta.Inclusion_Status=='Excluded').sum()} excluded",
      meta.Sample_ID.is_unique and meta.Sample_ID.isin(sm_ids).all() and (meta.Inclusion_Status != "Excluded").all())
    s2 = pd.read_csv(DER / "ds02_arg_summary.csv", index_col=0); ph = pd.read_csv(DER / "ds02_prokaryote_phylum_counts.csv", index_col=0)
    chk("C. Data integrity", "DS02: 447 unique samples with both ARG and taxonomy data; identical sample sets", 447, f"{len(s2)} / {len(ph)}", len(s2) == 447 and set(s2.index) == set(ph.index) and s2.index.is_unique)
    us = pd.read_csv(DER / "rf_unit_summary.csv", index_col=0); rm = pd.read_csv(DER / "rf_run_to_unit.csv")
    chk("C. Data integrity", "ResFinder: run lacking ARG rows (ERR1725953) absent from AMR units (missing not zero)", "absent",
      "absent" if not rm[rm.Run_ID == "ERR1725953"].unit_id.isin(us.index).any() else "present", not rm[rm.Run_ID == "ERR1725953"].unit_id.isin(us.index).any())
    chk("C. Data integrity", "ResFinder primary units exclude PRJEB13831 (probable duplicate physical samples)", 0, int((us.Project_Accession == "PRJEB13831").sum()), (us.Project_Accession == "PRJEB13831").sum() == 0)
    chk("C. Data integrity", "Pipelines kept separate: no DS02 sample among ResFinder units and no ResFinder unit among DS02 samples", "disjoint",
      f"{len(set(us.Sample_ID_first) & set(s2.index))} overlap", not (set(us.Sample_ID_first) & set(s2.index)))
    pil = us[us.Project_Accession == "PRJEB27054"]; acc = rm[rm.Project_Accession == "PRJEB27054"]
    chk("C. Data integrity", "PRJEB27054 accessions merged to biological samples (units < accessions)", f"{acc.Sample_ID.nunique()} accessions", f"{len(pil)} units", len(pil) < acc.Sample_ID.nunique())
    chk("C. Data integrity", "No zero-imputation: minimum DS02 total RPKM and ResFinder FPKM > 0 (only observed rows summed)", "> 0", f"{s2.total_RPKM.min():.3f} / {us.total_FPKM.min():.3f}", s2.total_RPKM.min() > 0 and us.total_FPKM.min() > 0)

def excel_checks(path, label):
    import openpyxl
    from build_workbook import SHEETS
    wb = openpyxl.load_workbook(path)
    chk(label, "Workbook opens without error", "opens", "opened", True)
    chk(label, "All 27 required worksheets exist", len(SHEETS), len([s for s in SHEETS if s in wb.sheetnames]), all(s in wb.sheetnames for s in SHEETS))
    pos = json.load(open(LOG / "workbook_table_positions.json"))
    bad, rng = [], random.Random(1)
    for f, (sh, r0, c0, r1, ce) in pos.items():
        src = (VER / f) if f.startswith("V0") else (TAB / f)
        df = pd.read_csv(src, low_memory=False); ws = wb[sh]
        if (r1 - r0) != max(len(df), 1): bad.append(f"{f}: rows {r1-r0} vs {len(df)}"); continue
        if ce - c0 + 1 != df.shape[1]: bad.append(f"{f}: cols"); continue
        for i in ([0, len(df) - 1] + [rng.randrange(len(df)) for _ in range(3)] if len(df) else []):
            for j in range(df.shape[1]):
                v = df.iat[i, j]; x = ws.cell(r0 + 1 + i, c0 + j).value
                if (v is None or (isinstance(v, float) and np.isnan(v))) and x is None: continue
                if isinstance(v, (int, float, np.number)) and not isinstance(v, bool):
                    try: okv = x is not None and close(float(x), float(v), 1e-9, 1e-12)
                    except (TypeError, ValueError): okv = False
                    if not okv: bad.append(f"{f}[{i},{j}] {v} vs {str(x)[:40]}")
                elif str(x).lstrip("'") != str(v)[:32000]: bad.append(f"{f}[{i},{j}]")
    chk(label, "Every table: row/column counts equal source CSV and spot-checked cells match (first, last, 3 random rows)", f"{len(pos)} tables", f"{len(bad)} mismatches" + (": " + "; ".join(bad[:5]) if bad else ""), not bad)
    nimg = sum(len(ws._images) for ws in wb.worksheets); nfig = len(load_catalog("figures"))
    chk(label, "Embedded images: every figure in the gallery plus topical sheets", f">= {nfig}", nimg, nimg >= nfig)
    gal = len(wb["Figure_Gallery"]._images); chk(label, "Figure_Gallery contains every figure", nfig, gal, gal == nfig)
    nch = sum(len(ws._charts) for ws in wb.worksheets); chk(label, "Native Excel charts present", ">= 3", nch, nch >= 3)
    chk(label, "No external links", 0, len(wb._external_links), len(wb._external_links) == 0)
    nform = sum(1 for ws in wb.worksheets for row in ws.iter_rows() for c in row if isinstance(c.value, str) and c.value.startswith("="))
    chk(label, "No formulas or broken references (all values static)", 0, nform, nform == 0)
    ntab = sum(len(ws.tables) for ws in wb.worksheets); chk(label, "Excel tables with filters created", f"{len([p for p in pos])}", ntab, ntab >= len(pos) - 5)

def science_checks():
    F = pd.read_csv(TAB / "R03_findings_register.csv")
    pat = re.compile(r"\b(causes?d?|caused by|proves?|demonstrates? transmission|leads? to|results? in increased|risk of infection)\b", re.I)
    flagged = [r.id for _, r in F.iterrows() if pat.search(str(r.finding))]
    chk("E. Scientific", "Findings contain no causal or transmission language", 0, flagged or 0, not flagged)
    noref = F[F.table_or_figure.isna() | (F.sample_size.isna())].id.tolist()
    chk("E. Scientific", "Every finding states a sample size and a table/figure reference", 0, noref or 0, not noref)
    ids = set(pd.read_csv(TAB / "R07_table_catalog.csv").file.str.replace(".csv", "", regex=False).str.split("_").str[0])
    refs = set(re.findall(r"\b([A-Z]{1,2}\d{2}[a-z]?)\b", " ".join(F.table_or_figure.astype(str))))
    missing = sorted(r for r in refs if r not in ids and not r.startswith(("F")))
    chk("E. Scientific", "Table IDs cited in findings exist in the table catalogue", "all exist", missing or "all exist", not missing)
    cap = load_catalog("figures")
    noneg = [int(r.number) for _, r in cap.iterrows() if r.number in (9, 11, 13, 14) and not re.search(r"not (establish|carriage|comparable|resistance|phenotyp)|does not|not be read", r.caption)]
    chk("E. Scientific", "AMR figure captions state that detection is not phenotype/clinical risk (Figs 9, 11, 13-14)", "all", noneg or "all", not noneg)
    s07 = pd.read_csv(TAB / "S07_spatial_feasibility.csv"); m07 = pd.read_csv(TAB / "M07_predictive_feasibility.csv")
    chk("E. Scientific", "No interpolated or hotspot maps produced; flagged as not performed", "Not performed", s07[s07.method.str.contains("hotspot|Interpolated")].decision.str[:13].tolist(), s07[s07.method.str.contains("hotspot|Interpolated")].decision.str.startswith("Not performed").all())
    chk("E. Scientific", "Disease/outbreak/food-business prediction declared not feasible", "Not feasible", m07[m07.target.str.contains("Disease|Food-safety")].decision.tolist(), m07[m07.target.str.contains("Disease|Food-safety")].decision.eq("Not feasible").all())

def report(fn):
    df = pd.DataFrame(checks); df.to_csv(VER / fn, index=False)
    lines = [f"# Verification report ({fn})", f"Run {dt.datetime.now().isoformat(timespec='seconds')}", "", f"**{(df.status=='PASS').sum()} of {len(df)} checks passed.**", ""]
    for cat, g in df.groupby("category"):
        lines += [f"## {cat}", "| Check | Expected | Observed | Status |", "|---|---|---|---|"] + [f"| {r.check} | {r.expected} | {r.observed} | {r.status} |" for _, r in g.iterrows()] + [""]
    (VER / fn.replace(".csv", ".md")).write_text("\n".join(lines), encoding="utf-8")
    print(df.status.value_counts().to_dict()); print(df[df.status == "FAIL"].to_string())

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    wbp = OUT / "AMR_Analysis_Results_Verified.xlsx"
    if mode == "full":
        code_checks(); stat_checks(); integrity_checks(); science_checks(); excel_checks(wbp, "D. Excel workbook (pre-final build)")
        report("V01_verification_checks.csv")
    else:
        excel_checks(Path(sys.argv[2]) if len(sys.argv) > 2 else wbp, "D. Excel workbook (final build)"); report("V02_final_workbook_checks.csv")
