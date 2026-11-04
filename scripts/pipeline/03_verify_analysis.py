"""Step 03 - Independent re-computation of headline supplementary results from the cleaned datasets.

Objective: confirm that the values reported in the supplementary tables can be reproduced from the
cleaned, publicly released data, using independently written code (not the original notebooks).

Analyses (methods chosen to match the verified analysis design; equations in
08_Methodology_and_Reports/Statistical_Analysis_Methodology.docx):
 A. Descriptive statistics of the evidence base (counts, medians)                -> Tables S1, S2, S92
 B. Stream 1 (DS02 Indian open drains, CARD/RGI):
      total ARG abundance  RPKM_s = sum_g reads_gs / (L_g/1e3) / (N_s/1e6)          -> Table S26
      site-level Kruskal-Wallis on 19 site medians, whole-site permutation p        -> Table S27
      linear mixed model log10 RPKM ~ city + season + year + (1 | site), REML       -> Table S28
 C. Stream 2 (ResFinder resource): median total FPKM by urban matrix              -> Table S37
 D. Stream 3 (uniform reprocessing, Dhaka): ranges, medians and exact permutation tests -> Tables S61, S67, S68
Outputs: 03_Data_Analysis/{Descriptive_Statistics, Statistical_Analysis, Comparative_Analysis, Analytical_Results}
"""
import itertools

import numpy as np
import pandas as pd
from scipy import stats

import shutil

from amr_paths import CLEAN_CSV, ANALYSIS, verified_table, rf_unit_summary, REPO_LAYOUT, log

SEED = 20261104
MISSING = {"Not reported", "Not yet verified", "Not available", "Not applicable"}  # documented missing-value codes
N_PERM = 9999
checks = []
SUBS = {k: (v if REPO_LAYOUT else k) for k, v in (("Descriptive_Statistics", "descriptive_statistics"),
        ("Statistical_Analysis", "statistical_analysis"), ("Comparative_Analysis", "comparative_analysis"),
        ("Analytical_Results", "results"))}


def chk(area, statistic, table, reported, recomputed, tol_abs=None, tol_rel=None):
    rep, rec = float(reported), float(recomputed)
    diff = abs(rep - rec)
    ok = (tol_abs is not None and diff <= tol_abs) or (tol_rel is not None and diff <= tol_rel * max(abs(rep), 1e-12)) \
        or (tol_abs is None and tol_rel is None and diff == 0)
    checks.append(dict(Area=area, Statistic=statistic, Supplement_Table=table, Reported=rep, Recomputed=rec,
                       Abs_Difference=diff, Tolerance=(f"abs {tol_abs}" if tol_abs is not None else
                                                       f"rel {tol_rel}" if tol_rel is not None else "exact"),
                       Match="YES" if ok else "NO"))


def season(m):
    return {12: "Winter (Dec-Feb)", 1: "Winter (Dec-Feb)", 2: "Winter (Dec-Feb)", 3: "Pre-monsoon (Mar-May)",
            4: "Pre-monsoon (Mar-May)", 5: "Pre-monsoon (Mar-May)", 10: "Post-monsoon (Oct-Nov)",
            11: "Post-monsoon (Oct-Nov)"}.get(m, "Monsoon (Jun-Sep)")


def site_kw(d, value, site, group, rng):
    s = d.groupby(site).agg(v=(value, "median"), g=(group, "first"))
    groups = [s.v[s.g == u].values for u in sorted(s.g.unique())]
    H, p = stats.kruskal(*groups)
    k, n = len(groups), len(s)
    labs, vals = s.g.values, s.v.values
    cnt = sum(stats.kruskal(*[vals[lp == u] for u in np.unique(lp)])[0] >= H - 1e-12
              for lp in (rng.permutation(labs) for _ in range(N_PERM)))
    return dict(H=H, df=k - 1, p_asymptotic=p, p_perm=(cnt + 1) / (N_PERM + 1), epsilon2=(H - k + 1) / (n - k),
                n_sites=n)


def exact_perm(x, y, stat=lambda a, b: np.median(a) - np.median(b)):
    """Two-sided exact permutation p-value over all relabellings."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    pooled = np.concatenate([x, y])
    obs = stat(x, y)
    cnt = tot = 0
    for idx in itertools.combinations(range(len(pooled)), len(x)):
        m = np.zeros(len(pooled), bool); m[list(idx)] = True
        cnt += abs(stat(pooled[m], pooled[~m])) >= abs(obs) - 1e-12
        tot += 1
    return obs, cnt / tot


def cliffs_delta(x, y):
    x, y = np.asarray(x), np.asarray(y)
    return (np.sign(x[:, None] - y[None, :])).mean()


def main():
    rng = np.random.default_rng(SEED)
    for sub in SUBS.values():
        (ANALYSIS / sub).mkdir(parents=True, exist_ok=True)
    sm = pd.read_csv(CLEAN_CSV / "DS_A_sample_master.csv", low_memory=False)
    inc = sm[sm.inclusion_status == "Included"].copy()

    # ---------------- A. descriptive
    chk("A. Evidence base", "Included samples", "Table S1", 8053, len(inc))
    chk("A. Evidence base", "Datasets with included samples", "Table S1", 112, inc.dataset_id.nunique())
    chk("A. Evidence base", "Countries reported", "Table S1", 106, inc.loc[~inc.country.isin(MISSING), "country"].nunique())
    cities = inc.loc[~inc.city.isin(MISSING) & ~inc.country.isin(MISSING), ["city", "country"]].drop_duplicates()
    chk("A. Evidence base", "Source-labelled City|Country pairs", "Table S1", 274, len(cities))
    chk("A. Evidence base", "Samples with coordinates", "Table S10", 6134, inc.latitude.notna().sum())
    chk("A. Evidence base", "Samples with a processed AMR table", "Table S8", 1072,
        inc.amr_data_available.str.startswith("Yes").sum())
    inc["gbp"] = inc.base_count_total / 1e9
    d02 = pd.read_csv(verified_table("D02_by_sample_type.csv")).set_index("Sample_Type")
    by = inc.groupby("sample_type").agg(samples=("sample_id", "size"), datasets=("dataset_id", "nunique"),
                                        with_coordinates=("latitude", lambda s: s.notna().sum()),
                                        median_Gbp=("gbp", "median"))
    by["pct_of_included"] = (100 * by.samples / len(inc)).round(1)
    by.to_csv(ANALYSIS / SUBS["Descriptive_Statistics"] / "Included_samples_by_sample_type_recomputed.csv")
    for st, r in by.iterrows():
        if st in d02.index:
            chk("A. Evidence base", f"Samples: {st}", "Table S2", d02.loc[st, "samples"], r.samples)
            chk("A. Evidence base", f"Datasets: {st}", "Table S2", d02.loc[st, "datasets"], r.datasets)
            chk("A. Evidence base", f"Median depth (Gbp): {st}", "Table S2", d02.loc[st, "median_Gbp"], r.median_Gbp,
                tol_rel=0.005)
    sd01 = pd.read_csv(verified_table("SD01_depth_threshold_denominators.csv")).set_index("Sample_Type").loc["ALL"]
    for lab, thr in ((">= 10 Mbp", 1e7), (">= 100 Mbp", 1e8), (">= 500 Mbp", 5e8), (">= 1 Gbp", 1e9)):
        chk("A. Evidence base", f"Samples with depth {lab}", "Table S92", sd01[lab], (inc.base_count_total >= thr).sum())
    desc = inc.groupby(["sample_type"]).gbp.describe()
    desc.to_csv(ANALYSIS / SUBS["Descriptive_Statistics"] / "Sequencing_depth_Gbp_by_sample_type.csv")

    # ---------------- B. Stream 1: DS02 total RPKM from cleaned gene table
    amr = pd.read_csv(CLEAN_CSV / "DS_C_amr_genes_long.csv.gz", low_memory=False)
    amr = amr[amr.dataset_id == "DS02"]
    runs = pd.read_csv(CLEAN_CSV / "DS_B_run_master.csv", low_memory=False)[["run_id", "read_count"]]
    amr["rpk"] = amr.detection_value / (amr.reference_length_bp / 1e3)
    tot = amr.groupby(["sample_id", "run_id"]).rpk.sum().reset_index().merge(runs, on="run_id")
    tot["total_RPKM"] = tot.rpk / (2 * tot.read_count / 1e6)   # paired reads: N = 2 x spots
    ds = tot.merge(sm[["sample_id", "city", "site_id", "collection_date"]], on="sample_id")
    ds["month"] = pd.to_numeric(ds.collection_date.str.slice(5, 7), errors="coerce")
    ds["year"] = ds.collection_date.str.slice(0, 4)
    ds["season"] = ds.month.map(season)
    ds["log10_total_RPKM"] = np.log10(ds.total_RPKM)
    ds.to_csv(ANALYSIS / SUBS["Statistical_Analysis"] / "DS02_total_ARG_RPKM_per_sample_recomputed.csv", index=False)
    a04 = pd.read_csv(verified_table("A04_ds02_arg_load_by_city.csv")).set_index("City")
    med = ds.groupby("city").total_RPKM.agg(["size", "median"])
    med.to_csv(ANALYSIS / SUBS["Descriptive_Statistics"] / "DS02_total_RPKM_by_city_recomputed.csv")
    for c, r in med.iterrows():
        chk("B. Stream 1 (DS02)", f"n samples: {c}", "Table S26", a04.loc[c, "n"], r["size"])
        chk("B. Stream 1 (DS02)", f"Median total RPKM: {c}", "Table S26", a04.loc[c, "median_total_RPKM"], r["median"],
            tol_abs=0.01)
    k = site_kw(ds, "log10_total_RPKM", "site_id", "city", rng)
    a05 = pd.read_csv(verified_table("A05s_ds02_arg_site_aware_tests.csv")).iloc[0]
    chk("B. Stream 1 (DS02)", "Site-level Kruskal-Wallis H (19 site medians)", "Table S27", a05.H, k["H"], tol_abs=0.01)
    chk("B. Stream 1 (DS02)", "Site-level epsilon-squared", "Table S27", a05.epsilon2, k["epsilon2"], tol_abs=0.001)
    chk("B. Stream 1 (DS02)", "Whole-site permutation p (9,999 perm.; Monte Carlo)", "Table S27", a05.p_perm, k["p_perm"],
        tol_abs=0.003)
    pd.DataFrame([k]).to_csv(ANALYSIS / SUBS["Statistical_Analysis"] / "DS02_site_level_kruskal_recomputed.csv", index=False)
    import statsmodels.formula.api as smf
    f = "log10_total_RPKM ~ C(city, Treatment('Delhi')) + C(season, Treatment('Winter (Dec-Feb)')) + C(year)"
    fit = None
    for method in ("lbfgs", "powell", "nm"):
        try:
            fit = smf.mixedlm(f, ds, groups=ds.site_id).fit(reml=True, method=method)
            break
        except Exception as exc:  # singular Hessian etc.
            log(f"MixedLM ({method}) failed: {exc}")
    a07 = pd.read_csv(verified_table("A07s_ds02_arg_mixed_site.csv")).set_index("term")
    if fit is not None:
        fe = pd.DataFrame(dict(estimate=fit.fe_params, ci_low=fit.conf_int().loc[fit.fe_params.index, 0],
                               ci_high=fit.conf_int().loc[fit.fe_params.index, 1]))
        fe["fold_change"] = 10 ** fe.estimate
        fe.to_csv(ANALYSIS / SUBS["Statistical_Analysis"] / "DS02_mixed_model_recomputed.csv")
        for city in ("Chennai", "Kolkata", "Mumbai"):
            term = f"C(city, Treatment('Delhi'))[T.{city}]"
            rep = a07.loc[f"C(City, Treatment('Delhi'))[T.{city}]", "fold_change"]
            chk("B. Stream 1 (DS02)", f"Mixed-model fold change vs Delhi: {city}", "Table S28", rep,
                fe.loc[term, "fold_change"], tol_abs=0.005)

    # ---------------- C. Stream 2: ResFinder resource unit medians (unit table from the corrected analysis build)
    src = rf_unit_summary()
    inputs = ANALYSIS / ("comparative_analysis" if REPO_LAYOUT else "Comparative_Analysis") / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    if src.parent != inputs:
        shutil.copyfile(src, inputs / src.name)   # ship the unit table so the check runs from the package alone
    unit = pd.read_csv(src)
    unit = unit[unit.primary & unit.has_ARG]
    a12 = pd.read_csv(verified_table("A12_rf_arg_load_by_sample_type.csv")).set_index("Sample_Type")
    um = unit.groupby("Sample_Type").agg(units=("unit_id", "size"), studies=("Project_Accession", "nunique"),
                                         median_total_FPKM=("total_FPKM", "median"),
                                         median_gene_richness=("gene_richness", "median"))
    um.to_csv(ANALYSIS / SUBS["Comparative_Analysis"] / "ResFinder_matrix_summary_recomputed.csv")
    for st, r in um.iterrows():
        if st in a12.index:
            chk("C. Stream 2 (ResFinder)", f"Units: {st}", "Table S37", a12.loc[st, "units"], r.units)
            chk("C. Stream 2 (ResFinder)", f"Median total FPKM: {st}", "Table S37", a12.loc[st, "median_total_FPKM"],
                r.median_total_FPKM, tol_rel=0.001)
            chk("C. Stream 2 (ResFinder)", f"Median gene richness: {st}", "Table S37", a12.loc[st, "median_gene_richness"],
                r.median_gene_richness)

    # ---------------- D. Stream 3: Dhaka uniform reprocessing
    cv = pd.read_csv(CLEAN_CSV / "DS_K_reprocessed_comparative_variables.csv")
    am = pd.read_csv(CLEAN_CSV / "DS_J_reprocessed_amr_summary.csv")
    dh = cv[cv.city_group == "Dhaka"].merge(am, on="run_accession")
    chk("D. Stream 3 (Dhaka)", "Dhaka metagenomes reprocessed", "Table S61", 16, len(dh))
    chk("D. Stream 3 (Dhaka)", "Median high-confidence ARG genes", "Table S61", 35.5, dh.genes_high_confidence.median())
    chk("D. Stream 3 (Dhaka)", "Minimum total FPKM-16S", "Table S61", 13248, dh.total_fpkm_16s.min(), tol_abs=1)
    chk("D. Stream 3 (Dhaka)", "Maximum total FPKM-16S", "Table S61", 1286405, dh.total_fpkm_16s.max(), tol_abs=1)
    ds84 = dh[dh.dataset_id == "DS84"]
    hosp = ds84[ds84.sample_type == "Hospital wastewater"].log10_total_fpkm_16s
    stp = ds84[ds84.sample_type == "Municipal wastewater"].log10_total_fpkm_16s
    obs, p = exact_perm(hosp, stp)
    chk("D. Stream 3 (Dhaka)", "DS84 hospital vs STP inlet: median difference log10 FPKM-16S", "Table S68", 0.29, obs,
        tol_abs=0.005)
    chk("D. Stream 3 (Dhaka)", "DS84 hospital vs STP inlet: exact permutation p", "Table S68", 0.057, p, tol_abs=0.0005)
    chk("D. Stream 3 (Dhaka)", "DS84 hospital vs STP inlet: Cliff's delta", "Table S68", 1.0, cliffs_delta(hosp, stp))
    kar = cv[cv.city_group == "Karachi"].merge(am, on="run_accession")
    dh_gss = dh[dh.dataset_id == "DS01"]
    obs2, p2 = exact_perm(dh_gss.genes_high_confidence, kar.genes_high_confidence)
    chk("D. Stream 3 (Dhaka)", "Dhaka vs Karachi (GSS protocol): median difference, high-confidence genes", "Table S68",
        -23.5, obs2)
    chk("D. Stream 3 (Dhaka)", "Dhaka vs Karachi: exact permutation p", "Table S68", 0.012, p2, tol_abs=0.0005)
    comp = pd.DataFrame([dict(comparison="DS84 hospital vs STP inlet (log10 FPKM-16S)", n1=len(hosp), n2=len(stp),
                              median_difference=obs, exact_p=p, cliffs_delta=cliffs_delta(hosp, stp)),
                         dict(comparison="Dhaka DS01 vs Karachi (high-confidence genes)", n1=len(dh_gss), n2=len(kar),
                              median_difference=obs2, exact_p=p2,
                              cliffs_delta=cliffs_delta(dh_gss.genes_high_confidence, kar.genes_high_confidence))])
    comp.to_csv(ANALYSIS / SUBS["Comparative_Analysis"] / "Dhaka_exact_permutation_tests_recomputed.csv", index=False)

    res = pd.DataFrame(checks)
    out = ANALYSIS / SUBS["Analytical_Results"] / "Independent_Verification_of_Supplementary_Results.xlsx"
    with pd.ExcelWriter(out, engine="openpyxl") as xw:
        res.groupby("Area").Match.value_counts().unstack(fill_value=0).reset_index().to_excel(
            xw, sheet_name="Summary", index=False)
        res.to_excel(xw, sheet_name="Checks", index=False)
    res.to_csv(out.with_suffix(".csv"), index=False)
    log(f"Verification: {len(res)} checks, {(res.Match == 'YES').sum()} match, {(res.Match == 'NO').sum()} mismatch")
    print(res[res.Match == "NO"].to_string())


if __name__ == "__main__":
    main()
