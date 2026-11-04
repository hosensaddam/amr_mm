"""Shared loaders and calculations for the uniformly reprocessed raw-read set (Dhaka + comparators).
Inputs live in Analysis_Results/derived/reprocessed/ (tidy tables collected from the per-run outputs of the
uniform pipeline: fastp 1.4.0 -> Kraken 2 2.17.1 / SILVA 138 16S + Bracken -> KMA 1.7.1 / ResFinder DB).
Nothing here is imputed: runs, genes or taxa that were not observed are absent, not zero-filled, except
inside a per-run gene/taxon matrix, where an absent feature means "not detected in that run"."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from amr_common import *  # noqa: F401,F403  (paths, palettes, statistics helpers)

REP = DER / "reprocessed"
DHAKA_DATASETS = ["DS01", "DS03", "DS84"]
GROUP_ORDER = ["Dhaka", "Karachi", "Kathmandu", "Kochi", "Colombo", "Delhi"]
CITY_COL = {"Dhaka": "#D55E00", "Karachi": "#0072B2", "Kathmandu": "#009E73", "Kochi": "#CC79A7", "Colombo": "#E69F00", "Delhi": "#56B4E9",
            "Calibration": "#7F7F7F"}
DATASET_COL = {"DS01": "#0072B2", "DS03": "#009E73", "DS84": "#D55E00", "DS02": "#56B4E9"}
TYPE_COL = {"Municipal wastewater": "#0072B2", "Hospital wastewater": "#D55E00", "Community wastewater": "#56B4E9", "Surface water": "#009E73"}
HC_IDENTITY, HC_COVERAGE = 90.0, 60.0   # ResFinder default thresholds (identity 90 %, coverage 60 %)

_RF47 = {"Beta-lactam": "Beta-lactam", "Aminoglycoside": "Aminoglycoside", "Folate pathway antagonist": "Folate pathway antagonist",
         "Tetracycline": "Tetracycline", "Quinolone": "Fluoroquinolone/quinolone", "Polymyxin": "Peptide/polymyxin", "Macrolide": "MLS",
         "Lincosamide": "MLS", "Streptogramin A": "MLS", "Streptogramin B": "MLS", "Pleuromutilin": "MLS", "Amphenicol": "Phenicol",
         "Fosfomycin": "Fosfomycin", "Glycopeptide": "Glycopeptide", "Oxazolidinone": "Oxazolidinone", "Rifamycin": "Rifamycin",
         "Nitroimidazole": "Other", "Steroid antibacterial": "Other", "Pseudomonic acid": "Other", "Ionophores": "Other"}
def rf47_class_group(c) -> str:
    """ResFinder DB (resfinder-4.7.2) phenotype-table class string -> harmonised display group (as amr_common.CLASS_ORDER)."""
    if pd.isna(c) or str(c).strip() == "": return "Unclassified"
    g = {_RF47.get(t.strip(), "Other") for t in str(c).split(",")}
    return g.pop() if len(g) == 1 else "Multi-class"
def rf47_mechanism(m) -> str:
    if pd.isna(m) or str(m).strip() == "": return "Not stated"
    m = str(m).strip(); ml = m.lower()
    if "efflux" in ml: return "Efflux"
    if ml.startswith("target protection"): return "Target protection"
    if ml.startswith("target modification"): return "Target modification"
    if ml.startswith("class a"): return "Enzymatic inactivation"
    return m[0].upper() + m[1:]

def load_run_meta() -> pd.DataFrame:
    """One row per processed run: ENA manifest + QC + harmonized sample metadata (updated master)."""
    qc = pd.read_csv(REP / "run_qc.csv")
    sm = load_csv("sample_master.csv")
    keep = ["Sample_ID", "Dataset_ID", "Project_Accession", "City", "Country", "Locality_As_Reported", "Latitude", "Longitude",
            "Latitude_Alternative", "Longitude_Alternative", "Coordinate_Precision", "Coordinate_Provenance_Note", "Sample_Type",
            "Sample_Type_Original_Harmonized", "Sample_Type_Correction", "Sample_Subtype", "Collection_Date", "Site_ID", "Biological_Sample_ID",
            "Inclusion_Status", "Sequencing_Platform"]
    m = qc.merge(sm[keep], left_on="sample_accession", right_on="Sample_ID", how="left")
    m["City_group"] = np.where(m.analysis_group == "Dhaka", "Dhaka", m.City)
    m.loc[m.analysis_group.str.startswith("Calibration"), "City_group"] = "Calibration"
    m["Year"] = pd.to_numeric(m.Collection_Date.astype(str).str[:4], errors="coerce")
    m["Period"] = np.where(m.Dataset_ID == "DS01", "2017-2018", np.where(m.Dataset_ID == "DS84", "2024", np.where(m.Dataset_ID == "DS03", "2025", m.Collection_Date.astype(str).str[:4])))
    # DS84 sample names (HWW/UWW/Control) and DS02 sites are carried as short labels
    m["Label"] = m.Sample_ID
    return m.set_index("run_accession")

def taxonomy_matrix(rank: str, runs, prokaryotes_only=True, source="kraken") -> pd.DataFrame:
    """Run x taxon fragment counts at one rank. Kraken 2 clade fragments for D/P/C/O/F; Bracken re-estimated
    fragments for genus (source='bracken'). Restricted to Bacteria + Archaea 16S fragments by default."""
    if rank == "G" and source == "bracken":
        b = pd.read_csv(REP / "bracken_genus_long.csv")
        k = pd.read_csv(REP / "taxonomy_kraken_long.csv")
        lin = k[k["rank"] == "G"].drop_duplicates("taxid").set_index("taxid").lineage
        b["domain"] = b.taxid.map(lin).fillna("").str.extract(r"D__([^;]+)")[0]
        if prokaryotes_only: b = b[b.domain.isin(["Bacteria", "Archaea"])]
        M = b.pivot_table(index="run_accession", columns="genus", values="bracken_est_fragments", aggfunc="sum", fill_value=0)
    else:
        k = pd.read_csv(REP / "taxonomy_kraken_long.csv")
        k = k[k["rank"] == rank].copy()
        k["domain"] = k.lineage.str.extract(r"D__([^;]+)")[0]
        if prokaryotes_only: k = k[k.domain.isin(["Bacteria", "Archaea"])]
        M = k.pivot_table(index="run_accession", columns="name", values="clade_fragments", aggfunc="sum", fill_value=0)
    return M.reindex(list(runs)).fillna(0)

def lineage_of(rank: str) -> pd.Series:
    k = pd.read_csv(REP / "taxonomy_kraken_long.csv")
    return k[k["rank"] == rank].drop_duplicates("name").set_index("name").lineage

def resfinder_table(meta: pd.DataFrame) -> pd.DataFrame:
    """Template-level ResFinder/KMA table with gene names, class groups, detection tiers and normalized abundance.
    Abundance uses fragmentCountAln (aligned fragments; the definition used by the Martiny et al. resource).
    FPKM_16S = fragments / (template kb) / (prokaryotic 16S fragments / 1e6)  [per million bacterial+archaeal 16S fragments]
    RPKM_depth = fragments / (template kb) / (QC read pairs / 1e6)            [per million quality-filtered read pairs]"""
    rf = pd.read_csv(REP / "resfinder_long.csv")
    rf = rf[rf.fragmentCountAln > 0].copy()
    L = pd.read_csv(REP / "resfinder_template_lengths.csv").set_index("template").length
    rf["length_bp"] = rf.template.map(L)
    assert rf.length_bp.notna().all(), "template length missing"
    rf["gene"] = rf.template.map(rf_gene_name)
    rf["allele"] = rf.template.str.split("_").str[0]
    rf["Class"] = rf.Class.fillna("Unclassified")
    rf["Class_Group"] = rf.Class.map(rf47_class_group)
    rf["Mechanism"] = rf["Mechanism of resistance"].map(rf47_mechanism)
    rf["Clinical"] = [clinical_label(a, "", "RF") for a in rf.template.str.split("_").str[0]]   # full allele name (fix 31 Oct 2026: collapsed gene names never matched OXA-48-like or MCR)
    rf["high_confidence"] = (rf.Template_Identity >= HC_IDENTITY) & (rf.Template_Coverage >= HC_COVERAGE)
    den16 = (meta.bacteria_16S_fragments + meta.archaea_16S_fragments)
    rf["FPKM_16S"] = rf.fragmentCountAln / (rf.length_bp / 1000) / (rf.run_accession.map(den16) / 1e6)
    rf["RPKM_depth"] = rf.fragmentCountAln / (rf.length_bp / 1000) / (rf.run_accession.map(meta.qc_pairs) / 1e6)
    return rf

def gene_matrix(rf, runs, value="FPKM_16S", tier="all") -> pd.DataFrame:
    d = rf if tier == "all" else rf[rf.high_confidence]
    return d.pivot_table(index="run_accession", columns="gene", values=value, aggfunc="sum", fill_value=0).reindex(list(runs)).fillna(0)

def rarefy_counts(M: pd.DataFrame, depth: int, seed=SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    out = [rng.multivariate_hypergeometric(M.loc[r].values.astype(np.int64), depth) if M.loc[r].sum() >= depth else np.full(M.shape[1], -1) for r in M.index]
    R = pd.DataFrame(out, index=M.index, columns=M.columns)
    return R.where(R >= 0)

def fold_change_ci(a, b, n_boot=4000, seed=SEED):
    """Ratio of medians a/b with percentile bootstrap 95% CI (resampling within groups)."""
    a = np.asarray(a, float); b = np.asarray(b, float); rng = np.random.default_rng(seed)
    bs = [np.median(rng.choice(a, len(a))) / np.median(rng.choice(b, len(b))) for _ in range(n_boot)]
    return float(np.median(a) / np.median(b)), *np.percentile(bs, [2.5, 97.5])

def exact_perm_test(a, b, stat=lambda x, y: np.median(x) - np.median(y), n_perm=20000, seed=SEED):
    """Two-sided permutation test (exact enumeration when feasible) for a difference between two small groups."""
    from itertools import combinations
    a = np.asarray(a, float); b = np.asarray(b, float); x = np.concatenate([a, b]); na = len(a)
    s0 = abs(stat(a, b)); from math import comb
    if comb(len(x), na) <= 50000:
        cnt = tot = 0
        for idx in combinations(range(len(x)), na):
            m = np.zeros(len(x), bool); m[list(idx)] = True
            cnt += abs(stat(x[m], x[~m])) >= s0 - 1e-12; tot += 1
        return cnt / tot, "exact"
    rng = np.random.default_rng(seed); cnt = 0
    for _ in range(n_perm):
        p = rng.permutation(x); cnt += abs(stat(p[:na], p[na:])) >= s0 - 1e-12
    return (cnt + 1) / (n_perm + 1), "Monte Carlo"

def cliffs_delta(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return float(np.mean(np.sign(a[:, None] - b[None, :])))
