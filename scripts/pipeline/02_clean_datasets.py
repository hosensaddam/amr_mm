"""Step 02 - Systematic, non-destructive cleaning of the chapter datasets.

Principles (see 08_Methodology_and_Reports/Data_Cleaning_Methodology.docx):
* No observation is deleted. Row counts before and after cleaning must be identical (asserted).
* Original values are never overwritten silently: numeric columns that contained missing-value
  codes get a typed numeric column plus a companion `<column>_missing_code` column that keeps the code.
* Outliers and suspicious values are FLAGGED in a `cleaning_flags` column, never removed or winsorised.
* Variable names are standardised to lower snake_case; the original name is kept in the dictionary.
* Every operation is counted and written to Cleaning_Operations_Log.csv.

Outputs: 02_Cleaned_Datasets/{Cleaned_CSV, Cleaned_Excel, Data_Dictionaries}, plus
03_Data_Analysis/Descriptive_Statistics/Dataset_Profiles.xlsx and 00_Project_Inventory/cleaning_qc.json
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from amr_paths import (SOURCE, HARMONIZED, HARMONIZED_UPDATED, REPROCESSED, CLEAN_CSV, CLEAN_XLSX,
                       DICT, ANALYSIS, INV, log)

MISSING_CODES = {"Not reported", "Not available", "Not yet verified", "Not applicable", "Restricted access",
                 "not applicable", "not collected", "missing", "NA", "N/A", "nan", ""}
EXCEL_MAX_ROWS = 200_000            # larger tables are delivered as (gzipped) CSV only
GZIP_ABOVE_MB = 40                  # CSVs larger than this are gzip-compressed (GitHub limit 100 MB)
TODAY = pd.Timestamp("2026-10-31")  # data-collection date; later collection dates are impossible

# --------------------------------------------------------------------------------------------
# Dataset configuration. desc: variable descriptions not already in the source data dictionary.
# --------------------------------------------------------------------------------------------
DATASETS = [
    dict(name="DS_A_sample_master", file=HARMONIZED_UPDATED / "sample_master.csv",
         title="Harmonized sample master (one row per repository sample; all screened samples)",
         key=["Sample_ID"], numeric=["N_Runs", "Latitude", "Longitude", "Base_Count_Total", "FASTQ_Bytes_Total",
                                     "Collection_Year", "Collection_Month", "Latitude_Alternative",
                                     "Longitude_Alternative"],
         date=["Collection_Date"], lat="Latitude", lon="Longitude", depth="Base_Count_Total",
         outlier=["Base_Count_Total", "FASTQ_Bytes_Total"],
         unit={"Latitude": "decimal degrees (WGS84)", "Longitude": "decimal degrees (WGS84)",
               "Base_Count_Total": "bases", "FASTQ_Bytes_Total": "bytes", "Collection_Month": "month (1-12)",
               "Collection_Year": "year"}),
    dict(name="DS_B_run_master", file=HARMONIZED / "run_master.csv",
         title="Harmonized run master (one row per sequencing run)",
         key=["Run_ID"], numeric=["Base_Count", "Read_Count", "FASTQ_Bytes", "SRC_lat", "SRC_lon"],
         outlier=["Base_Count", "Read_Count"], lat="SRC_lat", lon="SRC_lon",
         unit={"Base_Count": "bases", "Read_Count": "reads", "FASTQ_Bytes": "bytes"}),
    dict(name="DS_C_amr_genes_long", file=HARMONIZED / "amr_genes_long.csv",
         title="Producer-processed AMR gene detections (long format; CARD/RGI for DS02, ResFinder/KMA for the Martiny et al. resource)",
         key=["Dataset_ID", "Run_ID", "AMR_Gene", "Gene_Identifier"],
         numeric=["Detection_Value", "Completely_Mapped_Reads", "Average_Percent_Coverage", "Reference_Length_bp"],
         ranges={"Average_Percent_Coverage": (0, 100), "Detection_Value": (0, None)},
         outlier=["Detection_Value"], unit={"Detection_Value": "reads / fragments (see detection_unit)",
                                            "Average_Percent_Coverage": "%", "Reference_Length_bp": "bp"}),
    dict(name="DS_D_taxonomy_long", file=HARMONIZED / "taxonomy_long.csv",
         title="Producer-processed taxonomic profiles (DS02, Kraken 2/Bracken, phylum level)",
         key=["Run_ID", "Taxon_Name", "Taxonomic_Rank"],
         numeric=["Abundance_Value", "Relative_Abundance", "Kraken_Assigned_Reads"],
         ranges={"Relative_Abundance": (0, 1), "Abundance_Value": (0, None)},
         unit={"Abundance_Value": "Bracken-estimated reads", "Relative_Abundance": "fraction (0-1)"}),
    dict(name="DS_E_diversity_long", file=HARMONIZED / "diversity_long.csv",
         title="Producer-provided diversity indices (Martiny et al. ResFinder/rRNA resource)",
         key=["Run_ID", "Feature_Category"],
         numeric=["Richness_n", "Shannon", "Simpson", "Total_Fragments", "Category_Fragments"],
         ranges={"Simpson": (0, 1), "Shannon": (0, None), "Richness_n": (0, None)},
         outlier=["Total_Fragments"], unit={"Richness_n": "count", "Shannon": "nats", "Simpson": "unitless (0-1)",
                                            "Total_Fragments": "fragments"}),
    dict(name="DS_F_env_covariates_long", file=HARMONIZED / "env_covariates_long.csv",
         title="Environmental and design covariates reported by the sources (long format)",
         key=["Dataset_ID", "Sample_ID", "Covariate"]),
    dict(name="DS_G_reprocessed_run_qc", file=REPROCESSED / "run_qc.csv",
         title="Stream 3 uniform reprocessing: per-run read QC, 16S and ResFinder yields (51 runs)",
         key=["run_accession"],
         numeric=["requested_pairs", "acquired_pairs", "q30_rate_raw", "q30_rate_qc", "pct_reads_passed",
                  "duplication_rate", "bacteria_16S_fragments", "resfinder_aligned_fragments", "lat", "lon"],
         ranges={"q30_rate_raw": (0, 1), "q30_rate_qc": (0, 1), "pct_reads_passed": (0, 100)},
         lat="lat", lon="lon", outlier=["bacteria_16S_fragments", "resfinder_aligned_fragments"]),
    dict(name="DS_H_reprocessed_resfinder_annotated", file=REPROCESSED / "resfinder_annotated_all_runs.csv",
         title="Stream 3 uniform reprocessing: ResFinder/KMA template hits with annotation and normalisation",
         key=["run_accession", "template"],
         numeric=["fragmentCount", "Template_Identity", "Template_Coverage", "FPKM_16S", "RPKM_depth"],
         ranges={"Template_Identity": (0, 100), "Template_Coverage": (0, None)},
         unit={"FPKM_16S": "fragments per kb per million prokaryotic 16S fragments",
               "RPKM_depth": "reads per kb per million QC read pairs"}),
    dict(name="DS_I_reprocessed_taxonomy_kraken_long", file=REPROCESSED / "taxonomy_kraken_long.csv",
         title="Stream 3 uniform reprocessing: Kraken 2 (SILVA 138 16S) clade counts, all ranks",
         key=["run_accession", "rank", "taxid"], numeric=["clade_fragments", "direct_fragments"]),
    dict(name="DS_J_reprocessed_amr_summary", file=REPROCESSED / "amr_summary_all_runs.csv",
         title="Stream 3 uniform reprocessing: per-run ARG totals and diversity",
         key=["run_accession"],
         numeric=["total_FPKM_16S", "total_RPKM_depth", "genes_any_hit", "genes_high_confidence"],
         outlier=["total_FPKM_16S"]),
    dict(name="DS_K_reprocessed_comparative_variables", file=REPROCESSED / "comparative_variables_all_runs.csv",
         title="Stream 3 uniform reprocessing: analysis variables for Dhaka and South Asian comparators",
         key=["run_accession"], date=["Collection_Date"],
         numeric=["log10_total_FPKM_16S", "ARG_genes_high_conf", "genus_richness_rarefied", "genus_shannon"]),
    dict(name="DS_L_dhaka_geographic_verification", file=REPROCESSED / "dhaka_geographic_verification.csv",
         title="Dhaka-labelled samples: coordinate provenance and point-in-polygon verification",
         key=["Sample_ID", "coordinate_source"], numeric=["latitude", "longitude",
                                                           "distance_to_DNCC_DSCC_boundary_km"],
         lat="latitude", lon="longitude"),
]

# Descriptions for variables not covered by the source data dictionaries
EXTRA_DESC = {
    "Sample_Type_Original_Harmonized": "Sample type before the AUD-04 sample-level correction",
    "Sample_Type_Correction": "Audit correction applied to sample type (AUD-04), or 'No correction'",
    "Site_ID": "Open-drain site identifier (City | Site) recovered from ENA experiment titles (DS02 only)",
    "Site_ID_Source": "Provenance of Site_ID",
    "Latitude_Alternative": "Second source-reported latitude retained when two sources disagree (not averaged)",
    "Longitude_Alternative": "Second source-reported longitude retained when two sources disagree",
    "Coordinate_Provenance_Note": "Note on the origin and precision of the coordinates (AUD-18 verification)",
    "run_accession": "INSDC run accession", "processing_status": "Status of uniform reprocessing",
    "requested_pairs": "Read pairs requested (head-of-file subsample)", "acquired_pairs": "Read pairs acquired",
    "q30_rate_raw": "Fraction of bases with Phred >= 30 before QC", "q30_rate_qc": "Fraction of bases with Phred >= 30 after fastp",
    "pct_reads_passed": "Percentage of reads passing fastp filters", "duplication_rate": "fastp duplication rate",
    "bacteria_16S_fragments": "Read pairs classified as bacterial 16S by Kraken 2 (SILVA 138)",
    "resfinder_aligned_fragments": "Fragments aligned to ResFinder templates by KMA",
    "template": "ResFinder template (gene_allele_accession)", "fragmentCount": "KMA fragments assigned to template",
    "Template_Identity": "KMA template identity (%)", "Template_Coverage": "KMA template coverage (%)",
    "FPKM_16S": "Fragments per kilobase of template per million prokaryotic 16S fragments",
    "RPKM_depth": "Reads per kilobase per million QC read pairs", "high_confidence": "Identity >= 90% and coverage >= 60%",
    "Class_Group": "Antibiotic class group used in the chapter", "Clinical": "Clinically important determinant group (if any)",
    "rank": "Taxonomic rank code (D, P, C, O, F, G, S)", "taxid": "Kraken 2 SILVA taxon id", "name": "Taxon name",
    "clade_fragments": "Fragments assigned to the clade", "direct_fragments": "Fragments assigned directly to the taxon",
    "lineage": "Rank-prefixed lineage string",
    "total_FPKM_16S": "Total ARG abundance, FPKM-16S", "total_RPKM_depth": "Total ARG abundance, RPKM-depth",
    "genes_any_hit": "ResFinder genes with any alignment", "genes_high_confidence": "ResFinder genes at high confidence",
    "log10_total_FPKM_16S": "log10 of total ARG FPKM-16S", "ARG_genes_high_conf": "High-confidence ARG genes",
    "genus_richness_rarefied": "Genus richness rarefied to 850 genus-assigned fragments", "genus_shannon": "Genus Shannon index (nats)",
    "City_group": "Comparison group (city)", "Period": "Collection period group", "Unit": "Sampling unit used for permutation tests",
    "within_Dhaka_city_corporations_DNCC_DSCC": "Point inside Dhaka North/South City Corporation polygons",
    "distance_to_DNCC_DSCC_boundary_km": "Distance to nearest city-corporation boundary (km)",
}


def snake(name: str) -> str:
    s = re.sub(r"[^0-9A-Za-z]+", "_", str(name)).strip("_")
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s)
    s = re.sub(r"_+", "_", s).lower()
    return s.replace("16_s", "16s").replace("18_s", "18s")


def source_dictionary() -> dict:
    """Descriptions from the original workbooks' Data_Dictionary sheets (read-only)."""
    out = {}
    for wb, sheet in [("02_Sample_Metadata_Master.xlsx", "Data_Dictionary"),
                      ("03_Analysis_Ready_Metagenomic_Data.xlsx", "Data_Dictionary")]:
        try:
            d = pd.read_excel(SOURCE / wb, sheet_name=sheet)
            for _, r in d.iterrows():
                out.setdefault(str(r["Field"]), (str(r["Description"]), str(r.get("Type_or_Unit", "")),
                                                 "" if pd.isna(r.get("Allowed_Values_or_Codes")) else
                                                 str(r.get("Allowed_Values_or_Codes"))))
        except Exception as exc:
            log(f"WARNING could not read dictionary {wb}: {exc}")
    return out


def date_precision(s: str) -> str:
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", s): return "day"
    if re.fullmatch(r"\d{4}-\d{2}", s): return "month"
    if re.fullmatch(r"\d{4}", s): return "year"
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}/\d{4}-\d{2}-\d{2}", s): return "interval"
    return "unparsed"


def clean_one(cfg: dict, srcdict: dict, ops: list, profiles: list, dict_rows: list) -> pd.DataFrame:
    name = cfg["name"]
    df = pd.read_csv(cfg["file"], dtype=str, keep_default_na=False, low_memory=False)
    n0, c0 = df.shape
    flags = pd.Series([[] for _ in range(n0)], index=df.index)

    def op(step, detail, n):
        ops.append(dict(Dataset=name, Step=step, Detail=detail, Records_or_Cells_Affected=int(n)))

    op("load", f"Read {Path(cfg['file']).name} as text", n0)

    # 1. whitespace normalisation
    changed = 0
    for c in df.columns:
        new = df[c].str.strip().str.replace(r"\s{2,}", " ", regex=True)
        changed += int((new != df[c]).sum())
        df[c] = new
    op("whitespace", "Trimmed leading/trailing spaces and collapsed repeated spaces", changed)

    # 2. exact duplicate rows and key duplicates (flag, do not drop)
    dup_rows = df.duplicated(keep="first")
    op("duplicates", "Exact duplicate rows (flagged, retained)", dup_rows.sum())
    for i in df.index[dup_rows]:
        flags[i].append("exact_duplicate_row")
    key = cfg.get("key")
    if key:
        kd = df.duplicated(key, keep=False)
        op("duplicates", f"Rows sharing key {'+'.join(key)} (flagged, retained)", kd.sum())
        for i in df.index[kd & ~dup_rows]:
            flags[i].append("duplicate_key")

    # 3. typed numeric columns with missing-code companions
    for c in cfg.get("numeric", []):
        if c not in df.columns:
            continue
        raw = df[c]
        num = pd.to_numeric(raw.str.replace(",", "", regex=False), errors="coerce")
        nonnum = raw[num.isna() & (raw != "")]
        codes = nonnum[nonnum.isin(MISSING_CODES)]
        bad = nonnum[~nonnum.isin(MISSING_CODES)]
        if len(codes) or len(bad):
            df[c + "_missing_code"] = np.where(num.isna(), raw.replace("", "Blank"), "")
        df[c] = num
        op("type", f"{c}: converted to numeric; {len(codes)} missing-value codes moved to companion column", len(codes))
        if len(bad):
            op("type", f"{c}: {len(bad)} non-numeric values that are not documented codes (kept in companion column)",
               len(bad))
            for i in bad.index:
                flags[i].append(f"{snake(c)}_non_numeric")
        else:
            pass
        df[c] = df[c].replace({np.nan: None}).astype("Float64")

    # Placeholder: zero bases where no FASTQ exists (MG-RAST records) -> NA, flagged
    dcol = cfg.get("depth")
    if dcol and dcol in df.columns:
        z = df[dcol] == 0
        if z.any():
            df[dcol + "_missing_code"] = df.get(dcol + "_missing_code", pd.Series("", index=df.index))
            df.loc[z, dcol + "_missing_code"] = "Zero placeholder (no FASTQ base count reported)"
            df.loc[z, dcol] = pd.NA
            op("placeholder", f"{dcol}: 0 used as placeholder -> set to missing with code", int(z.sum()))

    # 4. ranges
    for c, (lo, hi) in cfg.get("ranges", {}).items():
        if c not in df.columns: continue
        v = df[c].astype(float)
        out = (v < lo) if lo is not None else pd.Series(False, index=df.index)
        if hi is not None:
            out |= v > hi
        out = out.fillna(False)
        op("range", f"{c}: values outside [{lo}, {hi}] (flagged)", out.sum())
        for i in df.index[out]:
            flags[i].append(f"{snake(c)}_out_of_range")
    lat, lon = cfg.get("lat"), cfg.get("lon")
    if lat and lat in df.columns:
        la, lo_ = df[lat].astype(float), df[lon].astype(float)
        bad = ((la.abs() > 90) | (lo_.abs() > 180)).fillna(False)
        zero = ((la == 0) & (lo_ == 0)).fillna(False)
        op("range", f"{lat}/{lon}: outside valid WGS84 range (flagged)", bad.sum())
        op("range", f"{lat}/{lon}: (0, 0) null-island coordinates (flagged)", zero.sum())
        for i in df.index[bad]: flags[i].append("coordinate_out_of_range")
        for i in df.index[zero]: flags[i].append("coordinate_zero_zero")
        half = (la.notna() ^ lo_.notna())
        op("consistency", f"{lat}/{lon}: only one of the pair reported", half.sum())
        for i in df.index[half]: flags[i].append("coordinate_pair_incomplete")

    # 5. dates
    for c in cfg.get("date", []):
        if c not in df.columns: continue
        raw = df[c]
        prec = raw.map(lambda s: "missing_code" if s in MISSING_CODES else date_precision(s))
        df[snake(c) + "_precision"] = prec
        parsed = pd.to_datetime(raw.str.slice(0, 10).where(prec.isin(["day", "interval"])), errors="coerce")
        df[snake(c) + "_iso"] = parsed.dt.strftime("%Y-%m-%d")
        op("date", f"{c}: precision classified (day/month/year/interval/missing); ISO date added where day-level",
           int(prec.isin(["day", "month", "year", "interval"]).sum()))
        unp = prec == "unparsed"
        op("date", f"{c}: values not matching an ISO 8601 pattern (flagged)", unp.sum())
        for i in df.index[unp]: flags[i].append("date_unparsed")
        ph = raw.str.startswith("1970-01-01")
        op("date", f"{c}: 1970-01-01 placeholder dates (flagged)", ph.sum())
        for i in df.index[ph]: flags[i].append("date_placeholder_1970")
        year = pd.to_numeric(raw.str.slice(0, 4), errors="coerce")
        fut = (parsed > TODAY).fillna(False) | (year > TODAY.year).fillna(False)
        op("date", f"{c}: dates after the data-collection date (flagged)", fut.sum())
        for i in df.index[fut]: flags[i].append("date_in_future")
        if "Collection_Year" in df.columns:
            cy = df["Collection_Year"].astype(float)
            mism = (year.notna() & cy.notna() & (year != cy))
            op("consistency", "Collection_Year differs from year of Collection_Date (flagged)", mism.sum())
            for i in df.index[mism]: flags[i].append("year_date_mismatch")
        if "Collection_Month" in df.columns:
            m = df["Collection_Month"].astype(float)
            badm = ((m < 1) | (m > 12)).fillna(False)
            op("range", "Collection_Month outside 1-12 (flagged)", badm.sum())

    # 6. outliers: Tukey far-out rule (3 x IQR) on log10 scale; flagged only
    for c in cfg.get("outlier", []):
        if c not in df.columns: continue
        v = df[c].astype(float)
        lv = np.log10(v.where(v > 0))
        q1, q3 = lv.quantile([0.25, 0.75])
        iqr = q3 - q1
        lo_, hi_ = q1 - 3 * iqr, q3 + 3 * iqr
        o = ((lv < lo_) | (lv > hi_)).fillna(False)
        op("outlier", f"{c}: log10 values beyond Q1-3*IQR / Q3+3*IQR [{10**lo_:.3g}, {10**hi_:.3g}] (flagged, retained)",
           o.sum())
        for i in df.index[o]: flags[i].append(f"{snake(c)}_extreme_value")

    # 7. categorical consistency: variants differing only in case/punctuation (reported, not merged)
    variants = []
    for c in df.columns:
        if df[c].dtype != object: continue
        u = df[c].dropna().unique()
        if 1 < len(u) <= 400:
            norm = pd.Series(u).str.lower().str.replace(r"[^0-9a-z]", "", regex=True)
            g = pd.Series(u).groupby(norm.values).apply(list)
            for vals in g[g.map(len) > 1]:
                variants.append(dict(Dataset=name, Variable=c, Variants=" | ".join(map(str, vals))))
    if variants:
        op("categorical", "Value variants differing only by case/punctuation (reported in Categorical_Variants; "
                          "not merged because they are source-reported labels)", len(variants))

    # Dataset-specific consistency checks
    if name == "DS_A_sample_master":
        iso = df.groupby("Country")["ISO_Country_Code"].nunique()
        op("consistency", "Countries mapped to more than one ISO code", (iso > 1).sum())
        inc = df["Inclusion_Status"] == "Included"
        notsg = inc & ~df["Sequencing_Type"].str.contains("Shotgun", case=False)
        op("consistency", "Included samples without a shotgun sequencing type (flagged)", notsg.sum())
        for i in df.index[notsg]: flags[i].append("included_not_shotgun")
        lowdepth = inc & (df["Base_Count_Total"].astype(float) < 1e8)
        op("depth", "Included samples with < 100 Mbp repository depth (flagged; analysed as sensitivity, Table S92)",
           lowdepth.sum())
        for i in df.index[lowdepth.fillna(False)]: flags[i].append("depth_below_100Mbp")

    df["cleaning_flags"] = flags.map(lambda l: ";".join(dict.fromkeys(l)))
    nflag = int((df["cleaning_flags"] != "").sum())
    op("flags", "Records with at least one cleaning flag", nflag)

    # 8. standardised names
    mapping = {c: snake(c) for c in df.columns}
    assert len(set(mapping.values())) == len(mapping), f"name collision in {name}"
    renamed = sum(1 for k, v in mapping.items() if k != v)
    df = df.rename(columns=mapping)
    op("rename", "Variable names standardised to lower snake_case", renamed)

    assert len(df) == n0, "row count changed during cleaning"
    op("validate", f"Row count preserved ({n0} -> {len(df)}); columns {c0} -> {df.shape[1]}", len(df))

    # Data dictionary
    inv = {v: k for k, v in mapping.items()}
    for c in df.columns:
        orig = inv.get(c, c)
        base = orig.replace("_missing_code", "")
        desc, typ, allowed = srcdict.get(orig, (EXTRA_DESC.get(orig, ""), "", ""))
        if orig.endswith("_missing_code"):
            desc = f"Missing-value code for {snake(base)} (original text where the numeric value is missing)"
        elif c.endswith("_precision"):
            desc = f"Precision of {c.replace('_precision', '')} (day, month, year, interval, missing_code, unparsed)"
        elif c.endswith("_iso"):
            desc = f"ISO 8601 date (YYYY-MM-DD) derived from {c.replace('_iso', '')} where day-level"
        elif c == "cleaning_flags":
            desc = "Semicolon-separated quality flags raised during cleaning (empty = none); records are never removed"
        elif not desc:
            desc = orig.replace("_", " ")
        s = df[c]
        nonmiss = s.notna() & (s.astype(str) != "")
        if pd.api.types.is_numeric_dtype(s):
            dtype = "numeric"
            sv = s.dropna().astype(float)
            possible = f"min {sv.min():.6g}; median {sv.median():.6g}; max {sv.max():.6g}" if len(sv) else "all missing"
            ncode = int(s.isna().sum())
        else:
            dtype = "text (categorical)" if s.nunique() <= 30 else "text"
            u = s[nonmiss].value_counts()
            possible = "; ".join(f"{k} ({v})" for k, v in u.head(12).items()) + (" ..." if len(u) > 12 else "")
            ncode = int(s.isin(MISSING_CODES).sum())
        dict_rows.append(dict(
            Dataset=name, Variable=c, Original_Name=orig, Description=desc,
            Data_Type=dtype, Unit=cfg.get("unit", {}).get(orig, typ if typ not in ("Text", "nan") else ""),
            Possible_Values_or_Range=possible[:500], Allowed_Codes=allowed,
            N_Records=len(s), N_Missing_or_Coded=ncode, Pct_Missing_or_Coded=round(100 * ncode / max(len(s), 1), 2),
            Transformation=("renamed" if c != orig else "") +
                           ("; typed numeric, codes moved to companion column" if orig in cfg.get("numeric", []) else "") +
                           ("; derived during cleaning" if orig not in inv.values() or orig.endswith(("_missing_code",)) or c.endswith(("_precision", "_iso")) or c == "cleaning_flags" else "")))
    profiles.append(dict(Dataset=name, Title=cfg["title"], Source_File=str(Path(cfg["file"]).relative_to(SOURCE)).replace("\\", "/"),
                         Rows_Original=n0, Columns_Original=c0, Rows_Cleaned=len(df), Columns_Cleaned=df.shape[1],
                         Exact_Duplicate_Rows=int(dup_rows.sum()), Records_Flagged=nflag,
                         Pct_Records_Flagged=round(100 * nflag / n0, 2)))
    return df, variants


def main():
    for d in (CLEAN_CSV, CLEAN_XLSX, DICT, ANALYSIS / "Descriptive_Statistics"):
        d.mkdir(parents=True, exist_ok=True)
    srcdict = source_dictionary()
    ops, profiles, dict_rows, all_variants = [], [], [], []
    cleaned = {}
    for cfg in DATASETS:
        df, variants = clean_one(cfg, srcdict, ops, profiles, dict_rows)
        all_variants += variants
        cleaned[cfg["name"]] = df
        csv_path = CLEAN_CSV / f"{cfg['name']}.csv"
        df.to_csv(csv_path, index=False)
        if csv_path.stat().st_size / 1e6 > GZIP_ABOVE_MB:
            df.to_csv(csv_path.with_suffix(".csv.gz"), index=False, compression="gzip")
            csv_path.unlink()
            log(f"  {cfg['name']}: CSV > {GZIP_ABOVE_MB} MB, written gzip-compressed")
        if len(df) <= EXCEL_MAX_ROWS:
            with pd.ExcelWriter(CLEAN_XLSX / f"{cfg['name']}.xlsx", engine="openpyxl") as xw:
                df.to_excel(xw, sheet_name="data", index=False)
                pd.DataFrame([r for r in dict_rows if r["Dataset"] == cfg["name"]]).to_excel(
                    xw, sheet_name="dictionary", index=False)
        log(f"Cleaned {cfg['name']}: {df.shape[0]} rows x {df.shape[1]} cols")

    # Cross-dataset referential integrity
    sm = cleaned["DS_A_sample_master"]
    samples = set(sm["sample_id"])
    runs = set(cleaned["DS_B_run_master"]["run_id"])
    for dname, col in [("DS_C_amr_genes_long", "sample_id"), ("DS_D_taxonomy_long", "sample_id"),
                       ("DS_E_diversity_long", "sample_id"), ("DS_F_env_covariates_long", "sample_id"),
                       ("DS_K_reprocessed_comparative_variables", "sample_id"),
                       ("DS_L_dhaka_geographic_verification", "sample_id")]:
        vals = set(cleaned[dname][col])
        miss = vals - samples
        ops.append(dict(Dataset=dname, Step="referential", Detail=f"{col} values not found in sample master",
                        Records_or_Cells_Affected=len(miss)))
    sm_runs = set(r for x in sm["run_id"] for r in str(x).split(";"))
    ops.append(dict(Dataset="DS_A_sample_master", Step="referential",
                    Detail="Run IDs listed in sample master but absent from run master",
                    Records_or_Cells_Affected=len(sm_runs - runs)))
    rq = set(cleaned["DS_G_reprocessed_run_qc"]["run_accession"])
    ops.append(dict(Dataset="DS_G_reprocessed_run_qc", Step="referential",
                    Detail="Reprocessed runs absent from run master", Records_or_Cells_Affected=len(rq - runs)))
    for dname in ["DS_H_reprocessed_resfinder_annotated", "DS_I_reprocessed_taxonomy_kraken_long",
                  "DS_J_reprocessed_amr_summary"]:
        miss = set(cleaned[dname]["run_accession"]) - rq
        ops.append(dict(Dataset=dname, Step="referential", Detail="run_accession absent from reprocessed run QC table",
                        Records_or_Cells_Affected=len(miss)))
    inc = sm[sm["inclusion_status"] == "Included"]
    ops.append(dict(Dataset="DS_A_sample_master", Step="validate",
                    Detail="Included samples (expected 8,053 as reported in Table S1)", Records_or_Cells_Affected=len(inc)))

    opsdf, prof, dd = pd.DataFrame(ops), pd.DataFrame(profiles), pd.DataFrame(dict_rows)
    var = pd.DataFrame(all_variants)
    opsdf.to_csv(CLEAN_CSV.parent / "Cleaning_Operations_Log.csv", index=False)
    dd.to_csv(DICT / "Data_Dictionary_All_Datasets.csv", index=False)
    with pd.ExcelWriter(DICT / "Data_Dictionary_All_Datasets.xlsx", engine="openpyxl") as xw:
        prof.to_excel(xw, sheet_name="Datasets", index=False)
        dd.to_excel(xw, sheet_name="Variables", index=False)
        for n in dd.Dataset.unique():
            dd[dd.Dataset == n].to_excel(xw, sheet_name=n[:31], index=False)
    with pd.ExcelWriter(ANALYSIS / "Descriptive_Statistics" / "Dataset_Profiles.xlsx", engine="openpyxl") as xw:
        prof.to_excel(xw, sheet_name="Dataset_Summary", index=False)
        opsdf.to_excel(xw, sheet_name="Cleaning_Operations", index=False)
        (var if len(var) else pd.DataFrame(columns=["Dataset", "Variable", "Variants"])).to_excel(
            xw, sheet_name="Categorical_Variants", index=False)
        flagcounts = []
        for n, df in cleaned.items():
            f = df["cleaning_flags"].str.split(";").explode()
            for k, v in f[f != ""].value_counts().items():
                flagcounts.append(dict(Dataset=n, Flag=k, Records=v))
        pd.DataFrame(flagcounts).to_excel(xw, sheet_name="Flag_Counts", index=False)
        miss = []
        for n, df in cleaned.items():
            for c in df.columns:
                s = df[c]
                m = s.isna() | s.astype(str).isin(MISSING_CODES)
                miss.append(dict(Dataset=n, Variable=c, Missing_or_Coded=int(m.sum()), Pct=round(100 * m.mean(), 2)))
        pd.DataFrame(miss).to_excel(xw, sheet_name="Completeness", index=False)
    json.dump(dict(profiles=profiles, n_ops=len(ops)), open(INV / "cleaning_qc.json", "w"), indent=1, default=int)
    log(f"Cleaning complete: {len(cleaned)} datasets; {len(dd)} dictionary entries; {len(var)} categorical variant groups")
    print(prof.to_string(index=False))


if __name__ == "__main__":
    main()
