"""Step 00 - Archive copies of the original source data and code used by this package.

Read-only copies (never modified) of the primary data-collection workbooks, the source log, the original
Supplementary Material.docx and the harmonized tables that step 02 cleans; CSVs larger than 40 MB are
stored gzip-compressed. A SHA-256 manifest records the checksum of every ORIGINAL file so the archive can
be verified against the source. Also copies the verified analysis code and executed notebooks.
"""
import gzip
import hashlib
import shutil
import stat

import pandas as pd

from amr_paths import SOURCE, RAW, PROJECT, HARMONIZED, HARMONIZED_UPDATED, REPROCESSED, UPDATED, BUILD_V2, log

ARCHIVE = {
    "Data_Collection_Workbooks": [SOURCE / f for f in (
        "00_Master_Data_Index.xlsx", "01_Public_Dataset_Inventory.xlsx", "02_Sample_Metadata_Master.xlsx",
        "03_Analysis_Ready_Metagenomic_Data.xlsx", "04_AMR_Reference_Database_Register.xlsx",
        "Public_Data_Source_Log.csv", "README_Data_Preparation.md")],
    "Original_Supplementary_Document": [SOURCE / "Supplementary Material.docx"],
    "Harmonized_Source_Tables": [HARMONIZED_UPDATED / "sample_master.csv", HARMONIZED_UPDATED / "sample_master_change_log.csv"]
                                + [HARMONIZED / f for f in ("run_master.csv", "amr_genes_long.csv", "taxonomy_long.csv",
                                                            "diversity_long.csv", "env_covariates_long.csv",
                                                            "transformation_log.csv")],
    "Reprocessed_Source_Tables": [REPROCESSED / f for f in (
        "run_qc.csv", "resfinder_annotated_all_runs.csv", "taxonomy_kraken_long.csv", "amr_summary_all_runs.csv",
        "comparative_variables_all_runs.csv", "dhaka_geographic_verification.csv", "versions.json")],
}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    rows = []
    for folder, files in ARCHIVE.items():
        out = RAW / folder
        out.mkdir(parents=True, exist_ok=True)
        for f in files:
            dest = out / f.name
            if f.suffix == ".csv" and f.stat().st_size > 40e6:
                dest = out / (f.name + ".gz")
                with open(f, "rb") as src, gzip.open(dest, "wb", compresslevel=9) as dst:
                    shutil.copyfileobj(src, dst)
            else:
                shutil.copyfile(f, dest)
            dest.chmod(stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)   # read-only archive copy
            rows.append(dict(Archive_File=str(dest.relative_to(PROJECT)).replace("\\", "/"),
                             Original_Path=str(f.relative_to(SOURCE)).replace("\\", "/"),
                             Original_Bytes=f.stat().st_size, Original_SHA256=sha(f),
                             Compressed=dest.suffix == ".gz"))
    pd.DataFrame(rows).to_csv(RAW / "Archive_Manifest_SHA256.csv", index=False)
    # verified analysis code and executed notebooks (archival; contain original absolute paths)
    code = PROJECT / "07_Code_and_Scripts" / "Original_Analysis_Code"
    nbs = PROJECT / "07_Code_and_Scripts" / "Jupyter_Notebooks"
    ign = shutil.ignore_patterns("__pycache__", "*.pyc")
    shutil.copytree(UPDATED / "02_Scripts", code / "analysis_scripts", dirs_exist_ok=True, ignore=ign)
    shutil.copytree(UPDATED / "03_Notebooks" / "notebook_sources", code / "notebook_sources", dirs_exist_ok=True, ignore=ign)
    shutil.copytree(UPDATED / "03_Notebooks" / "verification", code / "verification", dirs_exist_ok=True, ignore=ign)
    shutil.copytree(SOURCE / "AMR_Chapter_Reproducibility_Bundle" / "AMR_Chapter_Public_Data_Collection" / "scripts",
                    code / "data_collection_scripts", dirs_exist_ok=True, ignore=ign)
    (code / "supplement_build_v2").mkdir(exist_ok=True)
    for f in ("build_supp.py", "supp_figs.py", "overrides.json", "fig_interp_raw.json",
              "interp_g1.json", "interp_g2.json", "interp_g3.json", "interp_g4.json"):
        shutil.copyfile(BUILD_V2 / f, code / "supplement_build_v2" / f)
    nbs.mkdir(parents=True, exist_ok=True)
    for f in sorted((UPDATED / "03_Notebooks").glob("*.ipynb")):
        shutil.copyfile(f, nbs / f.name)
    shutil.copyfile(UPDATED / "06_Logs" / "environment.json", code / "environment_2026-10-31.json") \
        if (UPDATED / "06_Logs" / "environment.json").exists() else None
    log(f"Archive: {len(rows)} source files copied (read-only) with SHA-256 manifest; analysis code and "
        f"{len(list(nbs.glob('*.ipynb')))} notebooks copied")


if __name__ == "__main__":
    main()
