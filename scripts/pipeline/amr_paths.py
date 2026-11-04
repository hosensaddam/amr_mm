"""Shared paths for the supplementary-material workflow.

All paths are relative. PROJECT is the output package folder (two levels above
this file). SOURCE is the research directory that holds the original study
materials; by default it is the parent of PROJECT, and it can be overridden
with the environment variable AMR_SOURCE_DIR.
"""
import os
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
SOURCE = Path(os.environ.get("AMR_SOURCE_DIR", PROJECT.parent)).resolve()
# The same scripts run inside the GitHub repository (scripts/pipeline/), whose folders are named differently.
REPO_LAYOUT = (PROJECT / "data" / "cleaned").exists() and not (PROJECT / "02_Cleaned_Datasets").exists()

INV = PROJECT / "00_Project_Inventory"
RAW = PROJECT / "01_Raw_Data_Archive"
CLEAN_CSV = PROJECT / "02_Cleaned_Datasets" / "Cleaned_CSV"
CLEAN_XLSX = PROJECT / "02_Cleaned_Datasets" / "Cleaned_Excel"
DICT = PROJECT / "02_Cleaned_Datasets" / "Data_Dictionaries"
ANALYSIS = PROJECT / "03_Data_Analysis"
TABLES = PROJECT / "04_Tables"
FIGURES = PROJECT / "05_Figures"
DOCS = PROJECT / "06_Research_Documents"
REPORTS = PROJECT / "08_Methodology_and_Reports"
LOG = INV / "Processing_Log.txt"

# Verified analysis (31 Oct 2026) that the supplement is built from
UPDATED = SOURCE / "Updated_Research_Analysis"
HARMONIZED = SOURCE / "AMR_Chapter_Reproducibility_Bundle" / "AMR_Chapter_Public_Data_Collection" / "harmonized"
HARMONIZED_UPDATED = UPDATED / "01_Data" / "harmonized_updated"
REPROCESSED = UPDATED / "01_Data" / "reprocessed"
BUILD_V2 = SOURCE / "08_Final_Submission" / "Build_Sources" / "v2"
SUPP_DOCX_SOURCE = SOURCE / "Supplementary Material.docx"

RF_UNIT_SUMMARY_NAME = "rf_unit_summary.csv"   # ResFinder analytic-unit table (corrected analysis build)

if REPO_LAYOUT:
    INV = PROJECT / "documentation" / "quality_control"
    CLEAN_CSV = PROJECT / "data" / "cleaned"
    DICT = PROJECT / "data" / "data_dictionary"
    ANALYSIS = PROJECT / "analysis"
    TABLES = PROJECT / "tables"
    FIGURES = PROJECT / "figures"
    DOCS = PROJECT / "supplementary_materials"
    LOG = PROJECT / "analysis" / "verification_log.txt"


def verified_table(name: str) -> Path:
    """Verified analysis table by file name (e.g. 'D02_by_sample_type.csv'). Uses the original analysis folder
    when available, otherwise the exported supplementary table Table_Sxx_<ID>.csv of this package/repository."""
    p = UPDATED / "04_Tables" / name
    if p.exists():
        return p
    tid = name.split("_")[0]
    for d in (TABLES / "CSV_Tables", TABLES / "csv"):
        hits = sorted(d.glob(f"Table_S*_{tid}.csv"))
        if hits:
            return hits[0]
    raise FileNotFoundError(name)


def rf_unit_summary() -> Path:
    for p in (SOURCE / "Independent_Audit" / "03_Corrected_Analysis_Results" / "derived" / RF_UNIT_SUMMARY_NAME,
              ANALYSIS / "Comparative_Analysis" / "inputs" / RF_UNIT_SUMMARY_NAME,
              ANALYSIS / "comparative_analysis" / "inputs" / RF_UNIT_SUMMARY_NAME):
        if p.exists():
            return p
    raise FileNotFoundError(RF_UNIT_SUMMARY_NAME)


# Folders that are outputs of this workflow or a git clone; never inventoried as sources
EXCLUDE_DIRS = {PROJECT.name, "amr_mm"}


def log(msg: str) -> None:
    """Append a timestamped line to the processing log and echo it."""
    import datetime
    LOG.parent.mkdir(parents=True, exist_ok=True)
    line = f"{datetime.datetime.now().isoformat(timespec='seconds')}  {msg}"
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line)
