"""Step 10 - Stage the approved, public-release files into the local clone of the GitHub repository.

Copies only an explicit allow-list (below) into ../amr_mm (or AMR_REPO_DIR), mapping the package folders
onto the repository layout. Existing repository files are kept (the five data-collection workbooks and
the source log at the repository root are not touched). Nothing is committed or pushed by this script.
Withheld from publication (see documentation/Publication_Report.md): 01_Raw_Data_Archive (copies of
originals already public or superseded by the cleaned data), the full source-file inventory (internal),
text-review working records, cleaned Excel duplicates of the CSVs, and any chapter manuscript material.
"""
import os
import shutil
from pathlib import Path

from amr_paths import PROJECT, SOURCE, log

REPO = Path(os.environ.get("AMR_REPO_DIR", SOURCE / "amr_mm")).resolve()
P = PROJECT

MAP = [  # (source in package, destination in repository, glob or None for a single file)
    (P / "02_Cleaned_Datasets" / "Cleaned_CSV", "data/cleaned", "*"),
    (P / "02_Cleaned_Datasets" / "Data_Dictionaries", "data/data_dictionary", "*"),
    (P / "02_Cleaned_Datasets" / "Cleaning_Operations_Log.csv", "data/Cleaning_Operations_Log.csv", None),
    (P / "03_Data_Analysis" / "Descriptive_Statistics", "analysis/descriptive_statistics", "*"),
    (P / "03_Data_Analysis" / "Statistical_Analysis", "analysis/statistical_analysis", "*"),
    (P / "03_Data_Analysis" / "Comparative_Analysis", "analysis/comparative_analysis", "**/*"),
    (P / "03_Data_Analysis" / "Analytical_Results", "analysis/results", "*"),
    (P / "04_Tables" / "CSV_Tables", "tables/csv", "*"),
    (P / "04_Tables" / "Excel_Tables", "tables", "*"),
    (P / "04_Tables" / "Table_Index_S1-S96.csv", "tables/Table_Index_S1-S96.csv", None),
    (P / "05_Figures" / "PNG", "figures/png", "*"),
    (P / "05_Figures" / "TIFF", "figures/tiff", "*"),
    (P / "05_Figures" / "PDF", "figures/pdf", "*"),
    (P / "05_Figures" / "Editable_Figures", "figures/editable_plotting_code", "*"),
    (P / "05_Figures" / "Figure_Index_S1-S30.csv", "figures/Figure_Index_S1-S30.csv", None),
    (P / "06_Research_Documents" / "Supplementary_Material_Final.docx", "supplementary_materials/Supplementary_Material_Final.docx", None),
    (P / "06_Research_Documents" / "Supplementary_Material_Final.pdf", "supplementary_materials/Supplementary_Material_Final.pdf", None),
    (P / "06_Research_Documents" / "Supplementary_Appendices" / "Supplementary_Text_Corrections_Log.xlsx",
     "supplementary_materials/appendices/Supplementary_Text_Corrections_Log.xlsx", None),
    (P / "07_Code_and_Scripts" / "Python", "scripts/pipeline", "*.py"),
    (P / "07_Code_and_Scripts" / "Python", "scripts/pipeline", "*.ps1"),
    (P / "07_Code_and_Scripts" / "Original_Analysis_Code", "scripts/original_analysis", "**/*"),
    (P / "07_Code_and_Scripts" / "Jupyter_Notebooks", "notebooks", "*.ipynb"),
    (P / "07_Code_and_Scripts" / "Requirements.txt", "requirements.txt", None),
    (P / "08_Methodology_and_Reports" / "Data_Cleaning_Methodology.docx", "documentation/data_cleaning_report/Data_Cleaning_Methodology.docx", None),
    (P / "08_Methodology_and_Reports" / "Statistical_Analysis_Methodology.docx", "documentation/statistical_analysis_report/Statistical_Analysis_Methodology.docx", None),
    (P / "08_Methodology_and_Reports" / "Technical_Report.docx", "documentation/reproducibility/Technical_Report.docx", None),
    (P / "00_Project_Inventory" / "Data_Quality_Report.docx", "documentation/quality_control/Data_Quality_Report.docx", None),
    (P / "00_Project_Inventory" / "Quality_Control_Checklist.xlsx", "documentation/quality_control/Quality_Control_Checklist.xlsx", None),
    (P / "00_Project_Inventory" / "Publication_Report.md", "documentation/Publication_Report.md", None),
    (P / "README.md", "documentation/Supplementary_Package_README.md", None),
]
SKIP_NAMES = {"__pycache__", "10_stage_github.py"}


def main():
    assert (REPO / ".git").exists(), f"{REPO} is not a git clone"
    n = 0
    staged = []
    for src, dest, pattern in MAP:
        if pattern is None:
            files = [(src, REPO / dest)] if src.exists() else []
        else:
            files = [(f, REPO / dest / f.relative_to(src)) for f in sorted(src.glob(pattern))
                     if f.is_file() and not (set(f.parts) & SKIP_NAMES) and not f.name.startswith("~$")]
        for f, d in files:
            d.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(f, d)
            staged.append(str(d.relative_to(REPO)).replace("\\", "/"))
            n += 1
    (REPO / "documentation").mkdir(exist_ok=True)
    shutil.copyfile(P / "07_Code_and_Scripts" / "Python" / "10_stage_github.py", REPO / "scripts" / "pipeline" / "10_stage_github.py")
    log(f"GitHub staging: {n} files copied into the repository clone (not committed)")
    (PROJECT / "00_Project_Inventory" / "github_staged_files.txt").write_text("\n".join(staged), encoding="utf-8")


if __name__ == "__main__":
    main()
