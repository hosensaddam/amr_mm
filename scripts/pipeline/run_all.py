"""Master script: reproduce the complete supplementary-material workflow in order.

    python run_all.py            # full workflow (needs the original research directory as SOURCE)
    python run_all.py --verify   # recompute headline results from the cleaned data (works from the public repository)

SOURCE defaults to the parent of the package folder; set AMR_SOURCE_DIR to point elsewhere.
Step 07 needs Microsoft Word on Windows (field update + PDF export).
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FULL = ["00_archive_sources.py", "01_build_inventory.py", "02_clean_datasets.py", "03_verify_analysis.py",
        "06_build_supplementary_docx.py", "07_update_fields_and_export_pdf.ps1", "04_export_tables.py",
        "05_export_figures.py", "08_build_reports.py", "09_quality_control.py"]
VERIFY = ["03_verify_analysis.py"]


def run(step):
    if step.endswith(".ps1"):
        cmd = ["powershell", "-ExecutionPolicy", "Bypass", "-File", str(HERE / step)]
    else:
        cmd = [sys.executable, str(HERE / step)]
    print(f"\n=== {step} ===", flush=True)
    r = subprocess.run(cmd, cwd=HERE)
    if r.returncode != 0:
        sys.exit(f"Step {step} failed (exit code {r.returncode}); workflow stopped.")


if __name__ == "__main__":
    for s in (VERIFY if "--verify" in sys.argv else FULL):
        run(s)
    print("\nWorkflow complete.")
