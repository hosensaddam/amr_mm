"""Write README.md for Updated_Research_Analysis from the final tables (numbers are read, not typed)."""
import sys, json
from pathlib import Path
import pandas as pd
T = Path(sys.argv[1]); out = Path(sys.argv[2]); man = json.load(open(sys.argv[3]))
rd = lambda n: pd.read_csv(T / f"{n}.csv")
AC = rd("AC00_audit_issue_resolution_register"); F = rd("R03u_updated_findings_register"); R6 = rd("R06u_updated_figure_catalog"); qc = rd("DH05_sequencing_qc_all_runs")
lines = ["# Updated_Research_Analysis — audit corrections, Dhaka reprocessing and comparative analysis",
 "", "**Chapter:** Global Insights of the Urban Microbiome and AMR Landscape from Metagenomic Mapping (Academic Press, Elsevier)  ",
 "**Update date:** 31 October 2026 · **Basis:** Independent Research Audit (CONDITIONAL PASS) and its issue register  ",
 "**Original data and previous outputs:** unchanged — all 661 original files re-verified against their pre-update SHA-256 checksums at the end of the update (list and snapshot in `00_Backup_and_Checksums/`).", "",
 "## Main deliverables", "",
 "| File | Content |", "|---|---|",
 f"| `AMR_Updated_Final_Analysis.xlsx` (also copied to `..\\AMR_Updated_Final_Analysis.xlsx`) | Corrected global analysis + audit corrections + integration + Dhaka summary ({man['master']['sheets']} sheets, {man['master']['tables']} tables, {man['master']['charts']} native charts) |",
 f"| `Dhaka_City_Advanced_Analysis.xlsx` | Dedicated Dhaka workbook ({man['dhaka']['sheets']} sheets incl. README, {man['dhaka']['tables']} tables, Figures 31-45, {man['dhaka']['charts']} native charts) |",
 "| `Updated_Audit_Issue_Register.xlsx` | All 20 audit issues with update status, correction, before/after values and evidence |",
 "| `Updated_Statistical_Verification.xlsx` | Independent recalculation (V05), workbook checks (V06), raw-read reproducibility (V07), separate-agent recalculation (V08), internal checks (V04) |",
 "| `03_Notebooks/verification/V07_raw_read_reproducibility_check.csv` | End-to-end re-processing of one Dhaka run (SRR39037538) reproduced the stored outputs byte for byte |",
 "| `03_Notebooks/verification/V08_independent_agent_checks.csv` | Recalculation by a separate verification agent from raw per-run outputs, without access to the analysis code |",
 "| `07_Reports/Dhaka_Comparative_Analysis_Report.docx` | Dhaka report (methods, results, tables, figures, limitations) |",
 "| `07_Reports/Updated_Final_Verification_Report.docx` | Final verification record |", "",
 "## Folder layout", "",
 "| Folder | Content |", "|---|---|",
 "| `00_Backup_and_Checksums/` | SHA-256 of 661 original files; compressed snapshot of the previous outputs; copy of the harmonized inputs |",
 "| `01_Data/` | Updated harmonized sample master (+ cell-level change log: 16 cells changed), DS02 site map from ENA, ENA manifest of all processed runs, BioSample/BioProject records, boundary subset, reprocessed tidy tables |",
 "| `02_Scripts/` | Analysis helpers (`amr_common.py`, `dhaka_common.py`), notebook builder/runner, workbook builder, independent verifier, raw-read pipeline shell scripts, output collector, patch scripts |",
 "| `03_Notebooks/` | Executed notebooks 01-20 (fresh kernels) and their percent-format sources |",
 f"| `04_Tables/` | All analytical tables as CSV ({len(list(T.glob('*.csv')))} files) |",
 f"| `05_Figures/` | {len(R6)} figures: PNG 300 dpi (RGB) and `tiff/` LZW TIFF 600 dpi |",
 "| `06_Logs/` | Notebook run logs, environment, per-run processing logs in `reprocessing/` (fastp JSON and logs, Kraken 2 reports, Bracken tables, KMA .res/.mapstat), tool versions |",
 "| `07_Reports/` | Word reports |", "",
 "## Audit issues", "", "| Issue | Status | After |", "|---|---|---|"] + [f"| {r.issue} | {r.status} | {str(r.after)[:160]} |" for r in AC.itertuples()] + [
 "", "## Dhaka reprocessing (uniform workflow)", "",
 f"* Runs processed: {int((qc.processing_status == 'completed').sum())} (Dhaka {int(((qc.analysis_group == 'Dhaka') & (qc.processing_status == 'completed')).sum())}/16; comparators and calibration as listed in DH05).",
 "* Acquisition: first 4,000,000 read pairs of each ENA FASTQ pair (equal depth); fastp 1.4.0; Kraken 2.17.1 + SILVA 138 16S index (Bracken only as sensitivity); KMA 1.7.1 + ResFinder DB (tag resfinder-4.7.2).",
 "* Tools and databases were installed only in the isolated cloud processing workspace, not on the authors' computer; statistics, notebooks and workbooks ran in the local Python 3.12.10 environment.",
 "* Normalized ARG abundance = FPKM-16S (per million prokaryotic 16S fragments) or per million QC read pairs — relative measures, not absolute ARG load.", "",
 "## Findings (see R03u for denominators, evidence and limitations)", ""] + [f"* **{r.id}** — {r.finding}" for r in F.itertuples()] + [
 "", "## How to reproduce", "",
 "1. Windows, Python 3.12 with the packages in `06_Logs/environment.json`.",
 "2. Set `AMR_ROOT=C:\\Users\\USER\\Documents\\Research\\AMR` and `AMR_HARM=<path to 01_Data\\harmonized_updated>`; run `python 02_Scripts\\run_all.py` (rebuilds and executes notebooks 01-20 in fresh kernels). Windows Controlled Folder Access blocks Python writes inside Documents, so run from a working copy outside Documents and copy results back.",
 "3. `python build_updated_workbooks.py <out_dir>`, then `python verify_updated.py <out_dir>`, then rebuild the workbooks so that V05/V06 are embedded.",
 "4. Raw-read reprocessing: `02_Scripts/pipeline/process_sample_v2.sh RUN URL1 URL2 4000000` per run (Linux; fastp, Kraken 2, Bracken, KMA, ResFinder DB and the SILVA 16S index as in `06_Logs/reprocessing/versions.json`), then `collect_outputs.py`.", "",
 "Prepared with AI assistance (Claude) for author review; authors remain responsible for verification and interpretation."]
out.write_text("\n".join(lines), encoding="utf-8"); print("README written", len(lines))
