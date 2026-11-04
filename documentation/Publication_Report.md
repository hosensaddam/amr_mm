# GitHub Publication Report

Repository: https://github.com/hosensaddam/amr_mm (branch `master`). Prepared on 4 November 2026.

## Approach
* The existing repository content and history were preserved. The five data-collection workbooks and
  `Public_Data_Source_Log.csv` at the repository root were left unchanged.
* The README was rewritten to describe the whole repository. The content of the earlier README is
  retained in its section 8, with the broken backslash escaping repaired.
* New material was added in folders: `data/`, `analysis/`, `tables/`, `figures/`, `supplementary_materials/`,
  `scripts/`, `notebooks/`, `documentation/`, plus `requirements.txt` and `.gitignore`.
* There was no force-push, no history rewrite, no deletion and no change to repository visibility or licence.

## Public-release screen (before staging)
| Check | Result |
|---|---|
| Patient-level / personal data | None. The data are public INSDC/MG-RAST repository metadata and producer-processed AMR/taxonomy tables. |
| E-mail addresses | One institutional licensing address quoted in third-party reuse terms (`license@uw.edu`, inside an archived data-collection script). Not personal. |
| Credentials, API keys, tokens | None found. The only pattern match was a phrase ("multi-token") in a verification table, plus the QC script's own regular expression. |
| Local paths | The archived analysis code and notebooks contain the original working paths (`C:\Users\USER\Documents\Research\AMR`, `/home/claude`, `/mnt/user-data`). These are generic account or sandbox names, not personal names, and are documented as archival. The new workflow scripts use relative paths only. |
| File sizes | Largest file is 23 MB (`DS_D_taxonomy_long.csv`). The 145 MB AMR gene table is published gzip-compressed (8 MB). No file approaches GitHub's 100 MB limit. |
| Third-party copyright | No publisher PDFs or full texts are included, and no chapter manuscript text. Processed data come from CC BY 4.0 Zenodo records and are attributed in the README. |

## Withheld from publication (and why)
| Item | Reason |
|---|---|
| Chapter drafts, abstract and discussion drafts, main-text figures (`08_Final_Submission/`, root drafts) | Elsevier chapter content; publisher copyright. Never copied into the package. |
| Third-party raw downloads (`AMR_Chapter_Reproducibility_Bundle/.../raw_downloads`, including publication full texts) | Third-party copyright; can be re-downloaded from the sources listed in `Public_Data_Source_Log.csv`. |
| `01_Raw_Data_Archive/` | Read-only copies of the originals. They are already public at the repository root, or are superseded by the cleaned datasets. |
| `Source_File_Inventory.xlsx`, `Processing_Log.txt` | Internal file listing of the author's working folder. |
| Text-review working records | Internal; the resulting corrections log is published. |
| `Cleaned_Excel/` | Duplicates of the published CSV files. |
| `Supplementary_Material_Final.pdf` | Not yet generated: Word's PDF output stalled under automation on the build computer. Add it after saving from Word. |

## Licence
The repository has no licence file. None was added, because the licence is the authors' decision. Without one,
default copyright applies. Consider CC BY 4.0 for data and documents and MIT for code.
