# %% [markdown]
# # 12 · Dhaka City raw-read acquisition and quality control
# **Research objective:** make the 16 included Dhaka metagenomes analytically usable (Objectives 1-2, Dhaka focus) and document the comparator runs processed with the same workflow.
#
# **Input datasets.** (i) ENA read-run reports for every processed run (`derived/reprocessed/run_qc.csv`, retrieved 31 Oct 2026 from the ENA Portal API); (ii) the updated harmonized sample master (sample-level sample-type corrections AUD-04, coordinate provenance AUD-18); (iii) per-run fastp 1.4.0 JSON reports and Kraken 2 classification summaries collected from the processing workspace.
#
# **Processing workflow (identical for every run).** The first 4,000,000 read pairs of each public paired-end FASTQ pair were streamed from ENA (equal analysed depth for all runs; the smallest Dhaka library contains 4.03 M pairs); fastp 1.4.0 trimmed adapters (auto-detected for paired-end), poly-G tails and low-quality ends (sliding mean Q20), filtered reads with > 40 % bases below Q15 or shorter than 50 nt; Kraken 2.17.1 classified read pairs against the SILVA 138 16S index (confidence 0.1) with Bracken genus re-estimation; KMA 1.7.1 aligned read pairs to the ResFinder database (tag resfinder-4.7.2), one template per fragment (`-1t1`). Host (human) reads were not removed: neither the 16S index nor ResFinder contains human sequence, so human reads cannot be counted as bacterial taxa or ARGs; they only reduce the effective prokaryotic depth, which is measured directly (16S fragments).
#
# **Statistical assumptions.** Descriptive only; equal requested depth makes depth-dependent quantities (richness, detection) comparable across runs.
#
# **Limitations.** A head-of-file subsample uses the first clusters written by the instrument; this is not a random subsample of the whole library. Full-file MD5 checksums (ENA) are recorded but could not be verified for partially streamed files; read-pair counts and pairing were verified for every run.

# %%
#BOOTSTRAP
NB = "NB12"
nb_header(NB)
from dhaka_common import *
meta = load_run_meta()
ver = json.loads((REP / "versions.json").read_text())
print(meta.analysis_group.value_counts().to_dict()); print(json.dumps(ver, indent=1))

# %% [markdown]
# ## 1. Dhaka sample and run inventory (ENA/INSDC records)

# %%
dh = meta[meta.analysis_group == "Dhaka"].copy()
assert len(dh) == 16 and dh.Sample_ID.nunique() == 16, "expected 16 Dhaka runs/samples"
assert (dh.processing_status == "completed").all(), dh.processing_status.value_counts()
inv = dh.reset_index()[["Sample_ID", "run_accession", "study_accession", "Dataset_ID", "Project_Accession", "Sample_Type", "Sample_Type_Original_Harmonized", "Sample_Subtype",
                        "Collection_Date", "instrument_platform", "instrument_model", "library_layout", "library_strategy", "read_count", "base_count", "read1_mean_length_raw",
                        "fastq_ftp", "fastq_bytes", "fastq_md5"]].rename(columns={"read_count": "ENA_read_pairs", "base_count": "ENA_total_bases", "read1_mean_length_raw": "read_length_nt_observed"})
inv["ENA_Gbp"] = inv.ENA_total_bases / 1e9
inv["download_URLs"] = "https://" + inv.fastq_ftp.str.replace(";", " ; https://")
inv["fastq_GB_total"] = inv.fastq_bytes.map(lambda s: sum(int(x) for x in str(s).split(";")) / 1e9)
inv = inv.drop(columns=["fastq_ftp", "fastq_bytes"]).sort_values(["Dataset_ID", "Collection_Date"])
save_table(inv.pipe(rnd, 4), "DH01_dhaka_run_inventory", "Dhaka: sample, run and study accessions, sample type, platform, layout, read length, ENA total bases, download URLs, file sizes and ENA MD5 checksums (16 samples)", NB, sheet="Dhaka_Dataset_Inventory")
dsum = inv.groupby("Dataset_ID").agg(samples=("Sample_ID", "size"), projects=("Project_Accession", lambda s: "; ".join(sorted(set(s)))), platform=("instrument_model", "first"),
                                     first_date=("Collection_Date", "min"), last_date=("Collection_Date", "max"), median_Gbp=("ENA_Gbp", "median"), min_Gbp=("ENA_Gbp", "min"), max_Gbp=("ENA_Gbp", "max"),
                                     sample_types=("Sample_Type", lambda s: "; ".join(f"{k} {v}" for k, v in s.value_counts().items()))).reset_index()
save_table(dsum.pipe(rnd, 4), "DH02_dhaka_dataset_summary", "Dhaka: dataset-level inventory (3 datasets, 16 samples)", NB, sheet="Dhaka_Dataset_Inventory")
dsum

# %% [markdown]
# ## 2. Sample-type classification (AUD-04) and metadata

# %%
st = dh.reset_index()[["Sample_ID", "Dataset_ID", "Sample_Type_Original_Harmonized", "Sample_Type", "Sample_Type_Correction", "Sample_Subtype", "Collection_Date", "Latitude", "Longitude",
                       "Coordinate_Precision", "Locality_As_Reported", "City"]]
save_table(st, "DH03_dhaka_sample_metadata", "Dhaka: sample-level metadata after the AUD-04 correction (sample type verified against BioSample isolation_source and BioProject descriptions)", NB, sheet="Dhaka_Sample_Metadata")
tt = pd.crosstab(dh.Sample_Type, dh.Dataset_ID, margins=True, margins_name="Total")
save_table(tt.reset_index(), "DH04_dhaka_sample_type_by_dataset", "Dhaka: verified sample types by dataset (after AUD-04)", NB, sheet="Dhaka_Sample_Type")
print(tt)

# %% [markdown]
# ## 3. Sequencing and read-QC metrics for all processed runs

# %%
qcc = ["acquired_pairs", "qc_pairs", "pct_reads_passed", "q30_rate_raw", "q30_rate_qc", "read1_mean_length_raw", "read1_mean_length_qc", "gc_content_raw", "duplication_rate", "insert_size_peak",
       "adapter_trimmed_reads", "low_quality_reads", "too_short_reads", "kraken_classified_pairs", "pct_pairs_rRNA_classified", "bacteria_16S_fragments", "archaea_16S_fragments", "eukaryota_18S_fragments",
       "resfinder_aligned_fragments", "resfinder_templates_any"]
allqc = meta.reset_index()[["run_accession", "Sample_ID", "analysis_group", "City_group", "Dataset_ID", "Sample_Type", "Collection_Date", "instrument_model", "base_count", "processing_status"] + qcc]
allqc["ENA_Gbp"] = allqc.base_count / 1e9; allqc["fraction_of_library_analysed_pct"] = 100 * allqc.acquired_pairs / meta.read_count.values
allqc["prokaryotic_16S_fragments"] = allqc.bacteria_16S_fragments + allqc.archaea_16S_fragments
allqc["QC_flags"] = [("; ".join(f for f, c in [("<90% reads passed", r.pct_reads_passed < 90), ("Q30 < 80%", r.q30_rate_qc < .8), ("< 2,000 prokaryotic 16S fragments", r.prokaryotic_16S_fragments < 2000),
                                            ("duplication > 30%", r.duplication_rate > .3)] if c) or "none") for r in allqc.itertuples()]
save_table(allqc.drop(columns=["base_count"]).pipe(rnd, 4), "DH05_sequencing_qc_all_runs", "Read acquisition and QC (fastp 1.4.0; Kraken 2 16S classification; KMA ResFinder alignment) for all processed runs: 16 Dhaka and comparator/calibration runs, 4,000,000 read pairs requested per run", NB, sheet="Dhaka_Sequencing_QC")
dq = allqc[allqc.analysis_group == "Dhaka"]
save_table(dq.drop(columns=["base_count"]).pipe(rnd, 4), "DH06_dhaka_sequencing_qc", "Dhaka: read-level QC and depth assessment per sample (4 M read pairs analysed per sample)", NB, sheet="Dhaka_Data_Quality")
chk = pd.DataFrame([dict(check="All 16 Dhaka runs processed", result=int((dq.processing_status == "completed").sum()), expected=16),
                    dict(check="Runs with exactly 4,000,000 pairs acquired (all processed runs)", result=int((allqc.acquired_pairs == 4_000_000).sum()), expected=int((allqc.processing_status == "completed").sum())),
                    dict(check="Dhaka runs with < 90 % reads passing QC", result=int((dq.pct_reads_passed < 90).sum()), expected=0),
                    dict(check="Dhaka runs with < 2,000 prokaryotic 16S fragments", result=int((dq.prokaryotic_16S_fragments < 2000).sum()), expected=0)])
chk["status"] = np.where(chk.result == chk.expected, "PASS", "CHECK")
save_table(chk, "DH07_dhaka_qc_checks", "Dhaka: acquisition and QC verification checks", NB, sheet="Dhaka_Data_Quality")
print(chk); dq[["run_accession", "Dataset_ID", "pct_reads_passed", "q30_rate_qc", "prokaryotic_16S_fragments", "QC_flags"]]

# %%
order = allqc[allqc.processing_status == "completed"].sort_values(["analysis_group", "City_group", "Dataset_ID", "Collection_Date"])
order = pd.concat([order[order.analysis_group == "Dhaka"], order[order.analysis_group != "Dhaka"]])
fig, axes = plt.subplots(4, 1, figsize=(6.9, 6.6), sharex=True)
x = np.arange(len(order)); cols = [CITY_COL.get(c, "#999999") for c in order.City_group]
for ax, (c, lab, sc) in zip(axes, [("pct_reads_passed", "Reads passing QC (%)", 1), ("q30_rate_qc", "Q30 after QC (%)", 100), ("prokaryotic_16S_fragments", "Prokaryotic 16S\nfragments", 1), ("resfinder_aligned_fragments", "ResFinder-aligned\nfragments", 1)]):
    ax.bar(x, order[c] * sc, color=cols, width=.8)
    ax.set_ylabel(lab, fontsize=6.5)
    if c in ("prokaryotic_16S_fragments", "resfinder_aligned_fragments"): ax.set_yscale("log")
axes[0].set_ylim(80, 100); axes[1].set_ylim(70, 100)
axes[-1].set_xticks(x); axes[-1].set_xticklabels([f"{r.City_group} {r.Dataset_ID}" for r in order.itertuples()], rotation=90, fontsize=6)
from matplotlib.patches import Patch
axes[0].legend(handles=[Patch(color=v, label=k) for k, v in CITY_COL.items() if k in set(order.City_group)], frameon=False, ncol=7, fontsize=6, loc="lower left", bbox_to_anchor=(0, 1.0))
fig.tight_layout(h_pad=.4)
save_fig(fig, 31, "reprocessed_read_qc", "Read quality control of the uniformly reprocessed runs (Dhaka and comparators)",
         f"Each bar is one run (4,000,000 read pairs analysed per run): percentage of reads passing fastp quality control, Q30 base rate after QC, prokaryotic (bacterial + archaeal) 16S rRNA fragments classified by Kraken 2 against SILVA 138, and fragments aligned by KMA to ResFinder (log scales). Colours: city group. Dhaka n = {len(dq)}; comparators n = {len(order) - len(dq)}.",
         NB, "Dhaka objective (data usability)", "tables/DH05_sequencing_qc_all_runs.csv")

# %%
record_finding(NB, "DF01", "Dhaka (data usability)", f"All 16 included Dhaka metagenomes (DS01 4, DS03 4, DS84 8) were retrieved from ENA and processed with one workflow at an equal depth of 4,000,000 read pairs; {dq.pct_reads_passed.min():.1f}-{dq.pct_reads_passed.max():.1f}% of reads passed QC and each sample yielded {int(dq.prokaryotic_16S_fragments.min()):,}-{int(dq.prokaryotic_16S_fragments.max()):,} prokaryotic 16S fragments and {int(dq.resfinder_aligned_fragments.min()):,}-{int(dq.resfinder_aligned_fragments.max()):,} ResFinder-aligned fragments.",
               16, "Dhaka samples (16 of 16 processed)", "Descriptive QC (fastp, Kraken 2, KMA)", "Tables DH01-DH07; Figure 31",
               f"Dhaka now has sample-level taxonomic and AMR profiles; analysed depth is {dq.fraction_of_library_analysed_pct.min():.0f}-{dq.fraction_of_library_analysed_pct.max():.0f}% of each library.",
               "Head-of-file subsample; 16S-based taxonomy (genus resolution); host reads not removed (cannot enter 16S or ResFinder counts).", kind="Empirical (new reprocessing)")
