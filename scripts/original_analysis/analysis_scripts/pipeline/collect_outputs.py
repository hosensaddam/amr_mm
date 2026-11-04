"""Collect per-run outputs of the uniform reprocessing pipeline into tidy tables.
Inputs: runs/<RUN>/{acquire.tsv, fastp.json, kraken2.report, bracken_genus.tsv, kma_resfinder.res, kma_resfinder.mapstat}
Outputs (reprocessed/): run_qc.csv, taxonomy_kraken_long.csv, bracken_genus_long.csv, resfinder_long.csv, versions.json"""
import json, re, sys
from pathlib import Path
import pandas as pd
W = Path("/home/claude/work"); R = W / "runs"; O = W / "reprocessed"; O.mkdir(exist_ok=True)
DB = Path("/home/claude/tools/db/resfinder_db")
man = pd.read_csv(W / "ena/manifest_ena_filereport.tsv", sep="\t")
ph = pd.read_csv(DB / "phenotypes.txt", sep="\t", dtype=str)
ph = ph.rename(columns={"Gene_accession no.": "template"})[["template", "Class", "Phenotype", "Mechanism of resistance"]]
qc_rows, tax_rows, br_rows, rf_rows = [], [], [], []
for run in man.run_accession:
    d = R / run
    if not (d / "DONE").exists():
        qc_rows.append(dict(run_accession=run, processing_status="not completed")); continue
    acq = dict(l.rstrip("\n").split("\t", 1) for l in open(d / "acquire.tsv") if "\t" in l)
    j = json.load(open(d / "fastp.json")); b, a = j["summary"]["before_filtering"], j["summary"]["after_filtering"]
    fr = j["filtering_result"]
    # Kraken 2 report (paired: counts are read pairs = fragments)
    rep = []
    stack = {}
    for line in open(d / "kraken2.report"):
        pct, clade, direct, rank, taxid, name = line.rstrip("\n").split("\t")
        depth = (len(name) - len(name.lstrip(" "))) // 2; nm = name.strip()
        stack[depth] = (rank, nm); stack = {k: v for k, v in stack.items() if k <= depth}
        lineage = ";".join(f"{v[0]}__{v[1]}" for k, v in sorted(stack.items()) if v[0] in ("D", "P", "C", "O", "F", "G", "S"))
        rep.append(dict(rank=rank, taxid=int(taxid), name=nm, clade_fragments=int(clade), direct_fragments=int(direct), lineage=lineage))
    rep = pd.DataFrame(rep)
    tot_pairs = int(a["total_reads"]) // 2
    unclass = int(rep.loc[rep["rank"] == "U", "clade_fragments"].sum())
    dom = rep[rep["rank"] == "D"].set_index("name").clade_fragments
    for _, r in rep[rep["rank"].isin(["D", "P", "C", "O", "F", "G", "S"])].iterrows():
        tax_rows.append(dict(run_accession=run, **r.to_dict()))
    br = pd.read_csv(d / "bracken_genus.tsv", sep="\t")
    for _, r in br.iterrows():
        br_rows.append(dict(run_accession=run, genus=r["name"], taxid=int(r.taxonomy_id), kraken_assigned=int(r.kraken_assigned_reads), bracken_est_fragments=int(r.new_est_reads), bracken_fraction=float(r.fraction_total_reads)))
    # KMA ResFinder
    ms = pd.read_csv(d / "kma_resfinder.mapstat", sep="\t", skiprows=7)
    ms = ms.rename(columns={"# refSequence": "template"})
    res = pd.read_csv(d / "kma_resfinder.res", sep="\t").rename(columns={"#Template": "template"})
    res["template"] = res.template.str.strip()
    ms = ms.merge(res[["template", "Template_length", "Template_Identity", "Template_Coverage", "Query_Identity", "Query_Coverage", "Depth", "p_value"]], on="template", how="left")
    for _, r in ms.iterrows():
        rf_rows.append(dict(run_accession=run, template=r.template, readCount=int(r.readCount), fragmentCount=int(r.fragmentCount), readCountAln=int(r.readCountAln), fragmentCountAln=int(r.fragmentCountAln), refCoveredPositions=int(r.refCoveredPositions),
                            bpTotal=int(r.bpTotal), depthMax=int(r.depthMax), snpSum=int(r.snpSum), in_kma_res=pd.notna(r.Template_Identity), Template_length=r.Template_length,
                            Template_Identity=r.Template_Identity, Template_Coverage=r.Template_Coverage, Query_Identity=r.Query_Identity, Depth=r.Depth, kma_p_value=r.p_value))
    kma_frag = int([l for l in open(d / "kma_resfinder.mapstat") if l.startswith("## fragmentCount")][0].split("\t")[1])
    qc_rows.append(dict(run_accession=run, processing_status="completed", start_utc=open(d / "start.txt").read().strip(), end_utc=open(d / "end.txt").read().strip(),
                        requested_pairs=int(acq["requested_pairs"]), acquired_pairs=int(acq["raw_pairs_acquired_R1"]),
                        raw_bases_acquired=int(b["total_bases"]), read1_mean_length_raw=b["read1_mean_length"], read2_mean_length_raw=b["read2_mean_length"],
                        q30_rate_raw=b["q30_rate"], gc_content_raw=b["gc_content"], qc_reads=int(a["total_reads"]), qc_pairs=tot_pairs, qc_bases=int(a["total_bases"]),
                        q30_rate_qc=a["q30_rate"], read1_mean_length_qc=a["read1_mean_length"], pct_reads_passed=100 * int(a["total_reads"]) / int(b["total_reads"]),
                        low_quality_reads=fr.get("low_quality_reads"), too_short_reads=fr.get("too_short_reads"), too_many_N_reads=fr.get("too_many_N_reads"),
                        adapter_trimmed_reads=j.get("adapter_cutting", {}).get("adapter_trimmed_reads"), duplication_rate=j["duplication"]["rate"],
                        insert_size_peak=j.get("insert_size", {}).get("peak"), kraken_classified_pairs=tot_pairs - unclass,
                        pct_pairs_rRNA_classified=100 * (tot_pairs - unclass) / tot_pairs,
                        bacteria_16S_fragments=int(dom.get("Bacteria", 0)), archaea_16S_fragments=int(dom.get("Archaea", 0)), eukaryota_18S_fragments=int(dom.get("Eukaryota", 0)),
                        kma_fragments_input=kma_frag, resfinder_aligned_fragments=int(ms.fragmentCountAln.sum()), resfinder_templates_any=int(len(ms))))
qc = pd.DataFrame(qc_rows).merge(man, on="run_accession", how="left")
qc.to_csv(O / "run_qc.csv", index=False)
pd.DataFrame(tax_rows).to_csv(O / "taxonomy_kraken_long.csv", index=False)
pd.DataFrame(br_rows).to_csv(O / "bracken_genus_long.csv", index=False)
rf = pd.DataFrame(rf_rows).merge(ph, on="template", how="left")
rf.to_csv(O / "resfinder_long.csv", index=False)
ver = dict(fastp="1.4.0 (opengene.org prebuilt binary)", kraken2="2.17.1 (git tag v2.17.1)", bracken="Bracken git 8732ef5 (est_abundance.py; -l G -t 10, 150-mer distribution)",
           kma="1.7.1 (bitbucket genomicepidemiology/kma, commit b85e02a)", resfinder_db="resfinder_db git tag resfinder-4.7.2 (VERSION 2.5.1, commit 852f367); all.fsa 3,180 templates; SHA-256 2b20a863...",
           kraken2_db="16S_SILVA138_k2db (genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz; SHA-256 edc7dcf6...)",
           acquisition="First 4,000,000 read pairs of each ENA FASTQ pair streamed over HTTPS (head-of-file subsample; identical depth for every run)",
           fastp_params="--detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 --qualified_quality_phred 15 --length_required 50 --trim_poly_g -w 2",
           kraken2_params="--paired --confidence 0.1", kma_params="-ipe -mem_mode -ef -1t1 -nf -t 2")
json.dump(ver, open(O / "versions.json", "w"), indent=1)
print(qc.processing_status.value_counts().to_dict(), len(tax_rows), len(br_rows), len(rf))
