#!/usr/bin/env bash
# Uniform metagenomic processing of one paired-end run (Dhaka update, 31 Oct 2026)
# Usage: process_sample.sh RUN_ACC URL_R1 URL_R2 NPAIRS
set -euo pipefail
RUN=$1; U1=$2; U2=$3; N=${4:-4000000}
T=/home/claude/tools; DB=$T/db
OUT=/home/claude/work/runs/$RUN; S=/home/claude/work/scratch/$RUN
mkdir -p $OUT $S
[ -f $OUT/DONE ] && exit 0
L=$((N*4))
date -u +%FT%TZ > $OUT/start.txt
# 1. Acquire the first N read pairs (head-of-file subsample; identical N for all runs)
( curl -sS --retry 5 --retry-delay 10 https://$U1 | gzip -dc 2>/dev/null | head -n $L > $S/r1.fq ) &
( curl -sS --retry 5 --retry-delay 10 https://$U2 | gzip -dc 2>/dev/null | head -n $L > $S/r2.fq ) &
wait
n1=$(($(wc -l < $S/r1.fq)/4)); n2=$(($(wc -l < $S/r2.fq)/4))
echo -e "run\t$RUN\nraw_pairs_acquired_R1\t$n1\nraw_pairs_acquired_R2\t$n2\nrequested_pairs\t$N" > $OUT/acquire.tsv
[ "$n1" -eq "$n2" ] || { echo "pair mismatch" >> $OUT/acquire.tsv; exit 2; }
# 2-3. Adapter trimming, quality trimming/filtering, poly-G trimming (fastp 1.4.0)
$T/fastp -i $S/r1.fq -I $S/r2.fq -o $S/t1.fq -O $S/t2.fq --detect_adapter_for_pe \
  --cut_front --cut_tail --cut_mean_quality 20 --qualified_quality_phred 15 --length_required 50 \
  --trim_poly_g -w 2 -j $OUT/fastp.json -h $OUT/fastp.html 2> $OUT/fastp.log
rm -f $S/r1.fq $S/r2.fq
# 4-6. Taxonomic profiling of rRNA reads: Kraken2 2.17.1 vs SILVA 138 16S index; Bracken genus re-estimation
$T/k2bin/kraken2 --db $DB/k2_silva/16S_SILVA138_k2db --paired --threads 2 --confidence 0.1 \
  --report $OUT/kraken2.report --output - $S/t1.fq $S/t2.fq 2> $OUT/kraken2.log
python3 $T/Bracken/src/est_abundance.py -i $OUT/kraken2.report -k $DB/k2_silva/16S_SILVA138_k2db/database150mers.kmer_distrib \
  -l G -t 10 -o $OUT/bracken_genus.tsv --out-report $OUT/bracken_genus.report > $OUT/bracken.log 2>&1 || true
# 8. AMR genes: KMA 1.7.1 vs ResFinder DB (resfinder-4.7.2 tag), one template per read pair
$T/kma/kma -ipe $S/t1.fq $S/t2.fq -o $OUT/kma_resfinder -t_db $DB/resfinder_kma -mem_mode -ef -1t1 -nf -t 2 > $OUT/kma.log 2>&1
rm -rf $S
date -u +%FT%TZ > $OUT/end.txt
touch $OUT/DONE
