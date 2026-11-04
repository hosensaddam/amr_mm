#!/usr/bin/env bash
cd /home/claude/work
while pgrep -f "xargs -P 3 -L 1 bash -c pipeline/process_sample.sh" > /dev/null || pgrep -f "xargs -P 2 -L 1" > /dev/null; do sleep 20; done
cat pipeline/jobs.txt pipeline/jobs_calib.txt | while read r a b n; do [ -f runs/$r/DONE ] || echo "$r $a $b $n"; done > pipeline/jobs_remaining.txt
echo "remaining: $(wc -l < pipeline/jobs_remaining.txt)" > runs/finish_status.txt
for r in $(cut -d' ' -f1 pipeline/jobs_remaining.txt); do rm -rf scratch/$r; done
xargs -P 3 -L 1 bash -c 'pipeline/process_sample_v2.sh "$0" "$1" "$2" "$3" > runs/log_$0.txt 2>&1 || echo "FAIL3 $0" >> runs/FAILED.txt' < pipeline/jobs_remaining.txt
cat pipeline/jobs.txt pipeline/jobs_calib.txt | while read r a b n; do [ -f runs/$r/DONE ] || echo "$r $a $b $n"; done > pipeline/jobs_remaining2.txt
xargs -P 3 -L 1 bash -c 'pipeline/process_sample_v2.sh "$0" "$1" "$2" "$3" > runs/log_$0.txt 2>&1 || echo "FAIL4 $0" >> runs/FAILED.txt' < pipeline/jobs_remaining2.txt
echo FINALDONE > runs/FINALDONE
