#!/usr/bin/env bash
cd /home/claude/work
xargs -P 3 -L 1 bash -c 'pipeline/process_sample.sh "$0" "$1" "$2" "$3" > runs/log_$0.txt 2>&1 || echo "FAIL $0" >> runs/FAILED.txt' < pipeline/jobs.txt
echo ALLDONE > runs/ALLDONE
