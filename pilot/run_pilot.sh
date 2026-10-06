#!/usr/bin/env bash
# Reproduces the exploratory pilot (Table 5): N = 2500, five replications per condition,
# 499 permutations and 499 bootstrap resamples, master seed 20260801 (L'Ecuyer-CMRG).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
rm -rf "$HERE/output"
Rscript "$HERE/run_pilot.R" \
  --registry="$HERE/condition_registry.csv" \
  --output-dir="$HERE/output" \
  --condition-ids=C0_CORE_P1_NULL,C1_CORE_P1_D40,C2_CORE_P1_D40 \
  --sample-sizes=2500 \
  --replications=5 \
  --B-permutation=499 \
  --B-bootstrap=499 \
  --gammas=0.6,0.7,0.8 \
  --K-extrapolation=199 \
  --master-seed=20260801 \
  --cores=1
"${PYTHON:-python3}" "$HERE/summarize_pilot.py"
