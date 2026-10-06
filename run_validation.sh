#!/usr/bin/env bash
# Reproduces the cross-language verification: 20 fixed datasets, 1235 R-Python comparisons.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-python3}"
# 1. Regenerate the 20 datasets from the master seed and check that they are byte-identical to the frozen copies.
Rscript "$HERE/validation/generate_fixed_datasets.R" "$HERE/validation/results/regenerated"
cmp "$HERE/validation/fixed_datasets_combined.csv" "$HERE/validation/results/regenerated/fixed_datasets_combined.csv"
echo "Regenerated datasets are byte-identical to validation/fixed_datasets_combined.csv"
# 2. Evaluate the datasets independently in R and in Python, then compare.
Rscript "$HERE/validation/evaluate_fixed_datasets_R.R"
"$PY" "$HERE/validation/evaluate_fixed_datasets_python.py"
"$PY" "$HERE/validation/compare_r_python.py" > /dev/null
"$PY" - "$HERE" <<'PYEOF'
import json, sys
from pathlib import Path
s = json.loads((Path(sys.argv[1]) / "validation/results/comparison_summary.json").read_text())
c = s["comparison"]
print(f"R-Python comparisons passed: {c['passed']}/{c['component_comparisons'] + c['pairwise_comparisons']} "
      f"(failed: {c['failed']}; max absolute difference: {c['max_abs_difference']:.3g})")
PYEOF
