from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("egc_reference_v1", ROOT.parent.parent / "python" / "egc_reference.py")
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = module
spec.loader.exec_module(module)

def run_example(example: str, input_file: str) -> list[dict[str, object]]:
    with (ROOT / input_file).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    burden = [float(row["z_burden"]) for row in rows]
    groups = [row["ses_group"] for row in rows]
    weights = [float(row["weight"]) for row in rows]
    ranks = {row["ses_group"]: float(row["ses_rank"]) for row in rows}
    result = module.estimate_components(burden, groups, ranks, weights)
    values = {
        "egc_s": result.egc_s,
        "egc_lm": result.egc_lm,
        "egc_ld": result.egc_ld,
        "egc_dm_u": result.egc_dm_u,
        "egc_dm_r": result.egc_dm_r,
        "egc_dd": result.egc_dd,
        "mean_burden": result.mean_burden,
    }
    return [
        {"example": example, "metric": metric, "value": value, "implementation": "Python"}
        for metric, value in values.items()
    ]


output = run_example("A", "minimal_example_input.csv") + run_example(
    "B", "minimal_example_B_input.csv"
)
with (ROOT / "minimal_example_results_python.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["example", "metric", "value", "implementation"])
    writer.writeheader()
    writer.writerows(output)

for row in output:
    value = row["value"]
    rendered = "undefined" if value is None else f"{value:.17g}"
    print(f"Example {row['example']} — {row['metric']}: {rendered}")
