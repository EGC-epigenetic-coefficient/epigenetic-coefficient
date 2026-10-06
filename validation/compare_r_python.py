"""Compare frozen R and Python EGC results and run independent checks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_RUN = HERE / "results"
CODE_FILES = [
    HERE / "generate_fixed_datasets.R",
    ROOT / "R/egc_reference.R",
    HERE / "evaluate_fixed_datasets_R.R",
    ROOT / "python/egc_reference.py",
    HERE / "evaluate_fixed_datasets_python.py",
    HERE / "compare_r_python.py",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def value(text: str) -> float | None:
    return None if text == "" else float(text)


def close(observed: float, expected: float, atol: float, rtol: float) -> bool:
    return abs(observed - expected) <= atol + rtol * abs(expected)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN)
    args = parser.parse_args()
    run = args.run_dir
    config = json.loads((HERE / "config_frozen.json").read_text(encoding="utf-8"))
    atol = float(config["absolute_tolerance"])
    rtol = float(config["relative_tolerance"])

    dataset_rows = read_rows(HERE / "fixed_datasets_combined.csv")
    manifest = read_rows(HERE / "dataset_manifest.csv")
    dataset_ids = list(dict.fromkeys(row["dataset_id"] for row in dataset_rows))
    row_keys = [(row["dataset_id"], row["row_id"]) for row in dataset_rows]
    group_rank_values: dict[tuple[str, str], set[float]] = {}
    for row in dataset_rows:
        group_rank_values.setdefault((row["dataset_id"], row["group"]), set()).add(float(row["rank"]))
    data_checks = {
        "dataset_count_is_20": len(dataset_ids) == config["dataset_count"] == 20,
        "manifest_count_is_20": len(manifest) == 20,
        "row_keys_unique": len(row_keys) == len(set(row_keys)),
        "burden_finite": all(math.isfinite(float(row["burden"])) for row in dataset_rows),
        "weights_positive_finite": all(math.isfinite(float(row["weight"])) and float(row["weight"]) > 0 for row in dataset_rows),
        "one_rank_per_group": all(len(ranks) == 1 for ranks in group_rank_values.values()),
        "at_least_two_ranks_per_dataset": all(
            len({float(row["rank"]) for row in dataset_rows if row["dataset_id"] == dataset_id}) >= 2
            for dataset_id in dataset_ids
        ),
    }

    r_components = {row["dataset_id"]: row for row in read_rows(run / "r_components.csv")}
    py_components = {row["dataset_id"]: row for row in read_rows(run / "python_components.csv")}
    comparison_rows: list[dict[str, object]] = []
    defined_status_matches = 0
    for dataset_id in dataset_ids:
        for metric in config["comparison_components"]:
            r_value = value(r_components[dataset_id][metric])
            py_value = value(py_components[dataset_id][metric])
            status_match = (r_value is None) == (py_value is None)
            defined_status_matches += int(status_match)
            if r_value is None or py_value is None:
                abs_difference = None
                rel_difference = None
                passed = status_match
            else:
                abs_difference = abs(r_value - py_value)
                rel_difference = abs_difference / max(abs(py_value), atol)
                passed = close(r_value, py_value, atol, rtol)
            comparison_rows.append(
                {
                    "level": "component",
                    "dataset_id": dataset_id,
                    "lower_group": "",
                    "higher_group": "",
                    "metric": metric,
                    "r_value": r_value,
                    "python_value": py_value,
                    "abs_difference": abs_difference,
                    "relative_difference": rel_difference,
                    "defined_status_match": status_match,
                    "passed": passed,
                }
            )

    pair_key = lambda row: (row["dataset_id"], row["lower_group"], row["higher_group"])
    r_pairwise = {pair_key(row): row for row in read_rows(run / "r_pairwise.csv")}
    py_pairwise = {pair_key(row): row for row in read_rows(run / "python_pairwise.csv")}
    pair_keys_match = set(r_pairwise) == set(py_pairwise)
    for key in sorted(set(r_pairwise) & set(py_pairwise)):
        for metric in config["comparison_pairwise"]:
            r_value = value(r_pairwise[key][metric])
            py_value = value(py_pairwise[key][metric])
            abs_difference = abs(float(r_value) - float(py_value))
            rel_difference = abs_difference / max(abs(float(py_value)), atol)
            comparison_rows.append(
                {
                    "level": "pairwise",
                    "dataset_id": key[0],
                    "lower_group": key[1],
                    "higher_group": key[2],
                    "metric": metric,
                    "r_value": r_value,
                    "python_value": py_value,
                    "abs_difference": abs_difference,
                    "relative_difference": rel_difference,
                    "defined_status_match": True,
                    "passed": close(float(r_value), float(py_value), atol, rtol),
                }
            )

    def py(dataset_id: str, metric: str) -> float | None:
        return value(py_components[dataset_id][metric])

    independent_checks: dict[str, bool] = {}
    expected_ds01 = {"egc_s": 4.0, "egc_lm": 1.0, "egc_ld": 1.0, "egc_dm_u": 1.0, "egc_dm_r": 4.0, "egc_dd": 1.0, "mean_burden": 0.0}
    for metric, expected in expected_ds01.items():
        independent_checks[f"DS01_{metric}_exact"] = py("DS01", metric) is not None and close(float(py("DS01", metric)), expected, atol, rtol)
    for metric in ("egc_s", "egc_lm", "egc_dm_u", "egc_dm_r", "mean_burden"):
        independent_checks[f"DS02_{metric}_zero"] = py("DS02", metric) is not None and abs(float(py("DS02", metric))) <= atol
    independent_checks["DS02_egc_ld_undefined"] = py("DS02", "egc_ld") is None
    independent_checks["DS02_egc_dd_undefined"] = py("DS02", "egc_dd") is None
    for metric in ("egc_s", "egc_lm", "egc_ld", "egc_dm_u", "egc_dm_r", "egc_dd"):
        independent_checks[f"DS11_translation_{metric}"] = (
            py("DS03", metric) is None and py("DS11", metric) is None
        ) or (
            py("DS03", metric) is not None
            and py("DS11", metric) is not None
            and close(float(py("DS11", metric)), float(py("DS03", metric)), atol, rtol)
        )
    independent_checks["DS11_mean_shift_minus_5"] = close(float(py("DS11", "mean_burden")), float(py("DS03", "mean_burden")) - 5.0, atol, rtol)
    for metric in ("egc_s", "egc_lm", "egc_dm_u", "egc_dm_r", "mean_burden"):
        independent_checks[f"DS12_scale_{metric}"] = close(float(py("DS12", metric)), 0.08 * float(py("DS03", metric)), atol, rtol)
    for metric in ("egc_ld", "egc_dd"):
        independent_checks[f"DS12_dimensionless_{metric}"] = close(float(py("DS12", metric)), float(py("DS03", metric)), atol, rtol)
    independent_checks["DS08_s_zero"] = abs(float(py("DS08", "egc_s"))) <= 1e-12
    independent_checks["DS08_lm_zero"] = abs(float(py("DS08", "egc_lm"))) <= 1e-12
    independent_checks["DS08_dm_u_positive"] = float(py("DS08", "egc_dm_u")) > 0
    independent_checks["DS08_dd_undefined"] = py("DS08", "egc_dd") is None

    comparison_columns = [
        "level", "dataset_id", "lower_group", "higher_group", "metric", "r_value", "python_value",
        "abs_difference", "relative_difference", "defined_status_match", "passed",
    ]
    with (run / "comparison_long.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=comparison_columns)
        writer.writeheader()
        for row in comparison_rows:
            writer.writerow({key: ("" if value_ is None else format(value_, ".17g") if isinstance(value_, float) else value_) for key, value_ in row.items()})

    defined_differences = [float(row["abs_difference"]) for row in comparison_rows if row["abs_difference"] is not None]
    failed = [row for row in comparison_rows if not row["passed"]]
    maximum = max((row for row in comparison_rows if row["abs_difference"] is not None), key=lambda row: float(row["abs_difference"]))
    summary = {
        "comparison": {
            "absolute_tolerance": atol,
            "relative_tolerance": rtol,
            "component_comparisons": sum(row["level"] == "component" for row in comparison_rows),
            "pairwise_comparisons": sum(row["level"] == "pairwise" for row in comparison_rows),
            "passed": len(comparison_rows) - len(failed),
            "failed": len(failed),
            "defined_status_matches": defined_status_matches,
            "pair_keys_match": pair_keys_match,
            "max_abs_difference": max(defined_differences),
            "max_abs_difference_location": {
                key: maximum[key] for key in ("level", "dataset_id", "lower_group", "higher_group", "metric")
            },
        },
        "data_checks": data_checks,
        "independent_checks": independent_checks,
        "overall_pass": not failed and pair_keys_match and all(data_checks.values()) and all(independent_checks.values()),
    }
    (run / "comparison_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    audit_files = [
        HERE / "config_frozen.json",
        HERE / "dataset_manifest.csv",
        HERE / "fixed_datasets_combined.csv",
        run / "r_components.csv",
        run / "r_pairwise.csv",
        run / "python_components.csv",
        run / "python_pairwise.csv",
        run / "comparison_long.csv",
        run / "comparison_summary.json",
    ]
    hashes = {
        "code": {path.name: sha256(path) for path in CODE_FILES},
        "outputs": {path.name: sha256(path) for path in audit_files},
    }
    (run / "hash_manifest.json").write_text(json.dumps(hashes, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not summary["overall_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
