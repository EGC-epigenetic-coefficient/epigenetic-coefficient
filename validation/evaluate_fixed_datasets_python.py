"""Evaluate the frozen cross-language EGC datasets with the Python reference."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "python"))

from egc_reference import (  # noqa: E402
    estimate_components,
    probability_superiority_difference,
    wasserstein_1_transport,
)


DEFAULT_INPUT = HERE / "fixed_datasets_combined.csv"
DEFAULT_RUN = HERE / "results"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            encoded = {}
            for column in columns:
                value = row.get(column)
                if value is None:
                    encoded[column] = ""
                elif isinstance(value, float):
                    encoded[column] = format(value, ".17g")
                else:
                    encoded[column] = value
            writer.writerow(encoded)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RUN)
    args = parser.parse_args()
    rows = read_rows(args.input)
    dataset_ids = list(dict.fromkeys(row["dataset_id"] for row in rows))
    components_output: list[dict[str, object]] = []
    pairwise_output: list[dict[str, object]] = []
    for dataset_id in dataset_ids:
        subset = [row for row in rows if row["dataset_id"] == dataset_id]
        burden = np.asarray([float(row["burden"]) for row in subset], dtype=float)
        groups = np.asarray([row["group"] for row in subset], dtype=object)
        weights = np.asarray([float(row["weight"]) for row in subset], dtype=float)
        rank_map: dict[str, float] = {}
        for row in subset:
            rank = float(row["rank"])
            if row["group"] in rank_map and rank_map[row["group"]] != rank:
                raise ValueError(f"inconsistent rank in {dataset_id}")
            rank_map[row["group"]] = rank
        result = estimate_components(burden, groups, rank_map, weights)
        components_output.append(
            {
                "dataset_id": dataset_id,
                "egc_s": result.egc_s,
                "egc_lm": result.egc_lm,
                "egc_ld": result.egc_ld,
                "egc_dm_u": result.egc_dm_u,
                "egc_dm_r": result.egc_dm_r,
                "egc_dd": result.egc_dd,
                "mean_burden": result.mean_burden,
            }
        )
        total_weight = float(np.sum(weights))
        ordered = sorted(rank_map, key=rank_map.get)
        for lower_index, lower in enumerate(ordered[:-1]):
            lower_mask = groups == lower
            for higher in ordered[lower_index + 1 :]:
                higher_mask = groups == higher
                pairwise_output.append(
                    {
                        "dataset_id": dataset_id,
                        "lower_group": lower,
                        "higher_group": higher,
                        "p_lower": float(np.sum(weights[lower_mask]) / total_weight),
                        "p_higher": float(np.sum(weights[higher_mask]) / total_weight),
                        "rank_distance": rank_map[higher] - rank_map[lower],
                        "w1": wasserstein_1_transport(
                            burden[lower_mask], burden[higher_mask], weights[lower_mask], weights[higher_mask]
                        ),
                        "delta": probability_superiority_difference(
                            burden[lower_mask], burden[higher_mask], weights[lower_mask], weights[higher_mask]
                        ),
                    }
                )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    component_columns = ["dataset_id", "egc_s", "egc_lm", "egc_ld", "egc_dm_u", "egc_dm_r", "egc_dd", "mean_burden"]
    pairwise_columns = ["dataset_id", "lower_group", "higher_group", "p_lower", "p_higher", "rank_distance", "w1", "delta"]
    write_rows(args.output_dir / "python_components.csv", components_output, component_columns)
    write_rows(args.output_dir / "python_pairwise.csv", pairwise_output, pairwise_columns)
    print(f"Evaluated {len(dataset_ids)} datasets in Python")


if __name__ == "__main__":
    main()
