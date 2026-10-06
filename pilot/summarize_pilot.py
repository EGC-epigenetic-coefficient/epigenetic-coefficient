from __future__ import annotations

import csv
import hashlib
import json
import statistics
from pathlib import Path


PILOT = Path(__file__).resolve().parent / "output"
RUN = PILOT
RESULTS = RUN / "replication_results.csv"
METHODS = RUN / "method_results.csv"
ERRORS = RUN / "errors.csv"
RUNTIME = RUN / "runtime_metadata.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def summary_row(condition_id: str, truth: str, rows: list[dict[str, str]]) -> dict[str, object]:
    egc_s = [f(row, "egc_s") for row in rows]
    egc_dm = [f(row, "egc_dm_u") for row in rows]
    orientation = [f(row, "orientation_A") for row in rows]
    p_dm = [f(row, "p_dm") for row in rows]
    expected_sign = 0 if truth == "absent" else (1 if truth == "positive" else -1)
    sign_matches = None if expected_sign == 0 else sum((value > 0) == (expected_sign > 0) for value in egc_s)
    orientation_matches = None if expected_sign == 0 else sum((value > 0) == (expected_sign > 0) for value in orientation)
    return {
        "condition_id": condition_id,
        "truth_class": truth,
        "replications": len(rows),
        "egc_s_mean": statistics.mean(egc_s),
        "egc_s_median": statistics.median(egc_s),
        "egc_s_min": min(egc_s),
        "egc_s_max": max(egc_s),
        "egc_dm_u_mean": statistics.mean(egc_dm),
        "egc_dm_u_median": statistics.median(egc_dm),
        "egc_dm_u_min": min(egc_dm),
        "egc_dm_u_max": max(egc_dm),
        "p_dm_lt_0_05": sum(value < 0.05 for value in p_dm),
        "orientation_gate_true": sum(row["orientation_gate"].upper() == "TRUE" for row in rows),
        "egc_s_sign_matches_truth": sign_matches,
        "orientation_A_sign_matches_truth": orientation_matches,
        "mean_directional_coherence": statistics.mean(f(row, "directional_coherence") for row in rows),
        "mean_sign_stability_A": statistics.mean(f(row, "sign_stability_A") for row in rows),
        "max_bootstrap_failure_fraction": max(f(row, "bootstrap_failure_fraction") for row in rows),
    }


def main() -> None:
    replication = read_csv(RESULTS)
    methods = read_csv(METHODS)
    errors = read_csv(ERRORS)
    runtime = {row["key"]: row["value"] for row in read_csv(RUNTIME)}
    truth_by_condition = {
        "C0_CORE_P1_NULL": "absent",
        "C1_CORE_P1_D40": "positive",
        "C2_CORE_P1_D40": "negative",
    }

    summaries = []
    for condition_id, truth in truth_by_condition.items():
        rows = [row for row in replication if row["condition_id"] == condition_id]
        summaries.append(summary_row(condition_id, truth, rows))

    null = summaries[0]
    positive = summaries[1]
    negative = summaries[2]
    elapsed = float(runtime["elapsed_seconds"])
    checks = {
        "replications_15_of_15": len(replication) == 15 and len({row["job_id"] for row in replication}) == 15,
        "five_replications_per_condition": all(row["replications"] == 5 for row in summaries),
        "errors_zero": len(errors) == 0,
        "method_rows_present": len(methods) == 390,
        "unique_rng_streams": len({row["stream_signature"] for row in replication}) == 15,
        "runtime_within_requested_window": 300 <= elapsed <= 450,
        "null_primary_test_not_significant_5_of_5": null["p_dm_lt_0_05"] == 0,
        "null_orientation_gate_closed_5_of_5": null["orientation_gate_true"] == 0,
        "positive_primary_test_significant_5_of_5": positive["p_dm_lt_0_05"] == 5,
        "negative_primary_test_significant_5_of_5": negative["p_dm_lt_0_05"] == 5,
        "positive_signed_egc_correct_5_of_5": positive["egc_s_sign_matches_truth"] == 5,
        "negative_signed_egc_correct_5_of_5": negative["egc_s_sign_matches_truth"] == 5,
        "positive_orientation_correct_5_of_5": positive["orientation_A_sign_matches_truth"] == 5,
        "negative_orientation_correct_5_of_5": negative["orientation_A_sign_matches_truth"] == 5,
        "alternative_magnitude_separated_from_null": (
            null["egc_dm_u_max"] < positive["egc_dm_u_min"]
            and null["egc_dm_u_max"] < negative["egc_dm_u_min"]
        ),
        "no_bootstrap_failures": all(row["max_bootstrap_failure_fraction"] == 0 for row in summaries),
        "no_critical_support_flags": all(row["support_critical"].upper() == "FALSE" for row in replication),
        "all_extrapolations_admissible": all(row["extrapolation_admissible"].upper() == "TRUE" for row in replication),
        "runner_exclusion_flag_preserved": all(row["dry_run_excluded"].upper() == "TRUE" for row in replication),
    }

    summary_path = PILOT / "pilot_condition_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    audit = {
        "status": "PASS_EXPLORATORY_SANITY_CHECK" if all(checks.values()) else "REVIEW_REQUIRED",
        "scope": "exploratory sample pilot; not definitive validation",
        "elapsed_seconds": elapsed,
        "elapsed_minutes": elapsed / 60,
        "jobs_completed": len(replication),
        "errors": len(errors),
        "summaries": summaries,
        "checks": checks,
        "interpretation": {
            "result": "The EGC showed coherent discrimination and direction in this small controlled sample.",
            "what_is_supported": "Proceeding to broader validation is justified.",
            "what_is_not_supported": "The pilot does not establish final power, type-I error, robustness across clocks/cohorts, or readiness as an SDG KPI.",
        },
        "file_sha256": {
            path.name: sha256(path)
            for path in sorted(RUN.iterdir())
            if path.is_file()
        },
    }
    audit_path = PILOT / "pilot_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise SystemExit("pilot audit failed: " + ", ".join(failed))
    print(json.dumps({"status": audit["status"], "elapsed_minutes": audit["elapsed_minutes"], "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
