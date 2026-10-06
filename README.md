# Epigenetic Coefficient (EGC) — v1.0.0

Reference implementation, cross-language verification suite and exploratory pilot for the
**Epigenetic Coefficient (EGC)**, a distribution-sensitive framework for measuring inequality in
epigenetic age acceleration (EAA) across ordered socioeconomic strata.

This repository accompanies the article:

> Marigliano P., Capunzo M., Aliberti S.M. *Beyond the Socioeconomic Gradient: The Epigenetic
> Coefficient (EGC) as a Distribution-Sensitive Framework for Measuring Inequality in Epigenetic
> Ageing.* Life (under review).

Release **v1.0.0** is the version cited in the article. No individual-level human data are included:
all datasets are controlled or simulated.

## The EGC profile

For *G* ordered socioeconomic strata with population shares *p_g*, the EGC is reported as a profile
**{EGC-M, EGC-O, EGC-S; L}**:

| Paper | Code name | Meaning |
|---|---|---|
| EGC-M | `egc_dm_u` | Distributional magnitude: 2 Σ_{g<h} p_g p_h W1(F_g, F_h) (weighted pairwise Wasserstein-1 distances) |
| EGC-MR | `egc_dm_r` | Rank-calibrated magnitude (sensitivity analysis) |
| EGC-O | `egc_dd` | Conditional orientation A/D; undefined (empty) when D = 0 |
| EGC-S | `egc_s` | Rank-gradient companion (−β from the burden-on-rank model) |
| L | `mean_burden` | Population mean burden |
| — | `egc_lm`, `egc_ld` | Auxiliary location diagnostics (dispersion of stratum means; negative burden–rank correlation) |

Sign convention: socioeconomic rank 0 = maximum disadvantage, rank 1 = maximum advantage; higher
burden is worse.

## Repository structure

```
R/egc_reference.R               point estimators (R)
R/egc_inference.R               permutation test, stratified bootstrap, m-out-of-n intervals (R)
python/egc_reference.py         independent point estimators (Python)
examples/minimal/               two-stratum worked example with hand-calculated expected values
validation/                     20 controlled datasets, R and Python evaluators, comparison script
validation/expected/            frozen outputs of the verification reported in the article
pilot/                          exploratory simulation pilot (Table 5 of the article)
pilot/expected/                 frozen outputs of the pilot reported in the article
run_validation.sh               reproduces the cross-language verification
MANIFEST.sha256                 SHA-256 checksums of every file in the release
```

## Requirements

* R ≥ 4.3 (tested with R 4.6.1). Base R only; the `parallel` package ships with R.
* Python ≥ 3.10 (tested with Python 3.12.2) and NumPy (tested with 1.26.4): `pip install -r requirements.txt`.

## Reproducing the results

**1. Minimal example**

```bash
Rscript examples/minimal/run_minimal_example_R.R
python3 examples/minimal/run_minimal_example_python.py
```

Both must match `examples/minimal/minimal_example_expected.csv` (e.g. EGC-M = 0.8889, EGC-S = 2).

**2. Cross-language verification (about 10 seconds)**

```bash
./run_validation.sh          # set PYTHON=/path/to/python3 if needed
```

The script regenerates the 20 datasets from master seed 20260730 (L'Ecuyer-CMRG), checks that they
are byte-identical to the frozen copies, evaluates them independently in R and Python and compares
the results. Expected output:

```
R-Python comparisons passed: 1235/1235 (failed: 0; max absolute difference: 3.29e-14)
```

The tolerance is 1e-10 (absolute and relative). Results are written to `validation/results/` and can
be compared with `validation/expected/`.

**3. Exploratory pilot (about 5 minutes on one core)**

```bash
./pilot/run_pilot.sh
```

Settings: N = 2500 per replication, five replications per condition (null, positive and negative
gradient), 499 permutations, 499 bootstrap resamples, master seed 20260801. The summary is written to
`pilot/output/pilot_condition_summary.csv`. Expected values (Table 5 of the article):

| Condition | Mean EGC-M | Mean EGC-S | Global p < 0.05 | Orientation gate | Correct sign |
|---|---|---|---|---|---|
| Exact null | 0.0697 | +0.0326 | 0/5 | 0/5 | — |
| Positive gradient | 0.1279 | +0.3668 | 5/5 | 5/5 | 5/5 |
| Negative gradient | 0.1515 | −0.4335 | 5/5 | 5/5 | 5/5 |

With five replications per condition the pilot checks descriptive coherence only; it does not
estimate operating characteristics.

## Integrity

Verify the release files with:

```bash
shasum -a 256 -c MANIFEST.sha256
```

## Scope and limitations

The code implements the frozen EGC v1.0 definition. It demonstrates mathematical specification and
numerical reproducibility across software. Inferential calibration, biological validity and
incremental value over established inequality measures have not been established and require
empirical cohort data.

## License

MIT — see `LICENSE`. If you use the code, please cite the article (see `CITATION.cff`).
