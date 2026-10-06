"""Reference implementation of the candidate Epigenetic Coefficient estimands.

This module implements population/sample plug-in point estimators only.  It is
deliberately independent of any cohort-specific data preparation and does not
implement inferential procedures.  The sign convention is fixed throughout:
SES rank 0 is maximum disadvantage, SES rank 1 is maximum advantage, and
higher biological burden is worse.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, Mapping, Sequence

import numpy as np


ArrayLike = Sequence[float] | np.ndarray


@dataclass(frozen=True)
class EGCComponents:
    """Point-estimate components returned by :func:`estimate_components`."""

    egc_s: float
    egc_lm: float
    egc_ld: float | None
    egc_dm_u: float
    egc_dm_r: float
    egc_dd: float | None
    mean_burden: float


def _as_1d(values: Sequence[object] | np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return array


def _validated_inputs(
    burden: ArrayLike,
    groups: Sequence[Hashable] | np.ndarray,
    group_ranks: Mapping[Hashable, float],
    weights: ArrayLike | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[Hashable, float]]:
    z = _as_1d(burden, "burden").astype(float)
    g = _as_1d(groups, "groups")
    if len(z) != len(g) or len(z) == 0:
        raise ValueError("burden and groups must have the same positive length")
    w = np.ones(len(z), dtype=float) if weights is None else _as_1d(weights, "weights").astype(float)
    if len(w) != len(z):
        raise ValueError("weights must have the same length as burden")
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(w)):
        raise ValueError("burden and weights must be finite")
    if np.any(w <= 0):
        raise ValueError("weights must be strictly positive")

    observed = list(dict.fromkeys(g.tolist()))
    ranks: dict[Hashable, float] = {}
    for label in observed:
        if label not in group_ranks:
            raise ValueError(f"missing SES rank for group {label!r}")
        rank = float(group_ranks[label])
        if not np.isfinite(rank) or not 0 <= rank <= 1:
            raise ValueError("every SES rank must be finite and in [0, 1]")
        ranks[label] = rank
    if len(set(ranks.values())) < 2:
        raise ValueError("at least two distinct SES ranks are required")
    return z, g, w, ranks


def weighted_mean(values: ArrayLike, weights: ArrayLike) -> float:
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    return float(np.sum(w * x) / np.sum(w))


def weighted_covariance(x: ArrayLike, y: ArrayLike, weights: ArrayLike) -> float:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    w = np.asarray(weights, dtype=float)
    x_bar = weighted_mean(x_arr, w)
    y_bar = weighted_mean(y_arr, w)
    return float(np.sum(w * (x_arr - x_bar) * (y_arr - y_bar)) / np.sum(w))


def weighted_correlation(x: ArrayLike, y: ArrayLike, weights: ArrayLike, tol: float = 1e-14) -> float | None:
    var_x = weighted_covariance(x, x, weights)
    var_y = weighted_covariance(y, y, weights)
    if var_x <= tol or var_y <= tol:
        return None
    return weighted_covariance(x, y, weights) / float(np.sqrt(var_x * var_y))


def wasserstein_1_transport(
    x: ArrayLike,
    y: ArrayLike,
    x_weights: ArrayLike | None = None,
    y_weights: ArrayLike | None = None,
) -> float:
    """Exact weighted one-dimensional W1 via monotone mass transport."""

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    wx = np.ones(len(x_arr), dtype=float) if x_weights is None else np.asarray(x_weights, dtype=float)
    wy = np.ones(len(y_arr), dtype=float) if y_weights is None else np.asarray(y_weights, dtype=float)
    if len(x_arr) == 0 or len(y_arr) == 0 or len(wx) != len(x_arr) or len(wy) != len(y_arr):
        raise ValueError("W1 inputs must be non-empty with matching weights")
    if np.any(wx <= 0) or np.any(wy <= 0):
        raise ValueError("W1 weights must be strictly positive")
    order_x = np.argsort(x_arr, kind="mergesort")
    order_y = np.argsort(y_arr, kind="mergesort")
    sx = x_arr[order_x]
    sy = y_arr[order_y]
    px = wx[order_x] / np.sum(wx)
    py = wy[order_y] / np.sum(wy)

    i = j = 0
    rem_x = float(px[0])
    rem_y = float(py[0])
    distance = 0.0
    tol = 1e-15
    while i < len(sx) and j < len(sy):
        mass = min(rem_x, rem_y)
        distance += mass * abs(float(sx[i] - sy[j]))
        rem_x -= mass
        rem_y -= mass
        if rem_x <= tol:
            i += 1
            if i < len(sx):
                rem_x = float(px[i])
        if rem_y <= tol:
            j += 1
            if j < len(sy):
                rem_y = float(py[j])
    return float(distance)


def wasserstein_1_cdf(
    x: ArrayLike,
    y: ArrayLike,
    x_weights: ArrayLike | None = None,
    y_weights: ArrayLike | None = None,
) -> float:
    """Independent exact weighted W1 calculation by integrating |F-G|."""

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    wx = np.ones(len(x_arr), dtype=float) if x_weights is None else np.asarray(x_weights, dtype=float)
    wy = np.ones(len(y_arr), dtype=float) if y_weights is None else np.asarray(y_weights, dtype=float)
    if len(x_arr) == 0 or len(y_arr) == 0 or len(wx) != len(x_arr) or len(wy) != len(y_arr):
        raise ValueError("W1 inputs must be non-empty with matching weights")
    if np.any(wx <= 0) or np.any(wy <= 0):
        raise ValueError("W1 weights must be strictly positive")
    wx = wx / np.sum(wx)
    wy = wy / np.sum(wy)
    support = np.unique(np.concatenate((x_arr, y_arr)))
    if len(support) == 1:
        return 0.0
    distance = 0.0
    cdf_x = 0.0
    cdf_y = 0.0
    for index in range(len(support) - 1):
        point = support[index]
        cdf_x += float(np.sum(wx[x_arr == point]))
        cdf_y += float(np.sum(wy[y_arr == point]))
        distance += abs(cdf_x - cdf_y) * float(support[index + 1] - point)
    return float(distance)


def probability_superiority_difference(
    x: ArrayLike,
    y: ArrayLike,
    x_weights: ArrayLike | None = None,
    y_weights: ArrayLike | None = None,
) -> float:
    """Return P(X>Y)-P(X<Y), with ties contributing zero."""

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    wx = np.ones(len(x_arr), dtype=float) if x_weights is None else np.asarray(x_weights, dtype=float)
    wy = np.ones(len(y_arr), dtype=float) if y_weights is None else np.asarray(y_weights, dtype=float)
    wx = wx / np.sum(wx)
    wy = wy / np.sum(wy)
    comparison = np.sign(x_arr[:, None] - y_arr[None, :])
    pair_weights = wx[:, None] * wy[None, :]
    return float(np.sum(pair_weights * comparison))


def estimate_components(
    burden: ArrayLike,
    groups: Sequence[Hashable] | np.ndarray,
    group_ranks: Mapping[Hashable, float],
    weights: ArrayLike | None = None,
) -> EGCComponents:
    """Estimate F1, F2, F3-U/F3-R and the orientation component."""

    z, g, w, ranks = _validated_inputs(burden, groups, group_ranks, weights)
    total_weight = float(np.sum(w))
    individual_ranks = np.asarray([ranks[label] for label in g], dtype=float)

    rank_variance = weighted_covariance(individual_ranks, individual_ranks, w)
    if rank_variance <= 1e-14:
        raise ValueError("SES rank variance must be positive")
    beta = weighted_covariance(individual_ranks, z, w) / rank_variance
    egc_s = -float(beta)

    ordered_groups = sorted(ranks, key=lambda label: ranks[label])
    p: dict[Hashable, float] = {}
    means: dict[Hashable, float] = {}
    values: dict[Hashable, np.ndarray] = {}
    group_weights: dict[Hashable, np.ndarray] = {}
    for label in ordered_groups:
        mask = g == label
        values[label] = z[mask]
        group_weights[label] = w[mask]
        p[label] = float(np.sum(w[mask]) / total_weight)
        means[label] = weighted_mean(z[mask], w[mask])

    mean_burden = weighted_mean(z, w)
    egc_lm = float(np.sqrt(sum(p[label] * (means[label] - mean_burden) ** 2 for label in ordered_groups)))
    group_rank_values = np.asarray([ranks[label] for label in ordered_groups], dtype=float)
    group_mean_values = np.asarray([means[label] for label in ordered_groups], dtype=float)
    group_probabilities = np.asarray([p[label] for label in ordered_groups], dtype=float)
    correlation = weighted_correlation(group_rank_values, group_mean_values, group_probabilities)
    egc_ld = None if correlation is None else -float(correlation)

    dm_u = 0.0
    dm_r_numerator = 0.0
    dm_r_denominator = 0.0
    dd_numerator = 0.0
    dd_denominator = 0.0
    for lower_index, lower in enumerate(ordered_groups[:-1]):
        for higher in ordered_groups[lower_index + 1 :]:
            pair_weight = p[lower] * p[higher]
            rank_distance = ranks[higher] - ranks[lower]
            w1 = wasserstein_1_transport(
                values[lower], values[higher], group_weights[lower], group_weights[higher]
            )
            dm_u += 2.0 * pair_weight * w1
            dm_r_numerator += pair_weight * rank_distance * w1
            dm_r_denominator += pair_weight * rank_distance**2

            delta = probability_superiority_difference(
                values[lower], values[higher], group_weights[lower], group_weights[higher]
            )
            orientation_weight = pair_weight * rank_distance
            dd_numerator += orientation_weight * delta
            dd_denominator += orientation_weight * abs(delta)

    if dm_r_denominator <= 1e-14:
        raise ValueError("rank-distance denominator must be positive")
    egc_dm_r = float(dm_r_numerator / dm_r_denominator)
    egc_dd = None if dd_denominator <= 1e-14 else float(dd_numerator / dd_denominator)
    return EGCComponents(
        egc_s=egc_s,
        egc_lm=egc_lm,
        egc_ld=egc_ld,
        egc_dm_u=float(dm_u),
        egc_dm_r=egc_dm_r,
        egc_dd=egc_dd,
        mean_burden=mean_burden,
    )
