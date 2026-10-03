"""Tests for the Phase 53 canonical bifurcation identifiability methods."""

from __future__ import annotations

import math

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.recovery_scaling_identifiability import (
    CONTROL_SETTINGS,
    TRUE_MARGINS,
    _evaluate_frozen_criteria,
    estimate_exponential,
    estimate_scaling,
    exact_trajectory,
    fit_canonical_trajectory,
    generate_observation,
    local_recovery_rate,
    stable_equilibrium,
)


@pytest.mark.parametrize(
    ("model", "margin", "equilibrium", "rate"),
    [
        ("fold", 0.16, 0.4, 0.8),
        ("transcritical", 0.16, 0.16, 0.16),
        ("supercritical_pitchfork", 0.16, 0.4, 0.32),
    ],
)
def test_stable_branches_and_local_rates(
    model: str, margin: float, equilibrium: float, rate: float
) -> None:
    assert stable_equilibrium(model, margin) == pytest.approx(equilibrium)
    assert local_recovery_rate(model, margin) == pytest.approx(rate)


@pytest.mark.parametrize("model", ["fold", "transcritical", "supercritical_pitchfork"])
def test_exact_trajectory_starts_above_stable_branch_and_relaxes(model: str) -> None:
    margin = 0.16
    equilibrium = stable_equilibrium(model, margin)
    times = np.linspace(0.0, 80.0, 801)
    trajectory = exact_trajectory(model, margin, 1.2 * equilibrium, times)

    assert trajectory[0] == pytest.approx(1.2 * equilibrium)
    assert trajectory[-1] == pytest.approx(equilibrium, rel=1e-6)
    differences = np.diff(trajectory)
    assert np.all(differences <= 0.0)
    assert np.any(differences < 0.0)


@pytest.mark.parametrize("model", ["fold", "transcritical", "supercritical_pitchfork"])
def test_model_informed_fit_estimates_margin_without_truth_argument(model: str) -> None:
    margin = 0.16
    times, observed = generate_observation(
        model, margin, noise_sd=1e-6, duration_s=8.0, seed=9831
    )
    estimate = fit_canonical_trajectory(model, times, observed, noise_sd=1e-6)

    assert estimate.converged
    assert estimate.estimated_margin == pytest.approx(margin, rel=0.02)
    assert estimate.estimated_equilibrium == pytest.approx(
        stable_equilibrium(model, margin), rel=0.02
    )


def test_model_agnostic_power_fit_recovers_signatures_and_unknown_boundary() -> None:
    settings = np.asarray(CONTROL_SETTINGS)
    margins = np.asarray(TRUE_MARGINS)
    equilibrium = np.sqrt(margins)
    rate = 2.0 * np.sqrt(margins)
    estimate = estimate_scaling(settings, equilibrium, rate)

    assert estimate.converged
    assert estimate.equilibrium_exponent == pytest.approx(0.5, abs=1e-3)
    assert estimate.rate_exponent == pytest.approx(0.5, abs=1e-3)
    assert estimate.estimated_boundary == pytest.approx(1.0, abs=1e-3)


def test_free_exponential_fit_estimates_equilibrium_and_rate() -> None:
    times = np.linspace(0.0, 8.0, 81)
    equilibrium = 0.4
    rate = 0.8
    observed = equilibrium + 0.08 * np.exp(-rate * times)
    estimate = estimate_exponential(times, observed, noise_sd=0.001)

    assert estimate.converged
    assert estimate.equilibrium == pytest.approx(equilibrium, abs=1e-5)
    assert estimate.recovery_rate == pytest.approx(rate, rel=1e-5)
    assert estimate.r_squared == pytest.approx(1.0)


def test_observation_generation_is_deterministic_for_a_fixed_seed() -> None:
    first = generate_observation("fold", 0.09, 0.001, 3.0, seed=4422)
    second = generate_observation("fold", 0.09, 0.001, 3.0, seed=4422)

    np.testing.assert_array_equal(first[0], second[0])
    np.testing.assert_array_equal(first[1], second[1])


@pytest.mark.parametrize("margin", [0.0, -0.1, math.nan, math.inf])
def test_stable_equilibrium_rejects_boundary_and_invalid_margins(margin: float) -> None:
    with pytest.raises(ScientificValidationError):
        stable_equilibrium("fold", margin)


def test_unknown_model_and_invalid_scaling_data_are_rejected() -> None:
    with pytest.raises(ScientificValidationError):
        stable_equilibrium("unknown", 0.1)
    with pytest.raises(ScientificValidationError):
        estimate_scaling(
            np.asarray([1.0, 1.1, 1.2, 1.3]),
            np.asarray([0.1, -0.2, 0.3, 0.4]),
            np.asarray([0.1, 0.2, 0.3, 0.4]),
        )


def test_frozen_exponent_acceptance_uses_error_of_mean_estimate() -> None:
    signatures = {
        "fold": (0.5, 0.5),
        "transcritical": (1.0, 1.0),
        "supercritical_pitchfork": (0.5, 1.0),
    }
    scaling_summary = [
        {
            "model": model,
            "noise_sd": 0.0005,
            "window_s": 8.0,
            "mean_equilibrium_exponent": q + 0.1,
            "true_equilibrium_exponent": q,
            "equilibrium_exponent_ci95_low": q - 0.05,
            "equilibrium_exponent_ci95_high": q + 0.15,
            "mean_rate_exponent": p - 0.1,
            "true_rate_exponent": p,
            "rate_exponent_ci95_low": p - 0.15,
            "rate_exponent_ci95_high": p + 0.05,
            "median_absolute_equilibrium_exponent_error": 0.9,
            "median_absolute_rate_exponent_error": 0.9,
        }
        for model, (q, p) in signatures.items()
    ]
    metric_rows = [
        {
            "method": "model_informed",
            "noise_sd": 0.0005,
            "window_s": 8.0,
            "distance": "all",
            "accuracy_non_abstained": 0.8,
            "abstention_rate": 0.1,
        },
        {
            "method": "model_agnostic",
            "noise_sd": 0.0005,
            "window_s": 8.0,
            "distance": "all",
            "balanced_accuracy_non_abstained": 0.8,
            "abstention_rate": 0.1,
            "fold_sensitivity_non_abstained": 0.8,
            "linear_specificity_non_abstained": 0.8,
        },
    ]
    parameter_rows = [
        {
            "model": model,
            "model_informed_margin_median_relative_error": 0.1,
            "model_informed_equilibrium_median_relative_error": 0.1,
        }
        for model in signatures
    ]

    criteria = _evaluate_frozen_criteria(
        [], metric_rows, scaling_summary, parameter_rows
    )
    assert criteria["scaling_identifiability_all_classes"]["passed"]

    scaling_summary[0]["mean_equilibrium_exponent"] = 0.71
    criteria = _evaluate_frozen_criteria(
        [], metric_rows, scaling_summary, parameter_rows
    )
    assert not criteria["scaling_identifiability_all_classes"]["passed"]
