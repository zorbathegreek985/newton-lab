"""Tests for Phase 49 stability-boundary rate analysis."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from newton_lab.critical_slowing_down import (
    _fit_log_rate,
    _heat_condition,
    _pendulum_condition,
    pendulum_linearized_spectrum,
    run_critical_slowing_down_study,
)
from newton_lab.exceptions import ScientificValidationError


def test_pendulum_linearized_eigenvalues_and_envelope_rate() -> None:
    first, second, rate, frequency = pendulum_linearized_spectrum(0.2)
    assert first.real == pytest.approx(-0.2)
    assert second.real == pytest.approx(-0.2)
    assert first.imag == pytest.approx(frequency)
    assert second.imag == pytest.approx(-frequency)
    assert frequency == pytest.approx(math.sqrt(1.0 - 0.2**2))
    assert rate == pytest.approx(0.2)


@pytest.mark.parametrize("distance", [0.0, -0.1, 1.0, float("nan")])
def test_pendulum_spectrum_rejects_boundary_or_outside_values(distance: float) -> None:
    with pytest.raises(ScientificValidationError):
        pendulum_linearized_spectrum(distance)


def test_heat_recovery_prediction_is_first_mode_eigenvalue() -> None:
    row = _heat_condition(0.1, 0.01)
    assert row.equilibrium.startswith("uniform steady profile")
    assert row.stability.startswith("asymptotically stable")
    assert row.predicted_rate_per_s == pytest.approx(0.1 * math.pi**2)
    assert row.measured_rate_per_s == pytest.approx(row.predicted_rate_per_s, rel=1e-10)
    assert row.solver_rtol is None


def test_pendulum_nonlinear_measurement_agrees_with_local_prediction() -> None:
    row = _pendulum_condition(0.1, 0.01, 1e-10, 1e-12)
    assert row.equilibrium.startswith("downward rest")
    assert row.stability.startswith("locally asymptotically stable")
    assert row.predicted_rate_per_s == pytest.approx(0.1)
    assert row.measurement_points >= 3
    assert row.local_approximation_accepted


def test_nearest_swept_damping_to_boundary_remains_measurable() -> None:
    row = _pendulum_condition(0.025, 0.01, 1e-10, 1e-12)
    assert row.measurement_window_s == pytest.approx(480.0)
    assert row.measurement_points >= 50
    assert math.isfinite(row.measured_rate_per_s)


def test_pendulum_perturbation_sensitivity_is_measurable() -> None:
    small = _pendulum_condition(0.1, 0.01, 1e-10, 1e-12)
    larger = _pendulum_condition(0.1, 0.3, 1e-10, 1e-12)
    relative_difference = abs(
        small.measured_rate_per_s - larger.measured_rate_per_s
    ) / (small.measured_rate_per_s)
    assert relative_difference < 0.02


def test_log_rate_fit_recovers_exponential_and_rejects_nondecay() -> None:
    times = np.linspace(0.0, 5.0, 101)
    rate, r_squared = _fit_log_rate(times, np.exp(-0.7 * times))
    assert rate == pytest.approx(0.7)
    assert r_squared == pytest.approx(1.0)
    with pytest.raises(ScientificValidationError):
        _fit_log_rate(times, np.ones_like(times))


def test_study_writes_auditable_tables_and_supported_criteria(tmp_path: Path) -> None:
    metadata = run_critical_slowing_down_study(tmp_path)
    assert metadata["conditions"] == 30
    assert all(
        summary["model_prediction_supported"] for summary in metadata["model_summaries"]
    )
    assert {path.name for path in tmp_path.iterdir()} == {
        "per_condition_results.csv",
        "model_summary.csv",
        "metadata.json",
        "recovery_rate_vs_boundary_distance.png",
        "predicted_vs_measured_rates.png",
        "recovery_rate_scaling.png",
    }
