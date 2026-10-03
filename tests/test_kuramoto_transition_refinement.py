import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.kuramoto_transition_refinement import (
    _estimate_transition,
    _nested_population_sample,
    truncated_lorentzian_critical_coupling,
    truncated_lorentzian_density,
    truncated_lorentzian_normalization_constant,
)

PROTOCOL_PATH = Path("reports/phase_61_kuramoto_transition_finite_size/protocol.json")


def test_truncated_lorentzian_density_normalizes_and_gives_correct_kc() -> None:
    gamma = 1.0
    cutoff = 5.0
    normalization = truncated_lorentzian_normalization_constant(gamma, cutoff)

    assert normalization == pytest.approx(
        math.pi / (2.0 * math.atan(cutoff / gamma)), rel=1e-14
    )
    area, integration_error = quad(
        lambda omega: float(truncated_lorentzian_density(omega, gamma, cutoff)),
        -cutoff,
        cutoff,
        epsabs=1e-12,
    )
    assert area == pytest.approx(1.0, abs=1e-12)
    assert integration_error < 1e-8

    density_at_zero = 1.0 / (2.0 * gamma * math.atan(cutoff / gamma))
    assert truncated_lorentzian_density(0.0, gamma, cutoff) == pytest.approx(
        density_at_zero
    )
    assert truncated_lorentzian_critical_coupling(gamma, cutoff) == pytest.approx(
        2.0 / (math.pi * density_at_zero), rel=1e-14
    )


def test_nested_population_samples_pair_sizes_and_are_reproducible() -> None:
    frequencies, phases = _nested_population_sample(61001, 1024)
    repeated_frequencies, repeated_phases = _nested_population_sample(61001, 1024)
    other_frequencies, _other_phases = _nested_population_sample(61002, 1024)

    np.testing.assert_array_equal(frequencies, repeated_frequencies)
    np.testing.assert_array_equal(phases, repeated_phases)
    assert np.array_equal(frequencies[:64], frequencies[:256][:64])
    assert np.array_equal(phases[:256], phases[:1024][:256])
    assert not np.array_equal(frequencies, other_frequencies)
    assert np.max(np.abs(frequencies)) <= 5.0


def test_transition_estimator_distinguishes_seed_interval_and_grid_bracket() -> None:
    coupling = np.array([1.0, 1.1, 1.2, 1.3, 1.4])
    curve = np.array([0.1, 0.12, 0.2, 0.55, 0.78])
    curves = np.vstack(
        [
            curve + offset
            for offset in (-0.02, -0.01, 0.0, 0.01, 0.02, 0.0, -0.005, 0.005)
        ]
    )

    estimate = _estimate_transition(curves, coupling, bootstrap_seed=616161)

    assert estimate["censoring"] == "none"
    assert estimate["grid_bracket_low_rad_per_s"] == 1.2
    assert estimate["grid_bracket_high_rad_per_s"] == 1.3
    assert estimate["transition_estimate_rad_per_s"] == 1.25
    assert estimate["bootstrap_low_rad_per_s"] <= estimate["bootstrap_high_rad_per_s"]
    assert estimate["bootstrap_interval_censored"] is False
    assert estimate["bootstrap_boundary_fraction"] == 0.0
    assert estimate["selected_grid_interval_width_rad_per_s"] == pytest.approx(0.1)
    assert (
        estimate["seed_interval_width_rad_per_s"]
        != estimate["selected_grid_interval_width_rad_per_s"]
    )


def test_boundary_transition_is_censored_instead_of_reported_as_point() -> None:
    coupling = np.array([1.0, 1.1, 1.2, 1.3])
    curve = np.array([0.1, 0.1, 0.1, 0.9])
    curves = np.repeat(curve[None, :], 8, axis=0)

    estimate = _estimate_transition(curves, coupling, bootstrap_seed=616162)

    assert estimate["transition_estimate_rad_per_s"] is None
    assert estimate["censoring"] == "right"
    assert estimate["grid_bracket_low_rad_per_s"] == 1.2
    assert estimate["grid_bracket_high_rad_per_s"] == 1.3
    assert estimate["bootstrap_interval_censored"] is True
    assert estimate["bootstrap_low_rad_per_s"] is None
    assert estimate["bootstrap_high_rad_per_s"] is None


def test_transition_without_positive_slope_is_unresolved() -> None:
    coupling = np.array([1.0, 1.1, 1.2, 1.3])
    curves = np.repeat(np.array([[0.4, 0.3, 0.2, 0.1]]), 8, axis=0)

    estimate = _estimate_transition(curves, coupling, bootstrap_seed=616163)

    assert estimate["transition_estimate_rad_per_s"] is None
    assert estimate["censoring"] == "no_positive_slope"
    assert estimate["bootstrap_boundary_fraction"] is None
    assert estimate["bootstrap_interval_censored"] is True


def test_protocol_records_frozen_larger_sizes_and_numerical_check() -> None:
    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))

    assert protocol["protocol_status"] == "frozen_before_final_evaluation"
    assert protocol["finite_population_sizes"] == [64, 256, 512, 1024]
    assert protocol["numerical_resolution_check"]["population_size"] == 1024
    assert protocol["phase32_status"] == "BLOCKED_AUTHORIZATION"
    assert protocol["phase52_status"] == "NO-GO"


def test_invalid_truncated_density_parameters_are_rejected() -> None:
    with pytest.raises(ScientificValidationError):
        truncated_lorentzian_normalization_constant(0.0, 5.0)
    with pytest.raises(ScientificValidationError):
        truncated_lorentzian_normalization_constant(1.0, -5.0)
