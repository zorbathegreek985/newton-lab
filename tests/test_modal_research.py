"""Deterministic Phase 22 modal calculation and sampling tests."""

from __future__ import annotations

import math

import pytest

from newton_lab.dynamics import simulate_damped_oscillator
from newton_lab.exceptions import ScientificValidationError
from newton_lab.heat_diffusion_case_study import (
    HeatDiffusionModel,
    SpatialMode,
    evaluate_heat_diffusion,
    steady_temperature_profile,
)
from newton_lab.modal_research import (
    ModalDynamicsResearchResult,
    _estimate_log_decay_rate,
    heat_mode_decay_rate,
    oscillator_modal_properties,
    render_modal_dynamics_report,
    run_modal_dynamics_research,
)


@pytest.fixture(scope="module")
def research() -> ModalDynamicsResearchResult:
    return run_modal_dynamics_research()


def test_heat_modes_follow_phase19_analytical_attenuation(
    research: ModalDynamicsResearchResult,
) -> None:
    for mode in research.heat_modes:
        expected = tuple(
            math.exp(-mode.decay_rate_per_s * time) for time in mode.times_s
        )
        assert mode.attenuation_factors == pytest.approx(expected, rel=1e-14)
        assert mode.modal_amplitudes_k == pytest.approx(
            tuple(mode.initial_amplitude_k * factor for factor in expected), rel=1e-14
        )
    assert research.heat_modes[0].decay_rate_per_s == pytest.approx(
        research.heat_slow_mode_fit_per_s, rel=1e-12
    )


def test_heat_decay_rate_scales_with_diffusivity_mode_and_length() -> None:
    base = heat_mode_decay_rate(1e-4, 1, 1.0)
    assert heat_mode_decay_rate(2e-4, 1, 1.0) == pytest.approx(2 * base)
    assert heat_mode_decay_rate(1e-4, 2, 1.0) == pytest.approx(4 * base)
    assert heat_mode_decay_rate(1e-4, 1, 2.0) == pytest.approx(base / 4)
    with pytest.raises(ScientificValidationError):
        heat_mode_decay_rate(1e-4, 0, 1.0)


def test_mode_fit_respects_resolution_floor_and_combined_mode_can_cross_zero(
    research: ModalDynamicsResearchResult,
) -> None:
    fast = research.heat_modes[-1]
    assert fast.modal_amplitudes_k[-1] < research.declared_observation_floor_k
    assert research.heat_fast_mode_fit_above_resolution_per_s == pytest.approx(
        fast.decay_rate_per_s, rel=1e-12
    )
    values = research.combined_signal_perturbation_k
    assert min(values) < 0.0 < max(values)
    with pytest.raises(ScientificValidationError, match="strictly positive"):
        _estimate_log_decay_rate((0.0, 1.0, 2.0), (1.0, 0.0, -1.0))


def test_initial_amplitudes_scale_linear_responses_without_changing_rates() -> None:
    heat = HeatDiffusionModel(
        length_m=1.0,
        left_temperature_k=300.0,
        right_temperature_k=400.0,
        diffusivity_m2_per_s=1e-4,
        modes=(SpatialMode(mode_number=1, amplitude=1.0),),
    )
    scaled_heat = heat.model_copy(
        update={"modes": (SpatialMode(mode_number=1, amplitude=3.0),)}
    )
    base = steady_temperature_profile(heat, (0.25, 0.75))
    first = evaluate_heat_diffusion(heat, (0.25, 0.75), (0.0, 100.0))
    scaled = evaluate_heat_diffusion(scaled_heat, (0.25, 0.75), (0.0, 100.0))
    for first_row, scaled_row in zip(
        first.temperature_k, scaled.temperature_k, strict=True
    ):
        observed = tuple(value - b for value, b in zip(scaled_row, base, strict=True))
        expected = tuple(
            3.0 * (value - b) for value, b in zip(first_row, base, strict=True)
        )
        assert observed == pytest.approx(expected)
    assert scaled.attenuation_factors == first.attenuation_factors

    reference = simulate_damped_oscillator(1.0, 0.1, 4.0, 0.1, 0.0, 3.0, num_points=101)
    doubled = simulate_damped_oscillator(1.0, 0.1, 4.0, 0.2, 0.0, 3.0, num_points=101)
    assert doubled.displacement_m == pytest.approx(2.0 * reference.displacement_m)
    assert doubled.velocity_m_per_s == pytest.approx(2.0 * reference.velocity_m_per_s)


def test_oscillator_envelope_and_frequency_depend_separately_on_parameters(
    research: ModalDynamicsResearchResult,
) -> None:
    baseline = research.oscillator_rates[1]
    assert baseline.damping_ratio < 1.0
    assert baseline.envelope_decay_rate_per_s == pytest.approx(0.05)
    assert baseline.damped_frequency_rad_per_s < baseline.natural_frequency_rad_per_s
    assert baseline.sampled_decay_rate_per_s == pytest.approx(
        baseline.envelope_decay_rate_per_s, abs=2e-6
    )
    assert baseline.maximum_analytical_displacement_error_m < 1e-8

    damping_rows = research.oscillator_rates[:3]
    assert tuple(
        row.envelope_decay_rate_per_s for row in damping_rows
    ) == pytest.approx((0.025, 0.05, 0.1))
    mass_rows = (research.oscillator_rates[3], baseline, research.oscillator_rates[4])
    assert tuple(row.mass_kg for row in mass_rows) == (0.5, 1.0, 2.0)
    assert tuple(row.envelope_decay_rate_per_s for row in mass_rows) == pytest.approx(
        (0.1, 0.05, 0.025)
    )
    assert all(row.damping_ratio < 1.0 for row in research.oscillator_rates)
    zeta, omega_n, gamma, omega_d = oscillator_modal_properties(1.0, 0.1, 4.0)
    assert zeta == pytest.approx(0.025)
    assert omega_n == pytest.approx(2.0)
    assert gamma == pytest.approx(0.05)
    assert omega_d == pytest.approx(baseline.damped_frequency_rad_per_s)
    assert oscillator_modal_properties(1.0, 0.0, 4.0) == pytest.approx(
        (0.0, 2.0, 0.0, 2.0)
    )
    with pytest.raises(ScientificValidationError, match="underdamped"):
        oscillator_modal_properties(1.0, 4.0, 4.0)


def test_finite_oscillator_sampling_bias_and_short_window(
    research: ModalDynamicsResearchResult,
) -> None:
    dense = research.oscillator_dense_sampling
    sparse = research.oscillator_sparse_sampling
    short = research.oscillator_short_window
    assert dense.fitted_decay_rate_per_s == pytest.approx(0.05, abs=2e-6)
    assert sparse.fitted_decay_rate_per_s is not None
    assert abs(sparse.fitted_decay_rate_per_s - 0.05) > 0.02
    assert short.absolute_peak_count < 3
    assert short.fitted_decay_rate_per_s is None


def test_computed_report_contains_results_and_model_distinctions(
    research: ModalDynamicsResearchResult,
) -> None:
    report = render_modal_dynamics_report(research)
    assert "Attenuation factors" in report
    assert "gamma" in report
    assert "short window" in report
    assert "combined perturbation" in report
