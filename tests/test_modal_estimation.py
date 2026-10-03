"""Tests for Phase 23 imperfect-observation modal estimation."""

from __future__ import annotations

from math import isclose

import pytest

from newton_lab import modal_estimation
from newton_lab.exceptions import ScientificValidationError
from newton_lab.modal_estimation import (
    EstimateStatus,
    SyntheticObservationCase,
    estimate_exponential_rate,
    estimate_oscillator_envelope,
    evaluate_estimate,
    generate_exponential_observations,
    generate_heat_mode_observations,
    generate_mixed_heat_observations,
    generate_oscillator_observations,
    run_phase23_research,
)


def _times(stop: float, step: float) -> tuple[float, ...]:
    count = round(stop / step)
    return tuple(index * step for index in range(count + 1))


def test_scalar_exponential_fits_recover_exact_rate_and_keep_truth_separate() -> None:
    case = generate_exponential_observations(
        scenario_id="exact",
        amplitude=2.0,
        rate_per_s=0.4,
        times_s=_times(8.0, 0.2),
        value_unit="K",
    )

    for method in ("log_linear", "nonlinear"):
        estimate = estimate_exponential_rate(case.observations, method=method)
        evaluation = evaluate_estimate(estimate, case.ground_truth)
        assert estimate.status == EstimateStatus.INTERPRETABLE
        assert estimate.fit_converged
        assert estimate.estimated_rate_per_s == pytest.approx(0.4, rel=1e-8)
        assert evaluation.absolute_error_per_s == pytest.approx(0.0, abs=1e-9)
    assert "target_rate_per_s" not in case.observations.__class__.model_fields


def test_heat_mode_rate_matches_phase_19_analytical_rate() -> None:
    case = generate_heat_mode_observations(
        scenario_id="heat_n1",
        diffusivity_m2_per_s=1e-4,
        mode_number=1,
        length_m=1.0,
        initial_amplitude_k=1.0,
        times_s=_times(1200.0, 60.0),
    )
    estimate = estimate_exponential_rate(case.observations)
    assert estimate.estimated_rate_per_s == pytest.approx(
        case.ground_truth.target_rate_per_s, rel=1e-12
    )
    assert estimate.observation_rmse == pytest.approx(0.0, abs=1e-14)


def test_fixed_seed_noise_is_reproducible_and_noise_floor_is_recorded() -> None:
    def generate_noisy_case() -> SyntheticObservationCase:
        return generate_exponential_observations(
            scenario_id="noisy",
            amplitude=1.0,
            rate_per_s=0.3,
            times_s=_times(12.0, 0.25),
            value_unit="signal",
            noise_standard_deviation=0.05,
            noise_seed=2330,
        )

    first = generate_noisy_case()
    second = generate_noisy_case()
    assert first.observations.values == second.observations.values
    estimate = estimate_exponential_rate(first.observations)
    assert estimate.excluded_observations
    assert all(
        "noise floor" in item.reason or "nonpositive" in item.reason
        for item in estimate.excluded_observations
    )
    assert estimate.local_rate_noise_sensitivity_per_s is not None
    assert estimate.status == EstimateStatus.UNRELIABLE


def test_small_rate_against_noise_is_not_identifiable_by_declared_screen() -> None:
    case = generate_exponential_observations(
        scenario_id="weak_rate",
        amplitude=1.0,
        rate_per_s=0.001,
        times_s=_times(1.0, 0.1),
        value_unit="u",
        noise_standard_deviation=0.1,
        noise_seed=2340,
    )
    estimate = estimate_exponential_rate(case.observations, method="nonlinear")
    assert estimate.fit_converged
    assert estimate.status == EstimateStatus.NOT_IDENTIFIABLE
    assert any("not a confidence interval" in item for item in estimate.diagnostics)


def test_invalid_scalar_generator_inputs_are_rejected() -> None:
    with pytest.raises(ScientificValidationError):
        generate_exponential_observations(
            scenario_id="bad",
            amplitude=0.0,
            rate_per_s=1.0,
            times_s=(0.0, 1.0),
            value_unit="u",
        )
    with pytest.raises(ScientificValidationError):
        generate_exponential_observations(
            scenario_id="bad",
            amplitude=1.0,
            rate_per_s=float("nan"),
            times_s=(0.0, 1.0),
            value_unit="u",
        )
    with pytest.raises(ScientificValidationError):
        generate_exponential_observations(
            scenario_id="bad",
            amplitude=1.0,
            rate_per_s=1.0,
            times_s=(0.0, 0.0),
            value_unit="u",
        )


def test_short_scalar_series_is_insufficient() -> None:
    case = generate_exponential_observations(
        scenario_id="short",
        amplitude=1.0,
        rate_per_s=0.2,
        times_s=(0.0, 0.1),
        value_unit="u",
    )
    assert (
        estimate_exponential_rate(case.observations).status
        == EstimateStatus.INSUFFICIENT
    )


def test_nonlinear_optimizer_failure_is_reported(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = generate_exponential_observations(
        scenario_id="solver_failure",
        amplitude=1.0,
        rate_per_s=0.2,
        times_s=_times(2.0, 0.1),
        value_unit="u",
    )

    def fail_fit(*args: object, **kwargs: object) -> None:
        raise RuntimeError("injected solver failure")

    monkeypatch.setattr(modal_estimation, "least_squares", fail_fit)
    estimate = estimate_exponential_rate(case.observations, method="nonlinear")
    assert estimate.status == EstimateStatus.FIT_FAILED
    assert not estimate.fit_converged
    assert "injected solver failure" in estimate.convergence_message


def test_oscillator_envelope_separates_duration_cadence_and_status() -> None:
    short = generate_oscillator_observations(
        scenario_id="short",
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.1,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        times_s=_times(1.5, 0.05),
    )
    short_estimate = estimate_oscillator_envelope(short.observations)
    assert short_estimate.status == EstimateStatus.INSUFFICIENT
    assert short_estimate.estimated_rate_per_s is None

    sparse = generate_oscillator_observations(
        scenario_id="sparse",
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.1,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        times_s=_times(10.0, 0.8333333333333),
    )
    sparse_estimate = estimate_oscillator_envelope(sparse.observations)
    assert sparse_estimate.estimated_rate_per_s is not None
    assert sparse_estimate.status == EstimateStatus.UNRELIABLE
    assert any("samples per estimated" in item for item in sparse_estimate.diagnostics)

    dense = generate_oscillator_observations(
        scenario_id="dense",
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.1,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        times_s=_times(10.0, 0.05),
    )
    dense_estimate = estimate_oscillator_envelope(dense.observations)
    assert dense_estimate.status == EstimateStatus.INTERPRETABLE
    assert dense_estimate.estimated_rate_per_s == pytest.approx(0.05, abs=2e-4)
    assert len(dense_estimate.excluded_observations) > 0


def test_oscillator_observation_output_is_not_an_ode_trajectory() -> None:
    case = generate_oscillator_observations(
        scenario_id="osc",
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.1,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        times_s=_times(10.0, 0.05),
    )
    estimate = estimate_oscillator_envelope(case.observations)
    assert case.observations.observable == "signed_oscillator_displacement"
    assert estimate.method == "absolute_sample_peak_log_linear"
    assert isclose(case.ground_truth.component_rates_per_s[0], 0.05)


def test_mixture_estimates_are_window_dependent_effective_rates() -> None:
    case = generate_mixed_heat_observations(
        scenario_id="mix",
        mode_numbers=(1, 4),
        effective_amplitudes_k=(0.01, 1.0),
        diffusivity_m2_per_s=1e-4,
        length_m=1.0,
        times_s=_times(1200.0, 20.0),
    )
    estimate = estimate_exponential_rate(case.observations)
    evaluation = evaluate_estimate(estimate, case.ground_truth)
    assert case.ground_truth.target_rate_per_s is None
    assert evaluation.absolute_error_per_s is None
    assert "effective mixture rate" in evaluation.interpretation
    assert estimate.estimated_rate_per_s is not None
    assert case.ground_truth.component_rates_per_s[0] < estimate.estimated_rate_per_s
    assert estimate.log_rmse is not None and estimate.log_rmse > 0.15
    assert estimate.status == EstimateStatus.UNRELIABLE


def test_fast_heat_mode_can_fall_below_declared_observation_floor() -> None:
    experiment = next(
        item
        for item in run_phase23_research().experiments
        if item.observations.scenario_id == "mix_fast_component_below_noise_floor"
    )
    full_fit, early_fit, late_fit = (item.estimate for item in experiment.outcomes)
    assert experiment.observations.noise_standard_deviation == 0.005
    assert experiment.observations.noise_seed == 2340
    assert full_fit.status == EstimateStatus.UNRELIABLE
    assert any(
        "three-sigma noise floor" in item.reason
        for item in full_fit.excluded_observations
    )
    assert early_fit.status == EstimateStatus.INTERPRETABLE
    assert late_fit.status == EstimateStatus.INSUFFICIENT


def test_phase23_runner_is_deterministic_and_contains_all_four_effects() -> None:
    first = run_phase23_research()
    second = run_phase23_research()
    assert first == second
    ids = {item.observations.scenario_id for item in first.experiments}
    assert any(item.startswith("scalar_cadence_") for item in ids)
    assert any(item.startswith("scalar_duration_") for item in ids)
    assert any(item.startswith("scalar_noise_") for item in ids)
    assert any(item.startswith("oscillator_cadence_") for item in ids)
    assert any(item.startswith("mix_") for item in ids)
