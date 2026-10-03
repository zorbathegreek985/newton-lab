"""Phase 24 tests for robustness and multi-mode identification."""

from __future__ import annotations

import numpy as np
import pytest

from newton_lab import modal_identification
from newton_lab.exceptions import ScientificValidationError
from newton_lab.modal_identification import (
    BIC_MODEL_THRESHOLD,
    IdentificationStatus,
    ModelSelection,
    ModeObservationCase,
    compare_single_and_two_mode_models,
    evaluate_mode_identification,
    fit_two_mode_exponential,
    generate_irregular_times,
    generate_mode_observations,
    run_phase24_research,
)


def _times(duration: float = 15.0, count: int = 61) -> tuple[float, ...]:
    step = duration / (count - 1)
    return tuple(index * step for index in range(count))


def _case(
    *,
    rates: tuple[float, ...] = (0.2, 0.8),
    amplitudes: tuple[float, ...] = (1.0, 0.6),
    times: tuple[float, ...] | None = None,
    sigma: float = 0.015,
    seed: int = 2401,
) -> ModeObservationCase:
    return generate_mode_observations(
        scenario_id="test_case",
        amplitudes=amplitudes,
        rates_per_s=rates,
        times_s=_times() if times is None else times,
        noise_standard_deviation=sigma,
        noise_seed=seed,
    )


def test_observation_truth_and_noise_free_curve_are_separate_and_reproducible() -> None:
    first = _case()
    second = _case()
    assert first.observations.values == second.observations.values
    assert first.noise_free_signal_values == second.noise_free_signal_values
    assert first.observations.values != first.noise_free_signal_values
    assert "component_rates_per_s" not in first.observations.__class__.model_fields
    assert first.observations.noise_seed == 2401


def test_two_modes_are_recovered_and_mode_matching_ignores_label_swap() -> None:
    case = _case()
    reversed_starts = (
        (0.6, 0.8, 1.0, 0.2),
        (1.0, 0.2, 0.6, 0.8),
        (0.7, 0.1, 0.9, 1.0),
    )
    comparison = compare_single_and_two_mode_models(
        case.observations, initial_guesses=reversed_starts
    )
    evaluation = evaluate_mode_identification(comparison, case.ground_truth)
    assert comparison.two_mode_fit.converged
    assert comparison.two_mode_fit.fitted_rates_per_s is not None
    assert comparison.two_mode_fit.fitted_rates_per_s == pytest.approx(
        (0.2, 0.8), rel=0.1
    )
    assert comparison.selected_model is ModelSelection.TWO
    assert evaluation.status is IdentificationStatus.SUPPORTED
    assert len(evaluation.matched_rate_errors) == 2
    assert evaluation.matched_rate_errors[0].absolute_error_per_s == pytest.approx(
        abs(comparison.two_mode_fit.fitted_rates_per_s[0] - 0.2)
    )


def test_single_mode_model_wins_bic_for_single_exponential() -> None:
    case = _case(rates=(0.2,), amplitudes=(1.0,), seed=2420)
    comparison = compare_single_and_two_mode_models(case.observations)
    evaluation = evaluate_mode_identification(comparison, case.ground_truth)
    assert comparison.selected_model is ModelSelection.SINGLE
    assert comparison.bic_difference_single_minus_two is not None
    assert comparison.bic_difference_single_minus_two <= -BIC_MODEL_THRESHOLD
    assert evaluation.status is IdentificationStatus.SUPPORTED
    assert evaluation.matched_rate_errors[0].relative_error < 0.1
    expected_single_bic = sum(
        (observed - fitted) ** 2
        for observed, fitted in zip(
            case.observations.values,
            comparison.single_mode_fit.fitted_values,
            strict=True,
        )
    ) / case.observations.noise_standard_deviation**2 + 2 * np.log(
        len(case.observations.values)
    )
    assert comparison.single_bic == pytest.approx(expected_single_bic)


def test_close_modes_are_not_forced_into_a_two_mode_identification() -> None:
    case = _case(rates=(0.2, 0.24), amplitudes=(1.0, 1.0), seed=2423)
    comparison = compare_single_and_two_mode_models(case.observations)
    evaluation = evaluate_mode_identification(comparison, case.ground_truth)
    assert comparison.two_mode_fit.converged
    assert comparison.selected_model in (
        ModelSelection.SINGLE,
        ModelSelection.INCONCLUSIVE,
    )
    assert evaluation.status is not IdentificationStatus.SUPPORTED


def test_noise_free_lower_residual_does_not_select_more_complex_model() -> None:
    case = _case(sigma=0.0)
    comparison = compare_single_and_two_mode_models(case.observations)
    assert comparison.two_mode_fit.converged
    assert comparison.two_mode_fit.residual_rmse is not None
    assert comparison.single_mode_fit.observation_rmse is not None
    assert comparison.two_mode_fit.residual_rmse <= (
        comparison.single_mode_fit.observation_rmse
    )
    assert comparison.selected_model is ModelSelection.UNAVAILABLE
    assert comparison.single_bic is None and comparison.two_bic is None
    assert "residual ranking alone" in comparison.comparison_assumption


def test_irregular_timestamps_are_fixed_seed_ordered_and_have_matched_counts() -> None:
    for design in ("mild_jitter", "strong_jitter", "clustered", "missing_middle"):
        first = generate_irregular_times(
            duration_s=15.0, count=61, design=design, seed=2441
        )
        second = generate_irregular_times(
            duration_s=15.0, count=61, design=design, seed=2441
        )
        assert first == second
        assert first[0] == 0.0 and first[-1] == 15.0
        assert all(right > left for left, right in zip(first, first[1:], strict=False))
    missing = generate_irregular_times(
        duration_s=15.0, count=61, design="missing_middle", seed=2441
    )
    assert len(missing) < 61
    regular_match = _times(15.0, len(missing))
    assert len(regular_match) == len(missing)


def test_initialization_results_are_recorded_and_sorted_for_evaluation() -> None:
    case = _case()
    starts = ((0.6, 0.8, 1.0, 0.2), (1.0, 0.2, 0.6, 0.8))
    fit = fit_two_mode_exponential(case.observations, initial_guesses=starts)
    assert len(fit.start_records) == 2
    assert all(record.converged for record in fit.start_records)
    sorted_results = [
        tuple(sorted(record.fitted_rates_per_s or ())) for record in fit.start_records
    ]
    assert sorted_results[0] == pytest.approx(sorted_results[1], rel=1e-5)


def test_too_few_samples_are_insufficient_not_a_failed_fit() -> None:
    case = _case(times=_times(duration=3.0, count=7))
    comparison = compare_single_and_two_mode_models(case.observations)
    evaluation = evaluate_mode_identification(comparison, case.ground_truth)
    assert not comparison.two_mode_fit.converged
    assert evaluation.status is IdentificationStatus.INSUFFICIENT


def test_optimizer_exception_is_reported_as_numerical_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _case()

    def fail_fit(*args: object, **kwargs: object) -> None:
        raise RuntimeError("injected multi-mode failure")

    monkeypatch.setattr(modal_identification, "least_squares", fail_fit)
    comparison = compare_single_and_two_mode_models(case.observations)
    evaluation = evaluate_mode_identification(comparison, case.ground_truth)
    assert not comparison.two_mode_fit.converged
    assert evaluation.status is IdentificationStatus.FIT_FAILED
    assert all(
        "injected multi-mode failure" in record.convergence_message
        for record in comparison.two_mode_fit.start_records
    )


def test_invalid_signals_times_and_initial_guesses_are_rejected() -> None:
    with pytest.raises(ScientificValidationError):
        _case(rates=(0.2, 0.2))
    with pytest.raises(ScientificValidationError):
        _case(rates=(0.2, float("nan")))
    with pytest.raises(ScientificValidationError):
        _case(times=(0.0, 1.0, 1.0))
    case = _case()
    with pytest.raises(ScientificValidationError):
        fit_two_mode_exponential(
            case.observations, initial_guesses=((1.0, 0.2, 0.0, 0.8),)
        )


def test_research_runner_repeats_only_noise_in_repeated_conditions() -> None:
    results = run_phase24_research()
    assert results == run_phase24_research()
    well_separated = [
        item
        for item in results.repeated_noise_runs
        if item.design_name.endswith("well_separated_comparable")
    ]
    assert len(well_separated) == 12
    reference_times = well_separated[0].case.observations.times
    reference_truth = well_separated[0].case.ground_truth
    assert all(
        item.case.observations.times == reference_times for item in well_separated
    )
    assert all(
        item.case.ground_truth.component_rates_per_s
        == reference_truth.component_rates_per_s
        and item.case.ground_truth.component_amplitudes
        == reference_truth.component_amplitudes
        for item in well_separated
    )
    assert len({item.case.observations.noise_seed for item in well_separated}) == 12
    assert results.oscillator_sampling_summary
    assert results.sampling_coverages
    assert any(
        item.evaluation.status is IdentificationStatus.INSUFFICIENT
        for item in results.mode_matrix_runs
    )
    floor_crossing = next(
        item
        for item in results.mode_matrix_runs
        if item.design_name == "fast_mode_crosses_noise_floor"
    )
    assert floor_crossing.case.ground_truth.component_amplitudes == (1.0, 0.1)
    assert floor_crossing.case.observations.noise_standard_deviation == 0.015
    crossing_time = np.log(0.1 / (3.0 * 0.015)) / 0.8
    assert crossing_time == pytest.approx(0.999, abs=0.002)
    assert floor_crossing.comparison.selected_model is ModelSelection.TWO
    assert floor_crossing.evaluation.status is IdentificationStatus.INCONSISTENT
