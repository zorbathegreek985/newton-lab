"""Tests for Newton Lab's bounded Phase 18 case study."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
import pytest

import newton_lab.dynamics as dynamics
from newton_lab.cross_domain import (
    DampingTransferCaseStudy,
    analytical_step_response,
    run_damping_transfer_case_study,
)
from newton_lab.dynamics import DampedOscillatorResult
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.experiments import (
    ExperimentRunStatus,
    make_oscillator_experiment,
    run_experiment,
)


def _small_case(**kwargs: Any) -> DampingTransferCaseStudy:
    return run_damping_transfer_case_study(output_point_count=1201, **kwargs)


def test_closed_forms_cover_under_critical_and_overdamped_responses() -> None:
    tau = 0.5
    zeta = 0.25
    half_period = np.pi * tau / np.sqrt(1.0 - zeta**2)
    under = analytical_step_response(
        (0.0, half_period), damping_ratio=zeta, time_constant_s=tau
    )
    assert under[0] == pytest.approx(1.0)
    assert under[1] == pytest.approx(-np.exp(-zeta * np.pi / np.sqrt(1.0 - zeta**2)))
    assert analytical_step_response(
        (0.0, tau), damping_ratio=1.0, time_constant_s=tau
    ) == pytest.approx((1.0, 2.0 / np.e))
    assert analytical_step_response(
        (0.0,), damping_ratio=2.0, time_constant_s=tau
    ) == pytest.approx((1.0,))


@pytest.mark.parametrize(
    ("times", "damping", "tau"),
    [
        ((-1.0,), 0.2, 1.0),
        ((0.0, float("nan")), 0.2, 1.0),
        ((1.0, 0.0), 0.2, 1.0),
        ((0.0,), -1.0, 1.0),
        ((0.0,), 0.2, 0.0),
    ],
)
def test_analytical_response_rejects_invalid_inputs(
    times: tuple[float, ...], damping: float, tau: float
) -> None:
    with pytest.raises(ScientificValidationError):
        analytical_step_response(times, damping_ratio=damping, time_constant_s=tau)


def test_source_values_mapping_order_duplicates_and_provenance_are_retained() -> None:
    source = run_experiment(
        make_oscillator_experiment(
            experiment_id="phase18_duplicate_source_values",
            sweep_values=(0.0, 0.1, 0.1, 0.2),
            metrics=("maximum_absolute_velocity_m_per_s",),
        )
    )
    before = source.model_dump(mode="python")
    ratios = (0.0, 0.025, 0.025, 0.05, 0.7, 1.0, 2.0)
    case = _small_case(
        source_record=source,
        damping_ratios=ratios,
        baseline_damping_ratio=0.025,
    )
    assert source.model_dump(mode="python") == before
    assert case.source_experiment.model_dump(mode="python") == before
    assert case.source_metric_parameter_values == (0.0, 0.1, 0.1, 0.2)
    assert case.source_damping_ratio_values == pytest.approx((0.0, 0.025, 0.025, 0.05))
    expected_metrics = tuple(
        next(
            metric.value
            for metric in run.metrics
            if metric.metric_id == case.source_metric_id
        )
        for run in source.runs
    )
    assert case.source_metric_values == pytest.approx(expected_metrics)
    assert case.source_normalized_peak_rate_values_per_s == pytest.approx(
        tuple(value / 0.1 for value in expected_metrics)
    )
    assert all(
        difference is not None and difference < 0.02
        for difference in case.source_target_peak_rate_absolute_differences_per_s
    )
    assert tuple(point.damping_ratio for point in case.sweep_points) == ratios
    assert case.source_metric_unit == "m/s"
    assert case.source_hypothesis.validated_application is False
    assert case.real_world_application_status == "unvalidated"


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"stiffness_n_per_m": 5.0}, "stiffness_n_per_m"),
        ({"initial_displacement_m": 0.2}, "initial_displacement_m"),
        ({"metrics": ("final_displacement_m",)}, "maximum_absolute_velocity"),
    ],
)
def test_incompatible_or_missing_source_evidence_is_rejected(
    overrides: dict[str, Any], message: str
) -> None:
    source = run_experiment(
        make_oscillator_experiment(
            experiment_id="phase18_incompatible_source",
            sweep_values=(0.0, 0.1, 0.2),
            **{
                "metrics": ("maximum_absolute_velocity_m_per_s",),
                **overrides,
            },
        )
    )
    with pytest.raises(ScientificValidationError, match=message):
        _small_case(source_record=source)


def test_source_failures_are_rejected_instead_of_replaced() -> None:
    source = run_experiment(
        make_oscillator_experiment(
            experiment_id="phase18_failed_source",
            sweep_values=(0.0, -0.1, 0.2),
            metrics=("maximum_absolute_velocity_m_per_s",),
        )
    )
    assert source.runs[1].status == ExperimentRunStatus.FAILED
    with pytest.raises(ScientificValidationError, match="source failures"):
        _small_case(source_record=source)


def test_missing_run_sweep_provenance_is_rejected() -> None:
    source = run_experiment(
        make_oscillator_experiment(
            experiment_id="phase18_missing_provenance",
            sweep_values=(0.0, 0.1, 0.2),
            metrics=("maximum_absolute_velocity_m_per_s",),
        )
    )
    changed_run = source.runs[1].model_copy(update={"swept_parameter_name": None})
    incomplete = source.model_copy(
        update={"runs": (source.runs[0], changed_run, source.runs[2])}
    )
    with pytest.raises(ScientificValidationError, match="run order/provenance"):
        _small_case(source_record=incomplete)


def test_missing_metric_on_declared_successful_run_is_rejected() -> None:
    source = run_experiment(
        make_oscillator_experiment(
            experiment_id="phase18_missing_metric_value",
            sweep_values=(0.0, 0.1, 0.2),
            metrics=("maximum_absolute_velocity_m_per_s",),
        )
    )
    changed_run = source.runs[1].model_copy(update={"metrics": ()})
    incomplete = source.model_copy(
        update={"runs": (source.runs[0], changed_run, source.runs[2])}
    )
    with pytest.raises(
        ScientificValidationError, match="missing its finite maximum velocity metric"
    ):
        _small_case(source_record=incomplete)


def test_target_results_match_exact_response_and_show_settling_tradeoff() -> None:
    case = _small_case()
    underdamped = case.sweep_points[1]
    critical = next(point for point in case.sweep_points if point.damping_ratio == 1.0)
    overdamped = next(
        point for point in case.sweep_points if point.damping_ratio == 2.0
    )
    assert underdamped.metrics is not None
    assert underdamped.metrics.maximum_analytical_absolute_error < 1e-8
    assert underdamped.metrics.target_crossing_count > 0
    assert critical.metrics is not None
    assert critical.metrics.target_crossing_count == 0
    assert overdamped.metrics is not None
    assert overdamped.metrics.target_crossing_count == 0
    assert critical.metrics.sampled_settling_time_s is not None
    assert overdamped.metrics.sampled_settling_time_s is not None
    assert (
        overdamped.metrics.sampled_settling_time_s
        > critical.metrics.sampled_settling_time_s
    )
    assert case.peak_rate_trend_on_tested_grid == "nonincreasing_on_tested_grid"
    assert case.overshoot_trend_on_tested_grid == "nonincreasing_on_tested_grid"
    assert (
        case.settling_time_trend_on_tested_grid
        == "unavailable_censored_no_settling_within_horizon"
    )


def test_target_failures_preserve_sweep_order_and_status(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original: Callable[..., DampedOscillatorResult] = (
        dynamics.simulate_damped_oscillator
    )

    def fail_one(*args: Any, **kwargs: Any) -> DampedOscillatorResult:
        if kwargs.get("damping_coefficient_kg_per_s") == pytest.approx(0.1):
            raise IntegrationError("controlled test failure")
        return original(*args, **kwargs)

    monkeypatch.setattr("newton_lab.cross_domain.simulate_damped_oscillator", fail_one)
    case = _small_case()
    assert tuple(point.damping_ratio for point in case.sweep_points) == (
        0.0,
        0.025,
        0.05,
        0.7,
        1.0,
        2.0,
    )
    failed = next(point for point in case.sweep_points if point.status == "failed")
    assert failed.damping_ratio == pytest.approx(0.025)
    assert failed.failure_message == "controlled test failure"
    assert case.peak_rate_trend_on_tested_grid == "unavailable_incomplete_sweep"


def test_report_contains_equations_source_mapping_results_and_safeguards() -> None:
    case = _small_case()
    report = case.report_markdown
    assert report.count("Target sweep:") == 1
    assert "m x'' + c x' + k x = 0" in report
    assert "zeta = c/(2 sqrt(m k))" in report
    assert "Physical damping c (kg/s)" in report
    assert "Peak |q'| (s^-1)" in report
    assert case.source_hypothesis.hypothesis_id in report
    assert "unvalidated" in report
    assert "chronological out-of-sample" in report
    assert "profitability" in report


def test_case_study_is_reproducible_and_never_validates_application() -> None:
    first = _small_case()
    second = _small_case()
    assert first == second
    assert first.report_markdown == second.report_markdown
    assert first.target_domain_interpretation == "structural_analogy_only"
    assert first.real_world_application_status == "unvalidated"
    assert (
        first.source_hypothesis.evidence_status.value
        == "proposed_application_hypothesis"
    )
