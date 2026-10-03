"""Tests for the focused Phase 29 diffusion-inspired smoothing study."""

from __future__ import annotations

import math
from unittest.mock import patch

import numpy as np
import pytest

from newton_lab.diffusion_signal_smoothing import (
    DEFAULT_EXPOSURES,
    SmoothingRun,
    _moving_average,
    run_diffusion_smoothing_study,
)
from newton_lab.exceptions import ScientificValidationError


def _case_run(case_id: str, exposure: float) -> SmoothingRun:
    study = run_diffusion_smoothing_study(exposures=(exposure,))
    return next(
        run for run in study.runs if run.case_id == case_id and run.exposure == exposure
    )


def test_two_frequency_modes_follow_analytical_diffusion_attenuation() -> None:
    exposure = 0.001
    run = _case_run("two_frequency", exposure)
    assert run.status == "succeeded"
    assert run.expected_slow_attenuation == pytest.approx(
        math.exp(-(math.pi**2) * exposure)
    )
    assert run.expected_fast_attenuation == pytest.approx(
        math.exp(-((4.0 * math.pi) ** 2) * exposure)
    )
    assert run.slow_mode_retention == pytest.approx(
        run.expected_slow_attenuation, abs=1e-12
    )
    assert run.fast_mode_retention == pytest.approx(
        run.expected_fast_attenuation, abs=1e-12
    )
    assert run.analytical_modal_error_max == pytest.approx(0.0, abs=1e-15)
    assert run.fast_mode_retention is not None
    assert run.slow_mode_retention is not None
    assert run.fast_mode_retention < run.slow_mode_retention


@pytest.mark.parametrize("case_id", ["two_frequency", "slow_plus_high_frequency"])
def test_zero_exposure_is_identity_for_exact_finite_mode_inputs(case_id: str) -> None:
    run = _case_run(case_id, 0.0)
    assert run.status == "succeeded"
    assert run.diffusion_rmse_to_input == pytest.approx(0.0, abs=1e-14)
    assert run.slow_mode_retention == pytest.approx(1.0, abs=1e-12)
    assert run.fast_mode_retention == pytest.approx(1.0, abs=1e-12)


def test_step_projection_error_is_separate_from_zero_exposuremoothing() -> None:
    run = _case_run("abrupt_pulse", 0.0)
    assert run.status == "succeeded"
    assert run.projection_rmse_to_input is not None
    assert run.projection_rmse_to_input > 0.0
    assert run.diffusion_rmse_to_input == pytest.approx(
        run.projection_rmse_to_input, abs=1e-14
    )
    assert run.pre_transition_max_influence is not None
    assert run.pre_transition_change_from_projection == pytest.approx(0.0, abs=1e-14)


def test_high_frequency_component_is_attenuated_and_truth_error_is_measured() -> None:
    run = _case_run("slow_plus_high_frequency", 0.001)
    assert run.status == "succeeded"
    expected = math.exp(-((12.0 * math.pi) ** 2) * 0.001)
    assert run.fast_mode_retention is not None
    assert run.fast_mode_retention == pytest.approx(expected, abs=1e-12)
    assert run.fast_mode_retention < 0.25
    assert run.raw_rmse_to_truth is not None
    assert run.diffusion_rmse_to_truth is not None
    assert run.moving_average_rmse_to_truth is not None
    assert run.moving_average_rmse_to_input is not None
    assert run.diffusion_rmse_to_input is not None
    assert run.diffusion_rmse_to_truth < run.raw_rmse_to_truth
    assert run.diffusion_rmse_to_input > 0.0


def test_moving_average_baseline_is_centered_and_edge_renormalized() -> None:
    study = run_diffusion_smoothing_study(exposures=(0.0,))
    for run in study.runs:
        assert run.status == "succeeded"
        assert len(run.moving_average_values) == len(study.normalized_signal_time)
    assert "centered" in study.report_markdown
    assert "clipped" in study.report_markdown
    assert _moving_average(np.full(25, 3.0), 21) == pytest.approx(np.full(25, 3.0))


def test_diffusion_preserves_fixed_endpoints_and_ma_does_not_change_with_exposure() -> (
    None
):
    study = run_diffusion_smoothing_study(exposures=(0.0, 0.001))
    for case_id in ("two_frequency", "slow_plus_high_frequency", "abrupt_pulse"):
        rows = [run for run in study.runs if run.case_id == case_id]
        assert rows[0].diffusion_values[0] == pytest.approx(0.0, abs=1e-14)
        assert rows[0].diffusion_values[-1] == pytest.approx(0.0, abs=1e-14)
        assert rows[0].moving_average_values == rows[1].moving_average_values


def test_step_response_reports_delay_width_plateau_and_noncausal_influence() -> None:
    run = _case_run("abrupt_pulse", 0.001)
    assert run.status == "succeeded"
    assert run.step_rising_delay is not None
    assert run.step_falling_delay is not None
    assert run.step_rising_10_90_width is not None
    assert run.step_falling_90_10_width is not None
    assert run.step_rising_10_90_width > 0.0
    assert run.step_falling_90_10_width > 0.0
    assert run.step_plateau_retention is not None
    assert run.pre_transition_max_influence is not None
    assert run.pre_transition_change_from_projection is not None
    assert run.pre_transition_change_from_projection > 0.0
    assert run.moving_average_pre_transition_influence is not None


def test_default_exposureweep_is_deterministic_and_preserves_order() -> None:
    first = run_diffusion_smoothing_study()
    second = run_diffusion_smoothing_study()
    assert first == second
    assert first.exposures == DEFAULT_EXPOSURES
    assert len(first.runs) == 3 * len(DEFAULT_EXPOSURES)
    for case_id in ("two_frequency", "slow_plus_high_frequency", "abrupt_pulse"):
        rows = [run for run in first.runs if run.case_id == case_id]
        assert [run.sweep_index for run in rows] == list(range(len(DEFAULT_EXPOSURES)))
        assert [run.exposure for run in rows] == list(DEFAULT_EXPOSURES)


def test_repeated_and_invalid_exposures_are_retained_as_ordered_rows() -> None:
    study = run_diffusion_smoothing_study(
        exposures=(0.001, 0.0, -0.1, 0.001, math.nan, math.inf, True)
    )
    for case_id in ("two_frequency", "slow_plus_high_frequency", "abrupt_pulse"):
        rows = [run for run in study.runs if run.case_id == case_id]
        assert [run.sweep_index for run in rows] == list(range(7))
        assert rows[0].exposure == rows[3].exposure == 0.001
        assert [run.status for run in rows] == [
            "succeeded",
            "succeeded",
            "failed",
            "succeeded",
            "failed",
            "failed",
            "failed",
        ]
        assert rows[2].diffusion_rmse_to_truth is None
        assert rows[4].failure_message is not None


def test_numerical_evaluation_failure_is_recorded_for_each_case() -> None:
    with patch(
        "newton_lab.diffusion_signal_smoothing.evaluate_heat_diffusion",
        side_effect=ScientificValidationError("controlled numerical failure"),
    ):
        study = run_diffusion_smoothing_study(exposures=(0.001,))
    assert [run.status for run in study.runs] == ["failed", "failed", "failed"]
    assert all(run.failure_message is not None for run in study.runs)
    assert all(run.diffusion_values == () for run in study.runs)


@pytest.mark.parametrize(
    ("sample_count", "modal_truncation", "moving_average_width", "exposures"),
    [
        (100, 80, 21, DEFAULT_EXPOSURES),
        (501, 500, 21, DEFAULT_EXPOSURES),
        (501, 80, 20, DEFAULT_EXPOSURES),
        (501, 80, 502, DEFAULT_EXPOSURES),
        (501, 80, 21, ()),
    ],
)
def test_invalid_study_configuration_is_rejected(
    sample_count: int,
    modal_truncation: int,
    moving_average_width: int,
    exposures: tuple[float, ...],
) -> None:
    with pytest.raises(ScientificValidationError):
        run_diffusion_smoothing_study(
            sample_count=sample_count,
            modal_truncation=modal_truncation,
            moving_average_width=moving_average_width,
            exposures=exposures,
        )


def test_report_separates_evidence_classes_and_disclaims_financial_validation() -> None:
    report = run_diffusion_smoothing_study(exposures=(0.0,)).report_markdown.lower()
    assert "established mathematical result" in report
    assert "experimentally supported result" in report
    assert "candidate cross-domain application" in report
    assert "unvalidated financial application" in report
    assert "noncausal" in report
    assert "profitability" in report
    assert "not a random-noise realization" in report
