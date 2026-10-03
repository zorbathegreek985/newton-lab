"""Tests for Phase 57's frozen synthetic noise-dose benchmark."""

from __future__ import annotations

import inspect
import json
import math
from pathlib import Path

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.noise_dose_estimator_comparison import (
    AMPLITUDES,
    BOOTSTRAP_REPLICATES,
    CALIBRATION_SEEDS,
    EVALUATION_SEEDS,
    NOISE_FRACTIONS,
    TIME_COUNT,
    TIME_STEP_S,
    TIME_STOP_S,
    _estimate_correct_form,
    _estimate_free_exponential,
    _estimate_tail_log_linear,
    _exponent_groups,
    _exponent_summaries,
    _fit_record,
    _seed_bootstrap_interval,
    _write_csv,
    design,
    observation,
    standardized_innovations,
    theoretical_exponents,
    theoretical_rate,
    validate_protocol,
)
from newton_lab.recovery_scaling_identifiability import exact_trajectory

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "reports/phase_57_noise_dose_estimator_comparison/protocol.json"


def test_theoretical_rates_and_exponents_follow_each_selected_branch() -> None:
    assert theoretical_rate("fold", 0.16) == pytest.approx(0.8)
    assert theoretical_rate("transcritical", 0.16) == pytest.approx(0.16)
    assert theoretical_rate("supercritical_pitchfork", 0.16) == pytest.approx(0.32)
    assert theoretical_exponents("fold") == (0.5, 0.5)
    assert theoretical_exponents("transcritical") == (1.0, 1.0)
    assert theoretical_exponents("supercritical_pitchfork") == (0.5, 1.0)
    with pytest.raises(ScientificValidationError):
        theoretical_exponents("unknown")


def test_frozen_design_cardinality_and_time_grid_are_explicit() -> None:
    assert len(design()) == 3 * 5 * len(AMPLITUDES) * len(NOISE_FRACTIONS)
    assert TIME_COUNT == 751
    assert TIME_STOP_S == pytest.approx((TIME_COUNT - 1) * TIME_STEP_S)
    assert set(CALIBRATION_SEEDS).isdisjoint(EVALUATION_SEEDS)
    assert len(CALIBRATION_SEEDS) == 8
    assert len(EVALUATION_SEEDS) == 50


def test_protocol_is_frozen_and_rejects_seed_overlap() -> None:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    validate_protocol(protocol)
    protocol["evaluation_seeds"] = [CALIBRATION_SEEDS[0]]
    with pytest.raises(ScientificValidationError, match="mismatch|overlap"):
        validate_protocol(protocol)


def test_noise_is_normalized_to_equilibrium_and_uses_paired_innovations() -> None:
    seed, model, margin, amplitude = 58051, "fold", 0.16, 0.0025
    times, clean, zero, zero_metadata = observation(model, margin, amplitude, 0.0, seed)
    _, _, low, low_metadata = observation(model, margin, amplitude, 0.0025, seed)
    _, _, reference, ref_metadata = observation(model, margin, amplitude, 0.01, seed)
    equilibrium = math.sqrt(margin)
    innovations = standardized_innovations(seed, model, margin)
    assert times.shape == clean.shape == (TIME_COUNT,)
    assert times[0] == 0.0 and times[-1] == TIME_STOP_S
    assert np.array_equal(zero, clean)
    assert np.allclose(low - clean, 0.0025 * equilibrium * innovations)
    assert np.allclose(reference - clean, 0.01 * equilibrium * innovations)
    assert zero_metadata["noise_sd"] == 0.0
    assert math.isinf(zero_metadata["snr_rms_excursion"])
    assert low_metadata["noise_sd"] == pytest.approx(0.0025 * equilibrium)
    assert ref_metadata["noise_sd"] == pytest.approx(0.01 * equilibrium)


def test_noise_generation_is_deterministic_for_a_fixed_seed() -> None:
    first = observation("supercritical_pitchfork", 0.25, 0.04, 0.01, 58053)
    second = observation("supercritical_pitchfork", 0.25, 0.04, 0.01, 58053)
    assert np.array_equal(first[0], second[0])
    assert np.array_equal(first[1], second[1])
    assert np.array_equal(first[2], second[2])


def test_zero_noise_observations_do_not_depend_on_replicate_seed() -> None:
    _, clean_a, observed_a, _ = observation("transcritical", 0.09, 0.04, 0.0, 58050)
    _, clean_b, observed_b, _ = observation("transcritical", 0.09, 0.04, 0.0, 58099)
    assert np.array_equal(clean_a, clean_b)
    assert np.array_equal(observed_a, observed_b)


def test_observation_rejects_unfrozen_noise_or_amplitude_values() -> None:
    with pytest.raises(ScientificValidationError, match="amplitude or noise"):
        observation("fold", 0.16, 0.03, 0.01, 58050)
    with pytest.raises(ScientificValidationError, match="amplitude or noise"):
        observation("fold", 0.16, 0.04, 0.02, 58050)
    with pytest.raises(ScientificValidationError, match="margin"):
        observation("fold", 0.0, 0.04, 0.01, 58050)


def test_estimators_have_no_true_state_or_margin_inputs() -> None:
    assert tuple(inspect.signature(_estimate_free_exponential).parameters) == (
        "times",
        "observed",
    )
    assert tuple(inspect.signature(_estimate_tail_log_linear).parameters) == (
        "times",
        "observed",
    )
    assert tuple(inspect.signature(_estimate_correct_form).parameters) == (
        "model",
        "times",
        "observed",
    )


def test_finite_tail_offset_exposes_expected_transient_offset_bias() -> None:
    times = np.linspace(0.0, 10.0, 501)
    equilibrium = 1.25
    rate = 0.7
    values = equilibrium + 0.1 * np.exp(-rate * times)
    exponential = _estimate_free_exponential(times, values)
    tail = _estimate_tail_log_linear(times, values)
    assert exponential["converged"]
    assert exponential["estimated_rate"] == pytest.approx(rate, rel=1e-3)
    assert tail["converged"]
    # The terminal 10% still contains signal; treating its mean as equilibrium
    # biases this intentionally simple alternative's rate upward.
    assert tail["estimated_rate"] > rate * 1.1
    assert tail["included_points"] >= 5


def test_correct_form_diagnostic_recovers_matching_noise_free_trajectory() -> None:
    times = np.arange(TIME_COUNT, dtype=np.float64) * TIME_STEP_S
    values = exact_trajectory("transcritical", 0.16, 0.16 * 1.04, times)
    fit = _estimate_correct_form("transcritical", times, values)
    assert fit["converged"]
    assert fit["estimated_margin"] == pytest.approx(0.16, rel=1e-5)
    assert fit["estimated_equilibrium"] == pytest.approx(0.16, rel=1e-5)
    assert fit["estimated_rate"] == pytest.approx(0.16, rel=1e-5)


def test_tail_estimator_preserves_failure_reason_for_flat_or_wrong_sign_data() -> None:
    times = np.arange(100, dtype=np.float64)
    flat = _estimate_tail_log_linear(times, np.ones(100))
    rising = _estimate_tail_log_linear(times, np.exp(times / 10.0))
    assert not flat["converged"]
    assert flat["failure_reason"] == "fewer_than_five_positive_residuals"
    assert not rising["converged"]
    assert rising["failure_reason"] == "nonnegative_or_undefined_log_slope"


def test_truth_is_used_after_fit_only_for_scoring() -> None:
    condition = {
        "model": "fold",
        "margin": 0.16,
        "control_setting": 1.16,
        "amplitude_fraction": 0.04,
        "noise_fraction": 0.0,
    }
    times = np.arange(101, dtype=np.float64) * 0.1
    values = 0.4 + 0.02 * np.exp(-0.8 * times)
    truth_a = {
        "true_equilibrium": 0.4,
        "true_rate": 0.8,
        "noise_sd": 0.0,
        "snr_rms_excursion": math.inf,
    }
    truth_b = {**truth_a, "true_equilibrium": 0.5, "true_rate": 0.7}
    fit_a = _fit_record("evaluation", 58050, condition, times, values, truth_a)
    fit_b = _fit_record("evaluation", 58050, condition, times, values, truth_b)
    for a, b in zip(fit_a, fit_b, strict=True):
        assert a["estimated_rate"] == b["estimated_rate"]
        assert a["estimated_equilibrium"] == b["estimated_equilibrium"]
    assert (
        fit_a[0]["rate_relative_signed_error"] != fit_b[0]["rate_relative_signed_error"]
    )


def test_fit_record_keeps_attempts_and_rate_error_formula() -> None:
    condition = {
        "model": "fold",
        "margin": 0.16,
        "control_setting": 1.16,
        "amplitude_fraction": 0.04,
        "noise_fraction": 0.01,
    }
    times = np.linspace(0.0, 75.0, 751)
    values = 0.4 + 0.02 * np.exp(-0.8 * times)
    truth = {
        "true_equilibrium": 0.4,
        "true_rate": 0.8,
        "noise_sd": 0.004,
        "snr_rms_excursion": 0.5,
    }
    records = _fit_record("evaluation", 58050, condition, times, values, truth)
    assert len(records) == 3
    assert {row["estimator"] for row in records} == {
        "free_offset_exponential",
        "tail_offset_log_linear",
        "correct_normal_form",
    }
    reference = records[0]
    assert reference["valid"]
    expected = (reference["estimated_rate"] - 0.8) / 0.8
    assert reference["rate_relative_signed_error"] == pytest.approx(expected)
    assert reference["noise_sd_state_units"] == pytest.approx(0.004)


def test_exponent_identifiability_keeps_validity_and_coverage_separate() -> None:
    group_rows = [
        {
            "model": "fold",
            "amplitude_fraction": 0.04,
            "noise_fraction": 0.01,
            "estimator": "free_offset_exponential",
            "seed": seed,
            "valid": True,
            "estimated_q": 0.5,
            "true_q": 0.5,
            "q_bias": 0.0,
            "q_absolute_error": 0.0,
            "estimated_p": 0.7,
            "true_p": 0.5,
            "p_bias": 0.2,
            "p_absolute_error": 0.2,
            "failure_reason": "",
        }
        for seed in range(50)
    ]
    summary = _exponent_summaries(group_rows)[0]
    assert summary["attempted_groups"] == 50
    assert summary["valid_groups"] == 50
    assert summary["p_interval_covers_theory"] is False
    assert summary["p_identifiability_criterion_met"] is False
    assert summary["p_mean_absolute_error"] == pytest.approx(0.2)
    assert BOOTSTRAP_REPLICATES == 2000


def test_exponent_group_uses_five_margins_and_preserves_invalidity() -> None:
    rows = []
    for margin in (0.04, 0.09, 0.16, 0.25, 0.36):
        rows.append(
            {
                "partition": "evaluation",
                "seed": 58050,
                "model": "fold",
                "amplitude_fraction": 0.04,
                "noise_fraction": 0.0025,
                "estimator": "free_offset_exponential",
                "margin": margin,
                "control_setting": 1.0 + margin,
                "valid": True,
                "estimated_equilibrium": math.sqrt(margin),
                "estimated_rate": 2.0 * math.sqrt(margin),
            }
        )
    result = _exponent_groups(rows)[0]
    assert result["attempted_margin_count"] == 5
    assert result["valid"]
    assert result["estimated_q"] == pytest.approx(0.5, abs=1e-5)
    assert result["estimated_p"] == pytest.approx(0.5, abs=1e-5)
    rows[0]["valid"] = False
    invalid = _exponent_groups(rows)[0]
    assert not invalid["valid"]
    assert invalid["failure_reason"] == "one_or_more_invalid_rates"


def test_csv_output_schema_is_stable_and_uses_declared_column_order(
    tmp_path: Path,
) -> None:
    output = tmp_path / "schema.csv"
    _write_csv(output, [{"model": "fold", "valid": True, "seed": 58050}])
    assert output.read_text(encoding="utf-8").splitlines() == [
        "model,valid,seed",
        "fold,True,58050",
    ]


def test_zero_noise_exponent_is_not_scored_and_bootstrap_is_repeatable() -> None:
    rows = [
        {
            "model": "fold",
            "amplitude_fraction": 0.04,
            "noise_fraction": 0.0,
            "estimator": "free_offset_exponential",
            "seed": seed,
            "valid": True,
            "estimated_q": 0.5,
            "true_q": 0.5,
            "q_bias": 0.0,
            "q_absolute_error": 0.0,
            "estimated_p": 0.5,
            "true_p": 0.5,
            "p_bias": 0.0,
            "p_absolute_error": 0.0,
            "failure_reason": "",
        }
        for seed in range(50)
    ]
    summary = _exponent_summaries(rows)[0]
    assert summary["p_ci_low"] is None
    assert summary["p_interval_covers_theory"] is None
    assert summary["p_identifiability_criterion_met"] is None
    assert (
        summary["p_identifiability_criterion_status"]
        == "not_assessed_zero_noise_diagnostic"
    )
    assert _seed_bootstrap_interval([1.0, 2.0, 3.0], "fixed-label") == (
        _seed_bootstrap_interval([1.0, 2.0, 3.0], "fixed-label")
    )
