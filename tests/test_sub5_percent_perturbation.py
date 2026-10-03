"""Scientific behavior tests for the Phase 55 small-perturbation study."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.recovery_scaling_identifiability import (
    MODEL_NAMES,
    TRUE_MARGINS,
    local_recovery_rate,
    stable_equilibrium,
)
from newton_lab.sub5_percent_perturbation import (
    AMPLITUDES,
    BOOTSTRAP_REPLICATES,
    CALIBRATION_SEEDS,
    EVALUATION_SEEDS,
    REFERENCE_AMPLITUDES,
    SAMPLE_COUNT,
    SAMPLE_INTERVAL_TAU,
    SIGNATURES,
    SUB5_AMPLITUDES,
    WINDOW_TAU,
    _bootstrap_interval,
    _deterministic_snr_slope,
    _evaluation_only,
    _exponent_summary,
    amplitude_specification,
    estimate_rate,
    fit_exponent_group,
    paired_innovations,
    paired_observation,
    signal_to_noise,
    validate_protocol,
)


def test_phase_53_model_rates_and_exponents_are_reused() -> None:
    assert SIGNATURES == {
        "fold": (0.5, 0.5),
        "transcritical": (1.0, 1.0),
        "supercritical_pitchfork": (0.5, 1.0),
    }
    for model in MODEL_NAMES:
        assert stable_equilibrium(model, 0.16) > 0.0
        assert local_recovery_rate(model, 0.16) > 0.0


def test_frozen_amplitudes_and_exact_displacements_cover_all_settings() -> None:
    assert SUB5_AMPLITUDES == (0.0025, 0.005, 0.01, 0.02, 0.04)
    assert REFERENCE_AMPLITUDES == (0.05, 0.20)
    assert AMPLITUDES == SUB5_AMPLITUDES + REFERENCE_AMPLITUDES
    rows = amplitude_specification()
    assert len(rows) == len(MODEL_NAMES) * len(TRUE_MARGINS) * len(AMPLITUDES)
    for row in rows:
        assert row["initial_displacement"] == pytest.approx(
            row["amplitude_fraction_of_equilibrium"] * row["true_equilibrium"]
        )
        assert row["initial_state"] > row["true_equilibrium"] > 0.0
        assert row["displacement_over_margin"] > 0.0
        assert row["is_sub_5_percent_primary"] == (
            row["amplitude_fraction_of_equilibrium"] in SUB5_AMPLITUDES
        )


def test_all_amplitude_model_margin_trajectories_stay_on_positive_recovery_side() -> (
    None
):
    for model in MODEL_NAMES:
        for margin in TRUE_MARGINS:
            for amplitude in AMPLITUDES:
                times, clean, observed, _ = paired_observation(
                    model, margin, amplitude, "noisy", 56000
                )
                assert len(times) == SAMPLE_COUNT == 31
                np.testing.assert_allclose(
                    np.diff(times),
                    SAMPLE_INTERVAL_TAU / local_recovery_rate(model, margin),
                )
                assert clean[0] > stable_equilibrium(model, margin) > 0.0
                assert np.all(clean > 0.0)
                assert clean[-1] < clean[0]
                assert observed.shape == clean.shape


def test_fixed_observation_condition_and_reference_information_are_frozen() -> None:
    protocol = json.loads(
        Path("reports/phase_55_sub5_percent_perturbation/protocol.json").read_text(
            encoding="utf-8"
        )
    )
    validate_protocol(protocol)
    assert protocol["protocol_version"] == "1.2"
    assert protocol["primary_observation"]["noise_sd_fraction_of_equilibrium"] == 0.01
    assert protocol["primary_observation"]["window_tau"] == WINDOW_TAU == 3.0
    assert protocol["primary_observation"]["sampling_interval_tau"] == 0.1
    assert protocol["primary_observation"]["information_condition"] == "reference"
    assert (
        protocol["diagnostic_arm"]["only_changed_setting"] == "observation_noise_sd=0"
    )
    protocol["amplitude"]["sub_5_percent_primary_levels"] = [0.04]
    with pytest.raises(ScientificValidationError, match="amplitude protocol"):
        validate_protocol(protocol)


def test_deterministic_snr_trend_is_exact_and_seed_replication_independent() -> None:
    evaluation_rows = [
        {
            "partition": "evaluation",
            "model": "fold",
            "noise_arm": "noisy",
            "amplitude_fraction_of_equilibrium": amplitude,
            "snr_rms_clean_excursion_over_noise_sd": float(amplitude * 10),
        }
        for amplitude in SUB5_AMPLITUDES
    ]
    replicated = evaluation_rows + [
        {**row, "replicate_seed": 57000 + index}
        for index, row in enumerate(evaluation_rows)
    ]
    assert _deterministic_snr_slope(replicated, "fold") == pytest.approx(
        _deterministic_snr_slope(evaluation_rows, "fold")
    )
    assert math.isfinite(_deterministic_snr_slope(replicated, "fold"))


def test_seed_partition_is_disjoint_and_innovations_pair_amplitudes() -> None:
    assert set(CALIBRATION_SEEDS).isdisjoint(EVALUATION_SEEDS)
    assert len(EVALUATION_SEEDS) == 50
    first = paired_innovations(56001, "fold", 0.09)
    again = paired_innovations(56001, "fold", 0.09)
    different = paired_innovations(56002, "fold", 0.09)
    np.testing.assert_array_equal(first, again)
    assert not np.array_equal(first, different)


def test_calibration_rows_are_excluded_from_evaluation_summaries() -> None:
    rows = [
        {"partition": "calibration", "replicate_seed": 56000},
        {"partition": "evaluation", "replicate_seed": 57000},
    ]
    assert _evaluation_only(rows) == [rows[1]]


def test_observations_are_paired_and_zero_noise_is_exact() -> None:
    t_small, clean_small, noisy_small, _ = paired_observation(
        "transcritical", 0.16, 0.005, "noisy", 56003
    )
    t_large, clean_large, noisy_large, _ = paired_observation(
        "transcritical", 0.16, 0.04, "noisy", 56003
    )
    _, _, exact_zero, _ = paired_observation(
        "transcritical", 0.16, 0.005, "zero_noise", 56003
    )
    np.testing.assert_array_equal(t_small, t_large)
    np.testing.assert_array_equal(noisy_small - clean_small, noisy_large - clean_large)
    np.testing.assert_array_equal(exact_zero, clean_small)


def test_predeclared_snr_definition_and_edge_cases() -> None:
    clean = np.asarray([2.0, 1.5, 1.0])
    expected = float(np.sqrt(np.mean(np.square(clean - 1.0))) / 0.25)
    assert signal_to_noise(clean, 1.0, 0.25) == pytest.approx(expected)
    assert math.isinf(signal_to_noise(clean, 1.0, 0.0))
    with pytest.raises(ScientificValidationError):
        signal_to_noise(clean, 1.0, -0.1)
    with pytest.raises(ScientificValidationError):
        signal_to_noise(np.asarray([math.nan]), 1.0, 0.1)


def test_rate_estimation_is_deterministic_and_reports_rate_error() -> None:
    a = estimate_rate("fold", 0.09, 0.01, "noisy", 56004)
    b = estimate_rate("fold", 0.09, 0.01, "noisy", 56004)
    assert a == b
    assert a["rate_valid"]
    assert a["rate_absolute_error"] == pytest.approx(
        abs(a["estimated_rate"] - a["true_rate"])
    )
    assert a["rate_relative_signed_error"] == pytest.approx(
        (a["estimated_rate"] - a["true_rate"]) / a["true_rate"]
    )
    assert a["snr_rms_clean_excursion_over_noise_sd"] > 0.0


def test_zero_noise_uses_unchanged_estimator_and_exact_data() -> None:
    result = estimate_rate("transcritical", 0.16, 0.0025, "zero_noise", 56005)
    assert result["rate_valid"]
    assert result["noise_sd_state_units"] == 0.0
    assert result["estimator_noise_scale_state_units"] == pytest.approx(0.01 * 0.16)
    assert math.isinf(result["snr_rms_clean_excursion_over_noise_sd"])


def test_failed_rate_fit_invalidates_exponent_group_and_stays_unscored() -> None:
    records = []
    for margin in TRUE_MARGINS:
        eq = stable_equilibrium("fold", margin)
        records.append(
            {
                "margin": margin,
                "control_setting": 1.0 + margin,
                "rate_valid": True,
                "estimated_equilibrium": eq,
                "estimated_rate": local_recovery_rate("fold", margin),
            }
        )
    result = fit_exponent_group("fold", "noisy", 0.01, 56006, records)
    assert result["valid"]
    records[0]["rate_valid"] = False
    failed = fit_exponent_group("fold", "noisy", 0.01, 56006, records)
    assert not failed["valid"]
    assert failed["rate_valid_fraction"] == pytest.approx(0.8)
    assert math.isnan(failed["p_absolute_error"])


def test_exponent_summary_keeps_invalid_groups_in_denominator() -> None:
    groups = []
    for index, seed in enumerate(EVALUATION_SEEDS):
        groups.append(
            {
                "model": "fold",
                "amplitude_fraction_of_equilibrium": 0.01,
                "noise_arm": "noisy",
                "replicate_seed": seed,
                "valid": index % 2 == 0,
                "q": 0.5 if index % 2 == 0 else math.nan,
                "p": 0.5 if index % 2 == 0 else math.nan,
                "q_absolute_error": 0.0 if index % 2 == 0 else math.nan,
                "p_absolute_error": 0.0 if index % 2 == 0 else math.nan,
                "estimated_boundary": 1.0,
                "rate_valid_fraction": 1.0,
            }
        )
    summary = next(
        row
        for row in _exponent_summary(groups)
        if row["model"] == "fold"
        and row["amplitude_fraction_of_equilibrium"] == 0.01
        and row["noise_arm"] == "noisy"
    )
    assert summary["attempted_groups"] == 50
    assert summary["valid_exponent_groups"] == 25
    assert summary["valid_exponent_group_fraction"] == pytest.approx(0.5)
    assert summary["p_ci_covers_truth"]
    assert summary["p_ci95_width"] == pytest.approx(0.0)
    assert not summary["p_specific_identifiability_criterion"]
    assert summary["identifiability_criterion_applicability"] == "primary"


def test_zero_noise_exponent_criterion_is_diagnostic_not_scored() -> None:
    groups = [
        {
            "model": "fold",
            "amplitude_fraction_of_equilibrium": 0.01,
            "noise_arm": "zero_noise",
            "replicate_seed": seed,
            "valid": True,
            "q": 0.5,
            "p": 0.5,
            "q_absolute_error": 0.0,
            "p_absolute_error": 0.0,
            "estimated_boundary": 1.0,
            "rate_valid_fraction": 1.0,
        }
        for seed in EVALUATION_SEEDS
    ]
    row = next(
        item
        for item in _exponent_summary(groups)
        if item["model"] == "fold"
        and item["amplitude_fraction_of_equilibrium"] == 0.01
        and item["noise_arm"] == "zero_noise"
    )
    assert row["valid_exponent_groups"] == 50
    assert (
        row["identifiability_criterion_applicability"] == "diagnostic_only_not_scored"
    )
    assert row["p_specific_identifiability_criterion"] is None
    assert not row["ci_coverage_interpretable"]


def test_seed_bootstrap_interval_is_deterministic_and_covers_constant() -> None:
    values = [1.0] * 20
    first = _bootstrap_interval(values, 101)
    second = _bootstrap_interval(values, 101)
    assert first == second == pytest.approx((1.0, 1.0))
    assert BOOTSTRAP_REPLICATES == 2000
