"""Scientific behavior tests for the Phase 54 synthetic ablation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.exponent_identifiability_ablation import (
    FACTOR_LEVELS,
    SIGNATURES,
    _bootstrap_interval,
    _estimate_rate,
    _fit_group,
    _manifest,
    _ols_power,
    _paired_observation,
    _write_csv,
    factor_grid,
    paired_noise_innovations,
    validate_protocol,
)
from newton_lab.recovery_scaling_identifiability import (
    MODEL_NAMES,
    TRUE_MARGINS,
    local_recovery_rate,
    stable_equilibrium,
)


@pytest.mark.parametrize(
    ("model", "q", "p", "margin", "equilibrium", "rate"),
    [
        ("fold", 0.5, 0.5, 0.16, 0.4, 0.8),
        ("transcritical", 1.0, 1.0, 0.16, 0.16, 0.16),
        ("supercritical_pitchfork", 0.5, 1.0, 0.16, 0.4, 0.32),
    ],
)
def test_phase_53_theoretical_exponents_and_rates_are_reused(
    model: str,
    q: float,
    p: float,
    margin: float,
    equilibrium: float,
    rate: float,
) -> None:
    assert SIGNATURES[model] == (q, p)
    assert stable_equilibrium(model, margin) == pytest.approx(equilibrium)
    assert local_recovery_rate(model, margin) == pytest.approx(rate)


def test_factorial_has_72_unique_combinations_and_frozen_levels() -> None:
    grid = factor_grid()
    assert len(grid) == 72
    assert len({tuple(row.values()) for row in grid}) == 72
    assert {
        name: tuple(dict.fromkeys(row[name] for row in grid)) for name in FACTOR_LEVELS
    } == FACTOR_LEVELS


def test_protocol_levels_and_seed_partition_are_validated() -> None:
    protocol_path = (
        Path("reports") / "phase_54_exponent_identifiability_ablation" / "protocol.json"
    )
    protocol = json.loads(protocol_path.read_text(encoding="utf-8-sig"))
    validate_protocol(protocol)
    protocol["factors"]["amplitude_fraction_of_equilibrium"] = [0.2]
    with pytest.raises(ScientificValidationError, match="factor levels"):
        validate_protocol(protocol)


def test_replicate_noise_streams_are_paired_and_deterministic() -> None:
    a = paired_noise_innovations(55000, "fold", 0.04)
    b = paired_noise_innovations(55000, "fold", 0.04)
    c = paired_noise_innovations(55000, "transcritical", 0.04)
    np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(a[:11], a[:31][:11])
    assert not np.array_equal(a, c)
    with pytest.raises(ScientificValidationError):
        paired_noise_innovations(55000, "unknown", 0.04)


def test_calibration_and_evaluation_seed_ranges_are_disjoint() -> None:
    calibration = set(range(54000, 54008))
    evaluation = set(range(55000, 55020))
    assert calibration.isdisjoint(evaluation)


def test_oracle_equilibrium_is_only_supplied_to_fixed_offset_fit() -> None:
    reference = _estimate_rate("fold", 0.09, 0.2, 0.01, 3.0, "reference", 57214)
    oracle = _estimate_rate("fold", 0.09, 0.2, 0.01, 3.0, "oracle_equilibrium", 57214)
    assert oracle["estimated_equilibrium"] == oracle["true_equilibrium"]
    assert reference["estimated_equilibrium"] != oracle["true_equilibrium"]
    assert oracle["rate_valid"] and reference["rate_valid"]


@pytest.mark.parametrize("information", ["reference", "oracle_equilibrium"])
def test_zero_noise_leaves_observations_exact_and_rate_estimable(
    information: str,
) -> None:
    _, clean, observed = _paired_observation("transcritical", 0.16, 0.05, 0.0, 3.0, 101)
    result = _estimate_rate("transcritical", 0.16, 0.05, 0.0, 3.0, information, 101)
    assert result["noise_sd_state"] == 0.0
    assert result["rate_valid"]
    np.testing.assert_array_equal(observed, clean)
    # A finite-amplitude nonlinear transient is not exactly one exponential.
    assert result["rate_relative_error"] > 0.0


def test_fixed_boundary_power_fit_recovers_known_signature() -> None:
    margins = np.asarray(TRUE_MARGINS)
    values = 2.5 * margins**0.5
    exponent, intercept = _ols_power(margins, values)
    assert exponent == pytest.approx(0.5, abs=1e-12)
    assert np.exp(intercept) == pytest.approx(2.5)


def test_group_validity_retains_rate_failure_and_does_not_score_it() -> None:
    records = []
    for margin in TRUE_MARGINS:
        equilibrium = stable_equilibrium("fold", margin)
        records.append(
            {
                "control_setting": 1.0 + margin,
                "margin": margin,
                "rate_valid": True,
                "rate_relative_error": 0.0,
                "estimated_equilibrium": equilibrium,
                "estimated_rate": local_recovery_rate("fold", margin),
                "true_equilibrium": equilibrium,
            }
        )
    result = _fit_group("fold", "reference", records)
    assert result["valid"]
    assert result["attempted_margin_count"] == 5
    records[0]["rate_valid"] = False
    result = _fit_group("fold", "reference", records)
    assert not result["valid"]
    assert np.isnan(result["q_abs_error"])
    assert result["rate_valid_fraction"] == pytest.approx(0.8)


def test_oracle_equilibrium_and_boundary_affect_only_the_declared_stage() -> None:
    records = []
    for margin in TRUE_MARGINS:
        equilibrium = stable_equilibrium("fold", margin)
        records.append(
            {
                "control_setting": 1.0 + margin,
                "margin": margin,
                "rate_valid": True,
                "rate_relative_error": 0.0,
                "estimated_equilibrium": 1.1 * equilibrium * margin**0.1,
                "estimated_rate": local_recovery_rate("fold", margin),
                "true_equilibrium": equilibrium,
            }
        )
    oracle_boundary = _fit_group("fold", "oracle_boundary", records)
    both = _fit_group("fold", "oracle_both", records)
    assert oracle_boundary["estimated_boundary"] == pytest.approx(1.0)
    assert oracle_boundary["q"] == pytest.approx(0.6, abs=1e-10)
    assert both["q"] == pytest.approx(0.5, abs=1e-10)
    assert oracle_boundary["p"] == pytest.approx(0.5, abs=1e-10)
    assert both["p"] == pytest.approx(0.5, abs=1e-10)


def test_seed_bootstrap_interval_is_deterministic_and_covers_constant() -> None:
    values = np.asarray([0.5] * 12)
    first = _bootstrap_interval(values, seed=91)
    second = _bootstrap_interval(values, seed=91)
    assert first == second == pytest.approx((0.5, 0.5))


def test_output_schema_and_manifest_are_stable(tmp_path: Path) -> None:
    output = tmp_path / "report.csv"
    _write_csv(output, [{"seed": 1, "valid": True}, {"seed": 2, "valid": False}])
    assert output.read_text(encoding="utf-8").splitlines()[0] == "seed,valid"
    manifest = _manifest(tmp_path)
    assert manifest == {"report.csv": manifest["report.csv"]}
    assert len(manifest["report.csv"]) == 64


def test_three_phase_53_models_are_the_only_models() -> None:
    assert tuple(SIGNATURES) == MODEL_NAMES
