import json
from pathlib import Path

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.kuramoto_synchronization_study import (
    _frequency_sample,
    _transition_estimate,
    critical_coupling,
    run_synchronization_study,
)

PROTOCOL_PATH = Path("reports/phase_60_kuramoto_synchronization/protocol.json")


def test_continuum_critical_couplings_match_density_at_zero() -> None:
    assert critical_coupling("gaussian_sigma_1") == pytest.approx(np.sqrt(8.0 / np.pi))
    assert critical_coupling("truncated_lorentzian_gamma_1_cutoff_5") == pytest.approx(
        4.0 * np.arctan(5.0) / np.pi
    )


def test_transition_estimator_returns_sweep_interval_and_paired_bootstrap() -> None:
    coupling = np.linspace(0.5, 1.5, 11)
    base = 0.08 + 0.7 / (1.0 + np.exp(-12.0 * (coupling - 1.0)))
    curves = np.vstack(
        [
            base + offset
            for offset in (-0.02, -0.01, 0.0, 0.01, 0.02, 0.005, -0.005, 0.0)
        ]
    )

    result = _transition_estimate(curves, coupling, bootstrap_seed=24680)
    repeated = _transition_estimate(curves, coupling, bootstrap_seed=24680)

    assert result == repeated
    assert coupling[0] <= result[0] <= coupling[-1]
    assert result[1] <= result[0] <= result[2]
    assert not result[3]


def test_protocol_is_frozen_and_records_population_and_seed_design() -> None:
    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))

    assert protocol["protocol_version"] == "1.1"
    assert protocol["protocol_status"] == "frozen_before_final_evaluation"
    assert protocol["finite_population_sizes"] == [64, 256]
    assert len(protocol["replicate_seeds"]) == 8
    assert len(protocol["coupling_sweep_factors_of_continuum_Kc"]) == 11
    assert protocol["frequency_populations"]["truncated_lorentzian_gamma_1_cutoff_5"][
        "continuum_critical_coupling_rad_per_second"
    ] == pytest.approx(critical_coupling("truncated_lorentzian_gamma_1_cutoff_5"))
    assert protocol["phase32_status"] == "BLOCKED_AUTHORIZATION"
    assert protocol["phase52_status"] == "NO-GO"


def test_frequency_samples_follow_the_declared_support_and_seed() -> None:
    population = "truncated_lorentzian_gamma_1_cutoff_5"
    first = _frequency_sample(population, 256, 60001, 1)
    repeated = _frequency_sample(population, 256, 60001, 1)

    np.testing.assert_array_equal(first[0], repeated[0])
    np.testing.assert_array_equal(first[1], repeated[1])
    assert np.max(np.abs(first[0])) <= 5.0


def test_runner_rejects_modified_frozen_protocol_without_overwriting(
    tmp_path: Path,
) -> None:
    protocol_path = tmp_path / "protocol.json"
    protocol_path.write_text('{"phase": 60}\\n', encoding="utf-8")
    with pytest.raises(ScientificValidationError, match="protocol hash mismatch"):
        run_synchronization_study(tmp_path)

    assert sorted(path.name for path in tmp_path.iterdir()) == ["protocol.json"]


def test_runner_refuses_existing_result_files(tmp_path: Path) -> None:
    protocol_text = PROTOCOL_PATH.read_text(encoding="utf-8")
    (tmp_path / "protocol.json").write_text(protocol_text, encoding="utf-8")
    existing = tmp_path / "user_notes.txt"
    existing.write_text("preserve", encoding="utf-8")

    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        run_synchronization_study(tmp_path)

    assert existing.read_text(encoding="utf-8") == "preserve"
