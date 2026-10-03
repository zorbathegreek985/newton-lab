import hashlib
import json
import shutil
from pathlib import Path
from typing import TypedDict

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.pitchfork_path_dependence import (
    PitchforkCycle,
    RampLeg,
    _loop_area,
    pitchfork_equilibrium_magnitude,
    run_pitchfork_path_study,
    simulate_pitchfork_cycle,
)
from newton_lab.simulation import ODESolverConfiguration

PROTOCOL = Path("reports/phase_62_pitchfork_path_dependence/protocol.json")


class _CycleKwargs(TypedDict):
    mu_min: float
    mu_max: float
    leg_duration: float
    output_interval: float
    solver: ODESolverConfiguration


def _solver() -> ODESolverConfiguration:
    return ODESolverConfiguration(
        method="DOP853",
        relative_tolerance=1e-9,
        absolute_tolerance=1e-11,
        maximum_step=0.05,
    )


def test_pitchfork_stable_equilibrium_magnitude_matches_branches() -> None:
    values = np.array([-0.2, 0.0, 0.04, 0.16])

    np.testing.assert_allclose(
        pitchfork_equilibrium_magnitude(values), [0.0, 0.0, 0.2, 0.4]
    )
    assert pitchfork_equilibrium_magnitude(0.09) == pytest.approx(0.3)
    with pytest.raises(ScientificValidationError):
        pitchfork_equilibrium_magnitude(float("nan"))


def test_up_down_ramp_is_reproducible_and_initial_sign_selects_branch() -> None:
    kwargs: _CycleKwargs = {
        "mu_min": -0.16,
        "mu_max": 0.16,
        "leg_duration": 20.0,
        "output_interval": 0.1,
        "solver": _solver(),
    }
    positive = simulate_pitchfork_cycle(initial_state=0.01, **kwargs)
    positive_again = simulate_pitchfork_cycle(initial_state=0.01, **kwargs)
    negative = simulate_pitchfork_cycle(initial_state=-0.01, **kwargs)

    np.testing.assert_array_equal(positive.upward.state, positive_again.upward.state)
    np.testing.assert_array_equal(
        positive.downward.state, positive_again.downward.state
    )
    assert positive.upward.state[-1] > 0.0
    assert negative.upward.state[-1] < 0.0
    np.testing.assert_allclose(
        positive.upward.state, -negative.upward.state, atol=1e-11
    )
    np.testing.assert_allclose(
        positive.downward.state, -negative.downward.state, atol=1e-11
    )
    assert positive.downward.mu[0] == pytest.approx(0.16)
    assert positive.downward.mu[-1] == pytest.approx(-0.16)


def test_exact_zero_is_an_invariant_initial_state_control() -> None:
    cycle = simulate_pitchfork_cycle(
        mu_min=-0.16,
        mu_max=0.16,
        initial_state=0.0,
        leg_duration=20.0,
        output_interval=0.1,
        solver=_solver(),
    )

    assert np.max(np.abs(cycle.upward.state)) == 0.0
    assert np.max(np.abs(cycle.downward.state)) == 0.0


def test_equilibrium_branch_has_zero_quasistatic_loop_area() -> None:
    mu = np.array([-0.16, 0.0, 0.16])
    state = np.asarray(pitchfork_equilibrium_magnitude(mu))
    upward = RampLeg("up", mu, np.arange(3.0), state, 0)
    downward = RampLeg("down", mu[::-1], np.arange(3.0), state[::-1], 0)
    cycle = PitchforkCycle(1.0, 0.1, upward, downward)

    assert _loop_area(cycle) == pytest.approx(0.0, abs=1e-15)


def test_frozen_protocol_and_end_to_end_artifact_manifest(tmp_path: Path) -> None:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    protocol_hash = hashlib.sha256(PROTOCOL.read_bytes()).hexdigest().upper()
    assert (
        protocol_hash
        == "A2D35E4AD3B756AA6071AC3C7B4F0141ACEC05BC9DDE9BFE09981997A66351E6"
    )
    assert protocol["protocol_status"] == "frozen_before_primary_evaluation"
    assert protocol["kuramoto_thread_status"] == (
        "CLOSED_INCONCLUSIVE_UNDER_DECLARED_PROTOCOL"
    )

    output = tmp_path / "phase62"
    output.mkdir()
    shutil.copyfile(PROTOCOL, output / "protocol.json")
    metadata = run_pitchfork_path_study(output)

    assert metadata["primary_integrations"] == 12
    assert metadata["numerical_resolution_integrations"] == 8
    assert metadata["total_successful_integrations"] == 20
    assert metadata["path_metric_rows"] == 6
    assert metadata["numerical_check_rows"] == 4
    assert metadata["phase32_status"] == "BLOCKED_AUTHORIZATION"
    assert metadata["phase52_status"] == "NO-GO"

    manifest = json.loads((output / "sha256_manifest.json").read_text(encoding="utf-8"))
    for filename, expected_hash in manifest.items():
        actual_hash = (
            hashlib.sha256((output / filename).read_bytes()).hexdigest().upper()
        )
        assert actual_hash == expected_hash

    findings = json.loads(
        (output / "hypothesis_results.json").read_text(encoding="utf-8")
    )
    allowed = {
        "Supported within the tested conditions",
        "Not supported within the tested conditions",
        "Inconclusive",
    }
    assert findings["H1_finite_rate_path_difference"]["classification"] in allowed
    assert findings["H2_slower_ramp_reduces_loop"]["classification"] in allowed
    assert (
        findings["H3_initial_sign_selects_pitchfork_branch"]["classification"]
        in allowed
    )
    assert (
        findings["zero_initial_state_invariant_control"]["observed_invariant"] is True
    )
    assert findings["quasi_static_equilibrium_hysteresis"]["classification"].startswith(
        "Not supported"
    )
