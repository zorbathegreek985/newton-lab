# Phase 64 - Non-Normal Transient Amplification

## Design gate and feasibility

### Repository inspection

The Phase 63 proposal selected a matched-spectrum comparison of stable normal and non-normal matrices and was explicitly marked unapproved at that point. This Phase 64 request authorizes the bounded study. The closest prior work includes Phase 48's fixed-end path/heat modal decay, Phase 49's heat and pendulum recovery rates, Phase 50's synthetic fold mapping, Phases 53-58 recovery-rate identifiability and closure, Phases 60-61's closed Kuramoto transition study, and Phase 62's pitchfork parameter-path study. None measures finite-time propagator norm growth for matched-spectrum normal and non-normal state matrices.

The study reuses `ODESolverConfiguration` in `src/newton_lab/simulation.py`, the existing SciPy `solve_ivp` conventions in model-specific ODE modules, NumPy/SciPy already declared in `pyproject.toml`, and the project's focused-test and report/manifest conventions. `run_experiment`'s default adapters are tied to supported model families; none handles arbitrary matrices. A narrow model-specific matrix/propagator/gain module is the minimal addition. No dependency, generic matrix framework, or change to existing API semantics was needed.

Earlier modal/recovery studies measure asymptotic or fitted modal decay; Phase 62 studies response under a changing parameter. This study measures the induced norm of a time-invariant propagator over finite time. That is a different observable and mechanism.

### State, norm, normalization, and coordinate gate

The state is an **abstract dimensionless vector** `x=(x1,x2)` in R^2. Each component is a coefficient on one declared unit basis direction. The basis directions are orthonormal and share the same dimensionless scale. They have no physical interpretation. The norm is not energy, displacement, or another dimensional quantity.

The declared norm is `||x||_2 = sqrt(x1^2 + x2^2)`. It is a standard mathematical metric for the explicitly declared abstract orthonormal coordinate space. The induced gain `G(t)=||P(t)||_2` is the largest singular value of the propagator; it is not the amplification of every initial condition. The reported maximum is a maximum over a declared finite time grid, not an exact continuous-time maximizer.

The gain is invariant under orthogonal basis changes `A -> Q.T @ A @ Q`. Arbitrary rescaling or nonorthogonal similarity transforms are excluded because they change the declared metric and can change the gain. No coordinate-independent physical amplification is claimed.

**Design gate: PASS for this abstract mathematical state space.** The shared dimensionless normalization and Euclidean metric are explicit conventions. If interpreted as a physical model, this gate would fail without a physical state scale or justified energy metric.

### Matrix family and stability

Dimensionless time is `tau=t/t_ref`; the entries below are dimensionless rates per unit `tau`:

```text
dx/dtau = A_kappa x
A0      = [[-a, 0], [0, -b]]
A_kappa = [[-a, kappa], [0, -b]]
```

The frozen values are `a=1`, `b=2`, and `kappa in {0, 2, 4, 8}`. Since the matrices are triangular, both have eigenvalues `-1` and `-2`, so both are asymptotically stable. `A0` is normal. For nonzero `kappa`, `A_kappa.T @ A_kappa` differs from `A_kappa @ A_kappa.T`, so it is non-normal. These properties were checked numerically and by tests.

## Frozen protocol

The abstract-state design gate was documented before protocol freeze. A separate analytic-only design exploration (`reports/phase_64_non_normal_transient_growth/design_exploration.json`) used the closed-form propagator and singular values; it ran no ODE integrations and was not primary evaluation. It informed the compact coupling grid by including zero coupling, a weak case below the effect threshold, and stronger cases on either side of visible transient gain. The exploratory calculation is kept as a separate artifact.

The frozen protocol is `reports/phase_64_non_normal_transient_growth/protocol.json`, SHA-256 `6AE93A4020FE12E19D1B26F5EF36DF7EE07F547ACBAB5C590A74B73B78B1D523`. It was frozen before the primary evaluation. The time horizon is `[0,4]`, four decay times of the slow mode. Primary gain sampling uses 401 points at spacing 0.01; resolution sampling uses 801 points at spacing 0.005. Peak gain and peak time are explicitly sampled-grid values. The maximum accepted gain change is `1e-4` and the maximum peak-time shift is 0.01. The declared minimum descriptive effect is an absolute gain increase of 0.01 (one percent relative to a unit-norm perturbation); it is an operational threshold, not a universal physical threshold.

Selected `solve_ivp` trajectories use DOP853, `rtol=1e-10`, `atol=1e-12`, maximum step 0.025, and 401 output times. Tighter checks for kappa 4 and 8 use `rtol=1e-12`, `atol=1e-14`, maximum step 0.005, and 801 points. Initial directions for operator-gain comparisons come from the leading right singular vector of the propagator at the primary sampled peak; its largest-magnitude coordinate is made positive for a deterministic sign convention. A coordinate-axis initial condition at kappa 8 is a secondary example, not an estimate of the operator norm. The trajectory acceptance bound is maximum Euclidean state error `1e-8` against the analytic solution.

## Analytic propagator and predictions

For `a != b`, the solution of the second state equation is `x2(tau)=exp(-b*tau)*x2(0)`. Substitution into the first equation gives the propagator:

```text
P_kappa(tau) = exp(A_kappa*tau)
             = [[exp(-a*tau), kappa*(exp(-a*tau)-exp(-b*tau))/(b-a)],
                [0,           exp(-b*tau)]]
```

This follows by integrating the forced first component:

```text
x1(tau) = exp(-a*tau)*x1(0)
        + kappa*x2(0)*integral_0^tau exp(-a*(tau-s))*exp(-b*s) ds
```

The formula gives `P(0)=I`; at `kappa=0` it reduces to the diagonal normal propagator. As `b` approaches `a`, the off-diagonal term tends to `kappa*tau*exp(-a*tau)`. The frozen family uses distinct `a` and `b`, so this limiting case is documented but not included in the evaluation.

The normal gain is `G0(tau)=max(exp(-a*tau), exp(-b*tau))=exp(-tau)` for the frozen parameters. It is 1 at time zero and strictly below 1 for positive time. Both non-normal matrices still decay asymptotically because their eigenvalues are negative; transient growth is finite-time behavior and does not contradict asymptotic stability.

## Results

The analytic closed form was evaluated on the frozen time grids. The maximum below includes `tau=0`; the positive-time maximum excludes that point.

| kappa | Primary sampled maximum | Primary peak time | Refined sampled maximum | Refined peak time | Positive-time primary maximum |
|---:|---:|---:|---:|---:|---:|
| 0 | 1.000000 | 0.000 | 1.000000 | 0.000 | 0.990050 |
| 2 | 1.000000 | 0.000 | 1.000000 | 0.000 | 0.996187 |
| 4 | 1.174353 | 0.480 | 1.174353 | 0.480 | 1.174353 |
| 8 | 2.080046 | 0.650 | 2.080092 | 0.645 | 2.080046 |

For the normal baseline, the full-horizon sampled maximum is 1 at the initial time and there is no positive-time gain above 1. The weak non-normal case `kappa=2` also has no sampled positive-time amplification above 1.01. At `kappa=4` and `8`, the non-normal family exceeds 1.01 and the matched normal maximum by more than the frozen 0.01 threshold. The peak difference between primary and refined grids is zero for `kappa=0,2,4` and `4.68e-5` for `kappa=8`; the latter peak-time shift is 0.005. All meet the frozen resolution criteria.

Operator gain is not a particular trajectory's amplification. At `kappa=8`, the maximizing unit initial direction was approximately `(0.248949, 0.968516)` and its numerical trajectory reached norm 2.080046 near `tau=0.65`. The coordinate direction `(0,1)` reached about 2.016054 near `tau=0.68`, less than the operator gain. This confirms why one trajectory cannot stand in for the induced norm.

## Numerical verification and hypothesis outcomes

The closed form was compared with the independent `scipy.linalg.expm` implementation at four fixed times for all four couplings (16 comparisons); the maximum matrix-entry difference was `8.44e-15`, below `1e-12`. Eight semigroup comparisons had maximum difference `4.44e-16`. The zero-coupling gain curve differed from the direct diagonal normal reference by `2.22e-16`, below `1e-14`. Tests also check identity, the defining matrix differential equation by centered finite differences, and semigroup behavior.

Six `solve_ivp` integrations were attempted and all succeeded: four primary trajectories and two tighter comparisons. The largest maximum-state error against the analytic trajectory was `1.55e-15`, below `1e-8`. Interpolating the tighter trajectories to the primary output grid gave maximum primary-versus-tight trajectory differences of `8.95e-16` (`kappa=4`) and `1.11e-15` (`kappa=8`). The trajectories verify selected initial conditions; the singular-value calculation supplies the operator gain. Tighter tolerances confirm numerical agreement only; they say nothing about physical validity.

- **H1, finite-time amplification:** supported within the frozen family. `kappa=4` and `8` exceed 1.01 on both grids; `kappa=0` and `2` do not. The result is not a statement about all stable systems.
- **H2, greater gain than matched normal:** supported within the frozen family. The non-normal sampled maxima at `kappa=4` and `8` exceed the normal maximum 1 by at least 0.01 on both grids. The `kappa=2` case does not.
- **H3, analytic/numerical agreement:** supported within frozen tolerances. Propagator, semigroup, zero-coupling, grid-resolution, and all six trajectory checks met their criteria.
- **Normal baseline check:** supported; its gain never exceeds 1 on the sampled grid.

These are an analytic result for a declared matrix family and numerical verification of its implementation. Numerical agreement does not prove a theorem. No physical model, engineering system, financial application, or cross-domain utility was tested. The work makes no claim that all stable systems amplify, that eigenvalues are useless, or that this repository-level question is novel in the wider literature.

## Artifacts and reproduction

Implementation: [non_normal_transient_growth.py](../../src/newton_lab/non_normal_transient_growth.py). Focused tests: [test_non_normal_transient_growth.py](../../tests/test_non_normal_transient_growth.py). The output directory contains the separate `design_exploration.json`, frozen [protocol.json](../../reports/phase_64_non_normal_transient_growth/protocol.json), `gain_curves.csv`, `gain_summary.csv`, `analytic_reference_checks.csv`, `matrix_structure_checks.csv`, `trajectory_checks.csv`, `trajectory_samples.csv`, `hypothesis_results.json`, `metadata.json`, and `sha256_manifest.json`. Tables are sufficient to inspect this bounded study; no figure was needed.

To reproduce without overwriting the canonical outputs, copy the frozen protocol and design-exploration files into a new empty directory and call `run_non_normal_transient_growth_study` there. The isolated rerun is validated by comparing every manifest entry and the manifest itself. The canonical runner refuses to overwrite an output directory containing result files.

## Verification record

- Focused tests: `.venv\Scripts\python.exe -m pytest --basetemp .pytest-temp tests\test_non_normal_transient_growth.py -q` - **14 passed**.
- Full suite: `.venv\Scripts\python.exe -m pytest --basetemp .pytest-temp -q` - **614 passed**.
- `.venv\Scripts\ruff.exe check .` - passed.
- `.venv\Scripts\ruff.exe format --check .` - passed; 173 files already formatted.
- `.venv\Scripts\mypy.exe src` - passed; 55 source files.
- `.venv\Scripts\mypy.exe .` - passed; 105 files.
- The SHA-256 manifest validated all **10/10** entries. An isolated rerun with the final source matched every output hash and produced a byte-identical manifest.
- The Phase 51 Atlas JSON parsed, required historical statuses were preserved, and Atlas/report local links resolved. No Phase 63 proposal content was changed.

## Limitations and stopping decision

The result is specific to a deterministic two-state triangular matrix family, `a=1`, `b=2`, four tested coupling values, horizon `[0,4]`, a finite gain grid, and the declared Euclidean norm. Rescaling state coordinates can change the gain; only orthogonal changes preserve this comparison. The family is an abstract construction, not a fitted or measured physical model. The study does not address noise, nonlinear saturation, uncertainty in matrix coefficients, higher dimensions, or application value. A sampled-grid maximum may miss a narrow peak, though the prescribed half-step refinement met its criterion for this grid.

Phase 64 is complete as a bounded analytic-numerical study. Stop here: do not expand the coupling grid to seek larger amplification, do not reopen Phases 60-62, and do not begin Phase 65 automatically. Phase 32 remains `BLOCKED_AUTHORIZATION`, Phase 52 remains `NO-GO`, the Kuramoto thread remains closed, and Phase 62 remains complete.
