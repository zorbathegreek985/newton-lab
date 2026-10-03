# Phase 48 — Pinned Path Modal Decay and Convergence Validation

## Research question and scope

Does continuous-time neighbor relaxation on a uniformly spaced path reproduce
fixed-end heat-mode rates, and do its exact discrete modal rates converge to
the continuum rates at second order as the path is refined? The target system
is a mathematical network model. Its states are dimensionless agent estimates;
they are not temperatures, financial prices, or market observations.

Phase 47 proposed this specific test. This implementation tests the path model
and its modal predictions only. It does not evaluate arbitrary networks,
asynchronous protocols, delayed communication, or a practical application.
Phase 32 remains `BLOCKED_AUTHORIZATION`.

## Equations and assumptions

There are `N` interior agents and two leaders held at zero. Let
`h = 1/(N+1)`, so the leader locations are at either end of a unit interval.
For `i=1,...,N`, the dimensionless continuous-time update law is

```text
dz_i/dτ = (z_(i-1) - 2 z_i + z_(i+1)) / h²,
z_0 = z_(N+1) = 0.
```

Equivalently, `dz/dτ = -A z`, where `A` has diagonal `2/h²` and adjacent
off-diagonal entries `-1/h²`. The leaders are boundary conditions, not
integrated state variables. The model assumes a uniform path, equal symmetric
neighbor coupling, synchronous linear updates, fixed leaders, and no delays,
noise, dropouts, or state-dependent weights.

For mode `n`, the initial state is `z_i(0)=sin(nπ i h)`. Its exact discrete
eigenvalue is

```text
λ_(N,n) = 4 h^(-2) sin²(nπh/2).
```

The independent continuum reference uses a unit interval, unit diffusivity,
and homogeneous fixed ends: `λ_n=(nπ)²`. This is Newton Lab's existing
`heat_mode_decay_rate(1, n, 1)` API. It is not used to rescale the path
operator. The rate ratio is exactly
`[sin(nπh/2)/(nπh/2)]²`; for fixed `n`, the leading relative error is
`1-λ_(N,n)/λ_n ≈ (nπh)²/12`.

The exact discrete trajectory is `z(τ)=z(0) exp(-λ_(N,n)τ)`. This analytic
trajectory is separate from the numerical solution and from the continuum
trajectory. The numerical solution uses SciPy `solve_ivp` with DOP853. Samples
cover `s=λ_n τ` from 0 to 5. Modal rates are estimated by least-squares fitting
the log of the projection onto the initial sine mode against `τ`.

## Frozen design and reproducibility

The implementation is in `src/newton_lab/pinned_path_network.py`; automated
tests are in `tests/test_pinned_path_network.py`. Run the study from the
workspace root with:

```powershell
.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.pinned_path_network import run_pinned_path_study; run_pinned_path_study(Path('reports/phase_48_pinned_path_convergence'))"
```

The sweep uses `N = 7, 15, 31, 63, 127`, modes `n = 1, 2, 3`, 501 and 2,001
output samples, and both `(rtol, atol)=(1e-9,1e-11)` and
`(1e-12,1e-14)`. It contains 60 deterministic integrations. There is no
randomness or external dataset. The reported order is ordinary least-squares
on `log(relative rate error)` against `log(h)` over the finest four grid sizes.
Before interpreting the result, the declared criterion required, separately
for each mode: all discrete rates below the continuum rate, error strictly
decreasing at every refinement, and fitted order in `[1.7, 2.3]`.

Generated files are under
`reports/phase_48_pinned_path_convergence/`:

- `per_condition_results.csv` records all solver configurations, numerical
  rate fits, errors against the exact discrete solution, and model references.
- `convergence_summary.csv` records each mode's five rate ratios and spatial
  errors, fitted order, and criterion outcomes.
- `metadata.json` records equations, grid, solver, tolerances, acceptance
  criteria, environment, and whether all modes met the prediction.
- `modal_rate_convergence.png` shows the exact discrete spatial-rate errors.
- `modal_trajectory_comparison.png` compares numerical and exact discrete
  mode trajectories for selected grids and mode 1.

## Results

The table reports the exact discrete rate from the formula, its ratio to the
continuum rate, and relative spatial rate error. Values are independent of
solver tolerance and output density.

| Mode | N | `λ_(N,n)` | `λ_(N,n)/λ_n` | Relative spatial error |
|---:|---:|---:|---:|---:|
| 1 | 7 | 9.743419839 | 0.987214831 | 0.012785169 |
| 1 | 15 | 9.837936434 | 0.996791364 | 0.003208636 |
| 1 | 31 | 9.861679775 | 0.999197068 | 0.000802932 |
| 1 | 63 | 9.867622767 | 0.999799219 | 0.000200781 |
| 1 | 127 | 9.869108963 | 0.999949802 | 0.000050198 |
| 2 | 7 | 37.49033201 | 0.949641204 | 0.050358796 |
| 2 | 15 | 38.97367935 | 0.987214831 | 0.012785169 |
| 2 | 31 | 39.35174573 | 0.996791365 | 0.003208635 |
| 2 | 63 | 39.44671910 | 0.999197068 | 0.000802932 |
| 2 | 127 | 39.47049107 | 0.999799219 | 0.000200781 |
| 3 | 7 | 79.01652066 | 0.889560822 | 0.110439178 |
| 3 | 15 | 86.28755850 | 0.971417507 | 0.028582493 |
| 3 | 31 | 88.18619242 | 0.992792155 | 0.007207845 |
| 3 | 63 | 88.66603037 | 0.998194127 | 0.001805873 |
| 3 | 127 | 88.78631555 | 0.999548287 | 0.000451713 |

The continuum rates are 9.869604401, 39.478417604, and 88.826439610 for
modes 1, 2, and 3, respectively. Every ratio is below one and each mode's
error falls at all four refinements. The fitted orders on the finest four
grids were 1.99942, 1.99767, and 1.99476. All three modes meet the predeclared
criteria, so the second-order convergence prediction is supported for this
specified path model and tested grid.

Numerical integration was assessed separately. Across the 60 runs, the
largest fitted-rate relative error against the exact discrete eigenvalue was
`4.14e-11`; the largest absolute state discrepancy against the exact discrete
trajectory was `8.15e-6` across all runs and `1.59e-7` for the tighter
`rtol=1e-12` runs. The maximum fitted-rate difference between the coarsest
solver/output setting and the tightest/densest setting for a matched grid and
mode was `4.13e-11` relative. These are observed finite-computation errors,
not the spatial discretization error. In particular, the absolute state
discrepancy is larger than the rate-fit discrepancy because small numerical
components outside the ideal sine mode can affect individual state entries
without materially changing the projected modal slope. The exact formula,
not a fitted numerical trajectory, determines the convergence-order test.

The trajectory plot is illustrative for selected cases. The complete numeric
record is in the CSV files. No finite-grid observation is presented as proof
about a different network or an applied system.

## Interpretation and limitations

**Direct mathematical result:** the discrete eigenvalue formula gives a rate
ratio below one for these modes, with errors decreasing quadratically on the
specified nested grids. The measured fitted order also lies close to two.

**Numerical observation:** DOP853 solutions project onto the expected mode
with fitted rates close to the exact discrete rates under both declared
tolerance and output settings. The integration is not exact, as shown by the
state discrepancies. Output samples are observations of an internal
continuous-time solve; increasing their count does not refine the spatial
grid.

**Scope:** the correspondence follows from deliberately using the same
nearest-neighbor second-difference stencil and fixed-end boundary conditions.
It does not imply that general consensus protocols are heat equations, that
arbitrary graph topologies share these rates, or that the result validates a
real-world or financial application. The network time `τ` is dimensionless
algorithmic relaxation time; continuum heat time is normalized by `α/L²`.
No market data, forecasts, trading rules, or profitability claims are involved.
