# Phase 50 — Cross-Domain Hypothesis Selection and Validation

## Executive summary

Phase 49 found that the local recovery rate decreased linearly with distance
to zero dissipation in a fixed-end heat mode and an underdamped nonlinear
pendulum near its stable rest state. Phase 50 tests whether that general
recovery-rate viewpoint transfers usefully to a distinct engineering-control
stability problem with a different stability transition.

The target is the canonical scalar saddle-node normal form
`dx/dt = mu - x^2`, interpreted as a normalized control plant approaching
loss of a stable operating point. Its stable equilibrium has eigenvalue
`-2*sqrt(mu)`, so the positive local recovery rate scales as `mu^(1/2)`.
With 100 independent observation-noise replicates per condition, the
fold-specific rate model was evaluated on five held-out margins against an
exponent-one baseline transferred from Phase 49. The fold model's held-out
mean absolute log error was `0.00401`; the baseline's was `0.91840`. The
paired bootstrap 95% interval for their score difference (fold minus baseline)
was `[-0.91493, -0.91386]`. It supports this protocol's falsifiable prediction
for this synthetic benchmark.

This validates a model-specific rate law in a canonical simulated control
system. It does not validate any particular controller, converter, grid,
industrial process, or field measurement. The fold scaling law is established
normal-form theory; Newton Lab's result is a finite numerical and noisy
observation benchmark, not a new discovery or proof of broad application.

## Research question and prior Newton Lab evidence

**Question:** Can a physically motivated method that tracks local recovery
rate near a stability boundary help predict recovery in a distinct target
system, and does the exponent observed in Phase 49 transfer to that system?

The repository review found:

- **Phases 1–17:** established the Python package, physical solver families,
  equation registry, provenance-aware experiment records, BVP and algebraic
  solvers, discovery workflow, analysis, and reporting foundations. README,
  source APIs, tests, and architecture notes were inspected; Phase 48 and 49
  outputs were also checked directly.
- **Phases 18–19 and 22–28:** compared oscillator, adjustment, pendulum,
  steady heat, transient heat modes, and modal estimation. Phase 18's
  oscillator-to-adjustment mapping is an exact toy equation mapping, not a
  validated financial system. Phases 23–24 show that finite windows, cadence,
  noise, and close modes limit recovery-rate estimation. Phases 25–26 show
  finite-angle pendulum behavior departs from linearization in
  observable-dependent ways; Phase 28 separates solver error from model
  discrepancy.
- **Phases 29–30:** tested diffusion-inspired and causal smoothing on synthetic
  signals, documenting lag, attenuation, and noncausal/causal differences.
  These are not real-market validation and were not repeated here.
- **Phases 31–41:** specify chronological research safeguards, audit NSE
  sources and supplied single-day files, and record unresolved data terms.
  The data are not a temporal response record and Phase 32 remains
  `BLOCKED_AUTHORIZATION`, so finance was not selected.
- **Phases 42–46:** used synthetic signal diagnostics to demonstrate
  ambiguity, false alarms, and noise sensitivity. Those results argue against
  treating a diagnostic statistic as identification of mechanism.
- **Phases 47–48:** proposed and then validated the narrow fixed-end heat/path
  Laplacian connection. Phase 48 is complete prior work and remains unchanged.
- **Phase 49:** used sampled exact heat modes and numerical nonlinear-pendulum
  peaks. Each model's recovery rate was proportional to its own normalized
  dissipation parameter over the frozen five-point grid (fitted exponents
  `1.000000` and `1.000029`). It explicitly did not establish a universal
  critical exponent.

These findings motivate testing a different stability boundary. A fold
transition can have a square-root rate law even though the Phase 49 examples
had exponent one.

## Candidate selection

| Candidate | Assessment | Decision |
|---|---|---|
| Phase 49 recovery-rate signature applied to a control system approaching a saddle-node loss of operating point | Distinct nonlinear stability mechanism; explicit stable branch and Jacobian eigenvalue; matched rate estimator and a useful exponent-one transfer baseline; can be tested with a canonical, falsifiable simulation. | **Selected**, with the scope limited to a canonical benchmark rather than a named device. |
| Diffusion-inspired smoothing applied to signals | Phase 29–30 already tested smoothing, causal filtering, lag, and synthetic reconstruction; another such run would repeat that line rather than test a new target mechanism. | Not selected. |
| Oscillator damping mapped to finance risk adjustment | Phase 18 already tests an exact toy coefficient mapping, while the NSE data authorization gate remains unresolved and the supplied report file is single-day. | Not selected. |
| Pinned path network relaxation | Phase 48 already studies its exact Laplacian/heat connection; extending it would remain within that narrow structural comparison. | Not selected. |

No external services or data were accessed. Raw files in `data/raw/` are not
used: one is a single-day cross-section, and their provenance and permitted use
remain unresolved. There is no suitable authorized measurement record of an
engineering plant, so the strongest available defensible test is a controlled
simulation of an explicitly specified target mechanism.

## Source-to-target mapping and hypothesis

### Source result

Phase 49 measures a positive recovery rate from local perturbation decay. For
the tested heat mode, `lambda = alpha*(pi/L)^2`; for the underdamped pendulum,
the two signed eigenvalues are `-gamma +/- i*omega_d`, and the positive
envelope rate is `gamma = -Re(lambda)`. In both cases the rate tends to zero
as the selected dissipation distance tends to zero, with an exponent near one
under those particular parameterizations.

### Target normal form

We use the established saddle-node normal form

```text
dx/dt = mu - x^2
```

where `mu` is a normalized operating margin, `x` is a normalized scalar
control-state deviation, and time is normalized. This abstraction represents
the local behavior of a one-dimensional engineering operating point near a
fold loss of equilibrium. It is not a calibrated model of a specific device.
For `mu > 0`, the equilibria are `x_plus=+sqrt(mu)` (stable) and
`x_minus=-sqrt(mu)` (unstable). Linearizing at the stable state gives the
signed eigenvalue

```text
lambda_stable = d(mu-x^2)/dx at x_plus = -2*sqrt(mu).
```

The positive recovery rate is `rho=-lambda_stable=2*sqrt(mu)`. Thus this
target predicts exponent `1/2`, not Phase 49's exponent one. At `mu=0` the
equilibria coalesce and the recovery rate vanishes; the boundary is not
simulated as an ordinary stable condition. For `mu<0` no equilibrium exists.

### Explicit mapping

- **Source structure:** estimate local positive recovery rate from the decay
  of a small deviation after perturbing a stable equilibrium.
- **Mathematical definition:** `rho=-Re(lambda)` for a simple stable mode, or
  `-Re(lambda)` for each member of the pendulum's conjugate pair (both give
  the same positive envelope rate).
- **Target quantity:** decay rate of `x(t)-sqrt(mu)` after a positive pulse on
  the stable branch of a normalized control plant.
- **Assumptions:** scalar smooth fold normal form; known positive margin;
  constant margin during recovery; local perturbation; no process noise;
  independent additive observation noise; accurately known equilibrium for
  computing the observed deviation.
- **Potential value:** choose a rate-versus-margin model that predicts
  lower-margin recovery from higher-margin calibration and distinguishes a
  fold's square-root approach from an exponent-one extrapolation.
- **Evidence against:** fold-specific prediction fails the held-out score
  criterion, its paired uncertainty interval includes no improvement, the
  measured local rates miss the Jacobian prediction, or numerical integration
  error is too large to resolve the rates.

The candidate is not a claim that heat diffusion causes controller failure.
The transferable idea is a stability-aware recovery observable; the target's
exponent comes from the target model's own Jacobian.

## Frozen experimental protocol

The protocol was written to
[`reports/phase_50_cross_domain_validation/protocol.json`](../../reports/phase_50_cross_domain_validation/protocol.json)
before the final evaluation. Its SHA-256 is recorded in `metadata.json`. It
fixes all margins, seeds, tolerances, perturbations, estimator, metrics, and
acceptance criteria.

- **Calibration margins:** `mu = 0.25, 0.36, 0.49, 0.64, 0.81, 1.0`.
- **Held-out margins:** `mu = 0.04, 0.0625, 0.09, 0.1225, 0.16`; all lie below
  the smallest calibration margin. Margin values, not time samples, define
  the train/test split.
- **Perturbations:** initial state `x(0)=sqrt(mu)*(1+rho)` for `rho=0.01`
  (primary) and `0.05` (sensitivity check).
- **Repeated measurements:** 100 independent, seeded Gaussian observation
  noise realizations per margin, `rho`, with fixed standard deviation
  `2e-5` normalized state units and no process noise. The same noise vector is
  used for both perturbation sizes within each margin/replicate. Replicates
  are independent measurement-noise draws around the same deterministic
  trajectory; they do not represent independent physical systems or process
  realizations.
- **Numerical integration:** DOP853, `rtol=1e-11`, `atol=1e-13`, 101 requested
  samples on `[0, log(5)/(2*sqrt(mu))]`. This spans a five-fold ideal linear
  decay while retaining signal above the declared noise floor.
- **Rate measurement:** ordinary least-squares slope of log observed
  deviation `x_observed-sqrt(mu)` against time; the positive rate is the
  negative signed slope. Fit `R^2` and rate errors are retained per replicate.
- **Proposed method:** fix exponent `p=1/2`, derive it from the fold Jacobian,
  and fit only prefactor `C` on the six calibration margins using mean log
  measured rates.
- **Baseline:** fix exponent `p=1`, as observed in Phase 49, and fit only its
  prefactor on the identical calibration data.
- **Fairness:** both predictions use the same primary-perturbation records,
  six training margins, five held-out margins, estimator, and test replicates;
  only the predeclared exponent changes. No test result tunes either model.
- **Primary score:** mean absolute log error of predicted versus measured
  held-out rates, averaged equally over margins and observation-noise
  replicates. Also report RMSE of relative rate error.
- **Uncertainty:** 10,000 paired bootstrap resamples using fixed seed `64050`.
  Each replicate index is resampled as a block across calibration and held-out
  margins; calibration prefactors are refit in every resample. This interval
  describes observation-noise uncertainty conditional on the one simulated
  normal-form mechanism.
- **Acceptance:** proposed/baseline held-out score ratio at most `0.8`; upper
  endpoint of the paired 95% interval for score difference below zero; maximum
  relative ODE trajectory error below `1e-8`; and primary and larger
  perturbation mean rates each within 2% of the local linearized rate at every
  margin. These thresholds were fixed before final artifacts were inspected.

The design does not include a validation partition because the exponents and
all analysis choices are fixed theoretically before evaluation; only a
prefactor is calibrated on the training margins.

Reproduce with the project environment from the repository root:

```powershell
.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.fold_control_validation import run_fold_control_validation; run_fold_control_validation(Path('reports/phase_50_cross_domain_validation'))"
```

## Results

### Recovery-rate measurements on held-out margins

The table gives the primary perturbation's mean measured rate across 100
noise replicates, its replicate SD, and a normal-approximation 95% confidence
interval for the mean. The local theoretical rate is shown independently.

| Held-out `mu` | Predicted `2*sqrt(mu)` | Measured mean | Replicate SD | 95% CI for mean | Relative mean error |
|---:|---:|---:|---:|---:|---:|
| 0.04 | 0.400000 | 0.400718 | 0.002606 | [0.400207, 0.401228] | 0.1794% |
| 0.0625 | 0.500000 | 0.501307 | 0.002908 | [0.500737, 0.501877] | 0.2614% |
| 0.09 | 0.600000 | 0.601343 | 0.002980 | [0.600759, 0.601928] | 0.2239% |
| 0.1225 | 0.700000 | 0.702337 | 0.003050 | [0.701739, 0.702935] | 0.3338% |
| 0.16 | 0.800000 | 0.802264 | 0.002545 | [0.801765, 0.802763] | 0.2830% |

Across all primary margins, the largest relative error of the mean measured
rate from the local Jacobian rate was `0.334%`. For the larger `rho=0.05`
perturbation, the maximum was `1.196%`; both remain below the frozen 2%
screen. Maximum relative discrepancy of the numerical ODE trajectory from
the independent closed-form fold trajectory was `9.78e-12`, below its
`1e-8` criterion. Each displayed primary 95% interval for the measured mean
lies slightly above the exact Jacobian rate. The interval measures observation
noise variation around the finite-pulse estimator; it does not account for
systematic estimator bias. The consistent positive offset is expected from
fitting a straight log-amplitude slope to a finite perturbation of the
nonlinear fold trajectory, even within the local 2% acceptance threshold.

### Held-out prediction comparison

| Method | Fixed exponent | Calibration prefactor | Held-out mean absolute log error (95% bootstrap interval) | Held-out relative-rate RMSE |
|---|---:|---:|---:|---:|
| Fold-theory method | 0.5 | 2.00490 | 0.004010 (`[0.003773, 0.004261]`) | 0.005086 |
| Phase 49 exponent-one baseline | 1.0 | 2.74684 | 0.918405 (`[0.917927, 0.918892]`) | 0.596945 |

The proposed-to-baseline mean absolute log-error ratio was `0.00437`, a
`99.56%` reduction under this metric and split. The paired bootstrap interval
for the fold-minus-baseline score difference was
`[-0.91493, -0.91386]`, entirely below zero. All frozen criteria passed.
The difference is large because every test margin is below the training
margins and the two fixed exponents diverge during that extrapolation; it is
not a claim of comparable improvement on other systems or ranges.

## Interpretation

**Established theory:** The saddle-node normal form has a stable branch
`sqrt(mu)` and local signed eigenvalue `-2*sqrt(mu)`. Phase 49's selected
models have their own exponent-one relationships under their own parameter
definitions. Those source and target results are not conflicting: the local
eigenvalue depends on the model structure and on how distance to the boundary
enters that structure.

**Implementation check:** DOP853 trajectories agree with the independent
closed-form nonlinear transient to a maximum relative discrepancy below
`1e-11`. The log-rate measurement recovers the Jacobian rate within the frozen
2% mean-error criterion at both pulse sizes. The measured fold system is
nonlinear away from equilibrium; pulse amplitude introduces a small
finite-perturbation bias, more visible for the larger pulse.

**Experiment result:** On this simulated control benchmark, Phase 49's
qualitative recovery-rate observable transfers, but its exponent-one scaling
does not. The fold-derived exponent-one-half model added measurable held-out
predictive value over the exponent-one baseline on the prespecified low-margin
range. This is consistent with known saddle-node theory. It is a validation
of the implementation and protocol, not a discovery of the fold law.

The uncertainty interval uses independent observation-noise replicates and
does not represent variability across plants, process noise, parameter
uncertainty, or model misspecification. The small benchmark and specified
noise floor make the method comparison controlled, not empirically validated.

## Limitations, threats, and unresolved question

- The target is a canonical normal form, not a named engineering device with
  measured parameters. Real control plants may have multiple states, delays,
  saturation, stochastic forcing, changing inputs, and non-fold stability
  boundaries.
- `mu` and equilibrium `sqrt(mu)` are treated as known. In practice, margin
  estimation and equilibrium drift may dominate recovery-rate error.
- There is no process noise. The 100 replications vary observation noise only
  and share one deterministic underlying trajectory per condition.
- The simulated measurement uses the exact equilibrium to form deviations.
  Estimating that equilibrium from noisy/finite data is not tested.
- The fitting window scales with the theoretical rate. A blind early-warning
  application would need to estimate a useful window without knowing the
  target rate.
- Results cover only five held-out margins, one scalar fold, two relative
  pulse sizes, and one observation-noise level. The reported bootstrap is
  conditional on those choices.
- The Phase 49 exponent-one baseline is a deliberate transfer baseline, not
  the only plausible engineering model. A fitted exponent using the same
  held-out data would be leakage and was not used.

**Most important unresolved question:** Does the same rate-estimation and
model-selection procedure retain value when the operating margin and
equilibrium must be inferred from noisy, partially observed data from a
specified engineering plant? No such authorized and sufficiently documented
plant dataset is currently in the workspace.

## What this enables

Newton Lab can now test a useful cross-domain workflow: carry a measurable
stability signature into a distinct target mechanism, derive that target's
eigenvalue law independently, and compare its out-of-sample predictions
against a source-derived baseline. Here that workflow exposed why qualitative
critical slowing can transfer while the exponent does not. It provides a
benchmark for later empirical validation, but it does not itself validate a
real controller or bypass Phase 32 authorization.

## Generated artifacts

- [`protocol.json`](../../reports/phase_50_cross_domain_validation/protocol.json)
- [`per_condition_results.csv`](../../reports/phase_50_cross_domain_validation/per_condition_results.csv)
- [`condition_summary.csv`](../../reports/phase_50_cross_domain_validation/condition_summary.csv)
- [`predictor_comparison.csv`](../../reports/phase_50_cross_domain_validation/predictor_comparison.csv)
- [`paired_comparison.csv`](../../reports/phase_50_cross_domain_validation/paired_comparison.csv)
- [`metadata.json`](../../reports/phase_50_cross_domain_validation/metadata.json)
- [`held_out_recovery_rate_predictions.png`](../../reports/phase_50_cross_domain_validation/held_out_recovery_rate_predictions.png)
- [`held_out_prediction_error.png`](../../reports/phase_50_cross_domain_validation/held_out_prediction_error.png)
