# Phase 49 — Critical Slowing Down and Cross-System Stability Signatures

## Research question

As dissipation approaches zero from the stable side, does local recovery slow
in both fixed-end heat diffusion and a damped nonlinear pendulum? If so, do
their measured rate laws share an exponent, and how sensitive are those
measurements to perturbation amplitude and numerical settings?

This is a bounded comparison of two existing, structurally distinct Newton
Lab models. It does not claim a new law, a universal stability signature, or
physical equivalence. Phase 32 remains `BLOCKED_AUTHORIZATION`; no market data
or external service was used.

## Model selection and theory

### Fixed-end heat-diffusion mode

The Phase 19 model is the one-dimensional transient heat equation

```text
u_t = alpha*u_xx,       u(0,t)=u(L,t)=0
```

for a perturbation `u=T-T_base` about a fixed-end steady profile. For a single
sine mode, `u_n(x,t)=A sin(n*pi*x/L) exp(-lambda_n*t)`, with
`lambda_n=alpha*(n*pi/L)^2`. The selected equilibrium is the uniform `300 K`
profile on `L=1 m`, with fixed endpoints at `300 K`; only mode `n=1` is
perturbed. For `alpha>0`, the perturbation decays and the equilibrium attracts
compatible states. At `alpha=0`, every compatible profile is stationary and
there is no recovery to a unique steady profile.

The dimensionless boundary distance is
`delta_heat=alpha/(L^2*omega_ref)`, where `omega_ref=sqrt(g/L)=1 s^-1`
is the pendulum's reference frequency for the chosen `g=L=1` values. Thus
`lambda_1/omega_ref=pi^2*delta_heat`: the predicted power exponent is 1.

### Nonlinear damped pendulum

The Phase 3 full-sine model is

```text
theta_ddot + b/(m*L^2)*theta_dot + (g/L)*sin(theta) = 0
```

with the downward rest equilibrium `(theta,theta_dot)=(0,0)`. We set
`m=1 kg`, `L=1 m`, `g=1 m/s^2`, and vary the viscous rotational damping
coefficient `b`. For `b>0`, the equilibrium is locally asymptotically stable.
At `b=0`, it is a conservative center: nearby states remain bounded but do
not recover asymptotically.

Linearization gives eigenvalues
`r=-gamma +/- i*omega_d`, where `gamma=b/(2*m*L^2)` and
`omega_d=sqrt(g/L-gamma^2)` in the underdamped regime. The signed eigenvalues
have negative real part. The positive envelope recovery rate is
`gamma=-Re(r)`. The sweep remains underdamped. Its dimensionless distance is
the damping ratio `delta_pendulum=gamma/omega_0`, with
`omega_0=sqrt(g/L)`, so `gamma/omega_0=delta_pendulum` and the predicted
exponent is 1. Both eigenvalues have the same real part; the oscillatory
frequency is a separate mode property and is retained in the records.

These models have different states and mechanisms: a spatial continuum field
with real modal attenuation versus a nonlinear second-order mechanical system
with an oscillatory eigenpair. They share a zero-dissipation stability
boundary. The dimensionless distances are each model's natural normalized
dissipation parameter; their numerical values do not imply identical
physical control parameters. The heat prefactor `pi^2` and pendulum prefactor
1 also differ under this declared normalization.

## Frozen experiment design and measurement

The five predeclared distances are `0.025, 0.05, 0.1, 0.2, 0.4`, fitted
together by ordinary least squares of `log(measured rate)` against
`log(distance)`. The prediction is supported for a model only if its fitted
exponent lies in `[0.95, 1.05]`, the exponent-fit `R^2` is at least `0.999`,
all five primary measured rates are within 2% of their linearized predictions,
the measured rate strictly increases at each adjacent distance, and amplitude
sensitivity meets its separate criterion below.
The two initial perturbations are `0.01` and `0.3` in each model's own units
(kelvin for the heat mode, radians for the pendulum). Within-model amplitude
sensitivity passes if the largest relative measured-rate difference is at
most 2%. These criteria, distances, windows, tolerances, and fitting range
are recorded in `reports/phase_49_critical_slowing_down/metadata.json`.

For heat, the existing exact finite-mode evaluator supplies modal amplitudes
at 401 evenly spaced times over 12 predicted e-folds. A log-amplitude fit
measures the rate. This checks finite-window rate recovery from sampled
outputs; it is not an independent numerical PDE integration, so there is no
heat solver error estimate in this experiment.

For the pendulum, the existing full-sine solver uses DOP853 at both
`rtol=1e-10, atol=1e-12` and `rtol=1e-12, atol=1e-14`. The output cadence is at
least 100 samples per predicted damped period, and the window spans 12
linearized e-folds. Positive angular peaks above `1e-6` of the initial angle
are fit against their occurrence times. Small initial angle `0.01 rad` is
the primary local measurement; `0.3 rad` checks amplitude sensitivity. The
small-angle eigenvalue prediction is calculated independently of the
nonlinear trajectory integration.

The study is deterministic and uses no random seeds. It writes 30 conditions
to `per_condition_results.csv`, two model rows to `model_summary.csv`, three
figures, and environment/design metadata. Solver tolerance differences,
rate-fit residuals, and deviations from linearized rates are recorded
separately. The study command is:

```powershell
.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.critical_slowing_down import run_critical_slowing_down_study; run_critical_slowing_down_study(Path('reports/phase_49_critical_slowing_down'))"
```

## Results

### Primary local measurements

| Distance | Heat predicted (s⁻¹) | Heat measured (s⁻¹) | Pendulum predicted (s⁻¹) | Pendulum measured (s⁻¹) | Pendulum relative error |
|---:|---:|---:|---:|---:|---:|
| 0.025 | 0.2467401100 | 0.2467401100 | 0.025000 | 0.0250000522 | 0.0002087% |
| 0.05 | 0.4934802201 | 0.4934802201 | 0.050000 | 0.0499998647 | 0.0002705% |
| 0.1 | 0.9869604401 | 0.9869604401 | 0.100000 | 0.0999980476 | 0.0019524% |
| 0.2 | 1.9739208802 | 1.9739208802 | 0.200000 | 0.2000102519 | 0.0051260% |
| 0.4 | 3.9478417604 | 3.9478417604 | 0.400000 | 0.4000304189 | 0.0076047% |

The heat sampled-fit rates equal the analytic rates at displayed precision.
This is expected from the exact evaluator and should not be interpreted as an
independent PDE-solver validation. The primary pendulum rate errors were all
below `7.61e-5` relative, substantially inside the predeclared 2% criterion.
Primary peak fits used 5 to 77 peaks as damping increased/decreased across
the sweep; minimum reported fit `R^2` was approximately `0.99999994`.

| Model | Fitted exponent | Log-log `R^2` | Largest primary rate error | Largest amplitude sensitivity |
|---|---:|---:|---:|---:|
| Fixed-end heat, mode 1 | 1.000000 | 1.000000000 | 0 | 0 |
| Nonlinear damped pendulum | 1.000029 | 0.9999999995 | `7.60e-5` relative | `9.36e-5` relative |

Both models show slower local recovery as their declared zero-dissipation
boundary is approached. Over the complete predeclared range, both fitted
exponents satisfy the exponent criterion. The pendulum's maximum relative
rate difference between solver tolerance settings was approximately
`6.45e-7`; the largest absolute difference was `3.57e-8 s^-1`. These are
solver-setting sensitivity checks, not formal global error bounds.

The `0.3` amplitude pendulum fit differs from its `0.01 rad` fit by at most
`9.36e-5` relative in this design. This shows low amplitude sensitivity for
these particular initial states and fitting window. It does not imply that
finite-amplitude nonlinear pendulum recovery is generally exponential or
independent of amplitude.

The strictest-versus-looser solver-setting difference was at most
`6.45e-7` relative (`1.61e-8 s^-1` absolute) across matched pendulum
conditions. The minimum log-peak-fit `R^2` over the pendulum conditions was
`0.9999999367`. Both checks pass their recorded numerical-resolution criteria;
they do not bound all solver error.

## Cross-system interpretation

**Established model theory:** For the selected heat eigenmode,
`lambda_1` is linear in diffusivity. For the underdamped pendulum near its
downward equilibrium, the envelope rate is linear in damping. The two
different characteristic spectra have different observables: monotone modal
decay versus an oscillatory trajectory with a decaying envelope.

**Numerical validation in this experiment:** Sampled heat modal amplitudes
recover their known rate; numerical full-sine pendulum peaks recover the
linearized envelope rate with small error over the frozen sweep. Both show
the predicted exponent-one trend, and the pendulum results remain stable
under the tighter solver tolerance and the tested perturbation-amplitude
change.

**Cross-system hypothesis:** The limited hypothesis of a shared qualitative
slowing trend and exponent one is supported over the five tested distances.
The matching exponent follows from the selected linear dependence on each
system's dissipation parameter. It does not establish a shared mechanism,
common prefactor, common state space, or universal critical law. No nonlinear
bifurcation with a vanishing oscillation frequency was tested: the pendulum
frequency remains nonzero at zero damping, while heat loses diffusion and
becomes non-attracting.

## Limitations and broader research relevance

- Only one heat mode, one pendulum equilibrium, fixed physical parameters,
  one-dimensional homogeneous heat diffusion, and a five-point parameter
  range were examined.
- The heat output is analytical; its exact sampled fit cannot estimate
  independent PDE integration error. A future numerical PDE solver would be
  needed for that numerical cross-validation.
- The pendulum window and peak method become less well conditioned as damping
  tends to zero because the required observation duration grows. The sweep
  stops at dimensionless distance `0.025`; it does not approximate the exact
  boundary as a stable case.
- The `0.3 rad` perturbation is a within-model sensitivity check only. The
  amplitude values in kelvin and radians are not comparable across models.
- Results are from ideal equations and synthetic states, with no measured
  physical systems, external disturbances, or financial data.

This result gives Newton Lab a reusable way to ask whether a stability
signature survives across distinct model structures: define the boundary,
derive each model's own recovery law, then compare measured exponents and
measurement limits without assuming the mechanisms are equivalent. It
provides a bounded research pattern, not a claim of universal discovery.

## Artifacts

- [`per_condition_results.csv`](../../reports/phase_49_critical_slowing_down/per_condition_results.csv)
- [`model_summary.csv`](../../reports/phase_49_critical_slowing_down/model_summary.csv)
- [`metadata.json`](../../reports/phase_49_critical_slowing_down/metadata.json)
- [`recovery_rate_vs_boundary_distance.png`](../../reports/phase_49_critical_slowing_down/recovery_rate_vs_boundary_distance.png)
- [`predicted_vs_measured_rates.png`](../../reports/phase_49_critical_slowing_down/predicted_vs_measured_rates.png)
- [`recovery_rate_scaling.png`](../../reports/phase_49_critical_slowing_down/recovery_rate_scaling.png)
