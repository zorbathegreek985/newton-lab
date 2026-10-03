# Phase 22 — Modal dynamics, decay, and structural similarity

## Research question

Which mathematical properties of modal evolution are shared by Newton Lab's
Phase 18 oscillator and Phase 19 heat-diffusion examples, and which depend on
their different physical models?

The controlled calculations are reproducible with
`newton_lab.modal_research.run_modal_dynamics_research()`; the computed tables
below are rendered by `render_modal_dynamics_report`. The heat values come from
the existing Phase 19 case-study API. Oscillator traces use the existing SI
oscillator integration API. No market data or external data are used.

## Hypotheses

| Hypothesis | Result within this study |
|---|---|
| H1: Under the specified linear, time-invariant assumptions, individual modes evolve exponentially. | Supported. The heat finite-mode attenuation matches its analytic exponential; the underdamped oscillator numerical displacement matches the closed form and its sampled absolute peaks fit the analytic envelope rate under dense sampling. |
| H2: Modal decay rate depends on physical parameters as specified. | Supported for the tested values. Heat rate is proportional to diffusivity and mode number squared, and inversely proportional to length squared. Oscillator envelope rate is `c/(2m)`; damped frequency is distinct and changes with `m`, `k`, and damping ratio. |
| H3: Finite windows and sampling can mislead rate estimates. | Demonstrated. A sparse oscillator grid biases the peak-rate fit, a short window has too few peaks, and the fastest heat mode falls below a declared analysis floor at the final time. |
| H4: Structural similarity does not establish application validity. | Retained as a methodological principle. This phase supplies no empirical finance evidence. |

These are deterministic checks of the stated models and chosen parameter grid,
not universal theorems or statistical estimates.

## Mathematical models and assumptions

### Heat diffusion

The Phase 19 transient example uses

\[
\frac{\partial T}{\partial t}=\alpha\frac{\partial^2T}{\partial x^2},
\qquad 0<x<L.
\]

It assumes a homogeneous one-dimensional medium, constant positive
diffusivity, linear diffusion, fixed endpoint temperatures, and no internal
source in the transient perturbation. Subtracting the linear steady profile
`T_base(x)` makes the perturbation `u=T-T_base` satisfy homogeneous Dirichlet
conditions `u(0,t)=u(L,t)=0`. The sine eigenmodes are
`sin(n*pi*x/L)`, with

\[
\lambda_n=\alpha\left(\frac{n\pi}{L}\right)^2,\qquad
a_n(t)=a_n(0)e^{-\lambda_n t}.
\]

Here `alpha` is in m²/s, `L` in m, `lambda_n` in s⁻¹, and each temperature
modal coefficient is in K. These expressions do not apply unchanged to
arbitrary boundary conditions, spatially varying material properties, source
terms, or nonlinear diffusion. The Phase 19 steady BVP is separate from this
transient exact finite-mode evaluation.

### Damped oscillator

The Phase 18 source experiment uses the linear free oscillator

\[
m\ddot{x}+c\dot{x}+kx=0,
\quad \omega_n=\sqrt{k/m},
\quad \zeta=\frac{c}{2\sqrt{mk}}.
\]

The source cases use `m=1 kg`, `k=4 N/m`, initial displacement `0.1 m`, zero
initial velocity, and damping sweep `c=(0, 0.1, 0.2) kg/s`. The implementation
assumptions are a linear spring, linear viscous damping, constant parameters,
and no external driving force. The Phase 22 underdamped checks use the same
`k`, initial conditions, `c=(0.05, 0.1, 0.2) kg/s`, and `m=1 kg`; mass checks
use `m=(0.5, 1, 2) kg` with `c=0.1 kg/s`. All are underdamped.

For `0 <= zeta < 1`, the characteristic roots are

\[
r=-\gamma\pm i\omega_d,
\quad \gamma=\frac{c}{2m},
\quad \omega_d=\omega_n\sqrt{1-\zeta^2}.
\]

`gamma` is the displacement-envelope decay rate in s⁻¹. `omega_d` is the
oscillation frequency in rad/s. The displacement is an oscillatory linear
combination multiplied by `exp(-gamma*t)`. Its signed samples are not positive
modal amplitudes and cross zero. Initial displacement and velocity set the
modal coefficients and phase, not the characteristic roots; a mode can be
absent if its coefficient is zero. At zero damping the envelope does not
decay. Critical and overdamped cases have real roots and are outside the
oscillatory formulas used for the Phase 22 oscillator fits. The nonlinear
pendulum API is not treated as this linear oscillator.
For `c>0`, oscillator transients approach zero at long times; for `c=0` the
undamped oscillation persists. Heat perturbation modes approach zero while the
fixed linear steady profile remains. A finite run can support neither limit by
itself; these long-time statements follow from the specified equations.

## Experimental design and numerical method

The heat study reuses Phase 19's exact finite sine-mode solution and its Phase
19 values at `t=(0, 60, 300, 1200) s`; spatial values are sampled on its finite
position grid. Modal attenuation factors are calculated directly from the
closed form, not estimated by the steady BVP solver.

The oscillator study integrates each specified parameter set with Newton Lab's
existing SciPy `solve_ivp` adapter using RK45, 20,001 uniform output points over
10 s, `rtol=1e-10`, and `atol=1e-12`. Sampled absolute local peaks are fit by a linear
regression of `log(peak magnitude)` against time. Signed displacements and zero
crossings are never passed to the logarithmic fit. The analytical displacement
reference is evaluated at the same output times; its maximum absolute
discrepancy is reported below.

For heat rate estimation, only positive modal magnitudes are used. The
fast-mode threshold `1e-8 K` is an explicitly chosen analysis floor, not
instrument noise or a measured resolution. No confidence intervals or
statistical significance are inferred.

## Numerical results

### Heat modal predictions and attenuation

The default Phase 19 settings are `L=1 m`, `alpha=1e-4 m²/s`, with initial
modal amplitudes `a1(0)=1 K` and `a4(0)=0.5 K`.

| Mode | Initial amplitude (K) | Predicted `lambda_n` (s⁻¹) | Attenuation at t=(0, 60, 300, 1200) s |
|---:|---:|---:|---|
| n=1 | 1.0 | 0.000986960440 | 1, 0.942502, 0.743722, 0.305944 |
| n=4 | 0.5 | 0.0157913670 | 1, 0.387716, 0.00876131, 5.89218e-9 |

The positive-magnitude log fits return `0.000986960440 s⁻¹` for n=1. For n=4,
using only values at or above `1e-8 K`, the fitted rate is
`0.0157913670 s⁻¹`. The n=4 analytical amplitude at 1200 s is approximately
`2.94609e-9 K`, below the selected floor. With no finite floor and exact
arithmetic, the formula continues to decay; the result shows why a finite
observation cannot establish a long-time rate once a mode is beneath resolution.
In the full heat solution, the transient perturbation vanishes at long time and
the linear steady profile remains as `T_base(x)`.

### Heat parameter dependence

The following one-factor calculations hold other inputs fixed and evaluate
`lambda_n=alpha*(n*pi/L)^2`:

| Varied quantity | Values | Rates (s⁻¹) |
|---|---|---|
| `alpha` (m²/s), n=1, L=1 m | 5e-5, 1e-4, 2e-4 | 0.000493480, 0.000986960, 0.001973921 |
| n, alpha=1e-4 m²/s, L=1 m | 1, 2, 4 | 0.000986960, 0.003947842, 0.015791367 |
| `L` (m), alpha=1e-4 m²/s, n=1 | 0.5, 1, 2 | 0.003947842, 0.000986960, 0.000246740 |

Thus doubling diffusivity doubles the rate, doubling mode number multiplies
the rate by four, and doubling length divides it by four.

### Oscillator rate and closed-form comparison

| m (kg) | c (kg/s) | zeta | `gamma` (s⁻¹) | `omega_d` (rad/s) | Dense sampled peak fit (s⁻¹) | Max displacement error (m) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.05 | 0.0125 | 0.025 | 1.99984374 | 0.0250004175 | 2.13e-11 |
| 1 | 0.1 | 0.025 | 0.05 | 1.99937490 | 0.0500011950 | 1.73e-11 |
| 1 | 0.2 | 0.05 | 0.1 | 1.99749844 | 0.1000008580 | 1.16e-11 |
| 0.5 | 0.1 | 0.0353553 | 0.1 | 2.82665881 | 0.100014321 | 1.60e-11 |
| 2 | 0.1 | 0.0176777 | 0.025 | 1.41399257 | 0.0250004358 | 1.48e-11 |

The dense fits agree with `c/(2m)` for these inputs. At fixed mass and
stiffness, damping changes the envelope rate while `omega_d` changes only
slightly over this low-damping grid. At fixed `c=0.1 kg/s`, changing mass from
0.5 to 2 kg reduces `gamma` from 0.1 to 0.025 s⁻¹; the frequency also changes.
The maximum numerical-to-analytical displacement discrepancies are small for
this solver configuration and this interval. They assess implementation
agreement with the chosen equation, not whether the physical model describes
any particular real oscillator.

The separate Phase 19 steady BVP was numerically accepted and its maximum
profile error against the linear reference was 0 K on its reported mesh. This
supports only that steady numerical comparison; it is not evidence about the
transient integration, whose modal values above come from the separate exact
finite-mode formula.

## Observation-window and combined-signal checks

| Baseline oscillator observation | Points | Interval (s) | Absolute peaks | Fitted rate (s⁻¹) |
|---|---:|---:|---:|---:|
| Dense, 10 s | 20,001 | 0.0005 | 7 | 0.050001195 |
| Sparse, 10 s | 13 | 0.833333 | 6 | 0.076479020 |
| Short, 1.5 s | 3,001 | 0.0005 | 2 | unavailable (fewer than 3 peaks) |

The sparse fit overestimates the known `0.05 s⁻¹` envelope rate by about 53%
because sampled peak locations are coarse. A short window contains too few
peaks for the defined fit even though it is densely sampled. The sampling
interval and observation duration limit different parts of the inference.

The research helper also evaluates a heat perturbation with `a1(0)=0.5 K` and
`a2(0)=-1 K` at `x=0.25 m`. Its combined perturbation at
`t=(0, 100, 300, 400, 600) s` is `(-0.646447, -0.353500, -0.042999,
0.032080, 0.101956) K`. It crosses zero while each individual modal
coefficient retains its initial sign and decreases in magnitude. The signed
sum cannot be log-fit as one decaying exponential, and its absolute value can
increase near cancellation. Similarly, signed oscillator displacement crosses
zero even as its envelope declines.

Initial amplitudes affect observed magnitude but not these rates under the
linear assumptions: multiplying all heat modal amplitudes or oscillator
initial conditions by a constant multiplies the corresponding linear response
without changing its eigenvalues. Initial conditions can also suppress a
particular oscillator component or alter the observed phase.

## Cross-system comparison

| Property | Damped oscillator | Heat diffusion |
|---|---|---|
| Mathematical structure | Two-dimensional state; characteristic roots | Spatial eigenmode spectrum |
| Mode interpretation | Complex-conjugate roots when underdamped | Real sine eigenfunctions for homogeneous fixed-end perturbations |
| Rate formula | Envelope `c/(2m)` (s⁻¹) | `alpha*(n*pi/L)^2` (s⁻¹) |
| Oscillation | Yes when `zeta<1`; signed displacement crosses zero | A single coefficient does not oscillate; mode sums may cross zero |
| Parameter dependence | `c` and `m` set envelope rate; `k` also affects frequency | `alpha`, `n²`, and `1/L²` set rate |
| Key limitations | Undersampled peaks, finite window, damping regime | Boundary/material assumptions, mixed modes, finite grid and resolution |

The commonality is the exponential evolution factor for individual modes under
linear, time-invariant assumptions. Heat modes have real negative eigenvalues
`-lambda_n`. An underdamped oscillator has complex-conjugate roots whose real
part is `-gamma`; the imaginary part produces oscillation. A spatial spectrum
is not the same state representation as a two-dimensional oscillator state.
Their quantities may both have rate units of s⁻¹, but their amplitudes,
parameters, boundary/initial conditions, and physical meanings are distinct.
The comparison is framed through the dimensionless exponent `u=rate*time` and
the factor `exp(-u)`, not by equating raw amplitudes. This shared algebraic form
does not identify corresponding physical modes.

This phase reuses Phase 20's `damped_harmonic_oscillator` and
`transient_heat_diffusion_1d` registry identities. The oscillator/diffusion
structure-sharing edge remains a review candidate; no new equation records or
established cross-physics relationships are introduced, and the historical
Phase 18/19 experiment records and persistent artifacts are not rewritten.

## What is established and what is not

Within the declared equations and tested parameter values, analytical rates
match the parameter formulas. Heat attenuation is an exact evaluation of the
finite-mode solution at specified sample times and positions. The oscillator
is numerically integrated and compared with an analytical solution. Finite
sampling and finite windows materially affect estimates from those observations.

This does **not** establish that the two systems are physically equivalent,
that every observable decays monotonically, or that the numerical solver
independently validates either physical model. It does not establish a
profitable signal, predictive market model, validated risk-control method,
investment recommendation, or tradable inefficiency. No market data were used;
the Phase 18 application hypothesis remains unvalidated.

## Unresolved questions and next experiments

The next mathematical experiment should vary observation cadence, observation
floor, and mode mixtures systematically to map where the defined rate estimator
becomes biased or unavailable. A separate physical study would need measured
state variables, boundary conditions, model-specific noise, and independent
validation. Any future finance question would first need a defensible state,
intervention, objective, costs, risk constraints, and chronological out-of-
sample protocol. This report offers no financial conclusion.
