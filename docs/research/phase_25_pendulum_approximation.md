# Phase 25 - Nonlinear pendulum and small-angle approximation

## Research question and scope

How does the small-angle approximation's accuracy depend on release amplitude,
measured quantity, and observation horizon? This controlled study uses the
existing full-sine pendulum solver and compares it with the exact solution of
the linearized equation. It is an investigation of the declared mathematical
models, not validation against a laboratory pendulum.

The reproducible calculation is
`newton_lab.pendulum_approximation.run_pendulum_approximation_research()`;
the computed tables below are rendered by
`render_pendulum_approximation_report()`.

## Hypotheses

| Hypothesis | Prediction | Test | Result within the tested conditions | Limitation |
|---|---|---|---|---|
| H1 - Amplitude-dependent period | The nonlinear period increases with release amplitude while the small-angle period remains fixed. | Compare the elliptic-integral reference over the fixed amplitude grid and against `T0`. | Supported: the nonlinear period increased from 2.006079 s at 0.01 rad to 2.665871 s at 2 rad; `T0` remained 2.006067 s. | The grid stops at 2 rad and does not characterize near-separatrix behavior. |
| H2 - Observable-dependent error | Period, phase, and angle trajectory errors need not grow at the same rate. | Report period errors, maximum sampled angle errors, and state-phase lag separately over multiple windows. | Supported: at 0.5 rad, relative period error is 1.56%, maximum angle error is 0.0375 rad after 1 `T0` and 0.1359 rad after 3 `T0`; phase lag is 0.0955 rad after 1 `T0` and 0.4804 rad after 5 `T0`. | Metrics use this ideal model, the stated phase convention, and finite samples. |
| H3 - Accumulated phase error | A modest frequency/period mismatch can accumulate into a larger phase difference across cycles. | Compare unwrapped state phase against the linear phase at 1, 3, 5, and 10 small-angle periods. | Supported: at 0.2 rad, phase lag grows from 0.0156 rad at 1 `T0` to 0.1564 rad at 10 `T0`; at 0.5 rad it reaches 0.9717 rad at 10 `T0`. | The state-coordinate phase is an explicit comparison convention, not a universal phase estimator for nonlinear systems. |
| H4 - Approximation validity is tolerance-dependent | A useful amplitude bound must name its observable, tolerance, and time horizon. | Apply declared grid checks: 1% period error, 0.05 rad maximum angle error over 3 `T0`, and 0.1 rad phase lag at 5 `T0`. | Supported: the largest tested passing amplitudes were 0.35 rad for period and 3-period angle error, and 0.2 rad for 5-period phase lag. The next tested grid points failed each check. | These are grid-limited results for chosen thresholds and parameters, not universal angle cutoffs. |

## Models, assumptions, and units

The existing model is

\[
\ddot{\theta} + \frac{b}{mL^2}\dot{\theta}
  + \frac{g}{L}\sin(\theta)=0.
\]

The main experiment sets viscous torque coefficient `b=0`, uses zero initial
angular velocity, and releases from positive initial angle `theta0`. It
therefore evaluates the undamped, unforced equation

\[
\ddot{\theta}+\frac{g}{L}\sin(\theta)=0.
\]

The small-angle model substitutes `sin(theta) ~= theta`:

\[
\ddot{\theta}+\omega_0^2\theta=0,\qquad
\omega_0=\sqrt{g/L},\qquad
T_0=2\pi\sqrt{L/g}.
\]

With `theta(0)=theta0` and `theta_dot(0)=0`, its exact linear solution is
`theta_linear(t)=theta0*cos(omega0*t)`. The nonlinear and linear runs use the
same `L`, `g`, initial conditions, duration, and output times. Mass cancels
from the ideal angular equation, but remains an explicit input to the existing
simulation API.

Parameters use SI units: `m=1 kg`, `L=1 m`, `g=9.81 m/s^2`, and `b=0
kg m^2/s`; angle is in rad, angular velocity in rad/s, and time in seconds.
The assumptions are a point-mass bob, massless rigid rod, fixed pivot and
length, uniform gravitational field, no damping, and no external torque. The
existing API also permits damped simulations, but the no-damping period
formula below must not be applied unchanged to those cases.

## Analytical nonlinear period reference

For an ideal undamped pendulum released from rest at amplitude `theta0`,

\[
T(\theta_0)=4\sqrt{L/g}\,K(m),\qquad
m=\sin^2(\theta_0/2).
\]

SciPy's `scipy.special.ellipk(m)` uses the parameter convention
\[
K(m)=\int_0^{\pi/2}(1-m\sin^2 u)^{-1/2}\,du,
\]
not the modulus convention. This was checked against the installed SciPy
function documentation (`ellipk(m, ...)`). The library and dependency were
already part of the project; no dependency was added. As `theta0` tends to
zero, the expression approaches `T0`.

The nonlinear analytical reference and the period estimated from the
numerically integrated trajectory are reported separately. The estimate
detects positive-to-negative angle zero crossings and linearly interpolates
between adjacent requested samples; it averages the intervals between
successive crossings of that direction. This reduces but does not eliminate
finite-grid period-estimation error. Sampled maxima likewise do not guarantee
that continuous-time extrema were captured.

## Controlled design and numerical method

The amplitude grid was fixed before interpreting results:
`0.01, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0 rad`
(about `0.57` through `114.59 degrees`). It spans a near-linear regime through
large oscillations while stopping 1.14 rad below the unstable upright
configuration at `pi` rad. It is not intended to resolve behavior close to the
separatrix.

Each main run reuses `simulate_damped_pendulum` with `DOP853`, `rtol=1e-10`,
`atol=1e-12`, 12,001 uniformly requested output points, and duration `12*T0`
(about 24.073 s). The solver chooses adaptive internal steps; the requested
output grid is not its internal integration grid. Solver function-evaluation
counts are included in the results. The small-angle angular frequency is
`3.132092 rad/s` and its fixed period is `2.006067 s` for this parameter set.

For each amplitude, maximum absolute angle difference is computed on samples
within windows ending at 1, 3, 5, and 10 `T0`:

\[
E_\theta(W)=\max_{t_i\le W}
  |\theta_{\mathrm{nonlinear}}(t_i)-\theta_{\mathrm{linear}}(t_i)|.
\]

It is measured in radians. The reported value is a sampled maximum, not a
continuous supremum. Since the uniform run has 1,000 intervals per `T0`, the
window endpoints align with requested samples.

Phase is calculated from the nonlinear state in normalized linear
coordinates:

\[
\phi_{nl}(t)=\operatorname{unwrap}\left[
  \operatorname{atan2}\left(-\dot\theta/(\omega_0\theta_0),\theta/\theta_0\right)
\right],\qquad
\Delta\phi(t)=\omega_0t-\phi_{nl}(t).
\]

The unwrap starts from zero phase at release. Positive `Delta phi` means the
linear reference has advanced farther. This convention avoids a wrapped
single-sample difference hiding accumulated cycle lag. It is a consistent
state-space phase comparison, not a claim that nonlinear motion has a
sinusoidal phase at every amplitude.

## Results

The following table compares the two model predictions, the numerical period
estimate, and the two trajectory observables. Relative period error is
`abs(T_nonlinear-T0)/T_nonlinear`; the adjacent column uses `T0` as denominator.
`Delta T` is the analytical nonlinear period minus `T0`. Window vectors list
1, 3, 5, and 10 `T0`, in that order. All computed values below are generated
by the research module.

## Computed amplitude-sweep results

| Initial angle (rad) | Initial angle (deg) | Solver nfev | T small-angle (s) | T nonlinear analytical (s) | T nonlinear sampled estimate (s) | Estimate abs error (s) | Delta T (s) | Relative error / T nonlinear | Relative error / T small-angle | Max angle error at 1, 3, 5, 10 T0 (rad) | Phase lag at 1, 3, 5, 10 T0 (rad) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 0.01 | 0.5730 | 3302 | 2.00606668 | 2.00607922 | 2.00607922 | 9.03e-11 | 1.2538e-05 | 6.25e-06 | 6.25004e-06 | 2.97458e-07, 1.08073e-06, 1.86579e-06, 3.82906e-06 | 3.92695e-05, 0.000117809, 0.000196348, 0.000392696 |
| 0.05 | 2.8648 | 3467 | 2.00606668 | 2.00638017 | 2.00638017 | 1.38e-10 | 0.000313493 | 0.000156248 | 0.000156272 | 3.71852e-05, 0.000135103, 0.000233244, 0.00047867 | 0.000981479, 0.00294444, 0.0049074, 0.00981479 |
| 0.1 | 5.7296 | 3512 | 2.00606668 | 2.00732119 | 2.00732119 | 5.42e-11 | 0.00125451 | 0.000624967 | 0.000625358 | 0.000297561, 0.00108111, 0.00186644, 0.0038302 | 0.0039227, 0.0117681, 0.0196135, 0.039227 |
| 0.2 | 11.4592 | 3512 | 2.00606668 | 2.01109337 | 2.01109337 | 2e-10 | 0.00502669 | 0.00249948 | 0.00250574 | 0.00238303, 0.00865772, 0.0149445, 0.030647 | 0.0156393, 0.0469184, 0.0781987, 0.15641 |
| 0.35 | 20.0535 | 3527 | 2.00606668 | 2.02153429 | 2.02153429 | 3.55e-11 | 0.0154676 | 0.00765142 | 0.00771042 | 0.0128087, 0.04651, 0.0801786, 0.163381 | 0.0474644, 0.142429, 0.237495, 0.475962 |
| 0.5 | 28.6479 | 4196 | 2.00606668 | 2.03786792 | 2.03786792 | 1.25e-10 | 0.0318012 | 0.0156051 | 0.0158525 | 0.0375096, 0.135917, 0.233198, 0.464336 | 0.095531, 0.287174, 0.480362, 0.971674 |
| 0.75 | 42.9718 | 5084 | 2.00606668 | 2.07895269 | 2.07895269 | 1.85e-10 | 0.072886 | 0.035059 | 0.0363328 | 0.127949, 0.458061, 0.765263, 1.33385 | 0.208083, 0.636357, 1.08479, 2.22404 |
| 1 | 57.2958 | 5186 | 2.00606668 | 2.1391376 | 2.1391376 | 1.6e-10 | 0.133071 | 0.0622077 | 0.0663342 | 0.307614, 1.06302, 1.64285, 1.99994 | 0.356181, 1.14656, 1.97943, 3.86504 |
| 1.25 | 71.6197 | 5483 | 2.00606668 | 2.2215417 | 2.2215417 | 6.47e-11 | 0.215475 | 0.0969935 | 0.107412 | 0.610896, 1.94439, 2.49705, 2.49705 | 0.540631, 1.85586, 3.06207, 6.12344 |
| 1.5 | 85.9437 | 5531 | 2.00606668 | 2.33098678 | 2.33098678 | 7.83e-11 | 0.32492 | 0.139392 | 0.161969 | 1.07392, 2.89336, 2.99764, 2.99764 | 0.778135, 2.72103, 4.32564, 8.86106 |
| 1.75 | 100.2676 | 5516 | 2.00606668 | 2.47500563 | 2.47500563 | 1.44e-10 | 0.468939 | 0.18947 | 0.23376 | 1.7295, 3.48688, 3.48688, 3.4953 | 1.10315, 3.4567, 6.04592, 12.0478 |
| 2 | 114.5916 | 5711 | 2.00606668 | 2.66587094 | 2.66587094 | 2.5e-10 | 0.659804 | 0.2475 | 0.328904 | 2.59238, 3.99946, 3.99946, 3.99946 | 1.54993, 4.64978, 7.74963, 15.6121 |

### Numerical setting comparisons

| Amplitude (rad) | Comparison | Samples (base/comparison) | Maximum shared-sample angle difference (rad) | Period estimate base (s) | Period estimate comparison (s) | Period difference (s) |
|---:|---|---:|---:|---:|---:|---:|
| 0.1 | looser_tolerance_same_sampling | 12001/12001 | 5.91993e-09 | 2.00732119 | 2.00732119 | 1.42e-09 |
| 1 | denser_sampling_same_tolerance | 12001/24001 | 0 | 2.1391376 | 2.1391376 | -1.75e-10 |
| 2 | coarser_sampling_same_tolerance | 12001/601 | 1.87211e-14 | 2.66587094 | 2.66587287 | 1.93e-06 |

### Tolerance-based grid limits

- Relative period error <= 1% (denominator: nonlinear analytical period): largest tested passing amplitude `0.35 rad` (20.0535 deg); next tested failure at `0.5 rad`.
- Maximum angle error <= 0.05 rad over 3 T0: largest tested passing amplitude `0.35 rad` (20.0535 deg); next tested failure at `0.5 rad`.
- Accumulated phase lag <= 0.1 rad at 5 T0: largest tested passing amplitude `0.2 rad` (11.4592 deg); next tested failure at `0.35 rad`.

### Numerical settings

- Main sweep uses `DOP853`, `rtol=1e-10`, `atol=1e-12`, 12001 uniform requested output points over 24.0728 s (12 small-angle periods).
- Solver function evaluations in the first amplitude case: 3302. Adaptive internal steps differ from the uniform output grid.

## Numerical reliability and interpretation

The targeted numerical checks compare (a) looser tolerances at 0.1 rad on the
same 12,001-point grid, (b) 24,001 versus 12,001 points at 1 rad, and (c) 601
versus 12,001 points at 2 rad. They report both maximum angle difference at
shared sample times and the difference between interpolated-crossing period
estimates. The looser-tolerance trajectory differed by about `5.92e-9 rad` at
the sampled points; its estimated period differed by about `1.42e-9 s`. For
the 2 rad coarse grid, the sampled angle values at shared points agree to
roundoff, while the period estimate differs by about `1.93e-6 s`. This
distinguishes solver tolerance effects from output-grid effects on this
particular estimator. It is a focused sensitivity check, not a general proof
of convergence or a bound on continuous-time error.

The exact elliptic period increases monotonically on the tested amplitude
grid while `T0` is fixed. The numerical period estimates agree with that
reference to the displayed precision under the dense base output grid. At the
same time, the finite-window metrics expose different practical limits. For
example, at 0.5 rad the nonlinear period is about 1.56% longer than `T0`, the
sampled maximum angle error is about 0.0375 rad after one `T0` and 0.1359 rad
after three, and phase lag is about 0.4804 rad after five `T0`. A period-only
assessment would not describe those trajectory and phase differences.

The predeclared tolerance checks find a largest tested passing amplitude of
0.35 rad (20.05 degrees) for 1% relative period error and 0.05 rad maximum
angle error over three `T0`. Both next fail at 0.5 rad. For 0.1 rad accumulated
phase lag after five `T0`, the largest tested passing grid point is 0.2 rad
(11.46 degrees), with 0.35 rad the next failure. These limits depend on the
observable, chosen tolerance, and horizon. They are statements about this
grid, not exact crossing angles or universal small-angle validity thresholds.

## What this establishes and what it does not

The calculations demonstrate consequences of the ideal equations: the
nonlinear period depends on release amplitude; the linearized period does
not; trajectory and phase errors change with both amplitude and observation
window; and a dense numerical estimate can reproduce the analytical period
reference for the tested conditions. They do not validate either model
against real pendulum measurements, quantify unmodeled friction or rod
effects, establish continuous-time error bounds, or prove a universal angle
cutoff. They also establish no trading signal, market model, risk-control
method, or financial application. Any future analogy involving accumulated
phase error would be a hypothesis requiring independent empirical validation.

## Reproduction and next question

From the project root with the existing environment:

```powershell
.venv\Scripts\python.exe -c "from newton_lab.pendulum_approximation import run_pendulum_approximation_research, render_pendulum_approximation_report; print(render_pendulum_approximation_report(run_pendulum_approximation_research()))"
```

A useful next research question is how the stated tolerance limits change
under a measured damping law or finite-amplitude initial angular velocity,
while keeping the period reference appropriate to that altered model.
