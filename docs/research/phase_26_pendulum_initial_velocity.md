# Phase 26 - Pendulum approximation with nonzero initial angular velocity

## Research question and hypotheses

How do initial displacement and initial angular velocity jointly affect the
small-angle approximation's period, sampled trajectory error, accumulated
phase difference, and tolerance-dependent limits? The study uses Newton Lab's
existing undamped full-sine pendulum solver. It compares the nonlinear model
with the exact solution of its small-angle linearization under matched
physical parameters and initial state.

| Hypothesis | Prediction and test | Result in the declared grid | Limitation |
|---|---|---|---|
| H1 - Energy sets the turning amplitude and period. | Use energy conservation to compute turning amplitude; compare the elliptic period with zero-velocity Phase 25 and numerical crossing estimates. | Supported: reversing velocity preserves energy, turning amplitude, and analytical period. Nonzero speed can increase turning amplitude substantially; e.g. `(theta0, omega0)=(0.5 rad, 3 rad/s)` turns at `1.13860 rad`. | Ideal undamped librations only; no experimental measurements. |
| H2 - Initial direction affects observables beyond period. | Compare positive and negative initial speeds at the same initial angle, using matched time windows and initial-state phase. | Supported: at `theta0=0.5 rad`, `omega0=+1` and `-1 rad/s` have the same exact period `2.051802 s`, but their maximum angle errors after one linear period are `0.08037` and `0.04688 rad`; phase differences are `0.13979` and `0.13831 rad`. | These finite-window results depend on the chosen initial state and phase convention. |
| H3 - Phase error accumulates. | Compare unwrapped linear-coordinate phase from both trajectories at 1, 3, 5, and 10 linear periods. | Supported: for `(0.5 rad, +3 rad/s)`, phase difference grows from `0.54957 rad` after one linear period to `5.11159 rad` after ten. | The coordinate phase is a declared comparison convention, not a universal nonlinear phase estimator. |
| H4 - Validity limits depend on the observable and initial velocity. | Apply separate period, sampled-angle, and accumulated-phase tolerances to velocity-conditional grid slices. | Supported: at rest the largest passing tested initial angle is `0.35 rad` for period and 3-period angle error, but `0.2 rad` for 5-period phase difference. At `omega0=+/-1 rad/s`, the corresponding period limit is `0.2 rad`, angle limit `0.1 rad`, and no tested angle meets the 5-period phase criterion. | These are discrete-grid results conditional on signed initial velocity and selected thresholds, not universal bounds. |

## Equations, units, and assumptions

The reused model is the undamped, unforced equation

\[
\ddot{\theta}+\frac{g}{L}\sin\theta=0,
\]

and the linearized comparison is

\[
\ddot{\theta}_{lin}+\Omega^2\theta_{lin}=0,
\qquad \Omega=\sqrt{g/L}.
\]

Both trajectories start from the same `(theta0, omega_initial)`:

\[
\theta_{lin}(t)=\theta_0\cos(\Omega t)
 +\frac{\omega_{initial}}{\Omega}\sin(\Omega t),
\qquad T_{lin}=2\pi/\Omega.
\]

Code and report use `initial_angular_velocity_rad_per_s` for the initial
angular velocity and `small_angle_natural_frequency_rad_per_s` for `Omega`.
The parameters are `m=1 kg`, `L=1 m`, `g=9.81 m/s^2`, damping coefficient
`b=0 kg m^2/s`; angle is rad, angular velocity rad/s, and time seconds. The
mass cancels from the ideal angular equation but is passed to the existing
solver. Assumptions are a point-mass bob on a fixed rigid massless rod, a
uniform gravitational field, no damping, and no external torque. The initial
angle uses the principal convention `[-pi, pi]`; it is not silently reduced
modulo `2*pi`.

## Energy-derived turning amplitude and state classification

For the specified ideal model, the conserved energy divided by `m L^2` is

\[
\mathcal E=\frac12\omega_{initial}^2
 +\frac{g}{L}(1-\cos\theta_0).
\]

Define the dimensionless energy level

\[
h=\frac{L}{g}\mathcal E
 =1-\cos\theta_0+\frac{L\omega_{initial}^2}{2g}.
\]

For numerical stability near zero, the implementation evaluates
`1-cos(theta0)` as `2*sin(theta0/2)^2`. At a libration turning point,
`omega=0`, so `1-cos(theta_max)=h` and

\[
\theta_{max}=2\arcsin\sqrt{h/2},\qquad 0<h<2.
\]

The returned `theta_max` is a nonnegative magnitude. It equals `abs(theta0)`
only for release from rest; in general initial angle and turning amplitude
differ. A dimensionless level `h=2` is the separatrix, and `h>2` is rotation.
Equilibrium `(0,0)` is separately classified because it has no observed
oscillation period. A roundoff-sized energy excursion within 64 machine
epsilons of the separatrix is classified as separatrix-level and is not
clamped to a libration. The near-separatrix label flags librations with
`2-h <= 0.02`; this is a diagnostic convention for this study, not a new
physical boundary.

For valid librations, the analytical nonlinear period is

\[
T_{nl}=4\sqrt{L/g}\,K(m),\qquad
m=\sin^2(\theta_{max}/2)=h/2.
\]

SciPy `ellipk(m)` takes the parameter `m`, not the modulus. Its installed
documentation defines `K(m)=integral_0^(pi/2) (1-m sin^2 u)^(-1/2) du`.
The implementation refuses to evaluate this finite libration period for
rotation, separatrix, or equilibrium. At zero initial velocity its value
matches the Phase 25 helper; as `theta_max` tends to zero it approaches
`Tlin`.

## Experimental grid and measurements

The ordered core grid crosses initial angles
`0.1, 0.2, 0.35, 0.5, 1.2, 2.0 rad` with initial angular velocities
`-3, -1, 0, +1, +3 rad/s`. Additional cases include high-energy librations
`(0.1, +/-4.3)`, negative initial displacement `(-0.5, +/-1)`, near-boundary
rest releases at `3.0` and `3.12 rad`, an exact separatrix control, a 7 rad/s
rotational control, stable equilibrium, and an out-of-principal-range invalid
input. Sweep index, case ID, initial values, solver settings, and status are
retained in order. The grid includes paired speed signs and deliberately
substantial energy additions: for `(0.5, +/-3)`, the turning amplitude is
`1.13860 rad`, more than `0.63 rad` above the initial displacement magnitude.

The default accuracy criteria are fixed in the API: relative period error
`abs(Tnl-Tlin)/Tnl <= 0.01`; sampled maximum angle error `<=0.05 rad` over 3
linear periods; and absolute accumulated phase difference `<=0.1 rad` over 5
linear periods. The relative-period denominator and signed period difference
`Tnl-Tlin` are recorded explicitly. Criteria are configurable through
`AccuracyCriteria`; supported trajectory windows are 1, 3, 5, and 10 linear
periods.

The nonlinear trajectory uses `simulate_damped_pendulum` with `DOP853`,
`rtol=1e-10`, `atol=1e-12`, and 24,001 uniform requested output samples. Its
duration is `max(12*Tlin, 4*Tnl)`, so the 10-linear-period error window exists
and several nonlinear cycles are available, including near the boundary.
Internal solver steps are adaptive. A numerical period estimate collects
both positive-going and negative-going zero crossings, selects a direction
with repeated crossings, linearly interpolates each crossing between output
samples, and averages same-direction intervals. It does not assume that the
first downward crossing occurs at a fixed fraction of the period. A time-zero
crossing is included when the initial angle is zero and velocity is nonzero.
Fewer than two same-direction crossings produce an explicit incomplete-data
status, not a substituted period. Solver failures, invalid inputs, equilibrium,
separatrix, and rotation also remain explicit result statuses.

Trajectory error is the maximum sampled absolute difference between nonlinear
and linear angle at matching output times, in radians. It is not a guaranteed
continuous-time maximum. Phase is

\[
\phi(t)=\operatorname{unwrap}\left[\operatorname{atan2}
  (-\dot\theta/\Omega,\theta)\right],\qquad
\Delta\phi=\phi_{lin}-\phi_{nl}.
\]

Both models share their initial state, so their phase branches are aligned at
release; the initial phase `atan2(-omega_initial/Omega, theta0)` captures the
initial direction. Positive `Delta phi` means the linear reference is ahead.
The difference is unwrapped and therefore records accumulated cycle lag,
distinct from instantaneous angle error.

## Computed results

The generated tables below include each ordered case, energy class, initial
and turning amplitudes, analytical and numerical periods, signed and relative
period errors, crossing diagnostics, angle/phase windows, criteria outcomes,
conditional grid limits, and numerical sensitivities. Nonperiodic and invalid
controls have no fabricated period values.

## Computed initial-condition sweep

| # | Case | Initial angle (rad) | Initial angular velocity (rad/s) | Initial linear phase (rad) | State class | Outcome | Energy h | Separatrix gap | Turning amplitude (rad) | Amplitude increase (rad) | T nonlinear exact (s) | T nonlinear sampled (s) | T linear (s) | Signed Delta T (s) | Relative error / T nonlinear | Period crossings | Angle error at 1, 3, 5, 10 Tlinear (rad) | Phase difference at 1, 3, 5, 10 Tlinear (rad) | Criteria P/A/phase |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | grid_theta_0.1_omega_-3 | 0.1 | -3 | 1.46677 | librating | analyzed | 0.4637114 | 1.5363 | 1.004763 | 0.90476 | 2.1404877 | 2.1404877 | 2.0060667 | 0.134421 | 0.0627992 | negative:12 | 0.38299, 0.96777, 1.5682, 1.966 | 0.41995, 1.2049, 1.9257, 3.9838 | fail/fail/fail |
| 1 | grid_theta_0.1_omega_-1 | 0.1 | -1 | 1.267265 | librating | analyzed | 0.05596423 | 1.944 | 0.3361374 | 0.23614 | 2.0203255 | 2.0203255 | 2.0060667 | 0.0142588 | 0.00705766 | positive:12 | 0.014155, 0.041755, 0.068228, 0.14089 | 0.044651, 0.13393, 0.22318, 0.44593 | pass/pass/fail |
| 2 | grid_theta_0.1_omega_+0 | 0.1 | 0 | -0 | librating | analyzed | 0.004995835 | 1.995 | 0.1 | 0 | 2.0073212 | 2.0073212 | 2.0060667 | 0.00125451 | 0.000624967 | positive:12 | 0.00029756, 0.0010811, 0.0018664, 0.0038302 | 0.0039227, 0.011768, 0.019613, 0.039227 | pass/pass/pass |
| 3 | grid_theta_0.1_omega_+1 | 0.1 | 1 | -1.267265 | librating | analyzed | 0.05596423 | 1.944 | 0.3361374 | 0.23614 | 2.0203255 | 2.0203255 | 2.0060667 | 0.0142588 | 0.00705766 | negative:12 | 0.014937, 0.044788, 0.074549, 0.14817 | 0.044655, 0.13397, 0.22329, 0.44659 | pass/pass/fail |
| 4 | grid_theta_0.1_omega_+3 | 0.1 | 3 | -1.46677 | librating | analyzed | 0.4637114 | 1.5363 | 1.004763 | 0.90476 | 2.1404877 | 2.1404877 | 2.0060667 | 0.134421 | 0.0627992 | positive:11 | 0.40021, 0.99334, 1.5682, 1.966 | 0.42092, 1.2303, 1.9505, 3.994 | fail/fail/fail |
| 5 | grid_theta_0.2_omega_-3 | 0.2 | -3 | 1.364948 | librating | analyzed | 0.478649 | 1.5214 | 1.022363 | 0.82236 | 2.1455467 | 2.1455467 | 2.0060667 | 0.13948 | 0.0650091 | negative:12 | 0.3872, 1.0167, 1.6341, 2.0007 | 0.43405, 1.2281, 1.9828, 4.1102 | fail/fail/fail |
| 6 | grid_theta_0.2_omega_-1 | 0.2 | -1 | 1.011178 | librating | analyzed | 0.07090182 | 1.9291 | 0.3788295 | 0.17883 | 2.0242094 | 2.0242094 | 2.0060667 | 0.0181427 | 0.00896287 | positive:12 | 0.017813, 0.053548, 0.096137, 0.20095 | 0.056697, 0.16994, 0.28291, 0.5638 | pass/fail/fail |
| 7 | grid_theta_0.2_omega_+0 | 0.2 | 0 | -0 | librating | analyzed | 0.01993342 | 1.9801 | 0.2 | 0 | 2.0110934 | 2.0110934 | 2.0060667 | 0.00502669 | 0.00249948 | positive:12 | 0.002383, 0.0086578, 0.014945, 0.030647 | 0.015639, 0.046918, 0.078199, 0.15641 | pass/pass/pass |
| 8 | grid_theta_0.2_omega_+1 | 0.2 | 1 | -1.011178 | librating | analyzed | 0.07090182 | 1.9291 | 0.3788295 | 0.17883 | 2.0242094 | 2.0242094 | 2.0060667 | 0.0181427 | 0.00896287 | negative:12 | 0.021317, 0.064049, 0.10657, 0.21109 | 0.056736, 0.17029, 0.28391, 0.56802 | pass/fail/fail |
| 9 | grid_theta_0.2_omega_+3 | 0.2 | 3 | -1.364948 | librating | analyzed | 0.478649 | 1.5214 | 1.022363 | 0.82236 | 2.1455467 | 2.1455467 | 2.0060667 | 0.13948 | 0.0650091 | positive:11 | 0.42402, 1.0833, 1.6341, 2.0007 | 0.43683, 1.2824, 2.0279, 4.1426 | fail/fail/fail |
| 10 | grid_theta_0.35_omega_-3 | 0.35 | -3 | 1.220459 | librating | analyzed | 0.5193429 | 1.4807 | 1.069392 | 0.71939 | 2.1596142 | 2.1596142 | 2.0060667 | 0.153548 | 0.0710995 | negative:12 | 0.40595, 1.154, 1.8103, 2.0866 | 0.47148, 1.3044, 2.1671, 4.4332 | fail/fail/fail |
| 11 | grid_theta_0.35_omega_-1 | 0.35 | -1 | 0.7395232 | librating | analyzed | 0.1115957 | 1.8884 | 0.4769389 | 0.12694 | 2.0349641 | 2.0349641 | 2.0060667 | 0.0288974 | 0.0142004 | positive:12 | 0.027457, 0.10783, 0.19235, 0.3953 | 0.089335, 0.26689, 0.44292, 0.87923 | fail/fail/fail |
| 12 | grid_theta_0.35_omega_+0 | 0.35 | 0 | -0 | librating | analyzed | 0.06062729 | 1.9394 | 0.35 | 0 | 2.0215343 | 2.0215343 | 2.0060667 | 0.0154676 | 0.00765142 | positive:12 | 0.012809, 0.04651, 0.080179, 0.16338 | 0.047464, 0.14243, 0.2375, 0.47596 | pass/pass/fail |
| 13 | grid_theta_0.35_omega_+1 | 0.35 | 1 | -0.7395232 | librating | analyzed | 0.1115957 | 1.8884 | 0.4769389 | 0.12694 | 2.0349641 | 2.0349641 | 2.0060667 | 0.0288974 | 0.0142004 | negative:12 | 0.041816, 0.12707, 0.21129, 0.41283 | 0.089663, 0.26975, 0.45047, 0.90296 | fail/fail/fail |
| 14 | grid_theta_0.35_omega_+3 | 0.35 | 3 | -1.220459 | librating | analyzed | 0.5193429 | 1.4807 | 1.069392 | 0.71939 | 2.1596142 | 2.1596142 | 2.0060667 | 0.153548 | 0.0710995 | positive:11 | 0.48541, 1.2646, 1.8101, 2.0866 | 0.48075, 1.4134, 2.2229, 4.541 | fail/fail/fail |
| 15 | grid_theta_0.5_omega_-3 | 0.5 | -3 | 1.089692 | librating | analyzed | 0.581133 | 1.4189 | 1.138599 | 0.6386 | 2.1818172 | 2.1818172 | 2.0060667 | 0.175751 | 0.0805524 | positive:11 | 0.43473, 1.373, 2.0621, 2.2189 | 0.52519, 1.4354, 2.478, 4.9667 | fail/fail/fail |
| 16 | grid_theta_0.5_omega_-1 | 0.5 | -1 | 0.5682845 | librating | analyzed | 0.1733858 | 1.8266 | 0.5977319 | 0.097732 | 2.0518016 | 2.0518016 | 2.0060667 | 0.0457349 | 0.0222901 | positive:12 | 0.046879, 0.21486, 0.37859, 0.7473 | 0.13831, 0.4109, 0.68085, 1.3706 | fail/fail/fail |
| 17 | grid_theta_0.5_omega_+0 | 0.5 | 0 | -0 | librating | analyzed | 0.1224174 | 1.8776 | 0.5 | 0 | 2.0378679 | 2.0378679 | 2.0060667 | 0.0318012 | 0.0156051 | positive:12 | 0.03751, 0.13592, 0.2332, 0.46434 | 0.095531, 0.28717, 0.48036, 0.97167 | fail/fail/fail |
| 18 | grid_theta_0.5_omega_+1 | 0.5 | 1 | -0.5682845 | librating | analyzed | 0.1733858 | 1.8266 | 0.5977319 | 0.097732 | 2.0518016 | 2.0518016 | 2.0060667 | 0.0457349 | 0.0222901 | negative:12 | 0.080369, 0.24783, 0.41029, 0.7731 | 0.13979, 0.42302, 0.70889, 1.425 | fail/fail/fail |
| 19 | grid_theta_0.5_omega_+3 | 0.5 | 3 | -1.089692 | librating | analyzed | 0.581133 | 1.4189 | 1.138599 | 0.6386 | 2.1818172 | 2.1818172 | 2.0060667 | 0.175751 | 0.0805524 | positive:11 | 0.58573, 1.5081, 2.0614, 2.2189 | 0.54957, 1.6085, 2.5152, 5.1116 | fail/fail/fail |
| 20 | grid_theta_1.2_omega_-3 | 1.2 | -3 | 0.6736354 | librating | analyzed | 1.096358 | 0.90364 | 1.667304 | 0.4673 | 2.4229145 | 2.4229145 | 2.0060667 | 0.416848 | 0.172044 | positive:10 | 1.0128, 3.1937, 3.1937, 3.2011 | 0.87907, 3.2443, 5.2536, 10.557 | fail/fail/fail |
| 21 | grid_theta_1.2_omega_-1 | 1.2 | -1 | 0.2600386 | librating | analyzed | 0.6886106 | 1.3114 | 1.254142 | 0.054142 | 2.2231182 | 2.2231182 | 2.0060667 | 0.217052 | 0.0976338 | positive:11 | 0.52172, 1.8918, 2.489, 2.489 | 0.52298, 1.7915, 3.0757, 6.1499 | fail/fail/fail |
| 22 | grid_theta_1.2_omega_+0 | 1.2 | 0 | -0 | librating | analyzed | 0.6376422 | 1.3624 | 1.2 | 0 | 2.2030826 | 2.2030826 | 2.0060667 | 0.197016 | 0.0894274 | positive:11 | 0.53862, 1.7521, 2.3652, 2.3992 | 0.50032, 1.6969, 2.8536, 5.6832 | fail/fail/fail |
| 23 | grid_theta_1.2_omega_+1 | 1.2 | 1 | -0.2600386 | librating | analyzed | 0.6886106 | 1.3114 | 1.254142 | 0.054142 | 2.2231182 | 2.2231182 | 2.0060667 | 0.217052 | 0.0976338 | positive:11 | 0.70285, 2.0065, 2.464, 2.4959 | 0.59478, 1.9437, 3.077, 6.155 | fail/fail/fail |
| 24 | grid_theta_1.2_omega_+3 | 1.2 | 3 | -0.6736354 | librating | analyzed | 1.096358 | 0.90364 | 1.667304 | 0.4673 | 2.4229145 | 2.4229145 | 2.0060667 | 0.416848 | 0.172044 | positive:10 | 1.7972, 3.186, 3.186, 3.2027 | 1.2739, 3.2501, 5.5546, 11.054 | fail/fail/fail |
| 25 | grid_theta_2_omega_-3 | 2 | -3 | 0.4466363 | librating | analyzed | 1.874862 | 0.12514 | 2.635948 | 0.63595 | 3.5767529 | 3.5767529 | 2.0060667 | 1.57069 | 0.439137 | positive:7 | 4.1887, 4.748, 4.7701, 4.8123 | 2.624, 7.573, 13.154, 27.23 | fail/fail/fail |
| 26 | grid_theta_2_omega_-1 | 2 | -1 | 0.158302 | librating | analyzed | 1.467115 | 0.53288 | 2.056822 | 0.056822 | 2.7177295 | 2.7177295 | 2.0060667 | 0.711663 | 0.261859 | positive:9 | 2.5251, 4.0821, 4.0821, 4.0821 | 1.4725, 4.8159, 8.1594, 16.174 | fail/fail/fail |
| 27 | grid_theta_2_omega_+0 | 2 | 0 | -0 | librating | analyzed | 1.416147 | 0.58385 | 2 | 0 | 2.6658709 | 2.6658709 | 2.0060667 | 0.659804 | 0.2475 | positive:9 | 2.5924, 3.9995, 3.9995, 3.9995 | 1.5499, 4.6498, 7.7496, 15.612 | fail/fail/fail |
| 28 | grid_theta_2_omega_+1 | 2 | 1 | -0.158302 | librating | analyzed | 1.467115 | 0.53288 | 2.056822 | 0.056822 | 2.7177295 | 2.7177295 | 2.0060667 | 0.711663 | 0.261859 | positive:9 | 3.0331, 4.0325, 4.0325, 4.0794 | 1.8693, 5.2051, 8.5279, 16.377 | fail/fail/fail |
| 29 | grid_theta_2_omega_+3 | 2 | 3 | -0.4466363 | librating | analyzed | 1.874862 | 0.12514 | 2.635948 | 0.63595 | 3.5767529 | 3.5767529 | 2.0060667 | 1.57069 | 0.439137 | negative:7 | 4.4644, 4.8254, 4.8254, 4.8512 | 2.8705, 8.8718, 14.542, 27.881 | fail/fail/fail |
| 30 | high_energy_libration_theta_0_1_omega_+4.3 | 0.1 | 4.3 | -1.498085 | librating | analyzed | 0.9474015 | 1.0526 | 1.518174 | 1.4182 | 2.3401772 | 2.3401772 | 2.0060667 | 0.334111 | 0.142771 | positive:10 | 1.2529, 2.7163, 2.8752, 2.8752 | 1.0135, 2.6201, 4.5586, 8.8991 | fail/fail/fail |
| 31 | high_energy_libration_theta_0_1_omega_-4.3 | 0.1 | -4.3 | 1.498085 | librating | analyzed | 0.9474015 | 1.0526 | 1.518174 | 1.4182 | 2.3401772 | 2.3401772 | 2.0060667 | 0.334111 | 0.142771 | negative:11 | 1.1653, 2.7163, 2.8752, 2.8752 | 0.99002, 2.6172, 4.5124, 8.8961 | fail/fail/fail |
| 32 | negative_angle_positive_velocity | -0.5 | 1 | -2.573308 | librating | analyzed | 0.1733858 | 1.8266 | 0.5977319 | 0.097732 | 2.0518016 | 2.0518016 | 2.0060667 | 0.0457349 | 0.0222901 | positive:12 | 0.046879, 0.21486, 0.37859, 0.7473 | 0.13831, 0.4109, 0.68085, 1.3706 | fail/fail/fail |
| 33 | negative_angle_negative_velocity | -0.5 | -1 | 2.573308 | librating | analyzed | 0.1733858 | 1.8266 | 0.5977319 | 0.097732 | 2.0518016 | 2.0518016 | 2.0060667 | 0.0457349 | 0.0222901 | positive:12 | 0.080369, 0.24783, 0.41029, 0.7731 | 0.13979, 0.42302, 0.70889, 1.425 | fail/fail/fail |
| 34 | near_boundary_3_0_rest | 3 | 0 | -0 | near_separatrix_libration | analyzed | 1.989992 | 0.010008 | 3 | 4.4409e-16 | 5.1580668 | 5.1580667 | 2.0060667 | 3.152 | 0.611082 | negative:5 | 5.7063, 5.9863, 5.9863, 6 | 3.2912, 12.136, 18.898, 37.848 | fail/fail/fail |
| 35 | near_boundary_3_12_rest | 3.12 | 0 | -0 | near_separatrix_libration | analyzed | 1.999767 | 0.00023311 | 3.12 | 8.8818e-16 | 7.5540647 | 7.5540643 | 2.0060667 | 5.548 | 0.734439 | positive:4 | 6.0208, 6.2339, 6.2339, 6.2396 | 4.345, 13.137, 22.233, 46.969 | fail/fail/fail |
| 36 | separatrix_control | 0 | 6.2642 | -1.570796 | separatrix | separatrix_period_not_finite | 2 | 0 | n/a | n/a | n/a | n/a | 2.0060667 | n/a | n/a | n/a | n/a, n/a, n/a, n/a | n/a, n/a, n/a, n/a | n/a/n/a/n/a |
  Reason for `separatrix_control`: energy is at the separatrix within floating-point roundoff.
| 37 | rotational_control | 0 | 7 | -1.570796 | rotating | rotation_excluded_from_libration_period | 2.497452 | -0.49745 | n/a | n/a | n/a | n/a | 2.0060667 | n/a | n/a | n/a | n/a, n/a, n/a, n/a | n/a, n/a, n/a, n/a | n/a/n/a/n/a |
  Reason for `rotational_control`: normalized energy exceeds the separatrix value 2.
| 38 | equilibrium_control | 0 | 0 | -0 | stable_equilibrium | equilibrium_has_no_observed_period | 0 | 2 | 0 | 0 | n/a | n/a | 2.0060667 | n/a | n/a | n/a | n/a, n/a, n/a, n/a | n/a, n/a, n/a, n/a | n/a/n/a/n/a |
  Reason for `equilibrium_control`: zero angle and zero angular velocity are stable equilibrium.
| 39 | invalid_principal_angle_control | 3.3 | 0 | n/a | invalid_initial_condition | invalid_input | n/a | n/a | n/a | n/a | n/a | n/a | 2.0060667 | n/a | n/a | n/a | n/a, n/a, n/a, n/a | n/a, n/a, n/a, n/a | n/a/n/a/n/a |
  Reason for `invalid_principal_angle_control`: initial angle must use the principal convention [-pi, pi] rad.

## Conditional grid limits

| Observable | Initial angular velocity (rad/s) | Largest passing absolute initial angle (rad) | Next higher tested angle (rad) | Next point failed? | Passing cases | Failing cases |
|---|---:|---:|---:|---|---|---|
| relative_period_error | -3 | n/a | n/a | n/a | none | grid_theta_0.1_omega_-3, grid_theta_0.2_omega_-3, grid_theta_0.35_omega_-3, grid_theta_0.5_omega_-3, grid_theta_1.2_omega_-3, grid_theta_2_omega_-3 |
| relative_period_error | -1 | 0.2 | 0.35 | true | grid_theta_0.1_omega_-1, grid_theta_0.2_omega_-1 | grid_theta_0.35_omega_-1, grid_theta_0.5_omega_-1, negative_angle_negative_velocity, grid_theta_1.2_omega_-1, grid_theta_2_omega_-1 |
| relative_period_error | 0 | 0.35 | 0.5 | true | grid_theta_0.1_omega_+0, grid_theta_0.2_omega_+0, grid_theta_0.35_omega_+0 | grid_theta_0.5_omega_+0, grid_theta_1.2_omega_+0, grid_theta_2_omega_+0, near_boundary_3_0_rest, near_boundary_3_12_rest |
| relative_period_error | 1 | 0.2 | 0.35 | true | grid_theta_0.1_omega_+1, grid_theta_0.2_omega_+1 | grid_theta_0.35_omega_+1, grid_theta_0.5_omega_+1, negative_angle_positive_velocity, grid_theta_1.2_omega_+1, grid_theta_2_omega_+1 |
| relative_period_error | 3 | n/a | n/a | n/a | none | grid_theta_0.1_omega_+3, grid_theta_0.2_omega_+3, grid_theta_0.35_omega_+3, grid_theta_0.5_omega_+3, grid_theta_1.2_omega_+3, grid_theta_2_omega_+3 |
| relative_period_error | 4.3 | n/a | n/a | n/a | none | high_energy_libration_theta_0_1_omega_+4.3 |
| relative_period_error | -4.3 | n/a | n/a | n/a | none | high_energy_libration_theta_0_1_omega_-4.3 |
| sampled_angle_error | -3 | n/a | n/a | n/a | none | grid_theta_0.1_omega_-3, grid_theta_0.2_omega_-3, grid_theta_0.35_omega_-3, grid_theta_0.5_omega_-3, grid_theta_1.2_omega_-3, grid_theta_2_omega_-3 |
| sampled_angle_error | -1 | 0.1 | 0.2 | true | grid_theta_0.1_omega_-1 | grid_theta_0.2_omega_-1, grid_theta_0.35_omega_-1, grid_theta_0.5_omega_-1, negative_angle_negative_velocity, grid_theta_1.2_omega_-1, grid_theta_2_omega_-1 |
| sampled_angle_error | 0 | 0.35 | 0.5 | true | grid_theta_0.1_omega_+0, grid_theta_0.2_omega_+0, grid_theta_0.35_omega_+0 | grid_theta_0.5_omega_+0, grid_theta_1.2_omega_+0, grid_theta_2_omega_+0, near_boundary_3_0_rest, near_boundary_3_12_rest |
| sampled_angle_error | 1 | 0.1 | 0.2 | true | grid_theta_0.1_omega_+1 | grid_theta_0.2_omega_+1, grid_theta_0.35_omega_+1, grid_theta_0.5_omega_+1, negative_angle_positive_velocity, grid_theta_1.2_omega_+1, grid_theta_2_omega_+1 |
| sampled_angle_error | 3 | n/a | n/a | n/a | none | grid_theta_0.1_omega_+3, grid_theta_0.2_omega_+3, grid_theta_0.35_omega_+3, grid_theta_0.5_omega_+3, grid_theta_1.2_omega_+3, grid_theta_2_omega_+3 |
| sampled_angle_error | 4.3 | n/a | n/a | n/a | none | high_energy_libration_theta_0_1_omega_+4.3 |
| sampled_angle_error | -4.3 | n/a | n/a | n/a | none | high_energy_libration_theta_0_1_omega_-4.3 |
| absolute_accumulated_phase_difference | -3 | n/a | n/a | n/a | none | grid_theta_0.1_omega_-3, grid_theta_0.2_omega_-3, grid_theta_0.35_omega_-3, grid_theta_0.5_omega_-3, grid_theta_1.2_omega_-3, grid_theta_2_omega_-3 |
| absolute_accumulated_phase_difference | -1 | n/a | n/a | n/a | none | grid_theta_0.1_omega_-1, grid_theta_0.2_omega_-1, grid_theta_0.35_omega_-1, grid_theta_0.5_omega_-1, negative_angle_negative_velocity, grid_theta_1.2_omega_-1, grid_theta_2_omega_-1 |
| absolute_accumulated_phase_difference | 0 | 0.2 | 0.35 | true | grid_theta_0.1_omega_+0, grid_theta_0.2_omega_+0 | grid_theta_0.35_omega_+0, grid_theta_0.5_omega_+0, grid_theta_1.2_omega_+0, grid_theta_2_omega_+0, near_boundary_3_0_rest, near_boundary_3_12_rest |
| absolute_accumulated_phase_difference | 1 | n/a | n/a | n/a | none | grid_theta_0.1_omega_+1, grid_theta_0.2_omega_+1, grid_theta_0.35_omega_+1, grid_theta_0.5_omega_+1, negative_angle_positive_velocity, grid_theta_1.2_omega_+1, grid_theta_2_omega_+1 |
| absolute_accumulated_phase_difference | 3 | n/a | n/a | n/a | none | grid_theta_0.1_omega_+3, grid_theta_0.2_omega_+3, grid_theta_0.35_omega_+3, grid_theta_0.5_omega_+3, grid_theta_1.2_omega_+3, grid_theta_2_omega_+3 |
| absolute_accumulated_phase_difference | 4.3 | n/a | n/a | n/a | none | high_energy_libration_theta_0_1_omega_+4.3 |
| absolute_accumulated_phase_difference | -4.3 | n/a | n/a | n/a | none | high_energy_libration_theta_0_1_omega_-4.3 |

## Numerical sensitivity

| Case | Comparison | Samples base/other | Max angle difference at shared samples (rad) | Period base (s) | Period other (s) | Period estimate difference (s) | Solver success base/other |
|---|---|---:|---:|---:|---:|---:|---|
| grid_theta_0.1_omega_+1 | looser_tolerances | 24001/24001 | 1.327891e-07 | 2.020325 | 2.020325 | 6.188718e-09 | True/True |
| grid_theta_0.1_omega_+1 | coarser_output_grid | 24001/601 | 3.726186e-15 | 2.020325 | 2.020326 | 1.863965e-07 | True/True |
| grid_theta_0.5_omega_+3 | looser_tolerances | 24001/24001 | 2.436348e-06 | 2.181817 | 2.181817 | -6.84224e-08 | True/True |
| grid_theta_0.5_omega_+3 | coarser_output_grid | 24001/601 | 1.201816e-14 | 2.181817 | 2.181818 | 9.569762e-07 | True/True |
| near_boundary_3_12_rest | looser_tolerances | 24001/24001 | 0.01075938 | 7.554064 | 7.553512 | -0.0005521561 | True/True |
| near_boundary_3_12_rest | coarser_output_grid | 24001/601 | 2.220446e-14 | 7.554064 | 7.554064 | 7.645422e-10 | True/True |

## Numerical settings and criteria

- Main solver: DOP853, rtol=1e-10, atol=1e-12, 24001 uniform requested samples; duration is max(12 linear periods, 4 exact nonlinear periods).
- Criteria: relative period error <= 0.01; sampled angle error <= 0.05 rad over 3 linear periods; absolute phase difference <= 0.1 rad over 5 linear periods.

## Findings and numerical sensitivity

For the zero-velocity control at `theta0=0.5 rad`, the energy-derived
turning amplitude is exactly `0.5 rad` and the nonlinear period is
`2.037868 s`, compared with `Tlin=2.006067 s`. At the same initial angle and
`omega_initial=+/-1 rad/s`, the turning amplitude rises to `0.597732 rad` and
the analytical period is `2.051802 s`. Reversing velocity preserves energy
and period but changes initial phase and the path through a finite comparison
window. At `(0.5 rad, +3 rad/s)`, the turning amplitude is `1.138599 rad`, the
period is `2.181817 s`, and phase difference reaches `5.111594 rad` after ten
linear periods.

At rest, the largest tested passing absolute initial angle is 0.35 rad for
period and 3-period angle criteria; the next tested angle 0.5 rad fails both.
For phase at five linear periods, 0.2 rad is the largest passing point and
0.35 rad fails. Conditional limits differ for `omega_initial=+/-1 rad/s`:
period passes through 0.2 rad, sampled angle error through 0.1 rad, and the
phase criterion has no passing grid point. At `+/-3 rad/s` and `+/-4.3 rad/s`,
none of the tested initial displacements passes the default criteria. These
are conditional findings on the discrete grid; the full passing and failing
case IDs and each next-higher tested point are shown above. They are not
interpolated threshold estimates and do not collapse the two-dimensional
initial-condition domain to one universal amplitude.

Sensitivity tests compare base and looser tolerances on the same grid, and
base and coarser output grids, for a small, moderate/high-energy, and
near-separatrix case. The near-boundary state has `h=1.9997669`, only
`0.0002331` below the separatrix, and an analytical period of `7.554065 s`.
Under looser tolerances (`rtol=1e-7`, `atol=1e-9`) its sampled angle path
differs by about `0.01076 rad` and estimated period by about `-0.000552 s`.
The small case's corresponding differences are about `1.33e-7 rad` and
`6.19e-9 s`; the moderate/high-energy case differs by about `2.44e-6 rad` and
`-6.84e-8 s`. This focused comparison shows increased numerical sensitivity
near the separatrix; it is not a convergence proof. The coarser output grid
uses 601 samples and separately reports period-estimate changes. Solver
tolerance error, crossing interpolation/sampling error, and nonlinear-versus-
linear model discrepancy remain distinct quantities.

The existing plotting helper accepts a single pendulum or oscillator
trajectory and does not accept paired Phase 26 model traces or the conditional
grid. The report therefore uses explicit tables rather than adding plotting
infrastructure that could obscure initial displacement versus turning
amplitude.

## Limitations and conclusion

The calculations establish consequences of the ideal equations for the
specified parameter and initial-condition grid. They do not validate a real
pendulum, include damping, describe rotations with a libration formula, or
treat the separatrix as a periodic orbit. Sampled extrema and crossing-based
period estimates retain finite-grid error. Near-separatrix periods become
long and the classification is sensitive to floating-point energy accuracy.
The report creates no duplicate registry entry and does not rewrite Phase 25
results.

No physical equivalence, financial application, trading signal, predictive
market model, or risk-management method is established. A next research
question is how these conditional limits change under a measured damping law,
using an appropriate damped nonlinear period or decay reference rather than
the undamped elliptic formula.

## Reproduction

Run from the project root with the existing virtual environment:

```powershell
.venv\Scripts\python.exe -c "from newton_lab.pendulum_initial_velocity import run_initial_velocity_research, render_initial_velocity_report; print(render_initial_velocity_report(run_initial_velocity_research()))"
```
