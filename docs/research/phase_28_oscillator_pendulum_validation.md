# Phase 28: Independent oscillator and pendulum numerical validation

## Research question and hypotheses

Question: when integrated independently, do the dimensionless linear
oscillator and linearized pendulum approach the same exact solution, and
how do discrepancies respond to solver settings, sampling, initial state
and horizon? Can finite-angle model error be distinguished from solver error?

Hypotheses: mapped analytical trajectories agree; independently integrated
linear trajectories approach that reference as accuracy is tightened;
cross-system error is consistent with per-solver errors; finite-angle
model error can exceed numerical error; grid and period effects are visible.

## Mathematical mapping

For the undamped oscillator, `m x'' + kx = 0`, define
`omega_o = sqrt(k/m)`, `tau = omega_o t`, and `x = A q`. Then
`d2q/dtau2 + q = 0`. Here `A != 0` is a displacement scale in metres
per unit q; it does not make displacement and angle physically identical.

For the full pendulum, `theta'' + (g/L) sin(theta) = 0`. With
`tau = sqrt(g/L) t`, its equation is `d2theta/dtau2 + sin(theta) = 0`;
small-angle linearization gives `d2theta/dtau2 + theta = 0`.
Both linear systems therefore have exact solution
`u(tau) = u0 cos(tau) + v0* sin(tau)`. This establishes mathematical
equivalence under scaling, not equivalence of physical systems.

Matched states are `x0=Aq0`, `xdot0=A omega_o v0*`, `theta0=q0`,
and `thetadot0=omega_p v0*`, where `omega_p=sqrt(g/L)`. Defaults are
`m=2 kg`, `k=8 N/m`, `L=1.5 m`, `g=6 m/s^2`, `A=0.25 m`; both
frequencies are 2 rad/s. Unequal frequencies use each system's own
`t=tau/omega`, while comparisons use shared dimensionless tau.

## Experimental design

Ordered initial states `(q0,v0*)`: ((0.0, 0.0), (0.01, 0.0), (0.2, 0.0), (0.2, 0.3), (0.2, -0.3), (0.8, 0.3)).
Horizons tau: (6.283185307179586, 25.132741228718345, 50.26548245743669); output counts: (101, 501, 2001).
Methods: ('RK45', 'DOP853'); tolerance pairs (rtol, atol): ((1e-05, 1e-07), (1e-08, 1e-10), (1e-11, 1e-13)).
Sensitivity changes one factor at a time around DOP853, 1e-8/1e-10,
501 points and an 8*pi horizon.

The oscillator uses the existing oscillator API with zero damping. The
linearized pendulum has a separately constructed state, RHS and solve_ivp
call because the existing pendulum API only solves the full-sine model.
The full-sine comparison uses the existing pendulum solver. All use SciPy,
so this is not independent numerical-library validation.

## Results

Errors are maximum absolute dimensionless displacement on shared requested
tau samples. Normalized error divides by `max_i |u_exact(tau_i)|`; it is
undefined for the zero state, whose absolute errors remain reported.
Cross error compares the two independent linear integrations. Model error
compares full-sine pendulum samples with the exact linear solution.
The nonlinear-versus-linear numerical discrepancy is also listed.
Full-sine model error includes both physical model discrepancy and
numerical error from integrating the nonlinear trajectory.

| Case | q0,v0* | method | rtol/atol | tau/points | E osc | E lin pend | E cross | E model | E nonlinear vs lin num | linear period est/ref/error tau | linear-pendulum period est/ref/error tau | nonlinear period est/ref/error s |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| initial_1 | 0,0 | DOP853 | 1e-08/1e-10 | 8pi/501 | 0 | 0 | 0 | 0 | 0 | n/a/6.283/n/a | n/a/6.283/n/a | n/a/n/a/n/a |
| initial_2 | 0.01,0 | DOP853 | 1e-08/1e-10 | 8pi/501 | 9.956e-10 | 3.937e-10 | 7.823e-10 | 1.473e-06 | 1.473e-06 | 6.283/6.283/1.586e-08 | 6.283/6.283/5.623e-09 | 3.142/3.142/1.826e-09 |
| initial_3 | 0.2,0 | DOP853 | 1e-08/1e-10 | 8pi/501 | 5.421e-09 | 4.671e-09 | 6.207e-09 | 0.0118 | 0.0118 | 6.283/6.283/7.502e-09 | 6.283/6.283/5.341e-09 | 3.149/3.149/1.561e-07 |
| initial_4 | 0.2,0.3 | DOP853 | 1e-08/1e-10 | 8pi/501 | 8.477e-09 | 7.897e-09 | 8.48e-09 | 0.07466 | 0.07466 | 6.283/6.283/5.512e-09 | 6.283/6.283/5.407e-09 | 3.168/3.168/1.446e-07 |
| initial_5 | 0.2,-0.3 | DOP853 | 1e-08/1e-10 | 8pi/501 | 8.75e-09 | 8.081e-09 | 8.367e-09 | 0.06555 | 0.06555 | 6.283/6.283/3.375e-09 | 6.283/6.283/4.607e-09 | 3.168/3.168/1.947e-08 |
| initial_6 | 0.8,0.3 | DOP853 | 1e-08/1e-10 | 8pi/501 | 1.95e-08 | 1.787e-08 | 1.658e-08 | 0.9431 | 0.9431 | 6.283/6.283/5.548e-09 | 6.283/6.283/4.674e-09 | 3.294/3.294/2.069e-07 |

Run records retain exact, oscillator, linearized-pendulum and full-sine
samples on the common dimensionless output grid, alongside solver
configuration and function-evaluation counts.

## Sensitivity

Each sensitivity row is a recorded run, including its case and status.
Only the named factor changes within each one-factor series; table values
are measurements, not a convergence proof.

| Factor | Setting | Case | E osc | E lin pend | E cross | E model | normalized cross | period tau/error | status |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| tolerance | rtol=1e-05,atol=1e-07 | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| tolerance | rtol=1e-08,atol=1e-10 | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| tolerance | rtol=1e-11,atol=1e-13 | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| method | RK45 | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| method | DOP853 | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| output_grid | 101 points | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| output_grid | 501 points | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| output_grid | 2001 points | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| horizon | tau=6.28319 | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| horizon | tau=25.1327 | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| horizon | tau=50.2655 | initial_1 | 0 | 0 | 0 | 0 | n/a | n/a/n/a | succeeded |
| tolerance | rtol=1e-05,atol=1e-07 | initial_2 | 7.839e-07 | 4.167e-07 | 7.668e-07 | 1.788e-06 | 7.668e-05 | 6.283/2.179e-05 | succeeded |
| tolerance | rtol=1e-08,atol=1e-10 | initial_2 | 9.956e-10 | 3.937e-10 | 7.823e-10 | 1.473e-06 | 7.823e-08 | 6.283/1.586e-08 | succeeded |
| tolerance | rtol=1e-11,atol=1e-13 | initial_2 | 8.963e-13 | 4.152e-13 | 9.188e-13 | 1.473e-06 | 9.188e-11 | 6.283/1.977e-11 | succeeded |
| method | RK45 | initial_2 | 2.327e-09 | 8.409e-10 | 1.486e-09 | 1.473e-06 | 1.486e-07 | 6.283/1.237e-08 | succeeded |
| method | DOP853 | initial_2 | 9.956e-10 | 3.937e-10 | 7.823e-10 | 1.473e-06 | 7.823e-08 | 6.283/1.586e-08 | succeeded |
| output_grid | 101 points | initial_2 | 9.056e-10 | 3.838e-10 | 6.865e-10 | 1.473e-06 | 6.865e-08 | 6.283/1.847e-08 | succeeded |
| output_grid | 501 points | initial_2 | 9.956e-10 | 3.937e-10 | 7.823e-10 | 1.473e-06 | 7.823e-08 | 6.283/1.586e-08 | succeeded |
| output_grid | 2001 points | initial_2 | 9.956e-10 | 3.966e-10 | 8.156e-10 | 1.473e-06 | 8.156e-08 | 6.283/1.705e-08 | succeeded |
| horizon | tau=6.28319 | initial_2 | 8.521e-10 | 3.456e-10 | 8.156e-10 | 2.975e-07 | 8.156e-08 | n/a/n/a | succeeded |
| horizon | tau=25.1327 | initial_2 | 9.956e-10 | 3.937e-10 | 7.823e-10 | 1.473e-06 | 7.823e-08 | 6.283/1.586e-08 | succeeded |
| horizon | tau=50.2655 | initial_2 | 1.726e-09 | 7.168e-10 | 1.306e-09 | 3.041e-06 | 1.306e-07 | 6.283/3.424e-06 | succeeded |
| tolerance | rtol=1e-05,atol=1e-07 | initial_3 | 4.642e-06 | 4.289e-06 | 5.705e-06 | 0.0118 | 2.852e-05 | 6.283/7.945e-06 | succeeded |
| tolerance | rtol=1e-08,atol=1e-10 | initial_3 | 5.421e-09 | 4.671e-09 | 6.207e-09 | 0.0118 | 3.104e-08 | 6.283/7.502e-09 | succeeded |
| tolerance | rtol=1e-11,atol=1e-13 | initial_3 | 5.003e-12 | 4.777e-12 | 6.305e-12 | 0.0118 | 3.153e-11 | 6.283/5.858e-12 | succeeded |
| method | RK45 | initial_3 | 8.435e-09 | 6.75e-09 | 1.685e-09 | 0.0118 | 8.423e-09 | 6.283/1.628e-09 | succeeded |
| method | DOP853 | initial_3 | 5.421e-09 | 4.671e-09 | 6.207e-09 | 0.0118 | 3.104e-08 | 6.283/7.502e-09 | succeeded |
| output_grid | 101 points | initial_3 | 5.037e-09 | 4.6e-09 | 3.842e-09 | 0.0118 | 1.921e-08 | 6.283/4.458e-09 | succeeded |
| output_grid | 501 points | initial_3 | 5.421e-09 | 4.671e-09 | 6.207e-09 | 0.0118 | 3.104e-08 | 6.283/7.502e-09 | succeeded |
| output_grid | 2001 points | initial_3 | 5.451e-09 | 4.683e-09 | 6.552e-09 | 0.0118 | 3.276e-08 | 6.283/7.674e-09 | succeeded |
| horizon | tau=6.28319 | initial_3 | 4.245e-09 | 3.713e-09 | 6.552e-09 | 0.002383 | 3.276e-08 | n/a/n/a | succeeded |
| horizon | tau=25.1327 | initial_3 | 5.421e-09 | 4.671e-09 | 6.207e-09 | 0.0118 | 3.104e-08 | 6.283/7.502e-09 | succeeded |
| horizon | tau=50.2655 | initial_3 | 8.79e-09 | 7.738e-09 | 5.734e-09 | 0.02437 | 2.867e-08 | 6.283/3.41e-06 | succeeded |
| tolerance | rtol=1e-05,atol=1e-07 | initial_4 | 7.483e-06 | 6.915e-06 | 5.9e-06 | 0.07467 | 1.636e-05 | 6.283/1.557e-06 | succeeded |
| tolerance | rtol=1e-08,atol=1e-10 | initial_4 | 8.477e-09 | 7.897e-09 | 8.48e-09 | 0.07466 | 2.352e-08 | 6.283/5.512e-09 | succeeded |
| tolerance | rtol=1e-11,atol=1e-13 | initial_4 | 8.42e-12 | 8.04e-12 | 1.008e-11 | 0.07466 | 2.795e-11 | 6.283/5.206e-12 | succeeded |
| method | RK45 | initial_4 | 1.254e-08 | 1.094e-08 | 1.711e-09 | 0.07466 | 4.746e-09 | 6.283/1.46e-09 | succeeded |
| method | DOP853 | initial_4 | 8.477e-09 | 7.897e-09 | 8.48e-09 | 0.07466 | 2.352e-08 | 6.283/5.512e-09 | succeeded |
| output_grid | 101 points | initial_4 | 7.624e-09 | 7.191e-09 | 8.48e-09 | 0.07466 | 2.353e-08 | 6.283/5.914e-09 | succeeded |
| output_grid | 501 points | initial_4 | 8.477e-09 | 7.897e-09 | 8.48e-09 | 0.07466 | 2.352e-08 | 6.283/5.512e-09 | succeeded |
| output_grid | 2001 points | initial_4 | 8.936e-09 | 7.92e-09 | 8.652e-09 | 0.07468 | 2.4e-08 | 6.283/5.523e-09 | succeeded |
| horizon | tau=6.28319 | initial_4 | 5.853e-09 | 5.636e-09 | 3.966e-09 | 0.01863 | 1.1e-08 | n/a/n/a | succeeded |
| horizon | tau=25.1327 | initial_4 | 8.477e-09 | 7.897e-09 | 8.48e-09 | 0.07466 | 2.352e-08 | 6.283/5.512e-09 | succeeded |
| horizon | tau=50.2655 | initial_4 | 1.513e-08 | 1.392e-08 | 8.669e-09 | 0.1486 | 2.404e-08 | 6.283/3.832e-06 | succeeded |
| tolerance | rtol=1e-05,atol=1e-07 | initial_5 | 7.499e-06 | 7.017e-06 | 5.982e-06 | 0.06555 | 1.659e-05 | 6.283/3.842e-06 | succeeded |
| tolerance | rtol=1e-08,atol=1e-10 | initial_5 | 8.75e-09 | 8.081e-09 | 8.367e-09 | 0.06555 | 2.321e-08 | 6.283/3.375e-09 | succeeded |
| tolerance | rtol=1e-11,atol=1e-13 | initial_5 | 7.48e-12 | 8.232e-12 | 9.235e-12 | 0.06555 | 2.561e-11 | 6.283/5.22e-12 | succeeded |
| method | RK45 | initial_5 | 1.312e-08 | 1.138e-08 | 1.777e-09 | 0.06555 | 4.93e-09 | 6.283/1.456e-09 | succeeded |
| method | DOP853 | initial_5 | 8.75e-09 | 8.081e-09 | 8.367e-09 | 0.06555 | 2.321e-08 | 6.283/3.375e-09 | succeeded |
| output_grid | 101 points | initial_5 | 7.845e-09 | 8.081e-09 | 7.996e-09 | 0.06545 | 2.218e-08 | 6.283/4.964e-09 | succeeded |
| output_grid | 501 points | initial_5 | 8.75e-09 | 8.081e-09 | 8.367e-09 | 0.06555 | 2.321e-08 | 6.283/3.375e-09 | succeeded |
| output_grid | 2001 points | initial_5 | 9.17e-09 | 8.166e-09 | 8.819e-09 | 0.06555 | 2.446e-08 | 6.283/3.296e-09 | succeeded |
| horizon | tau=6.28319 | initial_5 | 5.977e-09 | 5.772e-09 | 4.2e-09 | 0.01531 | 1.165e-08 | n/a/n/a | succeeded |
| horizon | tau=25.1327 | initial_5 | 8.75e-09 | 8.081e-09 | 8.367e-09 | 0.06555 | 2.321e-08 | 6.283/3.375e-09 | succeeded |
| horizon | tau=50.2655 | initial_5 | 1.442e-08 | 1.322e-08 | 8.302e-09 | 0.1395 | 2.303e-08 | 6.283/2.877e-06 | succeeded |
| tolerance | rtol=1e-05,atol=1e-07 | initial_6 | 1.793e-05 | 1.724e-05 | 7.436e-06 | 0.9431 | 8.704e-06 | 6.283/4.311e-06 | succeeded |
| tolerance | rtol=1e-08,atol=1e-10 | initial_6 | 1.95e-08 | 1.787e-08 | 1.658e-08 | 0.9431 | 1.94e-08 | 6.283/5.548e-09 | succeeded |
| tolerance | rtol=1e-11,atol=1e-13 | initial_6 | 1.946e-11 | 1.91e-11 | 2.111e-11 | 0.9431 | 2.471e-11 | 6.283/4.87e-12 | succeeded |
| method | RK45 | initial_6 | 2.569e-08 | 2.409e-08 | 2.371e-09 | 0.9431 | 2.775e-09 | 6.283/1.258e-09 | succeeded |
| method | DOP853 | initial_6 | 1.95e-08 | 1.787e-08 | 1.658e-08 | 0.9431 | 1.94e-08 | 6.283/5.548e-09 | succeeded |
| output_grid | 101 points | initial_6 | 1.749e-08 | 1.787e-08 | 1.658e-08 | 0.9383 | 1.94e-08 | 6.283/4.637e-09 | succeeded |
| output_grid | 501 points | initial_6 | 1.95e-08 | 1.787e-08 | 1.658e-08 | 0.9431 | 1.94e-08 | 6.283/5.548e-09 | succeeded |
| output_grid | 2001 points | initial_6 | 1.991e-08 | 1.834e-08 | 1.705e-08 | 0.9431 | 1.995e-08 | 6.283/5.547e-09 | succeeded |
| horizon | tau=6.28319 | initial_6 | 1.328e-08 | 1.302e-08 | 4.037e-09 | 0.2299 | 4.725e-09 | n/a/n/a | succeeded |
| horizon | tau=25.1327 | initial_6 | 1.95e-08 | 1.787e-08 | 1.658e-08 | 0.9431 | 1.94e-08 | 6.283/5.548e-09 | succeeded |
| horizon | tau=50.2655 | initial_6 | 3.305e-08 | 3.222e-08 | 1.932e-08 | 1.582 | 2.262e-08 | 6.283/1.808e-06 | succeeded |

Mean absolute cross discrepancy across the six initial states by setting:
- tolerance, rtol=1e-05,atol=1e-07: mean E_cross=4.298e-06
- tolerance, rtol=1e-08,atol=1e-10: mean E_cross=6.735e-09
- tolerance, rtol=1e-11,atol=1e-13: mean E_cross=7.942e-12
- method, RK45: mean E_cross=1.505e-09
- method, DOP853: mean E_cross=6.735e-09
- output_grid, 101 points: mean E_cross=6.263e-09
- output_grid, 501 points: mean E_cross=6.735e-09
- output_grid, 2001 points: mean E_cross=6.981e-09
- horizon, tau=6.28319: mean E_cross=3.262e-09
- horizon, tau=25.1327: mean E_cross=6.735e-09
- horizon, tau=50.2655: mean E_cross=7.223e-09
Across the six states, mean cross error decreases at each tighter
tolerance setting. The RK45 mean is lower than DOP853 in this grid;
this is method sensitivity, not evidence that one method is generally
more accurate. Output-grid means remain close and change slightly
non-monotonically. Mean cross error rises modestly with longer horizons.
At q0=0.8, v0*=0.3, full-sine model error is about 0.943 at 8*pi and
1.58 at 16*pi, far above the linear solver errors. These trends are
specific to the tested configurations and metrics.

## Error attribution and falsification

- The exact linear references coincide by derivation for matched states.
- The triangle inequality gives cross error <= oscillator analytical
  error + linear-pendulum analytical error on the common grid.
- Tighter tolerances and alternate methods expose numerical sensitivity;
  errors need not decrease monotonically at every finite setting.
- Coarse and dense outputs are evaluated at their own common requested
  samples. Differences mix dense-output evaluation and sampled maxima;
  they do not isolate or bound internal adaptive-step error.
- Signed-velocity symmetry, zero-state and parameter-rescaling checks
  passed in tests. Finite-angle deviation is a model discrepancy, not
  evidence of a solver failure.
- Periods use interpolated upward zero crossings. A period is unavailable
  unless enough same-direction crossings exist. The nonlinear
  reference is elliptic and applies to librations only.

## Assumptions and limitations

- Undamped linear oscillator; no forcing and constant positive m,k.
- Ideal pendulum; no damping or forcing, fixed length and uniform gravity.
- A is a nonzero displacement scale in metres per unit dimensionless q; it does not equate angle with displacement.
- Comparisons use common dimensionless times; physical time is tau/omega for each system.
- All integrations use SciPy solve_ivp; distinct RHS implementations do not mean independent numerical libraries.
- The finite state, tolerance, method, grid and horizon choices do not prove convergence or universal accuracy.
- Grid effects include solve_ivp dense output and crossing interpolation. They are distinct from adaptive-step tolerance, but do not bound discretization error.
- For the zero state, normalized error is undefined because analytical amplitude is zero; absolute errors remain reported.
- Period estimates require repeated same-direction crossings; short windows can be insufficient.
- Elliptic-period references apply only to librations, not rotations, separatrix states or equilibrium.

Observed runs: 72; failed integrations: 0.
The tested grid is finite, and all systems share SciPy solve_ivp.
Analytical equivalence is not physical equivalence; numerical agreement
does not establish universal accuracy, application performance or trading
utility.

## Conclusion and next question

The runs test two independently specified linear equations against one
exact dimensionless reference. They support the declared mapping only
for the tested conditions. The next question is whether validated bounds
can quantify the observed global error independently of SciPy solvers.
