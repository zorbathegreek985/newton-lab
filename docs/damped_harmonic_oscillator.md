# Damped harmonic oscillator

Newton Lab's first dynamical system models a free oscillator with a linear
spring and linear viscous damping:

\[
m\ddot{x} + c\dot{x} + kx = 0.
\]

Here `m` is mass in kilograms, `c` is the damping coefficient in kilograms
per second, `k` is spring stiffness in newtons per metre, `x` is displacement
in metres, and time is in seconds. Initial velocity is in metres per second.
The model assumes constant parameters, a linear restoring force, damping
proportional to velocity, and no external driving force. Many real oscillators
have nonlinearities, changing parameters, friction, or external inputs, so this
model should not be treated as universally applicable.

## Numerical method and use

`simulate_damped_oscillator` rewrites the equation as a two-component
first-order system for displacement and velocity and integrates it with
SciPy's adaptive `solve_ivp`. It returns uniformly sampled NumPy arrays for
time, displacement, and velocity, together with the physical parameters,
solver settings, and number of derivative evaluations. A solver-reported
failure raises `IntegrationError`; invalid physical values or settings raise
`ScientificValidationError`.

```python
from newton_lab.dynamics import simulate_damped_oscillator

result = simulate_damped_oscillator(
    mass_kg=1.0,
    damping_coefficient_kg_per_s=0.2,
    stiffness_n_per_m=4.0,
    initial_displacement_m=0.1,
    initial_velocity_m_per_s=0.0,
    duration_s=10.0,
    num_points=1001,
)

print(result.time_s.shape, result.displacement_m[-1])
```

The numerical trajectory approximates the solution. For example, when damping
is zero the exact solution is available analytically; numerical integration
will agree within finite solver error, not necessarily bit-for-bit. `rtol` and
`atol` are represented by `relative_tolerance` and `absolute_tolerance` in the
API. They govern the adaptive solver's local error control, not a guarantee on
total trajectory error. Results can also vary slightly with solver choice,
software versions, and parameter scaling. Compare against analytical solutions
where available and perform tolerance or convergence checks for research uses.

`num_points` sets only the requested output sampling grid; it does not set the
integrator's internal step size. Dense-looking output does not imply the same
accuracy as a tighter integration tolerance.
