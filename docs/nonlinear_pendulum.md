# Damped nonlinear pendulum

The pendulum model is

\[
\ddot{\theta} + \frac{b}{mL^2}\dot{\theta}
  + \frac{g}{L}\sin(\theta) = 0,
\]

or, equivalently, the torque balance

\[
mL^2\ddot{\theta} + b\dot{\theta}
  + mgL\sin(\theta) = 0.
\]

Angle `theta` is in radians, angular velocity in radians per second, and angular
acceleration in radians per second squared. Mass `m` is in kilograms, length
`L` in metres, damping coefficient `b` in kg·m²/s, and gravitational
acceleration `g` in m/s².

## Physical assumptions

This model treats the bob as a point mass on a massless rigid rod, with a fixed
pivot and constant length in a uniform gravitational field. Damping is a
linear viscous torque proportional to angular velocity, and there is no
external driving torque. Rod flexibility and air-current effects are not
modeled separately. These assumptions describe an idealized system and may not
fit every physical pendulum.

## Nonlinearity and the small-angle approximation

The sine term retains the nonlinear restoring torque at every angle. Replacing
`sin(theta)` with `theta` gives the linear small-angle model, which is useful
when angular excursions are small and measured in radians. It is an
approximation: its error grows with amplitude, so it should not be silently
substituted for the full equation at larger angles.

## Numerical integration and energy

`simulate_damped_pendulum` expresses the second-order equation as a first-order
system in angle and angular velocity and integrates it with SciPy's adaptive
`solve_ivp`. The returned arrays are uniformly sampled at the requested output
times. `relative_tolerance` and `absolute_tolerance` govern the solver's local
error estimates and do not guarantee a fixed global trajectory error. Tighten
tolerances or compare solutions across settings when accuracy matters.

The mechanical energy for this model is

\[
E = \frac{1}{2}mL^2\dot{\theta}^2
    + mgL\left(1-\cos(\theta)\right).
\]

For the ideal undamped equation, this energy is conserved. With the specified
viscous damping, the model gives \(dE/dt=-b\dot{\theta}^2\leq0\), so energy
does not increase in the continuous model. Numerical trajectories can show
small deviations from these identities due to integration error and finite
sampling; they are checks on a computation, not exact numerical guarantees.

## Example

```python
from newton_lab.pendulum import simulate_damped_pendulum

result = simulate_damped_pendulum(
    mass_kg=0.2,
    length_m=0.75,
    damping_coefficient_kg_m2_per_s=0.01,
    gravity_m_per_s2=9.81,
    initial_angle_rad=0.8,
    initial_angular_velocity_rad_per_s=0.0,
    duration_s=8.0,
    num_points=801,
)

print(result.time_s.shape, result.angle_rad[-1])
```
