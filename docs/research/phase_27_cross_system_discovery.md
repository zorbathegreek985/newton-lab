# Phase 27 — Cross-system structural discovery and falsification

## Executive summary

The undamped linear oscillator and small-angle pendulum have identical dimensionless equations under matched scaling and initial states. That is mathematical, not physical, equivalence. Heat modes and oscillator envelopes share conditional exponential forms but not spectra or observables. Phase25–26 reject a universal small-angle threshold. The steady heat BVP is the transient model's equilibrium only under matched assumptions. Estimation results are methods evidence, not new physical laws. No application or financial performance is tested.

## Questions

H1: dimensionless oscillator/pendulum correspondence. H2: modal rates and inference. H3: nonlinear approximation limits. H4: steady versus transient heat. H5: estimation reliability.

## Source inventory

Phases 18–26 source modules, tests, and relevant checked-in reports were inspected. Phase18/19 models are available from source APIs; no Phase21 serialized results were found. The expected the expected Phase21 JSON artifact directory `.newton_lab/artifacts` directory is absent. Phase20 integration, registry, workflow, and Phase21 schemas (`knowledge/integration.py`, `knowledge/registry.py`, `knowledge/workflow.py`, `knowledge/artifacts.py`) were inspected; this module does not mutate them. Phase22–26 computed values below come from their checked-in generated reports. The focused H1 solver check is calculated here from existing solvers.

## Comparison matrix

| Phase / role | Model | Assumptions, states and units | Parameters / scales | Observables / evidence | Limits |
|---|---|---|---|---|---|
| Phase 18: Damped linear oscillator (physical model; `phase18_damping_transfer`) | m x¨ + c x˙ + k x = 0 | linear spring; linear viscous damping; constant parameters; no force; x: m; x˙: m/s | m: kg; c: kg/s; k: N/m; ωn=√(k/m); underdamped envelope time 1/γ, γ=c/(2m); ζ=c/(2√(mk)); τ=ωn t | displacement; velocity; sampled envelope. conditional structural analogy: Source report: m=1 kg, k=4 N/m, x0=0.1 m, v0=0; c=0,0.1,0.2 kg/s. | Envelope rate applies only when underdamped; zero damping does not decay. Phase 21 source JSON is absent; values are report-derived. No empirical application data. |
| Phase 26: Nonlinear pendulum (physical model; `phase26_nonlinear_pendulum`) | θ¨ + (g/L) sin θ = 0 | ideal point mass; fixed rigid length; uniform gravity; undamped, unforced; θ: rad; θ˙: rad/s | m: kg; L: m; g: m/s²; √(L/g); libration period also depends on energy; τ=√(g/L)t; energy h=θ˙²L/(2g)+1−cosθ | turning amplitude; period; angle error; state-coordinate phase. numerically supported relationship: Phase 26 joint initial-angle/velocity grid; report covers windows to 10 linear periods; rotations and separatrix excluded from libration-period formula. | Thresholds depend on state, observable, tolerance and window. Ideal equations only; no measurements. Phase 21 JSON record absent. |
| Phase 25: Small-angle pendulum (physical model; `phase25_small_angle`) | θ¨ + (g/L) θ = 0 | sin θ ≈ θ; fixed length; uniform gravity; undamped, unforced; θ: rad; θ˙: rad/s | L: m; g: m/s²; ω0=√(g/L); T0=2π√(L/g); τ=√(g/L)t | trajectory; period; phase. established mathematical relationship: Phase 25/26 compare DOP853 finite-grid trajectories with analytical references. | Linearization error accumulates and depends on metric and horizon. Linear model is not the full nonlinear pendulum. |
| Phase 22: Oscillator and heat modal dynamics (physical model; `phase22_modal_dynamics`) | Heat: a˙n=−α(nπ/L)²an; oscillator r=−γ±iωd | linear time-invariant equations; fixed-end heat perturbation; underdamped oscillator; heat coefficient: K; oscillator displacement: m | α: m²/s; L: m; m: kg; c: kg/s; k: N/m; heat τn=1/λn; oscillator envelope τ=1/γ and period=2π/ωd; heat λn t; oscillator ζ | mode coefficient; displacement envelope; sampled peaks. conditional structural analogy: Phase 22 uses Phase 18/19 APIs and reports analytical rates and sampled estimates. | Shared exponential factor does not imply same spectrum, observable or mechanism. Generated report is the available result record. |
| Phase 23: Single-rate estimation (inference procedure; `phase23_modal_estimation`) | Fit y(t)=A exp(−λt); estimator assumption, not physical law | synthetic truth; finite cadence/window; declared noise and method; signal: declared units; time: s | A: signal units; λ: s⁻¹; σ: signal units; window and cadence relative to 1/λ; λT and λΔt | rate estimate; fit error; status. inconclusive: Phase 23 reports cadence, duration, noise, mixed heat modes and oscillator envelope. | Fits can fail at noise floors, sparse cadence or short windows. Report-derived synthetic results do not guarantee behavior on measured data. |
| Phase 24: Single/two-mode identification (inference procedure; `phase24_mode_identification`) | Fit ΣAi exp(−λi t); compare models under declared criteria | synthetic positive modes; known-noise BIC where declared; multiple starts; observation: signal units; time: s | Ai: signal units; λi: s⁻¹; σ: signal units; window and rate separation jointly affect recovery; Δ=|λ2−λ1|/max(λ1,λ2) | model selection; rate recovery; residual; conditioning. numerically supported relationship: Phase 24 tests noise repetitions, irregular sampling, weak and close modes. | Extra parameters lower residual; low residual alone is not rate recovery. Finite synthetic design does not generalize to untested observations. |
| Phase 19: Steady heat conduction (physical model; `phase19_steady_heat`) | d/dx(k dT/dx)=0; constant k gives Tₓₓ=0 | steady; homogeneous; constant conductivity; no generation; fixed endpoints; T: K; x: m | k: W/(m K); endpoint temperatures: K; No time variable in this BVP.; not applicable | spatial profile; boundary/differential residual. established mathematical relationship: Report: 21 mesh points, accepted solve, zero boundary residual and 1.35e-15 differential residual. | BVP acceptance says nothing alone about transient rates. No heat equation registry identity; no Phase 21 JSON artifact. |
| Phase 19: Transient heat diffusion (physical model; `phase19_transient_heat`) | Tt=αTxx; u=T−Tbase obeys ut=αuxx | 1D homogeneous medium; constant α; no source; fixed endpoints; finite sine perturbation; T,u: K; x: m; t: s | α: m²/s; L: m; τn=1/[α(nπ/L)²]; λn=α(nπ/L)²; λnt | modal amplitudes; grid temperatures; analytic RMS gradient. numerically supported relationship: Phase 19: L=1 m, α=1e-4 m²/s, modes 1/4, four times and 101 positions. | Exact finite-series case, not a general PDE solver or denoising proof. No independent transient PDE solver cross-check. |

## Findings and falsification

### H1_dimensionless_structure: established mathematical relationship

Undamped linear oscillator and small-angle pendulum share X″+X=0 after scaling.

**Support:** Identical dimensionless equations plus matched normalized initial states give identical exact trajectories.

**Counterexample:** Damping, unmatched initial data or the retained finite-angle sine term breaks that correspondence.

**Assumptions:** c=0; τ=ωn t=√(g/L)t; matched coordinate and derivative scales

**Tested range:** Exact for all mapped states; numerical check uses θ0=0.2 rad, θ˙0=0.3 rad/s, τ=0…4π.

**Evidence:**

- Numerical oscillator vs exact linear pendulum normalized max error: 1.75e-10.
- Full-sine model deviation on same finite grid: 0.103.

**Untested:**

- Damped oscillator mapping and other solver implementations.

**Sources:** `src/newton_lab/dynamics.py`, `src/newton_lab/pendulum.py`, `docs/research/phase_25_pendulum_approximation.md`

### H2_modal_decay: conditional structural analogy

Heat modes and underdamped oscillator envelopes have exponential factors but different rates and observables.

**Support:** For fixed-end heat, λn=α(nπ/L)²; for the oscillator envelope, γ=c/(2m).

**Counterexample:** The Phase 22 signed two-mode sum crosses zero although each component decays, so aggregate magnitude need not decay monotonically.

**Assumptions:** linear time-invariant equations; modal amplitudes and observable specified

**Tested range:** Heat α=1e-4 m²/s,L=1 m,n=1,4; oscillator m=1 kg,c=0.05–0.2 kg/s,k=4 N/m.

**Evidence:**

- lambda1=0.00098696044/s, lambda4=0.015791367/s; gamma=0.05/s at c=0.1 kg/s,m=1 kg.
- Phase 22 Phase19-derived signed n=1,2 aggregate at x=L/4: -0.64645 K at 0 s, -0.04300 K at 300 s, +0.03208 K at 400 s, +0.10196 K at 600 s.
- Phase 23/24 report findings show short/noisy/sparse windows and weak or close modes limit estimation.

**Untested:**

- Untested real observations and arbitrary sampling/noise regimes.

**Sources:** `src/newton_lab/modal_research.py`, `src/newton_lab/modal_estimation.py`, `src/newton_lab/modal_identification.py`, `docs/research/phase_22_modal_dynamics.md`, `docs/research/phase_23_modal_estimation.md`, `docs/research/phase_24_mode_identification.md`

### H3_linearization_limits: numerically supported relationship

Small-angle validity depends on initial state, observable, tolerance and horizon.

**Support:** Phase25 grid limits differ by metric; Phase26 results vary across initial velocity slices.

**Counterexample:** At 0.35 rad the period and 3-period angle criteria pass while the 5-period phase criterion fails, refuting one shared cutoff.

**Assumptions:** undamped ideal pendulum; declared windows and tolerances; libration

**Tested range:** Phase25 rest grid 0.01–2 rad; Phase26 two-dimensional state grid and windows to 10 periods.

**Evidence:**

- Phase25 reports 0.35 rad period/3-period-angle limits and 0.2 rad phase limit.
- At 0.5 rad, period relative error 1.56%, angle error 0.1359 rad at 3 periods, phase lag 0.4804 rad at 5 periods.
- Phase26 rest controls report same 0.35/0.2 rad limits; (0.5,+3 rad/s) phase difference is 5.11159 rad at 10 periods.

**Untested:**

- No universal threshold; untested states and real apparatus.

**Sources:** `docs/research/phase_25_pendulum_approximation.md`, `docs/research/phase_26_pendulum_initial_velocity.md`

### H3b_universal_threshold: rejected or contradicted candidate

One initial-angle threshold predicts validity across observables and initial angular velocities.

**Support:** A shared threshold would require period, trajectory and phase criteria to change consistently across velocity slices.

**Counterexample:** Phase25 and Phase26 report distinct phase, period and trajectory limits; nonzero signed velocity changes the conditional limits.

**Assumptions:** same ideal undamped model; same declared accuracy criteria; threshold applies across observables and initial states

**Tested range:** Phase25 release from rest over 0.01–2 rad; Phase26 includes velocity slices at 0 and ±1 rad/s.

**Evidence:**

- At rest, period/angle largest passing tested angle is 0.35 rad; phase limit is 0.2 rad.
- At ±1 rad/s, period limit is 0.2 rad, angle limit is 0.1 rad, and no tested angle passes the 5-period phase criterion.

**Untested:**

- This rejects a universal threshold only for the tested grids and criteria; denser grids and other observables remain untested.

**Sources:** `docs/research/phase_25_pendulum_approximation.md`, `docs/research/phase_26_pendulum_initial_velocity.md`

### H4_steady_transient_heat: established mathematical relationship

The fixed-end steady heat line is the equilibrium baseline of the documented transient equation.

**Support:** Subtracting Tbase leaves homogeneous diffusion; each finite sine mode decays as exp(-λn t).

**Counterexample:** Changing boundaries, adding a source, or varying material properties requires a new derivation.

**Assumptions:** same endpoints; constant properties; no source; homogeneous Dirichlet perturbation

**Tested range:** Phase19 finite modes, L=1 m, α=1e-4 m²/s, t=0…1200 s, n=1,4.

**Evidence:**

- Exact modal factors and analytic RMS gradient; temperatures sampled at 101 grid positions.

**Untested:**

- Independent transient PDE solver and arbitrary forcing remain untested.

**Sources:** `docs/heat_diffusion_smoothing_case_study.md`, `src/newton_lab/heat_diffusion_case_study.py`

### H5_inference_reliability: inconclusive

A low residual is insufficient evidence that a fit recovered physical modes.

**Support:** Phase23/24 compare fits with held-out synthetic truth under varied designs.

**Counterexample:** Wrong model selection, ill-conditioning or rates outside the truth tolerance refute residual-only validation.

**Assumptions:** truth separate from fit; noise model declared for BIC

**Tested range:** Finite synthetic designs, seeds and noise assumptions in Phase23/24 reports.

**Evidence:**

- Inference code is not a physical model; weak/close modes and sparse data affect identifiability.

**Untested:**

- Structured run artifacts are absent; this synthesis uses report records and inspected source.

**Sources:** `src/newton_lab/modal_estimation.py`, `src/newton_lab/modal_identification.py`, `docs/research/phase_23_modal_estimation.md`, `docs/research/phase_24_mode_identification.md`

## Dimensionless trajectory experiment

Mapping: m=1 kg, k=1 N/m, c=0; L=1 m, g=1 m/s²; x/(1 m) maps numerically to θ/(1 rad). Thus both natural frequencies are 1 s⁻¹ and τ=t. This is an assumed coordinate mapping, not shared physical meaning.
Initial state: theta0=0.2 rad, theta_dot0=0.3 rad/s; tau=0 to 12.5664; 2001 samples; DOP853 rtol=1e-10, atol=1e-12.
Oscillator vs exact linearized pendulum normalized max error: 1.7452699e-10. Solver error only.
Oscillator vs full-sine pendulum normalized max deviation: 0.10305996; scale=0.36234932. Finite-angle difference, not physical validation.

## Reproducibility and limitations

Run `build_cross_system_report()` to compute the focused ODE comparison and render the report. `render_cross_system_report()` itself performs no simulations. Report-only historical values are not presented as structured artifacts. Phase19 uses an analytic finite-mode solution without an independent transient PDE solver; Phase23/24 run-level output is not persisted. Finite grids do not establish universal behavior. Smoothing is not denoising or prediction.

## Next research questions

1. Independently integrate the linearized pendulum over multiple matched states and solver settings.
2. Cross-check Phase19 finite-mode diffusion against a refined transient PDE discretization.
3. Preserve structured Phase23/24 result records to query recovery and failure cases directly.

No real-world control, predictive, financial, or trading claim is established.
