# Phase 63 — Research Direction Selection and Novelty Audit

## Status and scope

**This document is a repository-based research proposal, not an executed study. No new scientific results were generated, and the proposed question has not been approved for implementation.** “Distinct” below means distinct from the research recorded in Newton Lab; this audit makes no claim about novelty in the wider scientific literature.

The audit used existing reports, protocols, machine-readable records, manifested outputs, code, and documentation. It ran no simulations, parameter sweeps, scientific evaluations, or tests. No Git operations, external services, external data, or dependency changes were used.

## Research record audit

| Thread | Central question and structure | Evidence and recorded conclusion | Status; follow-ups that would repeat it |
|---|---|---|---|
| **[Phase 48: pinned-path convergence](phase_48_pinned_path_convergence.md)** | Does a uniform, fixed-end, nearest-neighbor path have discrete modal decay rates that approach fixed-end heat rates at second order? | The discrete rate is `4 h^(-2) sin^2(n pi h/2)`. For modes 1–3, the finest-grid fitted orders were 1.99942, 1.99767, and 1.99476; numerical fitted-rate error against the exact discrete reference was at most 4.14e-11 across 60 integrations. This is a deliberately shared stencil/boundary construction, not evidence about arbitrary graphs. | Complete. More path sizes or modes under the same uniform stencil would extend the convergence sweep, not test a distinct mechanism. |
| **[Phase 49: critical slowing down](phase_49_critical_slowing_down.md)** | Do distinct heat and nonlinear-pendulum models show their own predicted local recovery-rate laws near a zero-dissipation boundary? | Heat mode and underdamped pendulum each showed an exponent-one trend over five declared distances; fitted exponents were 1.000000 and 1.000029. Heat used an exact finite-mode evaluator; pendulum rates came from numerical peaks. There was no measured target or shared physical mechanism claim. | Complete. Another damping-distance sweep or a third model selected only to repeat a recovery exponent would be adjacent work. |
| **[Phase 50: fold cross-domain benchmark](phase_50_cross_domain_validation.md)** | Does the fold's square-root local rate outperform an exponent-one predictor on a synthetic control-style target? | Under the declared mapping and known source/target equations, the fold predictor's held-out mean absolute log error was 0.00401 versus 0.91840 for the baseline. This tests a synthetic mapping whose structure is specified; it does not qualify a real controller or plant. | Complete, synthetic. Repeating the same fold-to-scalar-target mapping with nearby parameters would add little. Phase 52 later found no qualified engineering target. |
| **[Phase 53](phase_53_recovery_scaling_identifiability.md) / [54](phase_54_exponent_identifiability_ablation.md) / [55](phase_55_sub5_percent_perturbation.md): recovery-scaling identifiability and robustness** | Can fold, transcritical, and pitchfork classes and exponents be distinguished from finite noisy trajectories, and what observation factors affect the estimates? | Phase 53's primary-condition classification criteria passed while all three model-specific exponent criteria failed. Phase 54 localized sensitivity to its tested information, amplitude, noise, and window factors. Phase 55 found evidence more consistent with a substantial low-SNR contribution plus a smaller deterministic finite-amplitude/model-fit signature; it did not identify unique causal shares, and no sub-5% exponent criterion passed. These are synthetic, equation-informed tests. | Complete. Further amplitude/noise/window sweeps on these same generators repeat the thread unless a genuinely new estimation mechanism is tested. |
| **[Phase 57: noise dose and estimator form](phase_57_noise_dose_estimator_comparison.md)** | How do noise dose and observation-only estimator form interact for the existing normal forms? | Effects depended on model and amplitude; no estimator was a universal winner. Rate estimation and exponent identification remained separate, and only three nonzero-noise exponent conditions passed their criterion. Exact-form fitting was a model-informed diagnostic. | Complete. Phase 58 explicitly closed the adjacent synthetic recovery-scaling thread. Re-running the dose grid or swapping a nearby estimator would reopen it without a distinct question. |
| **[Phase 58: closure](phase_58_research_thread_closure.md)** | What remains supported after consolidating the recovery-rate evidence, and what should follow? | It concluded that the synthetic scaling sequence was bounded, no independent target was qualified, and further adjacent variants had diminishing value. It recommended consolidation/pause rather than forcing another phase. | Closure record; not an active experiment. Repeating the consolidation or reopening the same synthetic thread is not recommended. |
| **[Phase 60](phase_60_kuramoto_synchronization.md) / [61](phase_61_kuramoto_transition_finite_size.md): Kuramoto synchronization and finite size** | How do frequency distribution, finite population, and resolution affect synchronization and the operational transition estimate? | Phase 60's 352 integrations showed finite-run coherence increasing over the declared sweeps but mixed transition estimates, including a boundary estimate. Phase 61's 592 integrations found high-size intervals that included the continuum value but were too wide under frozen criteria; high-size point estimates straddled the reference. Tighter solver settings shifted one restricted-grid estimate by zero grid intervals but did not resolve seed uncertainty. The conclusion is inconclusive. | **Closed as inconclusive under the declared protocol.** The Atlas and reports explicitly prohibit simply adding sizes, seeds, or coupling points to force a conclusion. |
| **[Phase 62: pitchfork path dependence](phase_62_pitchfork_path_dependence.md)** | Does an up/down control ramp produce a finite-rate response loop, and how do ramp duration and initial branch selection affect it? | The model has a single-valued stable equilibrium magnitude and thus no quasi-static magnitude hysteresis. The frozen deterministic protocol observed loop areas 0.00364812 and 0.02491512 for leg durations 20 and 80. Finite-rate path difference and initial-sign branch selection were supported within tested conditions; slower-ramp reduction was inconclusive because the slower case had a larger loop. | Complete bounded pilot. More ramp durations or nearby μ ranges in this same protocol would be an adjacent sweep, not a new mechanism. |

### Related prior work and infrastructure boundary

The proposed question is not wholly without context. [Phase 22](phase_22_modal_dynamics.md) and [Phase 24](phase_24_mode_identification.md) cover modal dynamics and mode identification; [Phase 28](phase_28_oscillator_pendulum_validation.md) compares oscillator and pendulum trajectories; the [heat/diffusion case study](../heat_diffusion_smoothing_case_study.md) covers heat-mode attenuation and smoothing; [Phase 26](phase_26_pendulum_initial_velocity.md) uses energy conservation for an ideal nonlinear pendulum. These studies concern modal decay/estimation, model comparison, diffusion, or conservative single-system motion. The audited record contains no study of **non-normal state-space operators, matched-spectrum systems, propagator singular values, or maximum finite-time norm amplification**.

Existing capabilities include typed ODE solver settings in `src/newton_lab/simulation.py` (`ODESolverConfiguration`), ODE integrations with SciPy `solve_ivp` in model-specific modules such as `src/newton_lab/dynamics.py` and `src/newton_lab/pendulum.py`, reusable result/experiment contracts in `simulation.py` and `experiments.py`, and established synthetic-test/report/manifest conventions. `run_experiment`'s default adapters support the existing oscillator, pendulum, algebraic, and BVP families; it is not a generic arbitrary-matrix adapter. No current model or analysis helper for transient growth/non-normality was identified. NumPy and SciPy are already declared dependencies; their existing matrix/ODE capabilities are sufficient for a narrowly scoped implementation if later approved. A general linear-systems framework is not required.

The Atlas continues to record Phase 32 as `BLOCKED_AUTHORIZATION`, Phase 52 as `NO-GO`, the Kuramoto thread as `CLOSED_INCONCLUSIVE_UNDER_DECLARED_PROTOCOL`, and Phase 62 as complete. The Phase 60, 61, and 62 protocols and outputs were inspected. Static SHA-256 checks confirmed 8/8 Phase 53, 7/7 Phase 54, 10/10 Phase 55, 9/9 Phase 57, 10/10 Phase 60, 9/9 Phase 61, and 8/8 Phase 62 manifest entries matched their recorded hashes.

## Candidate directions

### Candidate A — Finite-time amplification in stable non-normal dynamics

**Question.** For a declared state norm and finite horizon, when does eigenvalue-only stability information fail to characterize the maximum transient amplification of an asymptotically stable linear system, and how does eigenvector geometry account for the difference?

**Mechanism and distinction.** A non-normal matrix can have eigenvalues strictly in the left half-plane yet exhibit initial growth in state norm before eventual decay. This is not recovery-rate scaling near a bifurcation, collective synchronization, or a parameter-ramp loop. It asks how the propagator acts at finite time, a distinct stability concept absent from the audited research record.

**Existing capability and minimum addition.** Use the existing `ODESolverConfiguration` and SciPy `solve_ivp` conventions for trajectories, with SciPy's matrix exponential or a closed-form propagator as an independent reference. A small model-specific module could evaluate a two-state family
\[
A_\kappa=\begin{pmatrix}-a&\kappa\\0&-b\end{pmatrix},\quad a,b>0,
\]
against the normal control `A0 = diag(-a, -b)`. The eigenvalues remain `(-a, -b)` as κ varies. The primary finite-time gain could be `G(t) = ||exp(A_kappa t)||_2`, with peak gain and its time over a frozen horizon. No new general solver or dependency is proposed. The precise parameter envelope and whether the future work is analytic-only or includes trajectories remain decisions for a protocol, not this proposal.

**Smallest meaningful analysis.** Derive the propagator for the two-state family, compare peak gain against the matched-spectrum normal control, and independently verify selected state trajectories and gain values with `solve_ivp`. Include a unit-initial-state optimization reference from the leading right singular vector; check selected fixed initial directions as secondary diagnostics. Record both asymptotic eigenvalue decay and finite-time peak gain rather than collapsing them into one “stability” number.

**Hypotheses and outcomes.** (1) If matched-eigenvalue members with nonzero κ have peak gains distinguishable from the normal control under a predeclared norm/horizon, eigenvalues alone are insufficient for this finite-time metric; if not over the frozen envelope, the proposed envelope did not demonstrate that limitation. (2) If `G(t) > 1` occurs for a predeclared member despite stable eigenvalues and is reproduced by the independent numerical trajectories, transient amplification is observed under those coordinates and conditions; if no member exceeds unity, report that negative result without widening the range after inspection. (3) If analytic and numerical references disagree beyond a tolerance justified in the protocol, classify the numerical comparison as unresolved and repair only a diagnosed numerical issue before interpreting physical conclusions.

**Alternatives and risks.** A norm is coordinate-dependent: state scaling can create or hide Euclidean-norm amplification. Coordinates must therefore be explicitly dimensionless and frozen, or a justified energy/weighted norm must be derived and used. A finite horizon can miss a later peak; a coarse sampling grid can miss the true maximum; near-defective eigenvectors can make results sensitive to floating-point precision. This is a mathematical benchmark, not evidence that a particular real system is non-normal or that transient amplification causes damage.

**Reusable knowledge.** A qualified separation between asymptotic eigenvalue stability and finite-time propagator gain, plus an independent analytic/numerical verification pattern for linear systems.

### Candidate B — Energy exchange in coupled oscillators near resonance

**Question.** In a conservative two-oscillator system with fixed total energy, how do coupling and natural-frequency detuning govern the maximum energy transferred from one oscillator to the other?

This would add a coupled-system mechanism and quantify energy exchange rather than estimate modal decay. The existing oscillator solver can inform conventions and SciPy ODE integration can be reused, but no coupled mechanical model or energy-transfer analysis API was found; a new two-degree-of-freedom model and energy accounting would be needed. A bounded comparison could use the uncoupled, equal-frequency reference and predeclared detunings, verify total-energy conservation, and compare numerical exchange with a normal-mode solution. Risks include confusing beating with dissipation, interpreting coordinate energy as separately conserved when coupling energy exists, and overlap with Phase 22 modal dynamics and Phase 26 energy conservation. It is scientifically feasible but has a closer relationship to existing oscillation/modal work than Candidate A.

### Candidate C — Onset of spatial pattern formation in a reaction-diffusion system

**Question.** Under what conditions can diffusion destabilize a homogeneous equilibrium to spatial modes in a two-component reaction-diffusion model?

This would introduce a diffusion-driven instability rather than heat attenuation or synchronization. It is distinct in mechanism and could yield a useful dispersion-relation/grid-validation result. However, the current diffusion work uses exact finite heat modes and a steady BVP; no general time-dependent PDE or reaction-diffusion solver was identified. A faithful bounded study would need a new spatial discretization, nonlinear reaction terms, boundary handling, and grid/stability validation. That exceeds the desired narrow reuse of current capability and would confound the physical question with substantial numerical infrastructure. Defer until a specific PDE capability gap is assessed and approved.

## Qualitative feasibility and value comparison

| Candidate | Distinctness | Scientific clarity | Infrastructure fit | Interpretability | Boundedness | Knowledge gain |
|---|---|---|---|---|---|---|
| A — non-normal transient gain | High within this repository: finite-time operator growth is not the subject of Phases 48–62. | High: matched eigenvalues, declared norm, propagator gain, and normal baseline define a precise comparison. | Good: NumPy/SciPy and existing ODE settings suffice; a small matrix-specific analysis module is needed. | Good if coordinate scaling, norm, horizon, analytic reference, and grid error are handled explicitly. | High: two states, a fixed matrix family, one normal control, and a frozen parameter envelope. | Adds the reusable distinction between asymptotic stability and finite-time amplification. |
| B — coupled-oscillator energy exchange | Moderate-high, but adjacent to existing modal and pendulum-energy work. | High with a clearly defined coupled energy and detuning. | Moderate: no coupled model/energy API exists; one new model is needed. | Moderate-good; energy partition must include the coupling term. | Good for two degrees of freedom, though resonance, detuning, and energy bookkeeping widen the protocol. | Adds a coupled mode-transfer example, but with a closer relation to established oscillator work. |
| C — reaction-diffusion pattern onset | High. | High at the linear dispersion-relation level. | Low: no transient PDE/reaction-diffusion solver was found. | Lower until space/time convergence and boundary handling are validated. | Low-moderate because the numerical PDE method would itself be a major workstream. | Potentially useful, but requires infrastructure beyond a bounded next study. |

These are qualitative judgments, not scores. Candidate A is selected not just because it is easiest: it addresses a distinct stability mechanism with a precise matched-spectrum control and adds a reusable analytical distinction. Candidate B remains credible but is closer to existing modal/oscillator topics. Candidate C should wait for explicit PDE infrastructure planning.

## Proposed research brief — not approved for implementation

### Working title

**Finite-Time Amplification in Asymptotically Stable Non-Normal Linear Systems**

### Central question

For a frozen two-state, asymptotically stable matrix family and an explicitly declared state norm, does matching the eigenvalues leave materially different finite-horizon peak amplification as non-normal coupling changes, and can propagator analysis predict that difference accurately?

### Motivation and prior-work distinction

Newton Lab has tested asymptotic modal decay, recovery rates, finite noisy rate estimation, synchronization, and finite-rate path loops. Those outcomes do not measure the largest finite-time response of a stable non-normal propagator. This proposal would expand the stability vocabulary from eigenvalue decay to finite-time gain without reopening the closed Kuramoto or recovery-scaling threads. It would not establish that the structure applies to an engineering or financial system.

### Model and method

Use (\dot z=A_\kappa z), with `A_kappa` as Candidate A and fixed positive `a != b`; the normal baseline has the same eigenvalues and κ=0. The unequal-eigenvalue case avoids a degeneracy in the simple closed-form expression; whether to include the repeated-eigenvalue limiting case is an explicit unresolved scope choice, not a post-result extension. Freeze dimensionless coordinates, Euclidean or independently justified weighted norm, coupling values, horizon, output resolution, and initial-condition treatment before evaluation. Compare the exact matrix exponential/analytic propagator with `solve_ivp` at declared tolerances, and perform a time-grid refinement check for maxima. Include eigenvalues as the asymptotic reference and `A0` as the matched-spectrum normal baseline.

### Hypotheses

1. **Matched-spectrum distinction:** peak finite-time gain varies with the frozen non-normal coupling values despite invariant stable eigenvalues. A null/weak contrast would show that the selected envelope does not resolve a useful difference.
2. **Amplification condition:** at least one predeclared family member has `G_max > 1`; absence of such a member is a meaningful negative result for that envelope.
3. **Numerical agreement:** independently computed trajectories recover analytic propagator gains within a prespecified numerical tolerance; failure means the numerical evidence is unresolved, not that the matrix theory failed.

The exact parameter envelope, norm, meaningful-effect threshold, and numerical tolerances must be justified and frozen in a protocol before any evaluation. No numerical threshold is set here.

### Experimental outline

The primary estimands are peak gain (G_{max}=\max_{t\in[0,T]}\|e^{A_\kappa t}\|\), time of peak, and the long-time decay rate. Compare each κ with the normal κ=0 baseline having identical eigenvalues. The leading right singular vector supplies the unit initial condition attaining the instantaneous operator norm; add a small, predeclared set of coordinate-axis initial states to show direction dependence. For selected conditions, compare matrix-exponential trajectories with independently integrated `solve_ivp` trajectories. Controls are the normal matrix, zero initial state as a solver check only, and analytic norm/time-grid refinement. Report absolute and relative discrepancies, the declared norm and units/scaling, solver settings, and all attempted conditions. A full random-seed or held-out design is not necessary for this deterministic linear benchmark, but parameter choices and criteria must be frozen before calculating the results.

### Evidence and limitations

The study could establish a mathematical and numerical result for the declared matrix family: identical stable eigenvalues need not summarize finite-time gain in the selected norm, and a propagator-based metric may capture that difference. It cannot establish physical validity, robustness to model uncertainty, nonlinear saturation, stochastic response, a universal transient-growth law, or cross-domain utility. Since the benchmark is analytically specified, it is a method-capability and structural demonstration, not discovery of a law from data.

### Stopping rules

- **Complete:** all predeclared matrices and diagnostics are recorded; analytic references and numerical integrations agree within the frozen tolerance; the normal baseline and norm are explicit; the report states either the observed matched-spectrum distinction or its absence and its scope.
- **Inconclusive:** the peak is not resolved by the time grid, analytic/numerical disagreement exceeds tolerance, or a norm/coordinate choice changes interpretation without an independently justified metric.
- **Feasibility stop:** a defensible state normalization or energy metric cannot be specified, or the available matrix-exponential API proves unavailable; report the gap without adding a broad framework.
- **No parameter chasing:** do not extend coupling, horizon, matrix dimension, or norm after viewing outcomes merely to obtain `G_max > 1` or a desired contrast. Any materially different system requires a new rationale and approval.

### Reuse and future integration

Likely reuse: `ODESolverConfiguration`, NumPy/SciPy, SciPy `solve_ivp`, numerical-comparison conventions, and existing test/report/manifest practices. If later approved, keep implementation narrow under `src/newton_lab/`, focused synthetic tests under `tests/`, the protocol and deterministic results under `reports/phase_63_non_normal_transient_growth/`, and the final report under `docs/research/`. No such artifacts beyond this planning document are created here.

## Atlas update and decision

The Phase 51 machine-readable Atlas already has an `opportunities` list with explicit statuses for completed, candidate, recommended, and blocked directions. The proposal has therefore been added there as **`PROPOSED_UNAPPROVED`**, separate from completed research areas and active work. The Kuramoto status and Phase 62 completion record are unchanged. This planning document does not authorize implementation; wait for explicit approval.

## Validation and unresolved design decisions

- The Phase 51 JSON parsed; local repository path references were checked. The selected Atlas entries still record Phase 32 `BLOCKED_AUTHORIZATION`, Phase 52 `NO-GO`, the closed Kuramoto thread, and completed Phase 62.
- Existing manifests checked during this audit: Phase 53 8/8, 54 7/7, 55 10/10, 57 9/9, 60 10/10, 61 9/9, and 62 8/8 entries matched.
- This phase changed documentation only. No tests were run, no experiment or simulation was run, and no scientific source, test, protocol, result, figure, or manifest was modified.
- Before any future protocol could be frozen, the norm/state normalization, coupling envelope, finite horizon, peak-resolution rule, and effect threshold need explicit justification. If there is no defensible norm for the proposed states, pause rather than treating Euclidean coordinates as physically meaningful.
