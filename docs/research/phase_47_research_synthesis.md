# Phase 47 — Research Synthesis & Cross-Domain Hypothesis Design

## Executive summary

Phases 42–46 used reproducible synthetic time series to study descriptive diagnostics and abrupt-step rules. Across those studies, a statistic that detects a pattern did not reliably identify its generating process. Results depended on signal family, effect size, noise, threshold, and window. The studies support methodological conclusions within their declared synthetic designs; they do not establish real-world forecasting, financial value, causality, or trading profitability.

**Selected hypothesis:** For a uniformly spaced path of leader-pinned agents using local neighbor relaxation, the decay rate of each fixed low-frequency mode approaches the corresponding fixed-end heat-diffusion mode rate with second-order error as node spacing decreases.

This uses an existing physical heat model and modal-rate calculation, then tests a network algorithm with different state meaning and update mechanism. The project has not yet tested this transfer. Phase 32 remains `BLOCKED_AUTHORIZATION`; this proposal requires no market data.

## 1. Synthesis of Phases 42–46

| Phase | Question, methods, and design | Principal recorded findings | Limitations |
|---|---|---|---|
| 42 — [Synthetic time-series benchmark](phase_42_synthetic_timeseries_benchmark.md) | On 16 synthetic configurations, compare lag-one autocorrelation, a two-half mean statistic, a two-half variance ratio, and the existing causal exponential filter. Each configuration has 3 fixed seeds; 48 records total, 512 samples each. | The AR(1) flag appeared in 9/9 AR records, but also 9/9 damped oscillations, 3/9 trends, and 5/9 mean shifts. The two-half mean statistic flagged 9/9 mean shifts, 8/9 trends, and 4/9 AR records. The variance ratio flagged 5/9 variance shifts and 5/9 damped oscillations. At beta 0.8, filtering reduced synthetic-truth RMSE in 3/3 stationary-noise, 9/9 trend, 0/9 AR(1), 1/9 damped-oscillation, 9/9 variance-shift, and 9/9 mean-shift records. | Three seeds per setting provide little precision. Fixed midpoint and Gaussian constructions do not cover unknown change times or diverse observation processes. Aggregate RMSE can hide transition smearing or loss of genuine structure. |
| 43 — [Diagnostic robustness](phase_43_diagnostic_robustness.md) | Vary step magnitude, trend slope/noise, variance multiplier, damping/noise, and use 20 replicates per configuration. Reassess Phase 42’s fixed two-half mean and variance rules; 460 generated records are recorded. | Mean flags rose from 2/20 and 4/20 for steps 0.1 and 0.2 to 20/20 at 0.5; the same rule flagged 19/20 trends at slope 0.001 and 20/20 at 0.002. Variance flags rose from 0/20 at SD multiplier 1.0–1.1 to 20/20 at 2.0; damping 0.002 was flagged 19/20 at noise SD 0.2, but its flag rate changed from 20/20 at noise SD 0.05 to 0/20 at 0.5. | Fixed thresholds and 20 seeds characterize this grid only. Half-record statistics do not distinguish a step from drift or random variance change from deterministic envelope decay; no unknown-location search or independent threshold calibration was used. |
| 44 — [Robust step detection](phase_44_step_detection_robustness.md) | Preserve the Phase 42 baseline `abs(two_half_mean_z) >= 3`; compare it with global-slope, sinusoidal-envelope, and combined vetoes. Sixteen conditions use 40 calibration and 100 disjoint evaluation replicates each. Veto thresholds use calibration data only. | On 600 positives and 1,000 negatives, baseline TP/FP/FN/TN were 339/300/261/700 (FPR 30.0%, FNR 43.5%). Trend veto removed all baseline false alarms but also detected no positives. Envelope veto yielded 189/174/411/826 (17.4% FPR, 68.5% FNR). The combined veto likewise detected no positives. | The global checks can be influenced by the step itself; the envelope fit assumes the generated known frequency. The design is finite, synthetic, and oracle-midpoint. These are veto trade-offs, not general detector properties. |
| 45 — [Local level-shift detection](phase_45_local_step_detection.md) | Compare the preserved Phase 44 rules with raw and full-record-slope-adjusted midpoint differences using windows of 16, 32, and 64 samples per side. The same 16 conditions and 40/100 seed split are reused. Six local thresholds are selected on calibration data and frozen. | Overall, raw w=32 yielded 385/212/215/788 (FPR 21.2%, FNR 35.8%); adjusted w=32 yielded 418/288/182/712 (28.8%, 30.3%); adjusted w=64 yielded 347/91/253/909 (9.1%, 42.2%). The baseline detected 2/100 pure 0.1 steps and 96/100 0.2-step-plus-trend cases; local sensitivity varied considerably by window and mixture. | The midpoint is known, local comparisons are retrospective, and the full-record slope can be biased by a true step. No formal confidence intervals were reported. Noise levels were still coupled to families. |
| 46 — [Independent noise sensitivity](phase_46_noise_sensitivity.md) | Cross noise SDs 0.5, 1.0, and 1.5 with stationary nulls, a fixed trend, four pure-step magnitudes, and those steps plus a trend. Compare only the Phase 44 baseline and three Phase 45 rules, with their prior thresholds frozen. Use matched noise seeds across levels, disjoint calibration/evaluation seeds, and 40/100 replicates per cell. | On stationary nulls, FPR rose across noise levels for raw w=32 (5%, 43%, 57%), adjusted w=32 (20%, 52%, 67%), and adjusted w=64 (0%, 22%, 45%); baseline stayed at 2%. On pure trends, baseline FPR fell from 100% to 11%, while local rates generally rose at higher noise. Miss rates varied by method and step family; for pooled pure steps baseline FNR rose 22.0% → 53.2%, while adjusted w=32 was 27.5% → 22.0% under the frozen absolute threshold. The report and summary include pointwise Wilson intervals and all magnitude strata. | Rates are specific to additive Gaussian signals, fixed parameters, a known midpoint, and frozen thresholds calibrated on an earlier condition mixture. Wilson intervals are marginal and not adjusted for multiple comparisons. Increased local detection under higher noise can reflect more threshold-crossing excursions, not improved information. |

## 2. Evidence, interpretation, and limits

### Directly supported by these experiments

- On the generated records and declared seeds, diagnostics often flagged more than one signal family. In these tests, flags described statistic values; they did not uniquely identify causes.
- Effect size, noise, trend, damping, window size, and threshold changed measured detection and false-alarm rates. The direction of change was method- and condition-dependent.
- The Phase 44 global trend veto eliminated the tested baseline false positives and also suppressed all tested true steps. Phase 45 local methods retained some detections but traded false alarms against missed steps differently by window.
- Phase 46 independently varied Gaussian noise within four signal families. It documented interactions between noise, method, and step magnitude under the frozen thresholds.
- Fixed-seed repeated runs reproduced the recorded outputs for the phases that checked artifact hashes. Phase 46 included per-record data, subgroup denominators, and Wilson intervals.

These are results about the finite generated samples and the implemented statistics. They are not evidence about what generated an unobserved real series.

### Plausible interpretations

- Global statistics can respond to a broad trend or changing envelope even when the target event is absent; local statistics can respond to random local imbalance.
- With fixed absolute local thresholds, larger noise can increase false alarms and can also increase apparent detections among positive records. A larger detection fraction under that rule does not mean the signal became easier to distinguish.
- A veto that reacts to a property shared by confounders and true positives can reduce both false positives and true detections.

These interpretations are consistent with the observed results, but do not identify a unique underlying mechanism beyond each synthetic generator.

### Unresolved questions and claims not established

Unknown change-point location, non-Gaussian or dependent noise, missing or irregular observations, nonstationary noise, and broader signal families remain untested by Phases 42–46. Phase 46 did not test automatic threshold adaptation to noise. Finite replicate counts and the chosen condition mixtures do not establish universal operating characteristics.

None of these phases establishes real-world prediction, causal explanation, financial-market predictability, investment utility, or profitability. The synthetic outputs cannot validate an application to market data or bypass data authorization.

## 3. Reusable methodological lessons

- Change one declared factor at a time or use an explicit factorial crossing; document factors that remain confounded.
- Keep calibration and evaluation seeds disjoint. Freeze thresholds and configuration before evaluating held-out records.
- Use matched seeds for method comparisons where practical, and report TP/FP/FN/TN with the class denominator for each subgroup.
- Report results by signal family and effect size. A pooled score can hide opposite responses.
- Distinguish numerical reproducibility and solver accuracy from the validity of the model assumptions.
- Detection of a measured pattern is not identification of its generating process or cause.
- Describe uncertainty at the level the design supports; a small fixed-seed grid is not population validation.

## 4. Why move beyond the detector family

Phases 42–46 have now explored the principal ambiguity of fixed-window synthetic diagnostics through confounders, local windows, and an independent noise sweep. More variants of the same statistics would mostly extend the same methodological line. The next question should use Newton Lab’s existing physics mathematics to make a prediction about a distinct, explicitly defined system, with a baseline that can fail. The proposal below uses diffusion modes to predict a network algorithm’s relaxation rates; it does not reinterpret the detector results as an application.

## 5. Selected cross-domain hypothesis

### Source domain

One-dimensional transient heat conduction in a homogeneous medium with constant diffusivity and fixed endpoints. After subtracting the steady profile, the perturbation satisfies

`∂u/∂t = α ∂²u/∂x²`, with `u(0,t) = u(L,t) = 0`.

For a single sine mode, `u_n(x,t) = A_n sin(nπx/L) exp[-α(nπ/L)²t]`.

### Target domain

A leader-pinned distributed agreement protocol on a path of agents. Interior agent `i` changes its scalar estimate using only its two neighbors; two fixed boundary leaders hold the same reference value. This is a networked computation model: the state is an agent estimate and the independent variable is algorithmic relaxation time, not temperature or physical heat flow. Because the leaders are fixed, this is pinned/leader-following relaxation, not unconstrained average consensus.

### Shared mathematical structure and hypothesis

On a normalized interval of length 1, write `ξ=x/L` and `τ=αt/L²`. The heat mode rate is `λ_n=(nπ)²`. For `N` interior agents, let `h=1/(N+1)` and use the local update law

`dz_i/dτ = [z_(i-1) - 2z_i + z_(i+1)] / h²`,

with fixed leader states `z_0=z_(N+1)=0`. Its discrete mode `sin(nπih)` has rate

`λ_(N,n) = 4 h^(-2) sin²(nπh/2)`.

**Hypothesis:** for each fixed low mode `n ∈ {1,2,3}`, the normalized rate `λ_(N,n)/λ_n` approaches 1 from below with second-order spatial error as `h` is refined. Specifically, the expected ratio is `[sin(nπh/2)/(nπh/2)]²`, so `1 - λ_(N,n)/λ_n` is approximately `(nπh)²/12` at fine resolution.

This is an established Laplacian/diffusion mathematical structure. It is not a claim that an agent network conducts heat. The proposed experiment would test whether this precise mode-rate prediction is reproduced by Newton Lab’s continuum reference and a separately implemented network evolution, including their numerical error.

### Competing explanation

An apparent match could come only from deliberately choosing the same nearest-neighbor second-difference stencil and fixed boundaries. It would not show that arbitrary communication networks, delayed updates, weighted or directed graphs, or practical distributed algorithms behave like a continuous heat medium. A fit on a single low mode or a visually similar plot could conceal this narrow construction.

### Assumptions

- Uniform one-dimensional path; equal symmetric neighbor weights; two fixed leaders at the same reference value.
- Continuous-time, synchronous, linear updates with no delay, noise, dropout, or state-dependent weights.
- Constant diffusivity, fixed physical endpoints, and the corresponding homogeneous Dirichlet perturbation.
- Compare matched dimensionless state amplitudes, domain coordinate, and relaxation time. The physical heat variables have units; normalized agent estimates do not. Direct dimensional equality is not assumed.
- Only modes `n < N+1` are compared, and low-mode continuum behavior is expected only when the grid resolves the wavelength.

### Predictions and falsification

If supported, (1) the numerically estimated graph-mode rate agrees with the exact discrete eigenvalue at the declared solver tolerance; (2) the graph rate is below the continuum heat rate for each tested finite grid; (3) the relative error decreases as the node spacing shrinks; and (4) its refinement order approaches 2 for fixed modes.

The null is that the rate ratio does not show the predicted convergence over the declared path sizes, after numerical solver error is shown to be smaller than the discretization error. Predeclare path interior sizes `N ∈ {7, 15, 31, 63, 127}` and fixed modes `n ∈ {1,2,3}`. Estimate convergence order by linear regression of `log(relative rate error)` on `log(h)` over the finest four grids for each mode. Count the prediction as supported in this finite design only if each mode has rate ratio below 1 on those grids, errors decrease on each refinement, and fitted order lies in `[1.7, 2.3]`. Falsify the predicted finite-grid relation if a mode persistently has the wrong direction, fails to converge toward 1, or fails the predeclared order criterion after numerical convergence is established. If the ODE solver does not agree with the exact graph eigenmode to the numerical tolerance, classify the experiment as numerically unresolved rather than treating that as a failure of the mathematical hypothesis.

### Scope of a positive result

A positive result would support this low-mode continuum approximation for the specified pinned path protocol, parameterization, and grids. The discrete eigenvalue formula already predicts the relation analytically; the experiment checks the independently implemented network trajectory, extraction of modal rates, and convergence under this cross-domain interpretation. It would not establish physical equivalence, transfer to general graphs or discrete-time algorithms, empirical distributed-system performance, or a broad theory of consensus.

## 6. Alternatives considered

1. **Damped oscillator to second-order risk adjustment.** Newton Lab already implemented Phase 18’s exact coefficient mapping for a simplified adjustment model. Repeating that linear correspondence would add little beyond its existing report, and using market data would remain blocked by Phase 32 authorization.
2. **Damped oscillator to momentum optimization.** A heavy-ball optimizer near a quadratic minimum has oscillator-like linearized dynamics. It is not selected because the project contains no optimization-dynamics target implementation, and a purely quadratic test would reproduce the correspondence by construction. A meaningful nonlinear-objective study would require additional modeling choices beyond this bounded synthesis.

## 7. Reproducible experiment plan (design only)

### Models, factors, and initial/boundary conditions

- **Continuum source:** `HeatDiffusionModel` with one mode of unit normalized amplitude, homogeneous fixed endpoints, normalized `L=1`, `α=1`; `u_n(ξ,0)=sin(nπξ)`. This produces `a_n(τ)/a_n(0)=exp[-(nπ)²τ]`.
- **Network target:** `N` interior agent values with fixed zero leader states at indices 0 and `N+1`; initialize interior values as `z_i(0)=sin(nπih)`. Integrate `dz/dτ=-L_h z` using the path-graph Dirichlet Laplacian `L_h` with diagonal `2/h²` and adjacent entries `-1/h²`.
- **Vary:** `N ∈ {7,15,31,63,127}` and mode `n ∈ {1,2,3}`. Hold `α=1`, normalized length 1, update coefficient 1, boundaries, initial modal amplitude, and all solver/output settings fixed.
- **Independent variables:** grid spacing `h` and mode index `n`.
- **Dependent variables:** projected mode amplitude over time, estimated modal decay rate, rate ratio to the continuum heat prediction, and observed convergence order.
- **Time sampling:** use dimensionless modal exposure `s=λ_n τ` from 0 to 5 for each mode, so the continuum reference spans the same number of e-folding times.

### Baselines, metrics, and numerical checks

- **Primary null:** the graph rate does not converge at second order to the continuum heat rate for fixed low modes.
- **Reference baselines:** the exact continuum rate `λ_n=(nπ)²` and exact discrete graph rate `λ_(N,n)`. The latter separates ODE-solver error from finite-grid approximation error.
- Project the evolving network state onto its initialized sine eigenvector and estimate the log-amplitude slope. Compare this measured rate to both analytical rates.
- Report relative rate error, the analytical rate ratio, modal trajectory error against the exact graph exponential, and fitted refinement order. Retain the complete mode/grid table.
- Run the target integration at two declared tolerance levels, for example DOP853 `rtol=1e-9, atol=1e-11` and `rtol=1e-12, atol=1e-14`; require tighter integration to reduce or leave negligible the graph-solution discrepancy. Also compare at least two output-grid densities. Solver error must remain below the finite-grid continuum discrepancy on the coarser test grids.
- The model is deterministic, so random confidence intervals are not appropriate. Assess numerical uncertainty through solver tolerance, output-grid refinement, exact discrete-rate comparison, and the spread of convergence order across the three modes.

### How an apparent match could mislead

The same path-Laplacian stencil and boundary values are built into both sides, so agreement demonstrates this discretization’s spectral limit rather than independent causal discovery. Large `N`, low `n`, or a selected time window can make curves look similar while concealing rate bias. A target solver error could be mistaken for a physical discrepancy. The result would not include asynchronous messages, finite update steps, graph heterogeneity, stochastic disturbances, or measured network data.

### Minimum evidence and interpretation

Retain the hypothesis only for the tested setting if all three modes meet the frozen convergence criteria, tighter numerical settings confirm the rates, and the modal trajectory matches the exact graph evolution within declared numerical tolerance. If the analytic ratio is right but trajectory estimates do not converge numerically, repair or reject the implementation before interpreting the cross-domain comparison. If the graph rate fails the analytical ratio after solver convergence, reject the stated transfer for this implementation and investigate the boundary or scaling assumptions. Do not change grid sizes, modes, or criteria after viewing held-out-style outputs.

## 8. Existing code to reuse

- `src/newton_lab/heat_diffusion_case_study.py`: `HeatDiffusionModel`, `SpatialMode`, and `evaluate_heat_diffusion` provide the fixed-end continuum modal reference. Their existing boundary assumptions match the proposed pinned path endpoints.
- `src/newton_lab/modal_research.py`: `heat_mode_decay_rate(diffusivity_m2_per_s, mode_number, length_m)` supplies the analytical continuum rate `α(nπ/L)²`.
- `src/newton_lab/dynamics.py`: demonstrates the project’s validated `solve_ivp` use and numerical result conventions. Its two-state oscillator API should not be repurposed for a graph state vector.
- `src/newton_lab/simulation.py`: `SimulationModel` is the base for validated, serializable result records.
- `src/newton_lab/experiments.py`: provides typed experiment specifications, sweeps, run provenance, and adapters for existing model families. It has no path-network model today; do not claim an existing adapter can run this target without extension.
- `src/newton_lab/cross_system_discovery.py`: contains `ComparisonEntry`, `CrossSystemFinding`, and evidence classifications that can later record the scoped relationship without claiming physical equivalence.
- `tests/test_heat_diffusion_case_study.py` and `tests/test_cross_system_discovery.py`: establish model and evidence conventions. A later implementation will need focused tests for graph Laplacian eigenvalues, mean/reference boundary behavior, solver convergence, and rate extraction.

A later implementation will require a small network-relaxation model and result type plus a focused research runner; no such component is created in Phase 47.

## 9. Implementation checklist for the next phase

1. Implement the pinned path update law and validated grid/mode inputs without changing the existing heat or oscillator APIs.
2. Generate the declared single-mode initial conditions and retain exact condition, solver, and boundary metadata.
3. Produce the continuum heat and exact discrete graph analytical references.
4. Integrate each network condition at both solver tolerance settings and two output-grid densities.
5. Estimate modal amplitude/rate and calculate errors, refinement order, and criterion status exactly as predeclared.
6. Add deterministic tests for the update matrix, analytical eigenvalues, modal initialization/projection, seedless reproducibility, and numerical convergence.
7. Write the results and report only after the frozen design is implemented; label evidence as mathematical relation plus finite numerical validation, not a discovery of physical equivalence.
8. Keep Phase 32 at `BLOCKED_AUTHORIZATION`; do not use market data or infer financial performance.

## Project status

The Phase 42–46 detector studies are complete and remain unchanged. Phase 47 records a synthesis and a proposed experiment only. No new code, tests, datasets, plots, or generated output are part of this phase. Phase 32 remains `BLOCKED_AUTHORIZATION`.
