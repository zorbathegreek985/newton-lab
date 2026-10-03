# Phase 60 — Phenomenon-Driven Discovery Pilot: Kuramoto Synchronization

## Research question

How does increasing coupling produce collective synchronization in a population of phase oscillators, and how does the observed transition depend on the natural-frequency distribution and finite population size?

This is a bounded study of the Kuramoto phase model. The oscillator phases are the physical phenomenon being represented; this is not a simulation of full mechanical oscillator amplitudes. The study asks whether a macroscopic coherence measure changes with coupling in two frequency populations, and whether finite-population transition estimates agree approximately with the continuum prediction. It does not claim novelty for the model or validate synchronization in a measured device.

## Model, theory, and assumptions

For oscillator \(i=1,\ldots,N\), the simulated equation is

\[
\frac{d\theta_i}{dt}=\omega_i+\frac{K}{N}\sum_{j=1}^{N}\sin(\theta_j-\theta_i).
\]

Here \(\theta_i\) is phase (radians), \(\omega_i\) is its fixed natural angular frequency (radians per second), \(K\geq0\) is coupling (radians per second), and \(t\) is time (seconds). The \(1/N\) normalization keeps the mean-field coupling scale comparable across population sizes. The implemented vector field evaluates the interaction as \(K\,\mathrm{Im}(Z e^{-i\theta_i})\), where

\[
Z(t)=\frac{1}{N}\sum_j e^{i\theta_j(t)}=R(t)e^{i\Psi(t)}, \qquad R(t)=|Z(t)|.
\]

Thus \(R\in[0,1]\): zero is incoherence in the continuum idealization, and one is perfect phase alignment. For finite random populations, incoherent phases generally give nonzero \(R\) of order \(N^{-1/2}\); a nonzero order parameter alone is not evidence of a synchronized state.

The model is a phase reduction: it assumes oscillators remain near stable limit cycles, coupling is weak enough for phase-only dynamics, the interaction is identical and sinusoidal for all pairs, natural frequencies are fixed, and there is no forcing, delay, process noise, or adaptive network. The all-to-all mean-field form is deliberate. It isolates collective synchronization and does not reuse Newton Lab's existing single-oscillator ODE, diffusion, or pinned-path calculations as a substitute model.

For a centered symmetric unimodal continuous frequency density \(g(\omega)\), the standard continuum onset prediction for this coupling convention is

\[
K_c=\frac{2}{\pi g(0)}.
\]

It comes from loss of linear stability of the incoherent phase distribution at the zero-frequency mode. It applies to the infinite-population mean-field model under the stated density and coupling assumptions. It is not an exact threshold for a finite random sample, a universal threshold for every frequency law, or a promise that a finite-time sweep locates \(K_c\) precisely.

Two unit-scale populations were fixed before final evaluation:

* **Gaussian:** \(g(\omega)=\exp[-\omega^2/(2\sigma^2)]/(\sqrt{2\pi}\sigma)\), centered with \(\sigma=1\). Thus \(g(0)=1/\sqrt{2\pi}\) and \(K_c=\sqrt{8/\pi}\simeq1.59577\) rad/s.
* **Symmetrically truncated Lorentzian:** \(g(\omega)=\gamma/[2\arctan(W/\gamma)(\omega^2+\gamma^2)]\) for \(|\omega|\leq W\), zero outside, with \(\gamma=1\), \(W=5\) rad/s. The normalization is included, giving \(g(0)=1/[2\gamma\arctan(W/\gamma)]\) and \(K_c=4\gamma\arctan(W/\gamma)/\pi\simeq1.74867\) rad/s. This is not the unbounded Lorentzian; its cutoff is part of the declared model.

Both densities are symmetric and unimodal, so the continuum expression is applicable to these ideal populations. The formula uses the actual normalized \(g(0)\) for each density. Finite samples fluctuate around these densities.

## Frozen design and numerical method

The research hypothesis was that \(R\) should increase across a coupling-controlled transition; the transition location and sharpness should depend on \(g(0)\) and finite \(N\). The operational observable is the mean of \(R(t)\) over 80–120 seconds. The operational finite-grid transition is the midpoint of the adjacent \(K\) interval with the greatest positive finite-difference slope in ensemble-mean post-transient \(R\). No post-hoc coherence threshold is used.

The frozen protocol is [protocol.json](../../reports/phase_60_kuramoto_synchronization/protocol.json). It sets \(N=64,256\), eight replicate seeds (60001–60008), and 11 coupling factors from \(0.60K_c\) to \(1.40K_c\). Each replicate's sampled frequencies and initial phases are reused at every coupling for paired comparisons; frequency and phase draws use separate deterministic NumPy seed streams. Initial phases are independent uniform samples on \([-\pi,\pi)\). This is 2 populations × 2 sizes × 8 replicates × 11 couplings = **352 integrations**.

The ODE uses SciPy solve_ivp with DOP853, relative tolerance \(10^{-8}\), absolute tolerance \(10^{-10}\), maximum step 0.1 s, duration 120 s, and output every 0.1 s. To expose incomplete settling, replicate outputs also compare mean \(R\) over 80–100 s with 100–120 s. Transition uncertainty is a 95% percentile interval from 2,000 bootstrap resamples of the eight paired replicate curves. This interval reflects the small replicate ensemble and does not remove coupling-grid resolution or finite-time bias.

One exploratory feasibility attempt used an unbounded Lorentzian sample. Rare, very fast oscillators made adaptive integrations disproportionately expensive; the process was stopped before the runner wrote any results. That attempt is excluded. The final protocol explicitly records it and instead uses the normalized cutoff Lorentzian above. The final protocol hash is recorded in metadata and the output manifest.

## Results

All 352 final integrations completed. The coupling response is visible for both populations: ensemble-mean coherence generally rises as coupling increases, and the representative Gaussian trajectory moves from a diffuse, fluctuating state below predicted \(K_c\) to persistent high coherence above it. A final phase snapshot shows the same contrast for that one realization. The figures are [order parameter versus coupling](../../reports/phase_60_kuramoto_synchronization/order_parameter_vs_coupling.png), [representative coherence traces](../../reports/phase_60_kuramoto_synchronization/representative_synchronization_traces.png), and [representative phase distributions](../../reports/phase_60_kuramoto_synchronization/representative_phase_distribution.png).

At \(K/K_c=0.6,1.0,1.4\), mean \(R\) was respectively 0.203, 0.496, 0.834 for Gaussian \(N=64\); 0.094, 0.330, 0.808 for Gaussian \(N=256\); 0.172, 0.318, 0.607 for truncated-Lorentzian \(N=64\); and 0.089, 0.200, 0.590 for truncated-Lorentzian \(N=256\). These are finite-run observations, not population-independent constants. At the lower coupling, the larger populations' smaller baseline coherence is consistent with reduced finite-sample incoherent order; near and above onset, the curves and replicate spread depend on the distribution and size.

| Frequency population | \(N\) | Continuum \(K_c\) | Estimated transition (95% paired-seed bootstrap interval) | Sweep boundary? |
|---|---:|---:|---:|---|
| Gaussian, \(\sigma=1\) | 64 | 1.596 | 1.532 [1.404, 1.915] | No |
| Gaussian, \(\sigma=1\) | 256 | 1.596 | 1.787 [1.532, 1.915] | No |
| Truncated Lorentzian, \(\gamma=1,W=5\) | 64 | 1.749 | 2.378 [1.539, 2.378] | **Yes, upper edge** |
| Truncated Lorentzian, \(\gamma=1,W=5\) | 256 | 1.749 | 1.959 [1.819, 2.238] | No |

Rates and intervals are rad/s. Coupling-grid interval widths were 0.128 rad/s for Gaussian and 0.140 rad/s for truncated Lorentzian. Both Gaussian intervals include the continuum prediction. The truncated-Lorentzian \(N=64\) maximum-slope estimate is at the upper sweep boundary and therefore unresolved by the frozen rule; the \(N=256\) interval lies above the continuum value. Two sizes, eight frequency/phase realizations, a coarse grid, and finite-duration data do not establish convergence to \(K_c\). The result is mixed rather than a clean finite-size confirmation of the continuum threshold.

The early-to-late mean-\(R\) difference averaged near zero in each population/size stratum, but individual replicate-condition differences ranged from about \(-0.104\) to \(+0.094\). This spread warns against treating every trajectory as stationary. The truncated-Lorentzian \(N=64\) boundary estimate, broad bootstrap interval, and residual window variability make that case particularly inconclusive.

Machine-readable data, metadata, and integrity records are in [reports/phase_60_kuramoto_synchronization/](../../reports/phase_60_kuramoto_synchronization/). The replicate results file preserves seed-level observations; the synchronization summary gives ensemble means and standard errors; the transition estimates file records operational estimates and intervals. Representative time samples and phase snapshots are separately retained. The SHA-256 manifest covers the protocol, metadata, CSV files, and figures.

## Interpretation and reusable structure

**Established mathematical relationship (A).** For the declared continuum Kuramoto mean-field model and a centered symmetric unimodal density satisfying the stated assumptions, the incoherent-state threshold depends on the density through \(g(0)\) as \(2/[\pi g(0)]\). The truncated density needs its own normalization; substituting the unbounded Lorentzian value would be incorrect.

**Observed model behaviour (B/C).** In these finite integrations, higher coupling was associated with higher average coherence for both frequency populations. The response curves, finite-sample baseline, estimated onset, and uncertainty differed with the distribution and \(N\). The representative phase snapshots illustrate, but do not independently quantify, that collective organization.

**Candidate cross-domain connection.** A macroscopic order parameter can summarize collective coordination arising from many local pairwise interactions. This is a structural idea potentially testable in separately specified synchronization systems; no second system was tested here.

**Validated application (E).** None. These simulations do not validate a physical apparatus, biological population, power grid, market, or other real application. No financial data or external observations were used.

The experiment supports the qualitative coupling-to-coherence phenomenon in the implemented model, while the finite-size estimate comparison is only partly resolved. It does not establish universal critical behaviour or prove that finite \(N\) alone explains the offsets. The sweep has two sizes, eight seeds, and 11 coupling levels; it omits frequency sampling uncertainty beyond those draws, alternative initial-state ensembles, hysteresis tests, delays, network topology changes, forcing, and stochastic phase dynamics. The phase model also omits amplitude dynamics and may not describe strongly coupled or nonidentical physical oscillators outside phase-reduction conditions.

## Reproduction and conclusion

Use the project environment, create a fresh result directory, copy the frozen protocol, and run the same study function. Do not run in the canonical output directory because the runner intentionally refuses to overwrite existing artifacts:

    $out = '.pytest-temp/phase60-reproduction'
    New-Item -ItemType Directory -Force $out | Out-Null
    Copy-Item 'reports/phase_60_kuramoto_synchronization/protocol.json' "$out/protocol.json"
    .venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.kuramoto_synchronization_study import run_synchronization_study; print(run_synchronization_study(Path('.pytest-temp/phase60-reproduction')))"

The function verifies the frozen protocol SHA-256 before execution. Compare the resulting manifest with the canonical manifest to check byte-for-byte reproducibility. The recorded environment was Python 3.12.10, NumPy 2.5.3, SciPy 1.18.1, and Matplotlib 3.11.2. Matplotlib emitted a warning that it could not write its user-level font cache in this environment; it nevertheless completed the figures.

The pilot demonstrates how a collective coherence measure responds to coupling in two finite sampled frequency populations, while also showing why continuum theory and finite-time numerical transition estimates must be kept distinct. A justified next question is whether a denser, predeclared near-onset grid and larger populations reduce the observed transition uncertainty and resolve the truncated-Lorentzian offset. No Phase 61 work is started here.
