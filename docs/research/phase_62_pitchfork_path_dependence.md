# Phase 62 — Hysteresis and Path Dependence: A Bounded Research Pilot

## Question and feasibility decision

This pilot asks whether a controlled parameter cycle in a canonical system can produce path-dependent observations, and whether that effect should be called hysteresis. The feasible test uses the supercritical pitchfork normal form already implemented for Phase 53–55 recovery-scaling work. Those phases studied noisy trajectory estimation; Phase 62 instead varies the control parameter continuously and carries the state through an up/down cycle.

The model is

\[
\dot{x}=\mu(t)x-x^3.
\]

It is suitable for testing finite-rate lag and initialization-based branch selection. Its stable equilibrium magnitude is single-valued, so it does **not** have a quasi-static magnitude hysteresis loop. The pilot therefore tests finite-rate path dependence, not classical rate-independent hysteresis. This is a synthetic normal-form capability study, not evidence about a real system.

## Equilibria and stability

Equilibria satisfy \(x(\mu-x^2)=0\). The zero branch exists for every \(\mu\), with linear eigenvalue \(f_x(0)=\mu\): it is asymptotically stable for \(\mu<0\), nonhyperbolic but asymptotically stable at \(\mu=0\) (nearby solutions decay algebraically), and unstable for \(\mu>0\). For \(\mu>0\), the two branches \(x^*_\pm=\pm\sqrt{\mu}\) have eigenvalue \(f_x(x^*)=\mu-3(x^*)^2=-2\mu\), hence are stable. There is no interval in which the stable zero branch coexists with either stable nonzero branch. The stable equilibrium magnitude is \(0\) for \(\mu\leq0\) and \(\sqrt{\mu}\) for \(\mu>0\), and is single-valued. Consequently the quasi-static magnitude loop area is analytically zero. The signed nonzero branch is selected by the sign of the initial state; the exactly zero trajectory is invariant even after it becomes unstable.

## Frozen protocol and methods

The protocol was frozen before primary evaluation (`reports/phase_62_pitchfork_path_dependence/protocol.json`; SHA-256 `A2D35E4AD3B756AA6071AC3C7B4F0141ACEC05BC9DDE9BFE09981997A66351E6`). The continuous triangular control ramp spans \(\mu=-0.16\) to \(+0.16\), with leg durations 20 and 80 dimensionless time units. Each upward trajectory is continued directly into its downward leg. Initial states are \(+0.01\), \(-0.01\), and exactly zero. Thus the finite-rate result depends on the declared initial state and ramp protocol.

The implementation integrates the time-varying ODE with SciPy `solve_ivp` (DOP853; base `rtol=1e-9`, `atol=1e-11`, maximum step and output spacing 0.05). It separately checks output-grid refinement (spacing 0.025) and tighter solver settings (`rtol=1e-11`, `atol=1e-13`, maximum step 0.0125). The primary observable is the area between the up/down **absolute-state** curves, \(\int | |x_\uparrow|-|x_\downarrow| |\,d\mu\). It is normalized by the stable-branch reference area \((2/3)\mu_{max}^{3/2}=0.0426667\). This metric deliberately removes the sign difference between symmetry-related branches. The zero initial condition is retained as an invariant-control trajectory, not treated as a perturbed physical state.

Frozen acceptance requires numerical changes to remain below 1% of the branch-area reference. The three hypotheses and decision thresholds are recorded in the protocol; no thresholds were changed after evaluation. There were 12 primary integrations and 8 numerical-resolution integrations, all successful. The output includes 12,012 trajectory rows, six path-metric rows, and four numerical-check rows.

## Results

| Result | Fast ramp, leg duration 20 | Slow ramp, leg duration 80 |
|---|---:|---:|
| Mean loop area, nonzero initial signs | 0.00364812 | 0.02491512 |
| Normalized loop area | 0.08550 | 0.58395 |
| Maximum area change under output-grid refinement | \(1.26\times10^{-8}\) | \(2.01\times10^{-9}\) |

Tight-solver area changes were \(1.34\times10^{-17}\) (fast) and \(3.12\times10^{-17}\) (slow). All four changes are well below the frozen acceptance bound 0.000426667. Positive and negative initial states produced exactly symmetric absolute-state responses in the recorded outputs. The zero state remained exactly zero. At the upper turning point, the nonzero state magnitudes were about 0.00999 for both ramp durations, far below the instantaneous stable equilibrium magnitude 0.4. This indicates substantial finite-rate tracking lag over the tested protocol.

- **H1, finite-rate path difference:** supported within tested conditions. Nonzero-state loop areas exceed ten times the corresponding numerical area changes for both ramp speeds; the invariant zero control has zero loop area.
- **H2, slower ramp reduces loop:** inconclusive under its frozen rule. The measured slow-ramp area (0.0249151) is larger than the fast-ramp area (0.00364812). This is contrary to the simple expectation that more time per leg always reduces the loop. The result may reflect the chosen finite observation/ramp protocol and transient/initial-condition interaction; this pilot does not isolate those mechanisms or establish a general ramp-duration law.
- **H3, initial sign selects the pitchfork branch:** supported within tested conditions. The positive and negative initial states retain their respective signs, with symmetric absolute responses.
- **Quasi-static equilibrium magnitude hysteresis:** not supported by this model's equilibrium structure. The stable magnitude is single-valued; this analytical conclusion does not negate the observed finite-rate loop.

## Artifacts and reproduction

Implementation: `src/newton_lab/pitchfork_path_dependence.py`. Focused tests: `tests/test_pitchfork_path_dependence.py`. The frozen protocol and reproducible outputs are in `reports/phase_62_pitchfork_path_dependence/`: `ramp_trajectories.csv`, `path_metrics.csv`, `numerical_checks.csv`, `hypothesis_results.json`, `metadata.json`, two PNG figures, and `sha256_manifest.json`.

Run the experiment from the project root with:

```powershell
.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.pitchfork_path_dependence import run_pitchfork_path_study; run_pitchfork_path_study(Path('reports/phase_62_pitchfork_path_dependence'))"
```

The runner refuses to write over a populated output directory. An isolated rerun produced byte-identical files for all eight manifest entries; each matched the primary manifest, and the two manifest files were also byte-identical.

## Scope and stopping decision

The equations determine the equilibrium branches and their stability; the experiment does not discover those laws. The observed loop is a finite-rate lag/path effect under one bounded range, two ramp durations, three initial states, and deterministic dynamics. It is not evidence of rate-independent hysteresis, a universal rule that slower ramps enlarge loops, noise robustness, or behavior in a physical, engineering, financial, or other real system. No cross-domain connection was established. The exactly zero trajectory is mathematically invariant but is unstable on the positive-\(\mu\) segment and is sensitive to perturbations, so it is only a deterministic control. The pilot ends here; it does not initiate another study.
