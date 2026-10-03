# Phase 61 — Kuramoto Transition: Finite-Size and Resolution Study

## Research question and scope

Does the apparent offset between the measured synchronization transition and
the continuum critical coupling for a normalized truncated-Lorentzian frequency
population persist as population size increases and the coupling grid is
refined? This is a controlled simulation of the Phase 60 Kuramoto model. It is
not evidence about synchronization in a measured physical, social, biological,
or financial system.

The equations are

\[
\dot\theta_i=\omega_i+\frac{K}{N}\sum_{j=1}^{N}\sin(\theta_j-\theta_i),
\qquad i=1,\ldots,N.
\]

The continuum reference uses the symmetric truncated-Lorentzian density with
\(\gamma=1\) rad/s and \(W=5\) rad/s:

\[
g(\omega)=\frac{\gamma}{2\arctan(W/\gamma)(\omega^2+\gamma^2)}
\quad (|\omega|\le W), \qquad g(\omega)=0\quad (|\omega|>W).
\]

Its normalization follows from integrating over \([-W,W]\), and
\(g(0)=1/[2\gamma\arctan(W/\gamma)]\). Under the infinite-population
mean-field Kuramoto assumptions (centered, symmetric, unimodal normalized
frequency density; all-to-all sinusoidal coupling scaled by \(1/N\); fixed
frequencies and no forcing, delay, noise, or adaptive coupling),

\[
K_c=\frac{2}{\pi g(0)}=\frac{4\gamma\arctan(W/\gamma)}{\pi}
 =1.748668167243995\;\text{rad/s}.
\]

This is a continuum theoretical reference, not the exact transition of each
finite random population.

## Frozen design and estimator

Protocol version 1.0 was frozen before final evaluation. Its SHA-256 is
`ED307A5E8E7D6619F3D34CA1D2B5BC5B2EF04399DBB0F914EE0352C5FEC5DBF8`.
The primary study used \(N=64,256,512,1024\), eight replicate seeds
61001–61008, and 17 coupling values from \(0.75K_c\) to \(1.35K_c\), with a
denser grid near \(K_c\). For each seed, one 1024-member frequency sample and
one phase sample were generated from independent deterministic streams; smaller
populations use prefixes, so sizes are nested and paired. A seed's sample is
reused across coupling values.

Each trajectory ran from 0 to 120 s with 0.1 s output, DOP853, relative
tolerance \(10^{-8}\), absolute tolerance \(10^{-10}\), and maximum step 0.1 s.
The response is the time average of the order parameter \(R\) over 80–120 s.
At each size the operational transition is the midpoint of the adjacent
coupling interval with the steepest positive finite-difference slope in the
eight-seed mean response. The selected interval is retained as a grid bracket;
it is not the same quantity as the seed-bootstrap interval. A boundary maximum
is censored and is not reported as a reliable point estimate.

Uncertainty resamples the eight seed curves jointly over coupling values 2,000
times. Bootstrap draws whose selected interval touches a sweep edge are counted
as censored. If more than 2.5% are censored, the ordinary percentile interval
is omitted and seed uncertainty is marked unresolved. A separate 1024-member,
four-seed comparison used the same restricted 12-point coupling grid with
relative tolerance \(10^{-10}\), absolute tolerance \(10^{-12}\), and maximum
step 0.05 s. Frozen classification criteria are recorded in the protocol;
they require narrow, uncensored high-size intervals to call agreement or
stable same-side exclusion to call persistent disagreement. With four sizes,
no asymptotic convergence law is fitted.

## Results

The frozen run completed 592 integrations: 544 primary and 48 numerical
resolution comparisons. No integrations failed. The table separates the
selected finite-grid bracket from the seed-bootstrap interval.

| N | Transition midpoint (rad/s) | Grid bracket (rad/s) | 95% seed-bootstrap interval (rad/s) | Offset from continuum Kc (rad/s) |
|---:|---:|---:|---:|---:|
| 64 | 1.63938 | [1.61752, 1.66123] | [1.63938, 2.14212] | -0.10929 |
| 256 | 1.72681 | [1.70495, 1.74867] | [1.68309, 1.98911] | -0.02186 |
| 512 | 1.72681 | [1.70495, 1.74867] | [1.59566, 1.98911] | -0.02186 |
| 1024 | 1.81424 | [1.79238, 1.83610] | [1.59566, 1.85796] | +0.06558 |

All four primary seed-bootstrap intervals include the continuum value. The
N=512 and N=1024 intervals are both wider than the frozen maximum width
\(0.10K_c\) (about 0.17487 rad/s), so the agreement criterion fails. The two
highest-size point offsets also lie on opposite sides of the continuum value;
the persistent-disagreement criterion fails. The frozen conclusion is
**inconclusive**. The point estimates do not show a monotone approach to the
reference over these four sizes.

The tighter numerical settings changed the restricted-grid transition estimate
by 0 rad/s (0 fine-grid intervals). This supports numerical stability for that
particular four-seed, 12-coupling comparison; it does not resolve seed
uncertainty or finite-size variability.

Phase 60's independent N=256 truncated-Lorentzian estimate was 1.95851 rad/s
with interval [1.81861, 2.23830]. Phase 61's N=256 estimate is 1.72681 rad/s
with interval [1.68309, 1.98911]. The intervals overlap, but the realizations
are independent: this descriptive contrast cannot be attributed to population
size or interpreted as a paired replication.

The mean order parameter generally rose across the tested coupling sweeps, but
replicate late-minus-early window changes varied in sign. Those windows are
diagnostics, not convergence tests, and the 80–120 s average is an operational
finite-time response rather than proof of an asymptotic stationary state.

## Reproducible artifacts

The runner is `src/newton_lab/kuramoto_transition_refinement.py`; it reuses the
Phase 60 model in `src/newton_lab/kuramoto.py`. Focused tests are in
`tests/test_kuramoto_transition_refinement.py`. Run the study with:

```powershell
.\.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.kuramoto_transition_refinement import run_transition_refinement_study; print(run_transition_refinement_study(Path('reports/phase_61_kuramoto_transition_finite_size')))"
```

The output directory contains the frozen `protocol.json`, replicate-level and
summary CSVs, transition estimates, numerical-resolution comparison, three
figures, `metadata.json`, and `sha256_manifest.json`. Metadata records Python
3.12.10, NumPy 2.5.3, SciPy 1.18.1, and Matplotlib 3.11.2. Outputs are local
synthetic results; no external data or service was used. The runner refuses to
overwrite a populated output directory.

Figures:

- `reports/phase_61_kuramoto_transition_finite_size/coherence_by_coupling_and_size.png`
- `reports/phase_61_kuramoto_transition_finite_size/transition_vs_population_size.png`
- `reports/phase_61_kuramoto_transition_finite_size/numerical_resolution_check.png`

## Limitations and interpretation

- The continuum formula is conditional on its mean-field assumptions; finite
  samples, finite observation windows, and the operational max-slope rule can
  shift a measured transition.
- Eight seeds leave broad bootstrap intervals at larger N. The frozen width
  criterion therefore prevents a confident agreement claim even though both
  high-size intervals include \(K_c\).
- Nested prefixes create useful paired comparisons across N, but the four
  largest-size numerical-check seeds are only a subset and do not measure all
  solver or sampling uncertainty.
- The finite grid and finite-time averaging do not identify an asymptotic
  critical exponent or finite-size scaling law.
- Results demonstrate behavior of a declared synthetic model only. No
  experimental Kuramoto population, application system, or financial data was
  examined.
- Phase 32 remains `BLOCKED_AUTHORIZATION`; Phase 52 remains `NO-GO`.

## Status

Phase 61's frozen synthetic benchmark and report are complete. Its scientific
answer under the predeclared criteria is **inconclusive**: the observed
finite-size point estimates straddle the continuum reference, uncertainty is
too broad for the agreement criterion, and no stable same-side offset is
established. This does not establish that the Phase 60 offset persists or
vanishes asymptotically.

## Verification record

- Focused Phase 61 tests: 7 passed. The focused test file was rerun after the
  final test-only mypy annotation.
- Full suite: `.venv\Scripts\python.exe -m pytest --basetemp .pytest-temp` —
  595 passed in 39.05 s.
- Ruff: `.venv\Scripts\ruff.exe check .` — passed.
- Formatting: `.venv\Scripts\ruff.exe format --check .` — passed; 166 files
  already formatted.
- Source typing: `.venv\Scripts\mypy.exe src` — passed, 53 source files.
- Project typing: `.venv\Scripts\mypy.exe .` — passed, 101 files. The first
  whole-project run exposed missing SciPy stubs for `scipy.integrate.quad` in
  the synthetic normalization test; a specific `import-untyped` annotation
  was added to that import, then the focused test and whole-project mypy were
  rerun successfully.
- Output manifest: all 9 listed file hashes validated against the primary
  output directory.
- Reproducibility: an isolated rerun under `.pytest-temp` regenerated 9
  manifested outputs plus the manifest itself; all 10 file hashes matched the
  primary run. The temporary rerun directory was removed after comparison.
- Phase 60 protection: 16 recorded hashes (report, protocol, source, tests,
  and all output files) matched exactly. No Phase 48–50 or Phase 52 files were
  edited; Phase 51 Atlas/roadmap updates are intentional.
- Atlas JSON parsed successfully, and all referenced Phase 60–61 report,
  implementation, test, and artifact paths exist. Phase 32 remains
  `BLOCKED_AUTHORIZATION`; Phase 52 remains `NO-GO`.
