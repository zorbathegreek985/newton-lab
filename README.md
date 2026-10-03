# Newton Lab

**A Python research workspace for scientific simulation, numerical experiments, and careful comparison of mathematical structures.**

> **Current status:** Research is paused for consolidation. Narrow maintenance is permitted; no new research phase is authorized. See the [project checkpoint](docs/project_checkpoint.md).

Newton Lab investigates how selected mathematical models behave and when structures shared by different systems may support a testable hypothesis. It combines physical model implementations, numerical solvers, experiments, analysis, and evidence-aware reporting. Similar equations alone do not establish physical equivalence.

The project provides a set of connected capabilities, not one universal automatic pipeline. For supported models, the workflow is:

**Model and problem metadata -> solver adapter -> typed result or failure -> experiment record -> analysis, discovery descriptors, and report**

Study-specific research modules also run frozen protocols and produce their own outputs. The discovery API works on existing experiment records; it does not autonomously discover scientific laws.

## What is implemented

- **Scientific models and solvers:** damped harmonic oscillator and full-sine damped pendulum ODEs; a nonlinear algebraic root solver; and a separate collocation solver for boundary-value problems, including steady one-dimensional heat conduction.
- **Experiment workflows:** typed specifications and results, ordered single-parameter sweeps, retained point failures, analytical-reference comparisons, and central finite-difference sensitivity for supported examples.
- **Analysis and reporting:** read-only experiment summaries, sweep analysis, explicit metric-compatibility checks, optional visualizations, and Markdown reports.
- **Knowledge and discovery:** typed equation/source records, declared mathematical-structure metadata, feature comparisons, and `run_discovery` for selected behaviors and candidate hypotheses from existing experiment records.

Important boundaries: equation strings are descriptive metadata, not executable symbolic math. Structure matching compares declared features; it does not prove symbolic equivalence. Registry records distinguish candidate structural connections and proposed application hypotheses from empirically investigated or independently validated applications. ODE trajectories, algebraic candidates, and BVP profiles remain distinct result types. See the [technical project summary](docs/portfolio_project_summary.md) for the implementation map and API links.

## Selected research examples

- **Phase 48 - discrete modes and heat diffusion:** for a uniform nearest-neighbor path with fixed leaders, tested low-mode decay rates approached fixed-end heat-mode rates with approximately second-order grid refinement. This is a narrow result for the matched stencil and boundary setup, not a result for arbitrary graphs or a measured thermal system. [Report](docs/research/phase_48_pinned_path_convergence.md).
- **Phases 60-61 - Kuramoto finite-size study:** simulated phase-only all-to-all oscillators and compared transition estimates with a continuum reference. At the largest tested population sizes, intervals were too broad under the frozen criterion and point estimates straddled the reference. The thread is **CLOSED — INCONCLUSIVE UNDER THE DECLARED PROTOCOL**; this does not establish that convergence is absent. [Phase 60](docs/research/phase_60_kuramoto_synchronization.md) | [Phase 61](docs/research/phase_61_kuramoto_transition_finite_size.md).
- **Phase 62 - pitchfork path dependence:** a controlled parameter ramp produced finite-rate path differences in the tested supercritical pitchfork model. The slower-ramp-reduces-loop hypothesis remained inconclusive; the result does not establish rate-independent hysteresis or a general ramp law. [Report](docs/research/phase_62_pitchfork_path_dependence.md).
- **Phase 64 - non-normal transient gain:** a stable normal matrix and a matched-eigenvalue non-normal family were compared using a declared Euclidean norm and finite time grids. Sampled amplification was observed for selected couplings in that abstract matrix family; no physical application was tested, and sampled peaks are not asserted to be exact continuous-time maxima. [Report and outputs](docs/research/phase_64_non_normal_transient_growth.md).

The recovery-rate-scaling thread (Phases 48-57) is closed for now under the [Phase 58 decision](docs/research/phase_58_research_thread_closure.md). These studies demonstrate mathematical properties, numerical checks, or protocol-specific synthetic results; they do not establish cross-domain application effectiveness.

## Technology

Verified from `pyproject.toml`:

- Python 3.11 or newer; installable distribution `newton-lab`, import package `newton_lab`.
- Runtime: NumPy, SciPy, pandas, and Pydantic.
- Optional visualization: Matplotlib.
- Development: pytest, Ruff, and mypy; mypy is configured in strict mode.

## Get started on Windows

From the repository root, confirm Python is 3.11 or newer, create a virtual environment, and activate it:

```powershell
python --version
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
```

If PowerShell prevents activation, invoke the environment directly:

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
```

For optional plotting support, install the `visualization` extra in the activated environment:

```powershell
python -m pip install -e ".[dev,visualization]"
```

If activation is unavailable, use `.\.venv\Scripts\python.exe -m pip install -e ".[dev,visualization]"`.

Run the oscillator directly from Python:

```python
from newton_lab.dynamics import simulate_damped_oscillator

result = simulate_damped_oscillator(
    mass_kg=1.0,
    damping_coefficient_kg_per_s=0.2,
    stiffness_n_per_m=4.0,
    initial_displacement_m=0.1,
    initial_velocity_m_per_s=0.0,
    duration_s=10.0,
)
print(result.displacement_m[-1])
```

The example follows the existing public function signature; the model assumptions and numerical details are in the [oscillator notes](docs/damped_harmonic_oscillator.md).

## Checks

After installing the `dev` extra, run checks with the project environment:

```powershell
.\.venv\Scripts\python.exe -m pytest
.\.venv\Scripts\ruff.exe check .
.\.venv\Scripts\ruff.exe format --check .
.\.venv\Scripts\mypy.exe src
```

These are commands to run locally; this README refinement did not run them.

## Documentation

- [Technical project summary](docs/portfolio_project_summary.md) - architecture, selected findings, limitations, and demonstration ideas.
- [Project checkpoint](docs/project_checkpoint.md) - current pause, statuses, maintenance policy, and conditions for resuming work.
- [Research Atlas and roadmap](docs/research/phase_51_research_atlas_and_roadmap.md) - research index and preserved outcomes.
- [Research synthesis and capability-gap audit](docs/research/research_synthesis_capability_gap_audit.md) - evidence categories, API inventory, and capability boundaries.
- [Simulation architecture](docs/simulation_architecture.md) | [experiment guide](docs/experiments.md) | [analysis guide](docs/experiment_analysis.md) | [discovery guide](docs/discovery.md) | [visualization guide](docs/visualization.md).
- [Knowledge and equation registry](docs/physics_registry.md) | [mathematical structures](docs/mathematical_structure.md) | [algebraic solver](docs/nonlinear_algebraic_solver.md) | [boundary-value solver](docs/boundary_value_solver.md).

## Limitations and current status

Newton Lab is a scientific-computing and research platform, **not a validated financial trading system**. It makes no market-prediction or profitability claim. Cross-domain hypotheses require independent validation in the target domain, and no such engineering or financial validation is established by the cited research. Some findings are specific to synthetic models and declared protocols; the Kuramoto finite-size result and Phase 62 slower-ramp hypothesis remain inconclusive.

- Phase 32 financial-data readiness: `BLOCKED_AUTHORIZATION`.
- Phase 52 engineering-system qualification: `NO-GO`.
- Kuramoto thread: `CLOSED — INCONCLUSIVE UNDER THE DECLARED PROTOCOL`.
- Phase 62: complete; its slower-ramp hypothesis remains inconclusive.
- Phase 64: complete, limited to its declared matrix family, state norm, and finite-time sampling.

The project is **paused**. Narrow, justified maintenance is permitted; no new research phase is authorized. See the [project checkpoint](docs/project_checkpoint.md). No Phase 65 work is authorized.
