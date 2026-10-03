# Phase 20: Discovery knowledge integration

Phase 20 adds equation identities and query paths for the Phase 18 oscillator
case study and Phase 19 heat-diffusion case study. It reuses the equation
registry, `EquationRelation`, workflow relationship graph, and discovery
records. Equation text remains descriptive metadata; no symbolic parser or
second registry is introduced.

## Heat equation records

The example registry now contains three distinct records:

| ID | Expression | Scope |
|---|---|---|
| `steady_heat_conduction_divergence_1d` | `d/dx(k(x) dT/dx) = 0` | Source-free, steady, one-dimensional conduction; conductivity may vary with position. |
| `steady_heat_conduction_1d` | `d²T/dx² = 0` | The constant positive conductivity reduction of the divergence-form equation. |
| `transient_heat_diffusion_1d` | `∂T/∂t = α ∂²T/∂x²` | Source-free 1D diffusion with homogeneous constant properties. |

Both MIT OpenCourseWare references are recorded on the equations. The notes
*Introduction to Engineering Heat Transfer, Part 3* derive the 1D steady
conduction equation and its constant-property slab reduction; *1-D Thermal
Diffusion Equation and Solutions* gives the transient equation and
`α = k/(ρ cp)`. See the [steady conduction notes](https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/87d9f4544b7fd64a77201382500d057c_10_part3.pdf)
and the [transient diffusion notes](https://ocw.mit.edu/courses/3-185-transport-phenomena-in-materials-engineering-fall-2003/927ddb6b3dfc423de5570a0070614129_handout_htrans.pdf).

The registry stores two directed, mathematically derived relationships:

1. The constant-conductivity steady equation is a special case of the
   divergence-form equation: setting `k(x)=k>0` gives `k T_xx=0`, hence
   `T_xx=0`.
2. The transient equation has the steady equation as its stationary limit when
   `T_t=0`, `α` is constant and positive, material properties are constant, and
   the medium is 1D and source-free. With matching fixed endpoint temperatures,
   the selected steady profile is linear. Boundary conditions select that
   particular solution; the two equations are not identical.

These assumptions and derivations are recorded in each `EquationRelation`
evidence statement. `EquationRelation` has no separate structured `conditions`
field, so the conditions remain explicit prose alongside the typed relation.
No registry or knowledge-store schema migration was required.

The records classify the steady equations as linear ODEs and the transient
equation as a linear PDE. Thermal diffusivity is recorded in `m²/s`, thermal
conductivity in `W/(m K)`, density in `kg/m³`, and specific heat in
`J/(kg K)`. The built-in heat BVP still describes only its constant-property
steady problem; no transient simulation is inferred from it.

## Oscillator and diffusion structures

The existing `damped_harmonic_oscillator` record is reused. Its physical
relationship to Newton's second law remains in place. A separate workflow edge
compares it with `transient_heat_diffusion_1d` as a candidate shared structure
because both records explicitly classify their equations as linear and
deterministic. The candidate records their important differences: oscillator
responses may be oscillatory or non-oscillatory depending on damping, while
the diffusion modes decay without oscillation at rates proportional to spatial
frequency squared. It asserts no physical equivalence and remains
`candidate_requires_review`.

Phase 18's parameter map remains attached to its existing hypothesis and case
study: `τ=sqrt(m/k)`, `ζ=c/(2 sqrt(mk))`, and the displacement scale maps to a
dimensionless adjustment state. Its mathematical correspondence is limited to
the declared constant-target transformation. The target is a structural
analogy and the financial application remains unvalidated.

## Retrospective experiment associations

The source Phase 18 oscillator experiment and Phase 19 heat BVP had empty
`equation_record_ids`. The integration does not edit those specifications or
rewrite their historical discovery observations. Instead,
`ExperimentEquationAssociation` records:

- The source experiment and model IDs and the equation ID associated now.
- Whether the association was added retrospectively or was already declared.
- A compatibility rationale and limitations.
- The original equation IDs and the original discovery observation identities.
- Exact `ExperimentRun` snapshots, including run status, solver status,
  acceptance, metrics, and BVP mesh/residual diagnostics where applicable.

The association builder checks the oscillator model and SI parameter units, or
the heat BVP problem ID, units, constant-conductivity assumption, and absence
of internal generation. Unsupported model/equation pairs are rejected. The
Phase 19 observations continue to show that they originally had no equation
identity; the retrospective link does not make it appear otherwise.

## Query and report

Pass the already available Phase 18 and Phase 19 case-study objects to the
integration builder. It only reads the supplied objects and does not run a
solver. It returns an immutable `KnowledgeIntegrationIndex` containing the
equation registry, workflow graph, association records, and a deterministic
Markdown report.

```python
from newton_lab.knowledge.integration import (
    build_knowledge_integration,
    query_discovery_knowledge,
)

# oscillator_case and heat_case are existing returned case-study objects.
index = build_knowledge_integration(oscillator_case, heat_case)

heat_knowledge = query_discovery_knowledge(index, "transient_heat_diffusion_1d")
print(heat_knowledge.equation_records)
print(heat_knowledge.mathematical_structures)
print(heat_knowledge.registry_relationships)
print(heat_knowledge.workflow_relationships)
print(heat_knowledge.provenance)

# The complete integration result is serializable and can be restored without
# re-running the experiments:
saved = index.model_dump_json()
restored = type(index).model_validate_json(saved)
```

The same query accepts a case-study ID, such as
`phase19_heat_diffusion_smoothing`. Results separate registry-backed
mathematical derivations, review-required workflow relationships, and the
unvalidated Phase 18 application hypothesis. The report includes provenance,
relationship rationale, original observation identity, run evidence, and
limitations. It does not infer links from matching names or units.

## Evidence limits

- Registry mathematical relationships are scoped by their written
  assumptions; they do not assert physical equivalence outside those conditions.
- The oscillator/diffusion connection is a candidate structural comparison,
  not a stored mathematical equivalence.
- Retrospective association is metadata added after an experiment. It does not
  change what the original solver or discovery pipeline recorded.
- Numerical acceptance is not physical validation.
- The Phase 18 finance hypothesis remains unvalidated; Phase 19's synthetic
  smoothing example is not real market data or evidence of predictive value.
- The equation registry is a small starting set, not complete physics coverage.
