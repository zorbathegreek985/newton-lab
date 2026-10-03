# Mathematical structure and equation relationships

Newton Lab stores mathematical expressions as descriptive text. It does not
parse those strings, execute them, or prove symbolic equivalence. The typed
`MathematicalStructure` attached optionally to an `EquationRecord` captures
explicitly assessed metadata instead: ODE or PDE order, linearity,
autonomy/time-dependence, operators, variables, parameters, forcing,
conservation laws, symmetries, solution families, and numerical methods.

## Unknown versus assessed absence

For scalar fields, `None` means unknown or unassessed. For tuple fields,
`None` means unknown and `()` means assessed with no listed features. This
distinction prevents an omitted conservation-law assessment from being read as
proof that no conservation law exists. Existing equation classifications,
quantity definitions, units, initial conditions, boundary conditions, and
assumptions remain on `EquationRecord`; the structure model does not duplicate
them.

Add only features supported by analysis or a source. Keep operator names and
feature values stable and descriptive. The model is an extensible metadata
interface for a future symbolic backend; a future backend should have an
explicit input format/version and remain separate from this safe metadata
comparison layer.

## Comparisons

`compare_equations(left, right)` reports each declared feature as shared,
different, partially shared, or unknown. Set-valued features expose the
matching and record-specific values. It also compares recorded assumptions and
initial/boundary conditions, and can include graph relationships when given a
`KnowledgeBase`. It produces no score and always reports that physical
equivalence was not concluded.

The oscillator and pendulum records share second-order autonomous ODE structure
and derivative operators. Their classifications, linearity, states, and
parameters differ. The pendulum retains `sin(theta)`; the harmonic oscillator
has a linear restoring term. In the small-angle regime `|theta| << 1`, using
`sin(theta) approximately theta` and `x = L theta` yields an oscillator-shaped
linear equation when coefficients are matched. Newton Lab records that as a
directional approximation candidate with these conditions; it is not exact
equivalence and remains unreviewed in the example workspace.

## Relationships and discovery

`EquationRelationship` is a graph edge with a stable ID, source and target
equation IDs, typed relation, explanation, justification, optional conditions,
source IDs, directionality, and review status. Derivation, limit, special-case,
approximation, and transformation edges are directional. Shared structure,
operator, assumptions, symmetry, conservation, variables, solution method,
analogy, coupling, and equivalence-under-assumptions links are modeled as
undirected. Self-links are invalid. A relationship's identity is its ordered
endpoint pair and type when directional, or its unordered pair and type when
undirected; only one such edge is allowed. Different relation types between the
same equations remain distinct.

The workspace validates equation endpoints and source IDs against those
equations' source references. Candidate edges may have no source: automatic
feature matching is not a scientific source. `propose_relationships()` uses
only explicitly listed shared classifications and derivative operators. Every
suggestion has matched-feature reasoning and begins in
`candidate_requires_review`. Use `KnowledgeBase.review_relationship()` to
record an explicit accept/reject decision with reviewer and notes. Discovery
never promotes a candidate by itself.

```python
from newton_lab.knowledge import (
    RelationshipReviewStatus,
    build_example_knowledge_base,
    compare_equations,
    propose_relationships,
)

workspace = build_example_knowledge_base()
oscillator = workspace.get_entry("damped_harmonic_oscillator").record
pendulum = workspace.get_entry("damped_nonlinear_pendulum").record
comparison = compare_equations(oscillator, pendulum, knowledge_base=workspace)
suggestions = propose_relationships(workspace)

workspace = workspace.review_relationship(
    "small_angle_oscillator_approximation",
    RelationshipReviewStatus.ACCEPTED,
    reviewer="reviewer:physics",
    notes="Reviewed the small-angle and parameter-matching conditions.",
)
```

## Persistence and scientific limits

The new record fields and workspace relationship list are optional additions
with defaults, so Phase 5 schema-version-1 JSON continues to load without data
migration. New saves remain version 1 and persist structures, graph edges, and
review decisions. Unsupported versions continue to fail clearly. This is an
additive schema extension; a future incompatible change must increment the
version and add an explicit tested migration.

Metadata comparison is not symbolic reasoning. Shared form, operators,
variables, or solution methods do not prove physical equivalence, shared
mechanism, or validity in another domain. The example small-angle approximation
is valid only within its stated small-angle and coefficient-matching
conditions. Human review and suitable scientific evidence are still required.
