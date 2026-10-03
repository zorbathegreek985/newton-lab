# Physics knowledge and equation registry

The `newton_lab.knowledge` package stores equation metadata separately from
solver implementations. Its first taxonomy and sample records are deliberately
small; they are an extensible starting point, not an exhaustive catalogue.
The existing oscillator and pendulum solvers remain in their own modules and
are represented as tested implementations without being marked scientifically
verified.

## Taxonomy

`PhysicsTaxonomy` contains stable-ID `PhysicsBranch` entries, each with
`Subfield` entries. The initial set covers classical and analytical mechanics;
oscillations and waves; thermodynamics and statistical mechanics;
electromagnetism and optics; quantum mechanics and information; relativity;
gravitation and astrophysics; fluid and continuum mechanics; nonlinear dynamics,
chaos, and complex systems; condensed matter and materials; plasma; nuclear and
particle physics; geophysics; biophysics; and computational, mathematical, and
interdisciplinary physics.

Projects may construct a taxonomy with additional branches. IDs use lowercase
words separated by underscores and remain stable when display names change.
Subfield IDs are unique across the taxonomy, and equation records must point to
a subfield inside their selected branch.

## Equation record

`EquationRecord` is a Pydantic model. It requires an explicit stable
`equation_id`, canonical name, branch, and mathematical expression. Alternative
names, symbol and quantity definitions, SI units, dimensions, parameters, state
variables, assumptions, validity range, limitations, initial and boundary
conditions, classifications, relationships, interpretations, applications,
sources, verification state, implementation state, numerical methods, and tags
are optional.

The mathematical expression is descriptive text, not parsed executable math.
The schema can describe a formula, a system of equations, or a model; it does
not assert that the expression is solvable or ready to simulate. Classification
is a tuple so a record can, for example, be both a PDE and a field formulation,
or both nonlinear and deterministic.

## Evidence and provenance

Each `EvidenceClaim` carries one explicit category:

| Category | Meaning |
| --- | --- |
| `established_physical_result_or_application` | Established result or application, with source provenance required |
| `mathematically_derived_relationship` | A derivation under stated assumptions |
| `known_analogy_or_structural_similarity` | Similar mathematical structure, without asserting physical equivalence |
| `plausible_untested_application_hypothesis` | A proposal that has not been tested |
| `empirically_investigated_application` | Application investigated using evidence, with a source required |
| `independently_validated_result` | Result independently validated, with a source required |

`SourceReference` stores a stable local ID plus a title and available
bibliographic locator such as a DOI, URL, year, journal, or archive locator.
Claims refer to these local IDs, and the registry checks that every reference
resolves. The initial examples link Newton's *Principia* to the [Newton
Project catalogue](https://newtonproject.ox.ac.uk/catalogue/record/NATP00071),
Maxwell's 1865 paper to its [Royal Society DOI](https://doi.org/10.1098/rstl.1865.0008),
and Schrödinger's first 1926 communication to its [publisher DOI](https://doi.org/10.1002/andp.19263840404).
Modern notation in a record is identified as such; a source is not claimed to
reproduce that notation verbatim.

`VerificationStatus` describes scientific review of the record. The separate
`ImplementationStatus` describes software status. Passing code tests only marks
an implementation as tested; it does not establish that its physical
assumptions model a particular experiment.

## Equation relationships and application mappings

`EquationRelation` connects records using a named relation such as `derived_from`,
`limit_of`, `special_case_of`, shared assumptions, shared symmetries, or shared
solution methods. A relationship includes its own evidence claim. Registry
construction checks that target equation IDs exist.

`ApplicationMapping` identifies a domain and a proposed application, then
requires the evidence category to match its mapping level:

- `hypothesis` is an untested proposal.
- `structural_similarity` records an analogy only.
- `mathematical_mapping` records a mathematical derivation or mapping.
- `physically_justified_model` requires an established physical claim.
- `empirically_investigated` and `independently_validated` require their
  corresponding evidence categories and provenance.

These levels are not automatic steps. A formal similarity does not establish a
physical mechanism, and empirical investigation does not by itself establish
independent validation. Financial mappings in particular must not be presented
as predictive or profitable without appropriate empirical evidence.

For example, an exploratory cross-domain idea remains labeled as a hypothesis:

```python
from newton_lab.knowledge.registry import (
    ApplicationMapping,
    ApplicationMappingLevel,
    EvidenceCategory,
    EvidenceClaim,
)

ApplicationMapping(
    domain="complex systems",
    application="exploratory relaxation-model comparison",
    level=ApplicationMappingLevel.HYPOTHESIS,
    evidence=EvidenceClaim(
        category=EvidenceCategory.APPLICATION_HYPOTHESIS,
        statement="A proposed analogy requiring a domain-specific model and tests.",
    ),
)
```

## Using the initial records

Build the small example registry and inspect records:

```python
from newton_lab.knowledge.examples import build_example_registry

registry = build_example_registry()
newton_law = registry.get("newton_second_law")
quantum_records = registry.search("Schrödinger")

print(newton_law.mathematical_expression)
print(quantum_records[0].source_references[0].doi)
```

The examples include Newton's second law, Maxwell's equations in modern vacuum
SI differential form, the one-dimensional time-independent Schrödinger equation,
and registry entries for Newton Lab's existing oscillator and pendulum
implementations. The simulator entries are labeled `tested` in implementation
status but remain `unassessed` scientifically.

## Adding an equation

1. Choose a stable `equation_id` and an existing branch/subfield, or extend the
   taxonomy with a new branch and subfield.
2. Record the expression and its conventions, symbols, units, assumptions, and
   validity limits. Leave conditions or metadata empty when they do not apply.
3. Add claim categories and source references that support each claim. Do not
   use an analogy or a hypothesis as evidence of physical equivalence.
4. Add relationships only when their target records are included in the same
   registry build; unresolved references are rejected.
5. Set scientific verification independently from implementation status.
6. Validate and serialize with Pydantic:

```python
from newton_lab.knowledge.registry import EquationRegistry

registry = build_example_registry()
serialized = registry.model_dump_json(indent=2)
restored = EquationRegistry.model_validate_json(serialized)
```

`EquationRegistry.register(record)` returns a new validated registry and leaves
the original unchanged. To add several mutually related records, create a new
registry with their records together so all references can be checked at once.

## Future growth

The registry can provide metadata for future simulation adapters, source review
tools, mathematical-structure searches, and evidence tracking. Those systems
should consume this schema without moving physical equations into a generic
solver framework prematurely. Automated discovery may propose links, but such
links should enter as labeled hypotheses until reviewed and supported by
appropriate derivation or empirical work.
