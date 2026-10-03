"""Extensible branch and subfield taxonomy for physics knowledge records."""

from typing import Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

_IDENTIFIER_PATTERN = r"^[a-z0-9]+(?:_[a-z0-9]+)*$"


class Subfield(BaseModel):
    """A named physics subfield with a stable machine-readable identifier."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    id: str = Field(pattern=_IDENTIFIER_PATTERN)
    name: str = Field(min_length=1)
    description: str | None = None


class PhysicsBranch(BaseModel):
    """A high-level physics branch containing zero or more subfields."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    id: str = Field(pattern=_IDENTIFIER_PATTERN)
    name: str = Field(min_length=1)
    description: str | None = None
    subfields: tuple[Subfield, ...] = ()


class PhysicsTaxonomy(BaseModel):
    """An extensible taxonomy; callers may supply additional branches."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    branches: tuple[PhysicsBranch, ...]

    @model_validator(mode="after")
    def identifiers_are_unique(self) -> Self:
        branch_ids = [branch.id for branch in self.branches]
        if len(branch_ids) != len(set(branch_ids)):
            raise ValueError("physics branch identifiers must be unique")

        subfield_ids = [
            subfield.id for branch in self.branches for subfield in branch.subfields
        ]
        if len(subfield_ids) != len(set(subfield_ids)):
            raise ValueError("physics subfield identifiers must be globally unique")
        return self

    def branch(self, branch_id: str) -> PhysicsBranch | None:
        """Return a branch by stable ID, or ``None`` if it is absent."""
        return next(
            (branch for branch in self.branches if branch.id == branch_id), None
        )

    def contains_subfield(self, branch_id: str, subfield_id: str) -> bool:
        """Return whether a subfield belongs to the specified branch."""
        branch = self.branch(branch_id)
        return branch is not None and any(
            subfield.id == subfield_id for subfield in branch.subfields
        )


def default_physics_taxonomy() -> PhysicsTaxonomy:
    """Build an initial broad taxonomy, intended to grow as needed."""
    return PhysicsTaxonomy(
        branches=(
            PhysicsBranch(
                id="classical_mechanics",
                name="Classical mechanics",
                subfields=(
                    Subfield(id="newtonian_mechanics", name="Newtonian mechanics"),
                    Subfield(id="analytical_mechanics", name="Analytical mechanics"),
                ),
            ),
            PhysicsBranch(
                id="oscillations_and_waves",
                name="Oscillations and waves",
                subfields=(
                    Subfield(id="oscillations", name="Oscillations"),
                    Subfield(id="wave_physics", name="Wave physics"),
                    Subfield(id="acoustics", name="Acoustics"),
                ),
            ),
            PhysicsBranch(
                id="thermodynamics_and_statistical_mechanics",
                name="Thermodynamics and statistical mechanics",
                subfields=(
                    Subfield(id="thermodynamics", name="Thermodynamics"),
                    Subfield(id="statistical_mechanics", name="Statistical mechanics"),
                ),
            ),
            PhysicsBranch(
                id="electromagnetism_and_optics",
                name="Electromagnetism and optics",
                subfields=(
                    Subfield(id="electromagnetism", name="Electromagnetism"),
                    Subfield(id="optics", name="Optics"),
                ),
            ),
            PhysicsBranch(
                id="quantum_physics",
                name="Quantum physics",
                subfields=(
                    Subfield(id="quantum_mechanics", name="Quantum mechanics"),
                    Subfield(id="quantum_information", name="Quantum information"),
                ),
            ),
            PhysicsBranch(
                id="relativity",
                name="Relativity",
                subfields=(
                    Subfield(id="special_relativity", name="Special relativity"),
                    Subfield(id="general_relativity", name="General relativity"),
                ),
            ),
            PhysicsBranch(
                id="gravitation_and_astrophysics",
                name="Gravitation and astrophysics",
                subfields=(
                    Subfield(id="gravitation", name="Gravitation"),
                    Subfield(id="astrophysics", name="Astrophysics"),
                ),
            ),
            PhysicsBranch(
                id="fluids_and_continuum_mechanics",
                name="Fluid mechanics and continuum mechanics",
                subfields=(
                    Subfield(id="fluid_mechanics", name="Fluid mechanics"),
                    Subfield(id="continuum_mechanics", name="Continuum mechanics"),
                ),
            ),
            PhysicsBranch(
                id="nonlinear_and_complex_systems",
                name="Nonlinear dynamics, chaos, and complex systems",
                subfields=(
                    Subfield(id="nonlinear_dynamics", name="Nonlinear dynamics"),
                    Subfield(id="chaos", name="Chaos"),
                    Subfield(id="complex_systems", name="Complex systems"),
                ),
            ),
            PhysicsBranch(
                id="condensed_matter_and_materials",
                name="Condensed matter and materials physics",
                subfields=(
                    Subfield(id="condensed_matter", name="Condensed matter"),
                    Subfield(id="materials_physics", name="Materials physics"),
                ),
            ),
            PhysicsBranch(
                id="plasma_physics",
                name="Plasma physics",
                subfields=(Subfield(id="plasma_dynamics", name="Plasma dynamics"),),
            ),
            PhysicsBranch(
                id="nuclear_and_particle_physics",
                name="Nuclear and particle physics",
                subfields=(
                    Subfield(id="nuclear_physics", name="Nuclear physics"),
                    Subfield(id="particle_physics", name="Particle physics"),
                ),
            ),
            PhysicsBranch(
                id="geophysics",
                name="Geophysics",
                subfields=(
                    Subfield(id="solid_earth_physics", name="Solid Earth physics"),
                    Subfield(id="physical_oceanography", name="Physical oceanography"),
                ),
            ),
            PhysicsBranch(
                id="biophysics",
                name="Biophysics",
                subfields=(
                    Subfield(id="molecular_biophysics", name="Molecular biophysics"),
                    Subfield(id="systems_biophysics", name="Systems biophysics"),
                ),
            ),
            PhysicsBranch(
                id="computational_and_interdisciplinary_physics",
                name="Computational and interdisciplinary physics",
                subfields=(
                    Subfield(id="computational_physics", name="Computational physics"),
                    Subfield(id="mathematical_physics", name="Mathematical physics"),
                    Subfield(
                        id="interdisciplinary_physics",
                        name="Interdisciplinary physics",
                    ),
                ),
            ),
        )
    )
