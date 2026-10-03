"""Small, sourced examples for exercising the equation registry."""

from newton_lab.knowledge.registry import (
    EquationRecord,
    EquationRegistry,
    EquationRelation,
    EvidenceCategory,
    EvidenceClaim,
    ImplementationStatus,
    Linearity,
    MathematicalStructure,
    QuantityDefinition,
    RelationshipType,
    SourceReference,
    SourceType,
    VerificationStatus,
)
from newton_lab.knowledge.registry import (
    MathematicalClassification as Classification,
)
from newton_lab.knowledge.taxonomy import (
    PhysicsTaxonomy,
    default_physics_taxonomy,
)
from newton_lab.knowledge.workflow import (
    EquationRelationship,
    KnowledgeBase,
    RelationshipReviewStatus,
)


def example_equations() -> tuple[EquationRecord, ...]:
    """Return a few illustrative records, not an exhaustive physics catalogue."""
    newton_source = SourceReference(
        source_id="newton_principia_1687",
        title="Philosophiae Naturalis Principia Mathematica",
        authors=("Isaac Newton",),
        year=1687,
        publisher_or_journal="London",
        url="https://newtonproject.ox.ac.uk/catalogue/record/NATP00071",
        locator="Axiomata Sive Leges Motus",
    )
    maxwell_source = SourceReference(
        source_id="maxwell_electromagnetic_field_1865",
        title="A dynamical theory of the electromagnetic field",
        authors=("James Clerk Maxwell",),
        year=1865,
        publisher_or_journal=(
            "Philosophical Transactions of the Royal Society of London"
        ),
        doi="10.1098/rstl.1865.0008",
        url="https://doi.org/10.1098/rstl.1865.0008",
    )
    schrodinger_source = SourceReference(
        source_id="schrodinger_quantisierung_1926_first",
        title="Quantisierung als Eigenwertproblem (Erste Mitteilung)",
        authors=("Erwin Schrödinger",),
        year=1926,
        publisher_or_journal="Annalen der Physik, 384(4), 361–376",
        doi="10.1002/andp.19263840404",
        url="https://doi.org/10.1002/andp.19263840404",
    )
    steady_heat_source = SourceReference(
        source_id="mit_ocw_intro_engineering_heat_transfer_2002",
        title="Introduction to Engineering Heat Transfer, Part 3",
        responsible_organization="MIT OpenCourseWare",
        year=2002,
        source_type=SourceType.OTHER,
        url=(
            "https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/"
            "87d9f4544b7fd64a77201382500d057c_10_part3.pdf"
        ),
        locator="Section 2.1, equations 2.12b-2.17",
        reliability_notes=(
            "Official MIT OpenCourseWare course notes derive the one-dimensional "
            "steady conduction divergence equation, show the constant-property "
            "reduction for a plane slab, and explain the two fixed-temperature "
            "boundary conditions leading to a linear profile."
        ),
    )
    heat_transfer_source = SourceReference(
        source_id="mit_ocw_1d_thermal_diffusion_2003",
        title="1-D Thermal Diffusion Equation and Solutions",
        responsible_organization="MIT OpenCourseWare",
        year=2003,
        source_type=SourceType.OTHER,
        url=(
            "https://ocw.mit.edu/courses/3-185-transport-phenomena-in-"
            "materials-engineering-fall-2003/"
            "927ddb6b3dfc423de5570a0070614129_handout_htrans.pdf"
        ),
        locator=(
            "1-D Heat Conduction Solutions, sections 1(a)(i) and 2; "
            "constant-property Cartesian transient equation"
        ),
        reliability_notes=(
            "Official MIT OpenCourseWare course handout. It states the "
            "source-free 1D steady equation, the constant-property transient "
            "equation, and alpha = k/(rho cp)."
        ),
    )

    newton_law = EquationRecord(
        equation_id="newton_second_law",
        canonical_name="Newton's second law",
        alternative_names=("equation of motion",),
        branch_id="classical_mechanics",
        subfield_id="newtonian_mechanics",
        mathematical_expression="ΣF = dp/dt; for constant mass, ΣF = m a",
        description="Relates net force to the time rate of change of momentum.",
        symbol_definitions=(
            {
                "symbol": "ΣF",
                "name": "net force",
                "meaning": "Vector sum of forces on the body",
                "si_unit": "N",
                "dimensions": "kg·m·s^-2",
            },
            {
                "symbol": "p",
                "name": "linear momentum",
                "meaning": "Mass times velocity for a classical point mass",
                "si_unit": "kg·m/s",
                "dimensions": "kg·m·s^-1",
            },
        ),
        assumptions=("The reference frame is inertial.",),
        validity_domain=(
            "Classical mechanics; the constant-mass form assumes fixed mass."
        ),
        applicability_limits=(
            "Relativistic regimes require relativistic momentum and dynamics.",
        ),
        mathematical_classifications=(
            Classification.ODE,
            Classification.VECTOR_FORMULATION,
            Classification.DETERMINISTIC,
        ),
        physical_interpretations=(
            EvidenceClaim(
                category=EvidenceCategory.ESTABLISHED_PHYSICAL_RESULT,
                statement=(
                    "Newton's Principia presents the laws of motion underlying "
                    "classical inertial dynamics."
                ),
                source_reference_ids=(newton_source.source_id,),
            ),
        ),
        source_references=(newton_source,),
        verification_status=VerificationStatus.SOURCES_REVIEWED,
        tags=("force", "momentum", "classical mechanics"),
    )

    maxwell_equations = EquationRecord(
        equation_id="maxwell_equations_vacuum_si",
        canonical_name="Maxwell equations in vacuum (SI differential form)",
        alternative_names=("classical electromagnetic field equations",),
        branch_id="electromagnetism_and_optics",
        subfield_id="electromagnetism",
        mathematical_expression=(
            "∇·E = ρ/ε₀; ∇·B = 0; ∇×E = −∂B/∂t; ∇×B = μ₀J + μ₀ε₀∂E/∂t"
        ),
        description=(
            "A modern SI differential-form presentation of the classical "
            "electromagnetic field equations."
        ),
        symbol_definitions=(
            {
                "symbol": "E",
                "name": "electric field",
                "meaning": "Electric field vector",
                "si_unit": "V/m",
            },
            {
                "symbol": "B",
                "name": "magnetic flux density",
                "meaning": "Magnetic field vector in SI",
                "si_unit": "T",
            },
            {
                "symbol": "ρ",
                "name": "charge density",
                "meaning": "Electric charge per volume",
                "si_unit": "C/m^3",
            },
            {
                "symbol": "J",
                "name": "current density",
                "meaning": "Electric current per area",
                "si_unit": "A/m^2",
            },
        ),
        assumptions=(
            "Classical field description in vacuum with specified charge and "
            "current sources.",
            "The displayed equations use modern SI notation.",
        ),
        mathematical_classifications=(
            Classification.PDE,
            Classification.VECTOR_FORMULATION,
            Classification.FIELD_FORMULATION,
            Classification.LINEAR,
        ),
        physical_interpretations=(
            EvidenceClaim(
                category=EvidenceCategory.ESTABLISHED_PHYSICAL_RESULT,
                statement=(
                    "Maxwell's 1865 paper develops a dynamical theory of the "
                    "electromagnetic field; this record displays a modern SI form."
                ),
                source_reference_ids=(maxwell_source.source_id,),
            ),
        ),
        source_references=(maxwell_source,),
        verification_status=VerificationStatus.SOURCES_REVIEWED,
        tags=("electromagnetism", "field equations", "PDE"),
    )

    schrodinger_equation = EquationRecord(
        equation_id="schrodinger_equation_stationary_1d",
        canonical_name="Time-independent Schrödinger equation (one-dimensional)",
        alternative_names=("stationary Schrödinger equation",),
        branch_id="quantum_physics",
        subfield_id="quantum_mechanics",
        mathematical_expression=("[−(ℏ²/(2m)) d²/dx² + V(x)] ψ(x) = E ψ(x)"),
        description=(
            "Stationary states and energy eigenvalues for a one-dimensional "
            "nonrelativistic particle in a specified potential."
        ),
        symbol_definitions=(
            {
                "symbol": "ψ",
                "name": "wavefunction",
                "meaning": "Complex probability-amplitude field",
            },
            {
                "symbol": "ℏ",
                "name": "reduced Planck constant",
                "meaning": "Planck constant divided by 2π",
                "si_unit": "J·s",
            },
            {
                "symbol": "V",
                "name": "potential energy",
                "meaning": "Time-independent potential-energy function",
                "si_unit": "J",
            },
            {
                "symbol": "E",
                "name": "energy eigenvalue",
                "meaning": "Energy associated with a stationary state",
                "si_unit": "J",
            },
        ),
        parameters=(
            {
                "symbol": "m",
                "name": "particle mass",
                "meaning": "Mass parameter in the nonrelativistic kinetic operator",
                "si_unit": "kg",
            },
        ),
        state_variables=(
            {
                "symbol": "ψ(x)",
                "name": "wavefunction",
                "meaning": "Stationary quantum state in position space",
            },
        ),
        assumptions=(
            "Nonrelativistic quantum mechanics for a single particle.",
            "The potential V(x) and spatial boundary conditions must be specified.",
        ),
        applicability_limits=(
            "Relativistic particle dynamics and quantum field theory require "
            "broader models.",
        ),
        boundary_conditions=(
            "Boundary conditions and normalizability depend on the spatial domain.",
        ),
        mathematical_classifications=(
            Classification.ODE,
            Classification.EIGENVALUE_PROBLEM,
            Classification.OPERATOR_FORMULATION,
            Classification.LINEAR,
            Classification.DETERMINISTIC,
        ),
        physical_interpretations=(
            EvidenceClaim(
                category=EvidenceCategory.ESTABLISHED_PHYSICAL_RESULT,
                statement=(
                    "Schrödinger's 1926 first communication formulated the "
                    "wave-mechanics problem as an eigenvalue problem; this "
                    "record gives a modern one-dimensional differential form."
                ),
                source_reference_ids=(schrodinger_source.source_id,),
            ),
        ),
        source_references=(schrodinger_source,),
        verification_status=VerificationStatus.SOURCES_REVIEWED,
        tags=("quantum mechanics", "wavefunction", "eigenvalue problem"),
    )

    oscillator = EquationRecord(
        equation_id="damped_harmonic_oscillator",
        canonical_name="Damped harmonic oscillator",
        branch_id="oscillations_and_waves",
        subfield_id="oscillations",
        mathematical_expression="m x'' + c x' + k x = 0",
        mathematical_structure=MathematicalStructure(
            ode_order=2,
            linearity=Linearity.LINEAR,
            autonomous=True,
            explicit_time_dependence=False,
            derivative_operators=("first_time_derivative", "second_time_derivative"),
            state_variables=("displacement_x", "velocity_x_dot"),
            parameters=("mass_m", "damping_c", "stiffness_k"),
            forcing_terms=(),
            conservation_laws=(),
            symmetries=(),
            analytical_solution_families=("exponential_or_damped_trigonometric",),
            applicable_numerical_methods=("adaptive_initial_value_ODE_solver",),
        ),
        description="Free linear oscillator with viscous damping.",
        symbol_definitions=(
            {
                "symbol": "x",
                "name": "displacement",
                "meaning": "Position from equilibrium",
                "si_unit": "m",
            },
            {
                "symbol": "x'",
                "name": "velocity",
                "meaning": "Time derivative of displacement",
                "si_unit": "m/s",
            },
        ),
        parameters=(
            {
                "symbol": "m",
                "name": "mass",
                "meaning": "Oscillating mass",
                "si_unit": "kg",
            },
            {
                "symbol": "c",
                "name": "damping coefficient",
                "meaning": "Linear viscous damping coefficient",
                "si_unit": "kg/s",
            },
            {
                "symbol": "k",
                "name": "stiffness",
                "meaning": "Linear spring stiffness",
                "si_unit": "N/m",
            },
        ),
        assumptions=(
            "Linear spring force, linear viscous damping, and constant parameters.",
            "No external driving force.",
        ),
        applicability_limits=(
            "The model does not include nonlinear stiffness or non-viscous friction.",
        ),
        mathematical_classifications=(
            Classification.ODE,
            Classification.LINEAR,
            Classification.DETERMINISTIC,
        ),
        related_equations=(
            EquationRelation(
                target_equation_id="newton_second_law",
                relationship=RelationshipType.DERIVED_FROM,
                evidence=EvidenceClaim(
                    category=EvidenceCategory.MATHEMATICALLY_DERIVED,
                    statement=(
                        "The model follows from one-dimensional force balance with "
                        "linear restoring and viscous damping forces."
                    ),
                ),
            ),
        ),
        verification_status=VerificationStatus.UNASSESSED,
        implementation_status=ImplementationStatus.TESTED,
        numerical_methods=("SciPy solve_ivp",),
        tags=("oscillator", "damping", "implemented model"),
    )

    pendulum = EquationRecord(
        equation_id="damped_nonlinear_pendulum",
        canonical_name="Damped nonlinear pendulum",
        branch_id="nonlinear_and_complex_systems",
        mathematical_structure=MathematicalStructure(
            ode_order=2,
            linearity=Linearity.NONLINEAR,
            autonomous=True,
            explicit_time_dependence=False,
            derivative_operators=("first_time_derivative", "second_time_derivative"),
            state_variables=("angle_theta", "angular_velocity_theta_dot"),
            parameters=("bob_mass_m", "length_L", "damping_b", "gravity_g"),
            forcing_terms=(),
            conservation_laws=(),
            symmetries=(),
            analytical_solution_families=(),
            applicable_numerical_methods=("adaptive_initial_value_ODE_solver",),
        ),
        subfield_id="nonlinear_dynamics",
        mathematical_expression="θ'' + b/(m L²) θ' + (g/L) sin(θ) = 0",
        description=(
            "Unforced pendulum with nonlinear gravitational torque and viscous damping."
        ),
        symbol_definitions=(
            {
                "symbol": "θ",
                "name": "angular displacement",
                "meaning": "Angle from downward vertical",
                "si_unit": "rad",
            },
            {
                "symbol": "θ'",
                "name": "angular velocity",
                "meaning": "Time derivative of angle",
                "si_unit": "rad/s",
            },
        ),
        parameters=(
            {
                "symbol": "m",
                "name": "bob mass",
                "meaning": "Point-mass bob",
                "si_unit": "kg",
            },
            {
                "symbol": "L",
                "name": "length",
                "meaning": "Fixed pendulum length",
                "si_unit": "m",
            },
            {
                "symbol": "b",
                "name": "damping coefficient",
                "meaning": "Viscous torque coefficient",
                "si_unit": "kg·m²/s",
            },
            {
                "symbol": "g",
                "name": "gravitational acceleration",
                "meaning": "Uniform gravitational field strength",
                "si_unit": "m/s²",
            },
        ),
        assumptions=(
            "Point-mass bob, massless rigid rod, fixed pivot and length.",
            "Uniform gravity, linear viscous damping torque, and no external torque.",
        ),
        applicability_limits=(
            "Rod flexibility and damping beyond the specified torque law are "
            "not modeled.",
        ),
        initial_conditions=("Initial angle and angular velocity are required.",),
        mathematical_classifications=(
            Classification.ODE,
            Classification.NONLINEAR,
            Classification.DETERMINISTIC,
        ),
        related_equations=(
            EquationRelation(
                target_equation_id="newton_second_law",
                relationship=RelationshipType.DERIVED_FROM,
                evidence=EvidenceClaim(
                    category=EvidenceCategory.MATHEMATICALLY_DERIVED,
                    statement=(
                        "The rotational equation follows from torque balance for "
                        "the stated point-mass pendulum assumptions."
                    ),
                ),
            ),
        ),
        verification_status=VerificationStatus.UNASSESSED,
        implementation_status=ImplementationStatus.TESTED,
        numerical_methods=("SciPy solve_ivp",),
        tags=("pendulum", "nonlinear dynamics", "implemented model"),
    )

    steady_heat_divergence = EquationRecord(
        equation_id="steady_heat_conduction_divergence_1d",
        canonical_name="Steady one-dimensional heat conduction (divergence form)",
        alternative_names=("source-free steady heat conduction",),
        branch_id="thermodynamics_and_statistical_mechanics",
        subfield_id="thermodynamics",
        mathematical_expression="d/dx(k(x) dT/dx) = 0",
        mathematical_structure=MathematicalStructure(
            ode_order=2,
            linearity=Linearity.LINEAR,
            autonomous=True,
            explicit_time_dependence=False,
            derivative_operators=("second_spatial_derivative",),
            state_variables=("temperature_field_T",),
            parameters=("thermal_conductivity_field_k",),
            forcing_terms=(),
            analytical_solution_families=("integrated_flux_form",),
        ),
        description=(
            "One-dimensional source-free steady conduction in divergence form. "
            "The conductivity may vary with position; no constant-conductivity "
            "reduction is implied by this record."
        ),
        symbol_definitions=(
            QuantityDefinition(
                symbol="T(x)",
                name="temperature field",
                meaning="Temperature along the one-dimensional spatial domain",
                si_unit="K",
            ),
            QuantityDefinition(
                symbol="x",
                name="position",
                meaning="Cartesian coordinate along the conducting body",
                si_unit="m",
            ),
            QuantityDefinition(
                symbol="k(x)",
                name="thermal conductivity field",
                meaning="Local material thermal conductivity",
                si_unit="W/(m K)",
            ),
        ),
        parameters=(
            QuantityDefinition(
                symbol="k(x)",
                name="thermal conductivity field",
                meaning="Positive conductivity, potentially varying with position",
                si_unit="W/(m K)",
            ),
        ),
        assumptions=(
            "Steady state and one spatial Cartesian dimension.",
            "No internal heat generation.",
            "Continuum Fourier conduction; conductivity is positive and may vary "
            "with position.",
        ),
        validity_domain=(
            "A one-dimensional material segment with specified boundary conditions."
        ),
        applicability_limits=(
            "The equation does not specify boundary conditions or a particular "
            "temperature profile.",
            "Contact resistance, multidimensional effects, and non-Fourier "
            "transport are not represented.",
        ),
        mathematical_classifications=(
            Classification.ODE,
            Classification.LINEAR,
            Classification.DETERMINISTIC,
            Classification.CONSERVATION_LAW,
        ),
        physical_interpretations=(
            EvidenceClaim(
                category=EvidenceCategory.MATHEMATICALLY_DERIVED,
                statement=(
                    "For source-free 1D steady conduction, conservation of energy "
                    "and Fourier's law give d/dx(k(x) dT/dx)=0."
                ),
                source_reference_ids=(steady_heat_source.source_id,),
            ),
        ),
        source_references=(steady_heat_source,),
        verification_status=VerificationStatus.SOURCES_REVIEWED,
        implementation_status=ImplementationStatus.SPECIFIED,
        tags=("heat conduction", "steady state", "divergence form"),
    )

    steady_heat_constant = EquationRecord(
        equation_id="steady_heat_conduction_1d",
        canonical_name="Steady one-dimensional heat conduction (constant conductivity)",
        alternative_names=("linear steady heat profile equation",),
        branch_id="thermodynamics_and_statistical_mechanics",
        subfield_id="thermodynamics",
        mathematical_expression="d^2T/dx^2 = 0",
        mathematical_structure=MathematicalStructure(
            ode_order=2,
            linearity=Linearity.LINEAR,
            autonomous=True,
            explicit_time_dependence=False,
            derivative_operators=("second_spatial_derivative",),
            state_variables=("temperature_field_T",),
            parameters=("constant_thermal_conductivity_k",),
            forcing_terms=(),
            analytical_solution_families=("linear_profile",),
            applicable_numerical_methods=("two_point_boundary_value_solver",),
        ),
        description=(
            "Constant-conductivity reduction of source-free 1D steady heat "
            "conduction; it is distinct from the general divergence-form record."
        ),
        symbol_definitions=steady_heat_divergence.symbol_definitions,
        parameters=(
            QuantityDefinition(
                symbol="k",
                name="thermal conductivity",
                meaning="Positive constant thermal conductivity",
                si_unit="W/(m K)",
            ),
        ),
        state_variables=(steady_heat_divergence.symbol_definitions[0],),
        assumptions=(
            "Steady state and one spatial Cartesian dimension.",
            "No internal heat generation.",
            "Homogeneous material with constant positive thermal conductivity.",
        ),
        validity_domain=(
            "A one-dimensional homogeneous material segment with specified "
            "endpoint temperatures."
        ),
        applicability_limits=(
            "Boundary conditions select the particular profile.",
            "Variable conductivity, internal generation, and multidimensional "
            "effects require a different equation form.",
        ),
        boundary_conditions=(
            "For fixed endpoint temperatures T(0)=T_left and T(L)=T_right, "
            "the solution is linear in x.",
        ),
        mathematical_classifications=(
            Classification.ODE,
            Classification.LINEAR,
            Classification.DETERMINISTIC,
        ),
        related_equations=(
            EquationRelation(
                target_equation_id="steady_heat_conduction_divergence_1d",
                relationship=RelationshipType.SPECIAL_CASE_OF,
                evidence=EvidenceClaim(
                    category=EvidenceCategory.MATHEMATICALLY_DERIVED,
                    statement=(
                        "Setting k(x)=k for constant k in the divergence-form "
                        "equation yields k T_xx=0, hence T_xx=0 for k>0."
                    ),
                    source_reference_ids=(steady_heat_source.source_id,),
                ),
            ),
        ),
        physical_interpretations=(
            EvidenceClaim(
                category=EvidenceCategory.MATHEMATICALLY_DERIVED,
                statement=(
                    "With fixed endpoint temperatures, integrating T_xx=0 gives "
                    "the unique linear profile T(x)=T_left+(T_right-T_left)x/L."
                ),
                source_reference_ids=(steady_heat_source.source_id,),
            ),
        ),
        source_references=(steady_heat_source, heat_transfer_source),
        verification_status=VerificationStatus.SOURCES_REVIEWED,
        implementation_status=ImplementationStatus.TESTED,
        numerical_methods=("SciPy solve_bvp",),
        tags=("heat conduction", "steady state", "constant conductivity"),
    )

    transient_heat_diffusion = EquationRecord(
        equation_id="transient_heat_diffusion_1d",
        canonical_name="Transient one-dimensional heat diffusion",
        alternative_names=("one-dimensional heat equation",),
        branch_id="thermodynamics_and_statistical_mechanics",
        subfield_id="thermodynamics",
        mathematical_expression="partial T/partial t = alpha partial^2 T/partial x^2",
        mathematical_structure=MathematicalStructure(
            pde_order=2,
            linearity=Linearity.LINEAR,
            autonomous=True,
            explicit_time_dependence=False,
            derivative_operators=(
                "first_time_derivative",
                "second_spatial_derivative",
            ),
            state_variables=("temperature_field_T",),
            parameters=("thermal_diffusivity_alpha",),
            forcing_terms=(),
            analytical_solution_families=("separation_of_variables_sine_modes",),
        ),
        description=(
            "Source-free Cartesian heat diffusion in a homogeneous medium with "
            "constant thermal properties."
        ),
        symbol_definitions=(
            *steady_heat_divergence.symbol_definitions[:2],
            QuantityDefinition(
                symbol="t",
                name="time",
                meaning="Elapsed time",
                si_unit="s",
            ),
            QuantityDefinition(
                symbol="alpha",
                name="thermal diffusivity",
                meaning="Thermal conductivity divided by volumetric heat capacity",
                si_unit="m^2/s",
            ),
        ),
        parameters=(
            QuantityDefinition(
                symbol="alpha",
                name="thermal diffusivity",
                meaning="alpha = k/(rho cp) for constant material properties",
                si_unit="m^2/s",
                dimensions="L^2 T^-1",
            ),
            QuantityDefinition(
                symbol="k",
                name="thermal conductivity",
                meaning="Constant positive thermal conductivity",
                si_unit="W/(m K)",
            ),
            QuantityDefinition(
                symbol="rho",
                name="mass density",
                meaning="Constant material mass density",
                si_unit="kg/m^3",
            ),
            QuantityDefinition(
                symbol="cp",
                name="specific heat capacity",
                meaning="Constant specific heat capacity",
                si_unit="J/(kg K)",
            ),
        ),
        state_variables=(steady_heat_divergence.symbol_definitions[0],),
        assumptions=(
            "One Cartesian spatial dimension and a homogeneous continuum.",
            "Thermal conductivity, density, and specific heat are constant in "
            "space and time.",
            "No internal heat generation or bulk advection.",
        ),
        validity_domain=(
            "A one-dimensional body with a specified initial temperature field "
            "and appropriate boundary conditions."
        ),
        applicability_limits=(
            "Changing properties, internal sources, multidimensional transport, "
            "and complicated boundary histories are not represented.",
        ),
        initial_conditions=("An initial temperature field T(x,0) is required.",),
        boundary_conditions=(
            "Boundary conditions are required; the Phase 19 case uses fixed "
            "endpoint temperatures.",
        ),
        mathematical_classifications=(
            Classification.PDE,
            Classification.LINEAR,
            Classification.DETERMINISTIC,
        ),
        related_equations=(
            EquationRelation(
                target_equation_id="steady_heat_conduction_1d",
                relationship=RelationshipType.LIMIT_OF,
                evidence=EvidenceClaim(
                    category=EvidenceCategory.MATHEMATICALLY_DERIVED,
                    statement=(
                        "Setting partial T/partial t=0 gives T_xx=0 when alpha is "
                        "constant and positive; with alpha=k/(rho cp), this is the "
                        "source-free constant-conductivity steady conduction "
                        "equation. Fixed endpoint temperatures then select a linear "
                        "profile, but boundary conditions are needed to select it."
                    ),
                    source_reference_ids=(
                        steady_heat_source.source_id,
                        heat_transfer_source.source_id,
                    ),
                ),
            ),
        ),
        physical_interpretations=(
            EvidenceClaim(
                category=EvidenceCategory.MATHEMATICALLY_DERIVED,
                statement=(
                    "For a finite interval with fixed endpoint temperatures, "
                    "subtracting the linear steady profile produces homogeneous "
                    "Dirichlet conditions for the transient perturbation."
                ),
                source_reference_ids=(heat_transfer_source.source_id,),
            ),
        ),
        source_references=(steady_heat_source, heat_transfer_source),
        verification_status=VerificationStatus.SOURCES_REVIEWED,
        implementation_status=ImplementationStatus.TESTED,
        numerical_methods=("analytical finite sine-mode solution",),
        tags=("heat diffusion", "transient", "PDE", "thermal diffusivity"),
    )

    return (
        newton_law,
        maxwell_equations,
        schrodinger_equation,
        oscillator,
        pendulum,
        steady_heat_divergence,
        steady_heat_constant,
        transient_heat_diffusion,
    )


def build_example_registry(
    taxonomy: PhysicsTaxonomy | None = None,
) -> EquationRegistry:
    """Build a small demonstration registry with sourced, scoped records."""
    return EquationRegistry(
        taxonomy=taxonomy or default_physics_taxonomy(),
        equations=example_equations(),
    )


def build_example_knowledge_base() -> KnowledgeBase:
    """Build examples with an unreviewed small-angle approximation proposal."""
    registry = build_example_registry()
    entries = KnowledgeBase.from_equation_registry(registry).entries
    approximation = EquationRelationship(
        relationship_id="small_angle_oscillator_approximation",
        source_equation_id="damped_harmonic_oscillator",
        target_equation_id="damped_nonlinear_pendulum",
        relationship_type=RelationshipType.APPROXIMATION_OF,
        description=(
            "The linear oscillator approximates the nonlinear pendulum after "
            "mapping x = L theta."
        ),
        justification=(
            "Replacing sin(theta) with theta gives a linear second-order ODE "
            "when damping and restoring coefficients are matched."
        ),
        conditions=(
            "Absolute pendulum angle is much less than one radian.",
            "Use x = L theta and match oscillator coefficients to the pendulum.",
        ),
        directional=True,
        review_status=RelationshipReviewStatus.CANDIDATE,
    )
    return KnowledgeBase(
        taxonomy=registry.taxonomy,
        entries=entries,
        relationships=(approximation,),
    )
