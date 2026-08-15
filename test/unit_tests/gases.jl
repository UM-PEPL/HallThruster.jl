using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

@testset "Serialization" begin
    test_instances(het.Gas, (; het.Xenon, het.Krypton, het.Argon))

    for (name, gas) in pairs(het.propellants_v0_21_7)
        @test het.deserialize(het.Gas, name) == gas
    end
end

@testset "Gas and species" begin
    @test repr(het.Krypton) == "Kr"
    @test repr(het.Species(het.Xenon, 1)) == "Xe(+)"
    @test repr(het.Species(het.Xenon, 3)) == "Xe(3+)"
    @test repr(het.Species(het.Xenon, 0)) == "Xe"
    @test repr(het.Species(het.MolecularNitrogen, 1)) == "N2(+)"
    @test repr(het.Species(het.MolecularNitrogen, 0)) == "N2"
end

@testset "Molecules" begin
    # Basic parsing
    components = het.parse_chemical_formula("Ca(OH)2")
    info = het.molecule_info(components)

    expected_mass = het.ELEMENTS[:Ca].mass + 2 * (het.ELEMENTS[:O].mass + het.ELEMENTS[:H].mass)
    @test info.mass == expected_mass
    @test info.num_atoms == 5

    # no gamma specified
    @test_throws(ErrorException, het.Gas("Ca(OH)2"))

    # gamma specified
    CaOH2 = het.Gas("Ca(OH)2", γ = 1.5)
    @test CaOH2.M == expected_mass
    @test CaOH2.γ == 1.5

    co2_formula = het.parse_chemical_formula("CO2")
    info = het.molecule_info(co2_formula)
    @test info.mass > 44
    @test info.num_atoms == 3
    CO2 = het.Gas("CO2", γ = 1.28)
    @test CO2.M == het.CarbonDioxide.M

    # Complex parsing
    components = het.parse_chemical_formula("A4(B3(C2(D1)3)2)1")
    @test length(components) == 2
    comp1, comp2 = components
    @test comp1 == het.ElementTerm(:A, 4)

    comp2 = components[2]
    @test comp2 isa het.MoleculeTerm
    @test comp2.count == 1

    @test length(comp2.components) == 2
    comp21, comp22 = comp2.components
    @test comp21 == het.ElementTerm(:B, 3)

    @test comp22 isa het.MoleculeTerm
    @test comp22.count == 2
    @test length(comp22.components) == 2

    comp221, comp222 = comp22.components
    @test comp221 == het.ElementTerm(:C, 2)
    @test comp222 == het.ElementTerm(:D, 3)

    @test het.molecule_formula(components) == "A4B3(C2D3)2"

    M = 5.0
    Î³ = 1.0
    gas = het.Gas("Fake", "Fa"; Î³, M)
    @test repr(gas) == "Fake"
    @test gas.m == M / het.NA
    @test gas(2) == het.Species(gas, 2)
    # Redundancy checking
    components = het.parse_chemical_formula("(C)1((O002)1)")
    @test het.molecule_formula(components) == "CO2"
end

@testset "Excited state species" begin
    # Construction and display
    @test repr(het.Species(het.Xenon, 0, 1)) == "Xe(*)"
    @test repr(het.Species(het.Xenon, 0, 2)) == "Xe(2*)"
    @test repr(het.Species(het.MolecularNitrogen, 0, 1)) == "N2(*)"
    @test repr(het.Species(het.Xenon, 1, 1)) == "Xe(1+,1*)"

    # Convenience constructor
    @test het.Xenon(0, 1) == het.Species(het.Xenon, 0, 1)

    # Default excitation level is ground state, so old two-arg behavior is preserved
    @test het.Species(het.Xenon, 0).excited_level == 0
    @test het.Xenon(0) == het.Species(het.Xenon, 0, 0)
    @test repr(het.Species(het.Xenon, 0)) == "Xe"

    # Distinct levels are distinct species
    @test het.Xenon(0, 1) != het.Xenon(0)
    @test het.Xenon(0, 1) != het.Xenon(0, 2)

    # Helpers
    @test !het.is_excited(het.Xenon(0))
    @test het.is_excited(het.Xenon(0, 1))
    @test het.ground_state(het.Xenon(0, 2)) == het.Xenon(0)
    @test het.ground_state(het.Xenon(1, 1)) == het.Xenon(1)

    # Negative levels are invalid
    @test_throws ErrorException het.Species(het.Xenon, 0, -1)
end

@testset "Excited state allocation" begin
    ncells = 17

    propellant = het.Propellant(
        het.Xenon, 5.0e-6;
        max_charge = 2,
        velocity_m_s = 300.0,
        temperature_K = 500.0,
    )

    @testset "No excited levels (default, lumped behavior)" begin
        fluids = het.allocate_fluids(propellant, ncells)
        @test length(fluids.continuity) == 1
        @test length(fluids.isothermal) == 2
        @test isempty(het.excited_fluids(fluids))
        @test het.ground_neutral(fluids) === fluids.continuity[1]
        @test het.ground_neutral(fluids).species == het.Xenon(0)
    end

    @testset "Per-channel excitation" begin
        excited_levels = [2, 1]  # deliberately unsorted
        fluids = het.allocate_fluids(propellant, ncells; excited_levels)

        # One continuity fluid per excited level, plus the ground state
        @test length(fluids.continuity) == 3
        # Ion fluids unaffected
        @test length(fluids.isothermal) == 2
        @test all(f -> f.species.Z > 0, fluids.isothermal)

        # Ground state first, then excited levels in sorted order
        @test het.ground_neutral(fluids).species == het.Xenon(0)
        @test fluids.continuity[2].species == het.Xenon(0, 1)
        @test fluids.continuity[3].species == het.Xenon(0, 2)

        excited = het.excited_fluids(fluids)
        @test length(excited) == 2
        @test all(f -> het.is_excited(f.species), excited)

        ground = het.ground_neutral(fluids)
        for fluid in excited
            # Excited neutrals are continuity-only and advect with the
            # ground-state neutral background
            @test fluid.type == het._ContinuityOnly
            @test fluid.const_velocity == propellant.velocity_m_s
            @test fluid.const_velocity == ground.const_velocity
            @test fluid.sound_speed == ground.sound_speed

            # Same element and charge state as the ground neutral
            @test het.ground_state(fluid.species) == ground.species

            # Correctly-sized state arrays
            @test size(fluid.density) == (ncells + 2,)
            @test size(fluid.dens_ddt) == (ncells + 2,)
            @test size(fluid.flux_dens) == (ncells + 1,)
        end
    end

    @testset "Levels declared on Propellant" begin
        p = het.Propellant(het.Xenon, 5.0e-6; max_charge = 2, excited_levels = [2, 1])
        @test p.excited_levels == [1, 2]

        fluids = het.allocate_fluids(p, ncells)
        @test length(fluids.continuity) == 3
        @test fluids.continuity[2].species == het.Xenon(0, 1)
        @test fluids.continuity[3].species == het.Xenon(0, 2)

        # Default is no excited states
        @test het.Propellant(het.Xenon, 5.0e-6).excited_levels == Int[]

        # Level 0 is the ground state and is always present, so it is not a valid
        # entry; duplicates are rejected as for allowed_charges
        @test_throws ErrorException het.Propellant(het.Xenon, 5.0e-6; excited_levels = [0, 1])
        @test_throws ErrorException het.Propellant(het.Xenon, 5.0e-6; excited_levels = [1, 1])
    end
end

@testset "Excited state reaction parsing" begin
    # Neutral excitation uses a parenthesized level
    lhs, rhs = het._parse_reaction_equation("Xe + e -> Xe(*) + e")
    @test only(keys(lhs) |> collect |> t -> filter(x -> x.species != "e", t)).excited_level == 0
    @test only(keys(rhs) |> collect |> t -> filter(x -> x.species != "e", t)).excited_level == 1

    # Stepwise ionization out of an excited state balances charge
    lhs, rhs = het._parse_reaction_equation("Xe(*) + e -> Xe(+) + 2e")
    reactant = only(filter(x -> x.species != "e", collect(keys(lhs))))
    product = only(filter(x -> x.species != "e", collect(keys(rhs))))
    @test reactant.excited_level == 1
    @test reactant.charge == 0
    @test product.excited_level == 0
    @test product.charge == 1

    # Round-trip through the display string used for species lookup
    @test repr(reactant) == "Xe(*)"
    @test string(het.Xenon(0, 1)) == "Xe(*)"

    # Multiple levels
    lhs, _ = het._parse_reaction_equation("Xe(2*) + e -> Xe(*) + e")
    @test only(filter(x -> x.species != "e", collect(keys(lhs)))).excited_level == 2
end
