using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

@testset "Serialization" begin
    test_instances(het.Gas, (;het.Xenon, het.Krypton, het.Argon))

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
    CaOH2 = het.Gas("Ca(OH)2", γ=1.5)
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

    # Redundancy checking
    components = het.parse_chemical_formula("(C)1((O002)1)")
    @test het.molecule_formula(components) == "CO2"
end 

