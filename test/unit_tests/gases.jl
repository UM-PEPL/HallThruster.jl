using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

@testset "Serialization" begin
    test_instances(het.Gas, het.propellants)
end
@testset "Gas and species" begin
    @test repr(het.Krypton) == "Krypton"
    @test repr(het.Species(het.Xenon, 1)) == "Xe(+)"
    @test repr(het.Species(het.Xenon, 3)) == "Xe(3+)"
    @test repr(het.Species(het.Xenon, 0)) == "Xe"
    @test repr(het.Species(het.MolecularNitrogen, 1)) == "N2(+)"
    @test repr(het.Species(het.MolecularNitrogen, 0)) == "N2"

    M = 5.0
    γ = 1.0
    gas = het.Gas("Fake", "Fa"; γ, M)
    @test repr(gas) == "Fake"
    @test gas.m == M / het.NA
    @test gas(2) == het.Species(gas, 2)
end
