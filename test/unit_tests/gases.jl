using HallThruster: HallThruster as het

include("serialization_test_utils.jl")

function test_gases()
    @testset "Gases" begin
        @testset "Serialization" begin
            test_instances(het.Gas, het.propellants)
        end
        @testset "Gas and species" begin
            @test repr(het.Krypton) == "Krypton"
            @test repr(het.Electron) == "e-"
            @test repr(het.Species(het.Xenon, 1)) == "Xe+"
            @test repr(het.Species(het.Xenon, 3)) == "Xe3+"
            @test repr(het.Species(het.Xenon, 0)) == "Xe"

            M = 5.0
            γ = 1.0
            gas = het.Gas("Fake", "Fa"; γ, M)
            @test repr(gas) == "Fake"
            @test gas.m == M / het.NA
            @test gas.R == het.R0 / M
            @test gas.cp == γ / (γ - 1) * gas.R
            @test gas.cv == gas.cp - gas.R
            @test gas(2) == het.Species(gas, 2)
        end
    end
end