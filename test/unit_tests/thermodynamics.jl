using HallThruster: HallThruster as het

function test_property(property, laws)
    initval = 0.0
    for (i, (U, f)) in enumerate(laws)
        if i == 1
            initval = property(U, f)
        else
            if !isapprox(property(U, f), initval)
                return false
            end
        end
    end
    return true
end

function test_conservation_laws()
    @testset "Conservation laws" begin
        Xe_0 = het.Xenon(0)
        @test het.nvars(het.ContinuityOnly(Xe_0; u = 300, T = 300)) == 1
        @test het.nvars(het.IsothermalEuler(Xe_0; T = 300)) == 2
        @test het.nvars(het.EulerEquations(Xe_0)) == 3
    end
end

function test_properties()
    ρ = 1.0
    T = 300
    u = 300
    ϵ = het.Xenon.cv * T + 0.5 * u^2
    Xe_0 = het.Xenon(0)

    continuity_eq = het.Fluid(Xe_0; u, T)
    continuity_state = (ρ,)
    continuity = (continuity_state, continuity_eq)

    isothermal_eq = het.Fluid(Xe_0; T)
    isothermal_state = (ρ, ρ * u)
    isothermal = (isothermal_state, isothermal_eq)

    euler_eq = het.Fluid(Xe_0)
    euler_state = (ρ, ρ * u, ρ * ϵ)
    euler = (euler_state, euler_eq)

    laws = [continuity, isothermal, euler]
    laws_vector = [continuity, isothermal, euler]

    @testset "Thermodynamic property computation" begin
        # Check that thermodynamic property computations give identical
        # results for the different fluid types
        @test test_property(het.temperature, laws)
        @test test_property(het.pressure, laws)
        @test test_property(het.density, laws)
        @test test_property(het.number_density, laws)
        @test test_property(het.velocity, laws)
        @test test_property(het.sound_speed, laws)

        # check that versions which work for vectors instead of staticarrays work properly
        @test test_property(het.temperature, laws_vector)
        @test test_property(het.pressure, laws_vector)
        @test test_property(het.density, laws_vector)
        @test test_property(het.number_density, laws_vector)
        @test test_property(het.velocity, laws_vector)
        @test test_property(het.sound_speed, laws_vector)

        # Check that properties are being computed correctly
        @test het.temperature(continuity...) ≈ T
        @test het.velocity(continuity...) ≈ u
        @test het.number_density(continuity...) ≈ ρ / het.m(continuity_eq)
        @test het.density(continuity...) ≈ ρ
        @test het.pressure(continuity...) ≈ ρ * het.Xenon.R * T
        @test het.sound_speed(continuity...) ≈ √(het.γ(continuity_eq) * het.Xenon.R * T)
    end
end

function test_thermodynamics()
    @testset "Thermodynamics/Conservation laws" begin
        test_properties()
        test_conservation_laws()
    end
end
