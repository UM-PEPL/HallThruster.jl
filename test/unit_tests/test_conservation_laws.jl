@testset "Fluxes and conservation laws" begin
    let Xenon = HallThruster.Xenon,
        Fluid = HallThruster.Fluid,
        ContinuityOnly = HallThruster.ContinuityOnly,
        IsothermalEuler = HallThruster.IsothermalEuler,
        EulerEquations = HallThruster.EulerEquations,
        temperature = HallThruster.temperature,
        pressure = HallThruster.pressure,
        density = HallThruster.density,
        number_density = HallThruster.number_density,
        velocity = HallThruster.velocity,
        sound_speed = HallThruster.sound_speed,
        m = HallThruster.m,
        γ = HallThruster.γ,
        R = HallThruster.R,
        cp = HallThruster.cp,
        cv = HallThruster.cv,
        kB = HallThruster.kB,
        flux = HallThruster.flux,
        HLLE = HallThruster.HLLE,
        rusanov = HallThruster.rusanov,
        global_lax_friedrichs = HallThruster.global_lax_friedrichs,
        Xe_0 = HallThruster.Species(HallThruster.Xenon, 0),

        R = Xenon.R

        let Xe_0 = HallThruster.Species(HallThruster.Xenon, 0)
            @test HallThruster.ContinuityOnly(Xe_0; u = 300, T = 300) |> HallThruster.nvars == 1
            @test HallThruster.IsothermalEuler(Xe_0; T = 300) |> HallThruster.nvars == 2
            @test HallThruster.EulerEquations(Xe_0) |> HallThruster.nvars == 3
        end

        ρ = 1.0
        T = 300
        u = 300
        ϵ = Xenon.cv * T + 0.5 * u^2
        mXe = Xenon.m

        continuity_eq = Fluid(Xe_0; u, T)
        continuity_state = (ρ,)
        continuity = (continuity_state, continuity_eq)

        isothermal_eq = Fluid(Xe_0; T)
        isothermal_state = (ρ, ρ * u,)
        isothermal = (isothermal_state, isothermal_eq)

        euler_eq = Fluid(Xe_0)
        euler_state = (ρ, ρ * u, ρ * ϵ,)
        euler = (euler_state, euler_eq)

        laws = [continuity, isothermal, euler]
        laws_vector = [continuity, isothermal, euler]

        function test_property(property, laws)
            initval = 0.0
            for (i, (U, f)) in enumerate(laws)
                if i == 1
                    initval = property(U, f)
                else
                    if property(U, f) ≉ initval
                        return false
                    end
                end
            end
            return true
        end

        @testset "Thermodynamic property computation" begin
            # Check to make sure our property checking code works
            function fake_property(U, f::Fluid)
                if f.type == :EulerEquations
                    return 2
                else
                    return 1
                end
            end
            #@test !test_property(fake_property, laws)

            # Check that thermodynamic property computations give identical
            # results for the different fluid types
            @test test_property(temperature, laws)
            @test test_property(pressure, laws)
            @test test_property(density, laws)
            @test test_property(number_density, laws)
            @test test_property(velocity, laws)
            @test test_property(sound_speed, laws)

            # check that versions which work for vectors instead of staticarrays work properly
            @test test_property(temperature, laws_vector)
            @test test_property(pressure, laws_vector)
            @test test_property(density, laws_vector)
            @test test_property(number_density, laws_vector)
            @test test_property(velocity, laws_vector)
            @test test_property(sound_speed, laws_vector)

            # Check that properties are being computed correctly
            @test temperature(continuity...) ≈ T
            @test velocity(continuity...) ≈ u
            @test number_density(continuity...) ≈ ρ / m(continuity_eq)
            @test density(continuity...) ≈  ρ
            @test pressure(continuity...) ≈ ρ * HallThruster.Xenon.R * T
            @test sound_speed(continuity...) ≈ √(γ(continuity_eq) * R * T)
        end

        continuity_state_2 = continuity_state .* 2
        isothermal_state_2 = isothermal_state .* 2
        euler_state_2 = euler_state .* 2

        @testset "Flux computation" begin
            p = ρ * R * T
            f_euler = (ρ * u, ρ * u^2 + p, ρ * u * (ϵ + p / ρ),)
            @test flux(continuity...) == (f_euler[1],)
            @test flux(isothermal...) == (f_euler[1], f_euler[2],)
            @test flux(euler...) == f_euler

            isothermal_state_3 = (isothermal_state_2[1], -3 * isothermal_state[2],)
            euler_state_3 = (euler_state_2[1], -2 * euler_state_2[2], euler_state_2[3])

            # rusanov flux
            @test rusanov(continuity_state, continuity...) == flux(continuity...)
            @test rusanov(isothermal_state, isothermal...) == flux(isothermal...)
            @test rusanov(euler_state, euler...) == flux(euler...)

            #global_lax_friedrichs
            @test global_lax_friedrichs(continuity_state, continuity...) == flux(continuity...)
            @test global_lax_friedrichs(isothermal_state, isothermal...) == flux(isothermal...)
            @test global_lax_friedrichs(euler_state, euler...) == flux(euler...)

            # HLLE flux
            @test HLLE(continuity_state, continuity...) == flux(continuity...)
            @test HLLE(isothermal_state, isothermal...) == flux(isothermal...)
            @test HLLE(euler_state, euler...) == flux(euler...)
        end

        @testset "Reconstruction" begin
            # check that if the states actually lie along a line, we correctly reproduce the linear values at the inteface
            euler_state_0 = collect(euler_state)
            euler_state_L = 0.5 .* euler_state_0
            euler_state_R = 2.0 .* euler_state_0

            edge_L = zeros(length(euler_state_0), 2)
            edge_R = zeros(length(euler_state_0), 2)
            U_euler = hcat(euler_state_L, euler_state_0, euler_state_R)

            index = (ρi = [1], ρiui = [2])
            scheme = HallThruster.HyperbolicScheme(identity, HallThruster.no_limiter, false)
            config = (;scheme, ncharge = 1)
            params = (;index, config)

            HallThruster.compute_edge_states!(edge_L, edge_R, U_euler, params)
            @test edge_L[:, 1] == euler_state_L
            @test edge_R[:, end] == euler_state_R
            @test edge_L[:, 2] == euler_state_0
            @test edge_R[:, 1] == euler_state_0

            scheme = HallThruster.HyperbolicScheme(identity, HallThruster.no_limiter, true)
            HallThruster.compute_edge_states!(edge_L, edge_R, U_euler, params)
            @test edge_L[:, 1] == euler_state_L
            @test edge_R[:, end] == euler_state_R

            # check that if slopes have different signs, the avg slope resets to zero with any flux limiter
            euler_state_L2 = 2 .* euler_state_0
            U_euler_2 = hcat(euler_state_L2, euler_state_0, euler_state_R)

            limiters = [
                HallThruster.koren,
                HallThruster.minmod,
                HallThruster.osher(1.5),
                HallThruster.van_albada,
                HallThruster.van_leer
            ]

            for limiter in limiters
                scheme = HallThruster.HyperbolicScheme(identity, limiter, true)
                HallThruster.compute_edge_states!(edge_L, edge_R, U_euler_2, params)
                @test edge_L[:, 1] == euler_state_L2
                @test edge_R[:, end] == euler_state_R
                @test edge_L[:, 2] == euler_state_0
                @test edge_R[:, 1] == euler_state_0
            end

        end
    end
end
