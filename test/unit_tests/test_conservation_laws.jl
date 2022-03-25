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
        stagnation_energy = HallThruster.stagnation_energy,
        static_energy = HallThruster.static_energy,
        sound_speed = HallThruster.sound_speed,
        mach_number = HallThruster.mach_number,
        stagnation_enthalpy = HallThruster.stagnation_enthalpy,
        static_enthalpy = HallThruster.static_enthalpy,
        critical_sound_speed = HallThruster.critical_sound_speed,
        m = HallThruster.m,
        γ = HallThruster.γ,
        R = HallThruster.R,
        cp = HallThruster.cp,
        cv = HallThruster.cv,
        kB = HallThruster.kB,
        flux = HallThruster.flux,
        HLLE = HallThruster.HLLE,
        upwind = HallThruster.upwind,
        Xe_0 = HallThruster.Species(HallThruster.Xenon, 0),

        R = Xenon.R

        let Xe_0 = HallThruster.Species(HallThruster.Xenon, 0)
            @test HallThruster.Fluid(Xe_0, HallThruster.ContinuityOnly(u = 300, T = 300)) |> HallThruster.nvars == 1
            @test HallThruster.Fluid(Xe_0, HallThruster.IsothermalEuler(T = 300)) |> HallThruster.nvars == 2
            @test HallThruster.Fluid(Xe_0, HallThruster.EulerEquations()) |> HallThruster.nvars == 3
        end

        ρ = 1.0
        T = 300
        u = 300
        ϵ = Xenon.cv * T + 0.5 * u^2
        mXe = Xenon.m

        continuity_eq = Fluid(Xe_0, ContinuityOnly(; u, T))
        continuity_state = SA[ρ]
        continuity = (continuity_state, continuity_eq)

        isothermal_eq = Fluid(Xe_0, IsothermalEuler(T))
        isothermal_state = SA[ρ, ρ * u]
        isothermal = (isothermal_state, isothermal_eq)

        euler_eq = Fluid(Xe_0, EulerEquations())
        euler_state = SA[ρ, ρ * u, ρ * ϵ]
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
                if f.conservation_laws.type == :EulerEquations
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
            @test test_property(stagnation_energy, laws)
            @test test_property(static_energy, laws)
            @test test_property(sound_speed, laws)
            @test test_property(mach_number, laws)
            @test test_property(stagnation_enthalpy, laws)
            @test test_property(static_enthalpy, laws)
            @test test_property(critical_sound_speed, laws)

            # check that versions which work for vectors instead of staticarrays work properly
            @test test_property(temperature, laws_vector)
            @test test_property(pressure, laws_vector)
            @test test_property(density, laws_vector)
            @test test_property(number_density, laws_vector)
            @test test_property(velocity, laws_vector)
            @test test_property(stagnation_energy, laws_vector)
            @test test_property(static_energy, laws_vector)
            @test test_property(sound_speed, laws_vector)
            @test test_property(mach_number, laws_vector)
            @test test_property(stagnation_enthalpy, laws_vector)
            @test test_property(static_enthalpy, laws_vector)
            @test test_property(critical_sound_speed, laws_vector)

            # Check that properties are being computed correctly
            @test temperature(continuity...) ≈ T
            @test velocity(continuity...) ≈ u
            @test number_density(continuity...) ≈ ρ / m(continuity_eq)
            @test density(continuity...) ≈  ρ
            @test pressure(continuity...) ≈ ρ * HallThruster.Xenon.R * T
            @test static_energy(continuity...) ≈ cv(continuity_eq) * T
            @test stagnation_energy(continuity...) ≈ cv(continuity_eq) * T + 0.5 * u^2
            @test static_enthalpy(continuity...) ≈ cp(continuity_eq) * T
            @test stagnation_enthalpy(continuity...) ≈ cp(continuity_eq) * T + 0.5 * u^2
            @test sound_speed(continuity...) ≈ √(γ(continuity_eq) * R * T)
        end

        continuity_state_2 = continuity_state * 2
        isothermal_state_2 = isothermal_state * 2
        euler_state_2 = euler_state * 2

        @testset "Flux computation" begin
            p = ρ * R * T
            f_euler = SA[ρ * u, ρ * u^2 + p, ρ * u * (ϵ + p / ρ)]
            @test flux(continuity...) == SA[f_euler[1]]
            @test flux(isothermal...) == SA[f_euler[1], f_euler[2]]
            @test flux(euler...) == f_euler

            isothermal_state_3 = SA[isothermal_state_2[1], -3 * isothermal_state[2]]
            euler_state_3 = SA[euler_state_2[1], -2 * euler_state_2[2], euler_state_2[3]]

            # upwind flux
            @test upwind(continuity_state, continuity_state_2, continuity_eq) == flux(continuity...)
            @test upwind(isothermal_state, isothermal_state_2, isothermal_eq) == flux(isothermal...)
            @test upwind(isothermal_state, isothermal_state_3, isothermal_eq) == flux(isothermal_state_3, isothermal_eq)
            @test upwind(euler_state, euler_state_2, euler_eq) == flux(euler...)
            @test upwind(euler_state, euler_state_3, euler_eq) == flux(euler_state_3, euler_eq)

            # rusanov flux

            # HLLE flux
            @test HLLE(continuity_state, continuity...) == flux(continuity...)
            @test HLLE(isothermal_state, isothermal...) == flux(isothermal...)
            @test HLLE(euler_state, euler...) == flux(euler...)

            # electron pressure-coupled
            ne = 1e17
            Te = 6.0
            pe = HallThruster.e * ne * Te
            @test flux(continuity..., pe) ==  flux(continuity...)
            @test flux(isothermal..., pe) == flux(isothermal...) + SA[0.0, pe]
            @test flux(euler..., pe) == flux(euler...) + SA[0.0, pe, 0.0]
        end

        @testset "Reconstruction" begin
            # check that if the states actually lie along a line, we correctly reproduce the linear values at the inteface
            euler_state_L = 0.5 * euler_state
            euler_state_R = 2.0 * euler_state

            edge_L = zeros(length(euler_state), 2)
            edge_R = zeros(length(euler_state), 2)
            U_euler = hcat(euler_state_L, euler_state, euler_state_R)

            scheme = HallThruster.HyperbolicScheme(identity, HallThruster.no_limiter, false)

            HallThruster.compute_edge_states!(edge_L, edge_R, U_euler, scheme)
            @test edge_L[:, 1] == euler_state_L
            @test edge_R[:, end] == euler_state_R
            @test edge_L[:, 2] == euler_state
            @test edge_R[:, 1] == euler_state

            scheme = HallThruster.HyperbolicScheme(identity, HallThruster.no_limiter, true)
            HallThruster.compute_edge_states!(edge_L, edge_R, U_euler, scheme)
            @test edge_L[:, 1] == euler_state_L
            @test edge_R[:, end] == euler_state_R
            @test edge_L[:, 2] == 1.5 * euler_state
            @test edge_R[:, 1] == 0.75 * euler_state

            # check that if slopes have different signs, the avg slope resets to zero with any flux limiter
            euler_state_L2 = 2 * euler_state
            U_euler_2 = hcat(euler_state_L2, euler_state, euler_state_R)

            limiters = [
                HallThruster.koren,
                HallThruster.minmod,
                HallThruster.osher,
                HallThruster.superbee,
                HallThruster.van_albada,
                HallThruster.van_leer
            ]

            for limiter in limiters
                scheme = HallThruster.HyperbolicScheme(identity, limiter, true)
                HallThruster.compute_edge_states!(edge_L, edge_R, U_euler_2, scheme)
                @test edge_L[:, 1] == euler_state_L2
                @test edge_R[:, end] == euler_state_R
                @test edge_L[:, 2] == euler_state
                @test edge_R[:, 1] == euler_state
            end

        end
    end
end