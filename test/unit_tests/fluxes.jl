using HallThruster: HallThruster as het

function test_flux_serialization()
    @testset "Serialization" begin
        test_instances(het.FluxFunction, het.flux_functions)
    end
end

function test_flux_computation()
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

    @testset "Flux computation" begin
        p = ρ * het.Xenon.R * T
        f_euler = (ρ * u, ρ * u^2 + p, ρ * u * (ϵ + p / ρ))
        @test het.flux(continuity...) == (f_euler[1],)
        @test het.flux(isothermal...) == (f_euler[1], f_euler[2])
        @test het.flux(euler...) == f_euler

        # rusanov flux
        @test het.rusanov(continuity_state, continuity...) == het.flux(continuity...)
        @test het.rusanov(isothermal_state, isothermal...) == het.flux(isothermal...)
        @test het.rusanov(euler_state, euler...) == het.flux(euler...)

        #global_lax_friedrichs
        @test het.global_lax_friedrichs(continuity_state, continuity...) ==
              het.flux(continuity...)
        @test het.global_lax_friedrichs(isothermal_state, isothermal...) ==
              het.flux(isothermal...)
        @test het.global_lax_friedrichs(euler_state, euler...) == het.flux(euler...)

        # HLLE flux
        @test het.HLLE(continuity_state, continuity...) == het.flux(continuity...)
        @test het.HLLE(isothermal_state, isothermal...) == het.flux(isothermal...)
        @test het.HLLE(euler_state, euler...) == het.flux(euler...)
    end

    @testset "Reconstruction" begin
        # check that if the states actually lie along a line, we correctly
        # reproduce the linear values at the inteface
        euler_state_0 = collect(euler_state)
        euler_state_L = 0.5 .* euler_state_0
        euler_state_R = 2.0 .* euler_state_0

        edge_L = zeros(length(euler_state_0), 2)
        edge_R = zeros(length(euler_state_0), 2)
        U_euler = hcat(euler_state_L, euler_state_0, euler_state_R)

        index = (ρi = [1], ρiui = [2])
        scheme = het.HyperbolicScheme(het.FluxFunction(identity), het.no_limiter, false)
        config = het.Config(; discharge_voltage = 300.0, thruster = het.SPT_100,
            anode_mass_flow_rate = 1.0e-6,
            domain = (0, 1), scheme, ncharge = 1,)
        params = (; index, het.params_from_config(config)...,
            is_velocity_index = [false, true, false],)

        het.compute_edge_states!(
            edge_L, edge_R, U_euler, params, scheme.limiter, scheme.reconstruct,)
        @test edge_L[:, 1] == euler_state_L
        @test edge_R[:, end] == euler_state_R
        @test edge_L[:, 2] == euler_state_0
        @test edge_R[:, 1] == euler_state_0

        scheme = het.HyperbolicScheme(het.FluxFunction(identity), het.no_limiter, true)
        het.compute_edge_states!(
            edge_L, edge_R, U_euler, params, scheme.limiter, scheme.reconstruct,)
        @test edge_L[:, 1] == euler_state_L
        @test edge_R[:, end] == euler_state_R

        # check that if slopes have different signs, the avg slope resets to zero with any flux limiter
        euler_state_L2 = 2 .* euler_state_0
        U_euler_2 = hcat(euler_state_L2, euler_state_0, euler_state_R)

        for limiter in het.slope_limiters
            scheme = het.HyperbolicScheme(het.FluxFunction(identity), limiter, true)
            het.compute_edge_states!(
                edge_L, edge_R, U_euler_2, params, scheme.limiter, scheme.reconstruct,)
            @test edge_L[:, 1] == euler_state_L2
            @test edge_R[:, end] == euler_state_R
            @test edge_L[:, 2] == euler_state_0
            @test edge_R[:, 1] == euler_state_0
        end
    end
end

function test_fluxes()
    @testset "Fluxes" begin
        test_flux_serialization()
        test_flux_computation()
    end
end
