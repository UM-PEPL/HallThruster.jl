using HallThruster: HallThruster as het

function test_limiter()
    return @testset "Slope limiter" begin
        # List of slope limiters
        slope_limiters = [
            het.van_leer_limiter,
        ]

        # convert to flux limiters to test that they lie in 2nd order TVD region of sweby diagram
        flux_limiters = map(
            limiter -> (r -> limiter(r) * (r + 1) / 2),
            slope_limiters,
        )

        r1 = 0.0:0.1:1.0
        r2 = 1.0:0.1:5.0
        r3 = -2:0.1:0.0

        for limiter in flux_limiters
            # check 2nd order TVD properties
            ϵ = sqrt(eps(Float64))
            # ψ(0) = 0
            @test limiter(0) ≈ 0
            # ψ(r) = 0 ∀ r < 0
            @test all(@. limiter(r3) ≈ 0)

            # ψ(r) ≥ 0 ∀ r ∈ ℝ
            @test all(@. limiter(r1) ≥ 0 - ϵ)
            @test all(@. limiter(r2) ≥ 0 - ϵ)
            @test all(@. limiter(r3) ≥ 0 - ϵ)
            # ψ(1) = 1
            @test limiter(1) == 1
            # 1.0 ≤ ψ(r) ≤ 2.0 ∀ r ∈ [1, 2]
            @test limiter(100) ≤ 2.0 + ϵ
            @test limiter(100) ≥ 1.0 - ϵ
            @test all(@. limiter(r2) ≤ 2.0 + ϵ)
            @test all(@. limiter(r2) ≥ 1.0 - ϵ)
            # ψ(r) ≤ 2r ∀ r ∈ [0, 1]
            @test all(@. limiter(r1) ≤ 2r1 + ϵ)
            # ψ(r) ≤ r ∀ r ∈ [1, 2]
            @test all(@. limiter(r1) ≥ r1 - ϵ)
            # ψ(r) ≥ r ∀ r ∈ [1, 2]
            @test all(@. limiter(r2) ≤ r2 + ϵ)
        end
    end
end

function test_fluxes()
    # ρ = 1.0
    # T = 300
    # u = 300
    # ϵ = het.Xenon.cv * T + 0.5 * u^2
    # Xe_0 = het.Xenon(0)

    # continuity_eq = het.Fluid(Xe_0; u, T)
    # continuity_state = (ρ,)
    # continuity = (continuity_state, continuity_eq)

    # isothermal_eq = het.Fluid(Xe_0; T)
    # isothermal_state = (ρ, ρ * u)
    # isothermal = (isothermal_state, isothermal_eq)

    # @testset "Flux computation" begin
    #     p = ρ * het.Xenon.R * T
    #     flux_expected = (ρ * u, ρ * u^2 + p)
    #     @test het.flux(continuity...) == (flux_expected[1],)
    #     @test het.flux(isothermal...) == (flux_expected[1], flux_expected[2])
    #     @test het.flux(euler...) == flux_expected

    #     # rusanov flux
    #     @test het.rusanov(continuity_state, continuity...) == het.flux(continuity...)
    #     @test het.rusanov(isothermal_state, isothermal...) == het.flux(isothermal...)
    # end

    # @testset "Reconstruction" begin
    #     # check that if the states actually lie along a line, we correctly
    #     # reproduce the linear values at the inteface
    #     euler_state_0 = collect(euler_state)
    #     euler_state_L = 0.5 .* euler_state_0
    #     euler_state_R = 2.0 .* euler_state_0

    #     edge_L = zeros(length(euler_state_0), 2)
    #     edge_R = zeros(length(euler_state_0), 2)
    #     U_euler = hcat(euler_state_L, euler_state_0, euler_state_R)

    #     index = (ρi = [1], ρiui = [2])
    #     scheme = het.HyperbolicScheme(het.FluxFunction(identity), het.no_limiter, false)
    #     config = het.Config(; discharge_voltage = 300.0, thruster = het.SPT_100,
    #         anode_mass_flow_rate = 1.0e-6,
    #         domain = (0, 1), scheme, ncharge = 1,)
    #     params = (; index, het.params_from_config(config)...,
    #         is_velocity_index = [false, true, false],)

    #     het.compute_edge_states!(
    #         edge_L, edge_R, U_euler, params, scheme.limiter, scheme.reconstruct,)
    #     @test edge_L[:, 1] == euler_state_L
    #     @test edge_R[:, end] == euler_state_R
    #     @test edge_L[:, 2] == euler_state_0
    #     @test edge_R[:, 1] == euler_state_0

    #     scheme = het.HyperbolicScheme(het.FluxFunction(identity), het.no_limiter, true)
    #     het.compute_edge_states!(
    #         edge_L, edge_R, U_euler, params, scheme.limiter, scheme.reconstruct,)
    #     @test edge_L[:, 1] == euler_state_L
    #     @test edge_R[:, end] == euler_state_R

    #     # check that if slopes have different signs, the avg slope resets to zero with any flux limiter
    #     euler_state_L2 = 2 .* euler_state_0
    #     U_euler_2 = hcat(euler_state_L2, euler_state_0, euler_state_R)

    #     for limiter in het.slope_limiters
    #         scheme = het.HyperbolicScheme(het.FluxFunction(identity), limiter, true)
    #         het.compute_edge_states!(
    #             edge_L, edge_R, U_euler_2, params, scheme.limiter, scheme.reconstruct,)
    #         @test edge_L[:, 1] == euler_state_L2
    #         @test edge_R[:, end] == euler_state_R
    #         @test edge_L[:, 2] == euler_state_0
    #         @test edge_R[:, 1] == euler_state_0
    #     end
    # end
end

test_limiter()
test_fluxes()
