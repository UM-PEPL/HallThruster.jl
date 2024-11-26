using HallThruster: HallThruster as het

include("serialization_test_utils.jl")

function test_wall_loss_serialization()
    @testset "Serialization" begin
        @testset "NoWallLosses" begin
            test_subtype(het.WallLossModel, het.NoWallLosses())
        end

        @testset "ConstantSheathPotential" begin
            test_subtype(het.WallLossModel, het.ConstantSheathPotential(20.0, 0.1, 0.5))
        end

        @testset "WallSheath" begin
            test_subtype(het.WallLossModel, het.WallSheath(het.BNSiO2, 1.0))
        end

        @testset "Wall materials" begin
            test_instances(het.WallMaterial, het.wall_materials)
        end
    end
end

function test_electron_losses()
    @testset "Electron wall losses" begin
        no_losses = het.NoWallLosses()
        h = het.edge_to_center_density_ratio()
        mi = het.Xenon.m
        me = het.me

        αin = 0.3
        αout = 1.0

        landmark_losses = het.ConstantSheathPotential(20.0, αin, αout)

        geom = het.SPT_100.geometry
        L_ch = geom.channel_length
        A_ch = het.channel_area(het.SPT_100)

        Tev = 4.0
        ne = 1e18
        nϵ = 3 / 2 * ne * Tev
        cache = (;
            ne = [ne, ne, ne, ne],
            ni = [ne ne ne ne],
            Tev = [Tev, Tev, Tev, Tev],
            ϵ = 1.5 .* [Tev, Tev, Tev, Tev],
            nϵ = [nϵ, nϵ, nϵ, nϵ],
            Z_eff = [1.0, 1.0, 1.0, 1.0],
            γ_SEE = [0.0, 0.0, 0.0, 0.0],
            radial_loss_frequency = [0.0, 0.0, 0.0, 0.0],
            νew_momentum = [0.0, 0.0, 0.0, 0.0],
        )
        transition_length = 0.2 * L_ch
        config = het.Config(;
            thruster = het.SPT_100,
            anode_mass_flow_rate = 5e-6,
            domain = (0.0, 1.0),
            discharge_voltage = 300.0,
            transition_length,
            propellant = het.Xenon, ncharge = 1,
            electron_plume_loss_scale = 0.0,
        )
        index = (; ρi = [1], nϵ = 2)

        ncells = 2
        grid = het.generate_grid(geom, (0, 2 * L_ch), het.EvenGrid(2))

        mi_kr = het.Krypton.m
        γmax = 1 - 8.3 * sqrt(me / mi)
        γmax_kr = 1 - 8.3 * sqrt(me / mi_kr)

        params = (;
            het.params_from_config(config)...,
            cache, index, grid,
            γ_SEE_max = γmax,
        )

        arr = zeros(ncells + 2)
        het.wall_power_loss!(arr, no_losses, params)
        @test arr[2] == 0.0

        het.wall_power_loss!(arr, landmark_losses, params)
        @test arr[2] ≈ αin * 6.0 * 1e7 * exp(-20 / 6.0)
        @test arr[3] ≈ αout * 6.0 * 1e7 * exp(-20 / 6.0)

        γ1 = 0.5
        γ2 = 1.0
        @test het.sheath_potential(Tev, γ1, mi) ==
              Tev * log(0.5 * (1 - γ1) * sqrt(2 * mi / π / me))
        @test het.sheath_potential(Tev, γ2, mi) ==
              Tev * log(0.5 * (1 - γ2) * sqrt(2 * mi / π / me))

        BN = het.BoronNitride
        @test het.SEE_yield(BN, 10.0, γmax) ≈
              BN.σ₀ + 1.5 * 10.0 / BN.ϵ_star * (1 - BN.σ₀)
        # Check that SEE properly obeys space charge limit at high electron temps
        @test het.SEE_yield(BN, 9000.0, γmax) ≈ γmax

        @test het.SEE_yield(BN, 9000.0, γmax_kr) ≈ γmax_kr

        γ = het.SEE_yield(BN, Tev, γmax)

        Vs = het.sheath_potential(Tev, γ, mi)

        α = 1 / 4
        sheath_model = het.WallSheath(BN, α)
        Δr = geom.outer_radius - geom.inner_radius
        Δz = grid.edges[2] - grid.edges[1]
        V_cell = A_ch * Δz

        @test het.freq_electron_wall(no_losses, params, 2) == 0.0
        @test het.freq_electron_wall(no_losses, params, 3) == 0.0

        @test het.freq_electron_wall(landmark_losses, params, 2) == 1e7
        @test het.freq_electron_wall(landmark_losses, params, 3) *
              het.linear_transition(
            grid.cell_centers[3], L_ch, config.transition_length, 1.0, 0.0,) == 0.0e7

        γ = het.SEE_yield(BN, Tev, γmax)
        νiw = α * sqrt(het.e * Tev / mi) / Δr * h
        νew = νiw / (1 - γ)

        params.cache.γ_SEE .= γ
        params.cache.νew_momentum[1] = νew
        params.cache.νew_momentum[2] = νew
        params.cache.νew_momentum[3] = 0.0
        params.cache.νew_momentum[4] = 0.0
        params.cache.radial_loss_frequency[1:4] .= νew

        Iiw = het.wall_ion_current(sheath_model, params, 2, 1)
        Iew = het.wall_electron_current(sheath_model, params, 2)
        @test Iew ≈ Iiw / (1 - γ)
        @test Iiw ≈ νiw * het.e * V_cell * ne
        @test Iew ≈ νew * het.e * V_cell * ne

        @test het.freq_electron_wall(sheath_model, params, 2) *
              het.linear_transition(
            grid.cell_centers[2], L_ch, config.transition_length, 1.0, 0.0,) ≈ νew
        @test het.freq_electron_wall(sheath_model, params, 3) *
              het.linear_transition(
            grid.cell_centers[3], L_ch, config.transition_length, 1.0, 0.0,) ≈ 0.0

        het.wall_power_loss!(arr, sheath_model, params)
        @test arr[2] ≈ νew * (2 * Tev + (1 - γ) * Vs)
        @test arr[4] ≈ 0.0
    end
end

function test_ion_losses()
    @testset "Ion wall losses" begin
        config = (; thruster = het.SPT_100, propellant = het.Krypton,
            ncharge = 2, transition_length = 0.0, anode_mass_flow_rate = 5e-6,
            discharge_voltage = 300.0, domain = (0.0, 1.0),)
        geom = het.SPT_100.geometry
        Δr = geom.outer_radius - geom.inner_radius
        L_ch = geom.channel_length
        h = het.edge_to_center_density_ratio()

        Tn = 300.0
        Tev = 3.0

        α = 0.8
        mi = config.propellant.m

        γ_SEE_max = 1 - 8.3 * sqrt(het.me / mi)

        γ = het.SEE_yield(het.BoronNitride, Tev, γ_SEE_max)
        νiw = α * sqrt(het.e * Tev / mi) / Δr * h
        νew = νiw * γ / (1 - γ)

        grid = het.generate_grid(geom, (0, 2 * L_ch), het.EvenGrid(2))

        fluids = [
            het.Fluid(Xenon(0); u = 300.0, T = Tn),
            het.Fluid(Xenon(1); T = Tn),
            het.Fluid(Xenon(2); T = Tn),
        ]

        Tev = 4.0
        ne = 1.3e18
        ni_1 = 7e17
        ni_2 = 3e17

        cache = (;
            ne = [ne, ne, ne, ne], Tev = [Tev, Tev, Tev, Tev],
            Z_eff = [1.0, 1.0, 1.0, 1.0], ni = [ni_1 ni_1 ni_1 ni_1; ni_2 ni_2 ni_2 ni_2],
            γ_SEE = [0.0, 0.0, 0.0, 0.0],
            νew_momentum = [νew, νew, 0.0, 0.0],
        )

        index = (ρn = 1, ρi = [2, 4], ρiui = [3, 5], nϵ = 6)

        u_bohm = sqrt(het.e * Tev / mi)

        ρn = 10 * ne * mi

        ρi_1 = ni_1 * mi
        ρi_2 = ni_2 * mi
        ui_1 = 10.0
        ui_2 = sqrt(2) * 10.0

        nϵ = ne * Tev * 3 / 2

        U = [ρn ρn ρn ρn
             ρi_1 ρi_1 ρi_1 ρi_1
             ρi_1*ui_1 ρi_1*ui_1 ρi_1*ui_1 ρi_1*ui_1
             ρi_2 ρi_2 ρi_2 ρi_2
             ρi_2*ui_2 ρi_2*ui_2 ρi_2*ui_2 ρi_2*ui_2
             nϵ nϵ nϵ nϵ]

        dU = zeros(size(U))

        γ_SEE_max = 1 - 8.3 * sqrt(het.me / mi)
        base_params = (; cache, fluids, grid, index, γ_SEE_max)

        config_no_losses = het.Config(; config..., wall_loss_model = het.NoWallLosses())
        params_no_losses = (; base_params..., het.params_from_config(config_no_losses)...)

        for i in 1:4
            cache.Z_eff[i] = (ni_1 + 2 * ni_2) / (ni_1 + ni_2)
        end

        # Test 1: no wall losses
        het.apply_ion_wall_losses!(
            dU, U, params_no_losses, config_no_losses.wall_loss_model,)

        @test all(dU .≈ 0.0)

        # Test 2: LANDMARK wall losses
        αin, αout = 1.0, 1.0
        constant_sheath = het.ConstantSheathPotential(-20, αin, αout)
        config_constant_sheath = het.Config(; config..., wall_loss_model = constant_sheath)
        params_constant_sheath = (;
            base_params..., het.params_from_config(config_constant_sheath)...,)

        i = 2
        Iiw_1 = het.wall_ion_current(constant_sheath, params_constant_sheath, i, 1)
        Iiw_2 = het.wall_ion_current(constant_sheath, params_constant_sheath, i, 2)
        Iew = het.wall_electron_current(constant_sheath, params_constant_sheath, i)

        # Ion and electron wall currents are equivalent inside of channel
        @test Iiw_1 + Iiw_2 == Iew

        dU .= 0.0
        # Check that wall losses work correctly
        het.apply_ion_wall_losses!(
            dU, U, params_constant_sheath, config_constant_sheath.wall_loss_model,)

        # Neutrals should recombine at walls
        @test dU[index.ρn, i] ≈ -(dU[index.ρi[1], i] + dU[index.ρi[2], i])

        # Rate of ion loss is equal to ni νiw
        u_bohm = sqrt(het.e * Tev / mi)
        νiw = u_bohm / Δr * h

        @test dU[index.ρi[1], i] ≈ -U[index.ρi[1], i] * νiw
        @test dU[index.ρi[2], i] ≈ -U[index.ρi[2], i] * νiw * sqrt(2)

        # ion momentum loss is equal to Iiw / e / V_cell * ui
        @test dU[index.ρiui[1], i] ≈ -U[index.ρiui[1], i] * νiw
        @test dU[index.ρiui[2], i] ≈ -U[index.ρiui[2], i] * νiw * sqrt(2)

        # No wall current outside of channel
        dU .= 0.0
        i = 3
        @test het.wall_ion_current(constant_sheath, params_constant_sheath, i, 1) ==
              0.0
        @test het.wall_ion_current(constant_sheath, params_constant_sheath, i, 2) ==
              0.0
        @test het.wall_electron_current(constant_sheath, params_constant_sheath, i) ==
              0.0

        het.apply_ion_wall_losses!(
            dU, U, params_constant_sheath, config_constant_sheath.wall_loss_model,)
        @test all(dU[:, 3] .≈ 0.0)

        # Test 3: Self-consistent wall sheath
        dU .= 0.0
        material = het.BNSiO2
        wall_sheath = het.WallSheath(material, α)
        config_wall_sheath = het.Config(; config..., wall_loss_model = wall_sheath)
        params_wall_sheath = (;
            base_params..., het.params_from_config(config_wall_sheath)...,)

        γ = het.SEE_yield(material, Tev, γ_SEE_max)
        νiw = α * sqrt(het.e * Tev / mi) / Δr * h
        νew = νiw * γ / (1 - γ)

        params_wall_sheath.cache.νew_momentum[1:2] .= νew
        i = 2
        params_wall_sheath.cache.γ_SEE[i] = γ
        Iiw_1 = het.wall_ion_current(wall_sheath, params_wall_sheath, i, 1)
        Iiw_2 = het.wall_ion_current(wall_sheath, params_wall_sheath, i, 2)
        Iew = het.wall_electron_current(wall_sheath, params_wall_sheath, i)

        @test Iiw_1 ≈ Iew * ni_1 / ne * (1 - γ)
        @test Iiw_2 ≈ 2 * Iew * ni_2 / ne * (1 - γ)
        @test Iew ≈ inv(1 - γ) * (Iiw_1 + Iiw_2)

        het.apply_ion_wall_losses!(
            dU, U, params_wall_sheath, config_wall_sheath.wall_loss_model,)

        # Neutrals should recombine at walls
        @test dU[index.ρn, i] ≈ -(dU[index.ρi[1], i] + dU[index.ρi[2], i])

        # Rate of ion loss is equal to Iiw / e / V_cell
        @test dU[index.ρi[1], i] ≈ -νiw * U[index.ρi[1], i]
        @test dU[index.ρi[2], i] ≈ -νiw * U[index.ρi[2], i] * sqrt(2)

        # ion momentum loss is equal to Iiw / e / V_cell * ui
        @test dU[index.ρiui[1], i] ≈ -νiw * U[index.ρiui[1], i]
        @test dU[index.ρiui[2], i] ≈ -νiw * U[index.ρiui[2], i] * sqrt(2)

        # No wall losses in plume
        i = 3
        @test het.wall_ion_current(wall_sheath, params_wall_sheath, i, 1) == 0.0
        @test het.wall_ion_current(wall_sheath, params_wall_sheath, i, 2) == 0.0
        @test het.wall_electron_current(wall_sheath, params_wall_sheath, i) == 0.0

        het.apply_ion_wall_losses!(
            dU, U, params_wall_sheath, config_wall_sheath.wall_loss_model,)
        @test all(dU[:, 3] .≈ 0.0)
    end
end

function test_walls()
    test_wall_loss_serialization()
    test_electron_losses()
    test_ion_losses()
end
