using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

function test_wall_loss_serialization()
    return @testset "Serialization" begin
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
    return @testset "Electron wall losses" begin
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
        ne = 1.0e18
        nϵ = 3 / 2 * ne * Tev
        ncells = 2

        transition_length = 0.2 * L_ch
        config = het.Config(;
            thruster = het.SPT_100,
            anode_mass_flow_rate = 5.0e-6,
            domain = (0.0, 2 * L_ch),
            discharge_voltage = 300.0,
            transition_length,
            propellant = het.Xenon,
            ncharge = 1,
            electron_plume_loss_scale = 0.0,
        )
        simparams = het.SimParams(grid = het.EvenGrid(ncells))
        params = het.setup_simulation(config, simparams)
        (; cache, grid) = params

        # Initialize plasma state
        @. cache.ne = [ne, ne, ne, ne]
        # @. cache.ni = [ne ne ne ne]
        @. cache.Tev = [Tev, Tev, Tev, Tev]
        @. cache.ϵ = 1.5 .* [Tev, Tev, Tev, Tev]
        @. cache.nϵ = [nϵ, nϵ, nϵ, nϵ]

        mi_kr = het.Krypton.m
        γmax = 1 - 8.3 * sqrt(me / mi)
        γmax_kr = 1 - 8.3 * sqrt(me / mi_kr)

        arr = zeros(ncells + 2)
        het.wall_power_loss!(arr, no_losses, params)
        @test arr[2] == 0.0

        het.wall_power_loss!(arr, landmark_losses, params)
        @test arr[2] ≈ αin * 6.0 * 1.0e7 * exp(-20 / 6.0)
        @test arr[3] ≈ αout * 6.0 * 1.0e7 * exp(-20 / 6.0)

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

        @test het.freq_electron_wall(landmark_losses, params, 2) == 1.0e7
        @test het.freq_electron_wall(landmark_losses, params, 3) *
            het.linear_transition(
            grid.cell_centers[3], L_ch, config.transition_length, 1.0, 0.0,
        ) == 0.0e7

        γ = het.SEE_yield(BN, Tev, γmax)
        νiw = α * sqrt(het.e * Tev / mi) / Δr * h
        νew = νiw / (1 - γ)

        params.cache.γ_SEE .= γ
        params.cache.νew_momentum[1] = νew
        params.cache.νew_momentum[2] = νew
        params.cache.νew_momentum[3] = 0.0
        params.cache.νew_momentum[4] = 0.0
        params.cache.radial_loss_frequency[1:4] .= νew

        Iew = het.wall_electron_current(sheath_model, params, 2)
        @test Iew ≈ νew * het.e * V_cell * ne

        @test het.freq_electron_wall(sheath_model, params, 2) *
            het.linear_transition(
            grid.cell_centers[2], L_ch, config.transition_length, 1.0, 0.0,
        ) ≈ νew
        @test het.freq_electron_wall(sheath_model, params, 3) *
            het.linear_transition(
            grid.cell_centers[3], L_ch, config.transition_length, 1.0, 0.0,
        ) ≈ 0.0

        het.wall_power_loss!(arr, sheath_model, params)
        @test arr[2] ≈ νew * (2 * Tev + (1 - γ) * Vs)
        @test arr[4] ≈ 0.0
    end
end

function test_ion_losses()
    return @testset "Ion wall losses" begin
        propellant = het.Krypton
        Tn = 300.0
        Ti = 1000.0
        un = 150.0
        config = (;
            thruster = het.SPT_100, propellant = propellant,
            ncharge = 2, transition_length = 0.0, anode_mass_flow_rate = 5.0e-6,
            discharge_voltage = 300.0, domain = (0.0, 1.0),
        )
        geom = het.SPT_100.geometry
        Δr = geom.outer_radius - geom.inner_radius
        L_ch = geom.channel_length
        h = het.edge_to_center_density_ratio()

        Tev = 3.0

        α = 0.8
        mi = config.propellant.m

        γ_SEE_max = 1 - 8.3 * sqrt(het.me / mi)

        γ = het.SEE_yield(het.BoronNitride, Tev, γ_SEE_max)
        νiw = α * sqrt(het.e * Tev / mi) / Δr * h
        νew = νiw * γ / (1 - γ)

        grid = het.generate_grid(het.EvenGrid(2), geom, (0, 2 * L_ch))

        Tev = 4.0
        ne = 1.3e18
        ni_1 = 7.0e17
        ni_2 = 3.0e17

        cache = (;
            ne = [ne, ne, ne, ne], Tev = [Tev, Tev, Tev, Tev],
            Z_eff = [1.0, 1.0, 1.0, 1.0], ni = [ni_1 ni_1 ni_1 ni_1; ni_2 ni_2 ni_2 ni_2],
            γ_SEE = [0.0, 0.0, 0.0, 0.0],
            νew_momentum = [νew, νew, 0.0, 0.0],
        )

        ρn = 10 * ne * mi
        ρi_1 = ni_1 * mi
        ρi_2 = ni_2 * mi
        ui_1 = 10.0
        ui_2 = sqrt(2) * 10.0

        config_no_losses = het.Config(; config..., wall_loss_model = het.NoWallLosses())

        fluid_containers = het.allocate_fluids(config_no_losses.propellants[1], 2)
        (; continuity, isothermal) = fluid_containers
        @test length(continuity[1].density) == 4
        @. continuity[1].density = ρn
        @. isothermal[1].density = ρi_1
        @. isothermal[1].momentum = ρi_1 * ui_1
        @. isothermal[2].density = ρi_2
        @. isothermal[2].momentum = ρi_2 * ui_2
        fluids = [continuity..., isothermal...]

        γ_SEE_max = 1 - 8.3 * sqrt(het.me / mi)
        base_params = (; cache, grid, γ_SEE_max, fluid_containers = (; continuity, isothermal))

        params_no_losses = (; base_params..., het.params_from_config(config_no_losses)...)

        for i in 1:4
            cache.Z_eff[i] = (ni_1 + 2 * ni_2) / (ni_1 + ni_2)
        end

        # Test 1: no wall losses
        het.apply_ion_wall_losses!(fluid_containers, params_no_losses)
        for fluid in fluids
            @test all(fluid.dens_ddt .≈ 0.0)
            @test all(fluid.mom_ddt .≈ 0.0)
        end

        # Test 2: LANDMARK wall losses
        αin, αout = 1.0, 1.0
        constant_sheath = het.ConstantSheathPotential(-20, αin, αout)
        config_constant_sheath = het.Config(; config..., wall_loss_model = constant_sheath)
        params_constant_sheath = (;
            base_params..., het.params_from_config(config_constant_sheath)...,
        )

        in_index = 2
        out_index = 3

        for fluid in fluids
            fluid.mom_ddt .= 0.0
            fluid.dens_ddt .= 0.0
        end

        # Check that wall losses work correctly
        het.apply_ion_wall_losses!(fluid_containers, params_constant_sheath)

        u_bohm = sqrt(het.e * Tev / mi)
        νiw = u_bohm / Δr * h

        # Neutrals should recombine at walls
        @test continuity[1].dens_ddt[in_index] ≈ -(isothermal[1].dens_ddt[in_index] + isothermal[2].dens_ddt[in_index])
        @test continuity[1].dens_ddt[out_index] ≈ 0.0

        for ion in isothermal
            # Rate of ion loss is equal to ni νiw
            @test ion.dens_ddt[in_index] ≈ -ion.density[in_index] * νiw * sqrt(ion.species.Z)
            # ion momentum loss is equal to Iiw / e / V_cell * ui
            @test ion.mom_ddt[in_index] ≈ -ion.momentum[in_index] * νiw * sqrt(ion.species.Z)

            # No losses outside of channel
            @test ion.dens_ddt[out_index] ≈ 0.0
            @test ion.mom_ddt[out_index] ≈ 0.0
        end

        # No wall current outside of channel
        @test het.wall_electron_current(constant_sheath, params_constant_sheath, out_index) == 0.0

        # Test 3: Self-consistent wall sheath
        for fluid in fluids
            fluid.mom_ddt .= 0.0
            fluid.dens_ddt .= 0.0
        end

        material = het.BNSiO2
        wall_sheath = het.WallSheath(material, α)
        config_wall_sheath = het.Config(; config..., wall_loss_model = wall_sheath)
        params_wall_sheath = (;
            base_params..., het.params_from_config(config_wall_sheath)...,
        )

        γ = het.SEE_yield(material, Tev, γ_SEE_max)
        νiw = α * sqrt(het.e * Tev / mi) / Δr * h
        νew = νiw * γ / (1 - γ)

        params_wall_sheath.cache.νew_momentum[1:2] .= νew
        params_wall_sheath.cache.γ_SEE[in_index] = γ
        het.apply_ion_wall_losses!(fluid_containers, params_wall_sheath)

        # Neutrals should recombine at walls
        @test continuity[1].dens_ddt[in_index] ≈ -(isothermal[1].dens_ddt[in_index] + isothermal[2].dens_ddt[in_index])
        @test continuity[1].dens_ddt[out_index] ≈ 0.0

        for ion in isothermal
            # Rate of ion loss is equal to ni νiw
            @test ion.dens_ddt[in_index] ≈ -ion.density[in_index] * νiw * sqrt(ion.species.Z)
            # ion momentum loss is equal to Iiw / e / V_cell * ui
            @test ion.mom_ddt[in_index] ≈ -ion.momentum[in_index] * νiw * sqrt(ion.species.Z)

            # No losses outside of channel
            @test ion.dens_ddt[out_index] ≈ 0.0
            @test ion.mom_ddt[out_index] ≈ 0.0
        end
    end
end

test_wall_loss_serialization()
test_electron_losses()
test_ion_losses()
