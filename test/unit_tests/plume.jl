using Test
using HallThruster: HallThruster as het

@testset "Plume geometry starts at the physical exit" begin
    thruster = het.SPT_100
    geometry = thruster.geometry
    L = geometry.channel_length
    temperature_exit = 10.0 + 100.0 * L
    slope = 0.3
    velocity = sqrt((5 / 3) * het.e * temperature_exit / het.Xenon.m) / slope

    # Include a domain ending inside the channel, several exit alignments, and
    # a long domain where the annular hole closes.
    for edges in (
            range(0.0, 0.01; length = 8),
            range(0.0, 0.08; length = 14),
            range(0.0, 0.08; length = 18),
            [0.0, 0.005, 0.02, 0.03, 0.07, 0.15, 0.4],
        )
        grid = het.Grid1D(edges)
        cache = het.allocate_arrays(grid.num_cells, 0)
        @. cache.Tev = 10.0 + 100.0 * grid.cell_centers
        cache.m_eff .= het.Xenon.m
        cache.avg_ion_vel .= velocity
        params = (; cache, grid, thruster)
        het.update_plume_geometry!(params)
        for (i, z) in enumerate(grid.cell_centers)
            spread = slope * max(0.0, z - L)
            inner = max(0.0, geometry.inner_radius - spread)
            outer = geometry.outer_radius + spread
            @test cache.inner_radius[i] ≈ inner
            @test cache.outer_radius[i] ≈ outer
            @test cache.channel_area[i] ≈ π * (outer^2 - inner^2)
            derivative = z <= L ? 0.0 : 2π * slope * (inner + outer)
            @test cache.dA_dz[i] ≈ derivative
            @test cache.dlnA_dz[i] ≈ derivative / cache.channel_area[i]
            @test cache.tanδ[i] ≈ (z <= L ? 0.0 : slope)
        end
        # Repeated calls must rebuild the shape, without accumulating spread.
        area = copy(cache.channel_area)
        het.update_plume_geometry!(params)
        @test cache.channel_area == area
    end
end

@testset "Heavy-species stages conserve mass and momentum in a changing plume" begin
    for gas in (het.Xenon, het.Krypton)
        config = het.Config(
            thruster = het.SPT_100, domain = (0.0, 0.08), discharge_voltage = 300.0,
            propellants = [het.Propellant(gas, 5.0e-6; max_charge = 2, ion_temperature_K = 0.0)],
            solve_plume = true, wall_loss_model = het.NoWallLosses(),
            excitation_model = :None, anode_boundary_condition = :dirichlet,
        )
        sim = het.SimParams(
            grid = het.UnevenGrid(24), dt = 1.0e-8, duration = 1.0e-8,
            num_save = 2, verbose = false,
        )
        params = het.setup_simulation(config, sim)
        (; cache, grid, fluid_containers, fluid_array) = params
        for (j, fluid) in enumerate(fluid_array)
            @. fluid.density = gas.m * j * 1.0e18 * (1 + grid.cell_centers)
            if fluid.species.Z != 0
                @. fluid.momentum = fluid.density * (1.0e4 + 2.0e4 * grid.cell_centers)
            end
        end
        cache.Tev .= 20.0
        cache.∇ϕ .= 0.0
        cache.K .= 0.0
        het.update_heavy_species!(params)
        @. cache.nϵ = 1.5 * cache.ne * cache.Tev
        het.update_plume_geometry!(params)
        neutral = only(fluid_containers.continuity)
        ion = first(fluid_containers.isothermal)
        rate = 1.0e-14
        cost = 14.0
        reactions = [het.ElectronImpactReaction(cost, gas(0), [gas(1)], fill(rate, 256))]
        empty!(params.reaction_groups)
        append!(
            params.reaction_groups,
            het.build_electron_impact_groups(
                reactions, het.reactant_indices(reactions, fluid_array),
                het.product_indices(reactions, fluid_array), fluid_array,
            ),
        )

        cells = het.interior_cells(grid.cell_centers)
        total_mass() = sum(
            cache.channel_area[i] * grid.dz_cell[i] * sum(f.density[i] for f in fluid_array)
                for i in cells
        )
        total_momentum() = sum(
            cache.channel_area[i] * grid.dz_cell[i] * sum(f.momentum[i] for f in fluid_containers.isothermal)
                for i in cells
        )
        mass_before = total_mass()
        momentum_before = total_momentum()
        temperature_before = copy(cache.Tev)
        area_before = copy(cache.channel_area)
        mass_rates = Float64[]
        momentum_rates = Float64[]
        energy_losses = Float64[]
        stage_areas = Vector{Float64}[]
        acceleration = 1.0e10
        source = function (fluids, _)
            @test cache.Tev == temperature_before
            push!(stage_areas, copy(cache.channel_area))
            left_area = het.face_area(grid, cache.channel_area, 1)
            right_area = het.face_area(grid, cache.channel_area, length(grid.edges))
            push!(
                mass_rates, sum(
                    left_area * f.flux_dens[1] - right_area * f.flux_dens[end] for f in fluid_array
                )
            )
            momentum_rate = sum(
                left_area * f.flux_mom[1] - right_area * f.flux_mom[end] for f in fluids.isothermal
            )
            loss = 0.0
            for i in cells
                volume = cache.channel_area[i] * grid.dz_cell[i]
                ne = sum(f.species.Z * f.density[i] / gas.m for f in fluids.isothermal)
                mass_source = rate * ne * neutral.density[i]
                momentum_rate += volume * mass_source * neutral.vel_prim[i]
                loss += volume * mass_source / gas.m * cost
                for f in fluids.isothermal
                    f.mom_ddt[i] += acceleration * f.density[i]
                    momentum_rate += volume * acceleration * f.density[i]
                end
            end
            push!(momentum_rates, momentum_rate)
            push!(energy_losses, loss)
            return nothing
        end
        dt = sim.dt
        @test !het.step_heavy_species!(fluid_containers, params, source, dt; solve_plume = true)
        @test length(stage_areas) == 2
        @test maximum(abs, stage_areas[2] .- area_before) > 1.0e-8
        @test maximum(abs, cache.channel_area .- stage_areas[2]) > 1.0e-10
        @test total_mass() ≈ mass_before + 0.5 * dt * sum(mass_rates) rtol = 1.0e-13
        @test total_momentum() ≈ momentum_before + 0.5 * dt * sum(momentum_rates) rtol = 1.0e-13
        integrated_loss = sum(
            cache.channel_area[i] * grid.dz_cell[i] * cache.inelastic_losses[i] for i in cells
        )
        @test integrated_loss ≈ 0.5 * sum(energy_losses) rtol = 1.0e-13
        @test cache.channel_area_old == area_before
        @test cache.nϵ .* cache.channel_area ≈ cache.plume_energy_old .* area_before
        @test all(>(0.0), ion.density)
    end
end

@testset "Uniform pressure stays balanced in a tapered plume" begin
    grid = het.Grid1D([0.0, 0.01, 0.04, 0.08])
    propellant = het.Propellant(
        het.Xenon, 0.0; velocity_m_s = 0.0, ion_temperature_K = 1000.0,
    )
    fluids = het.allocate_fluids(propellant, grid)
    for fluid in [fluids.continuity; fluids.isothermal]
        fluid.density .= het.Xenon.m * 1.0e18
    end
    area = @. 0.004 * exp(10 * grid.cell_centers)
    het.update_convective_terms!(fluids, grid, true, area)
    for fluid in [fluids.continuity; fluids.isothermal]
        @test fluid.dens_ddt ≈ zeros(grid.num_cells) atol = 1.0e-14
    end
    @test only(fluids.isothermal).mom_ddt ≈ zeros(grid.num_cells) atol = 1.0e-14
end

@testset "Temporal expansion performs electron pressure work" begin
    cache = het.allocate_arrays(5, 0)
    for landmark in (false, true), factor in (0.5, 1.0, 2.0)
        cache.channel_area_old .= [1.0, 1.2, 1.5, 2.0, 3.0]
        cache.channel_area .= factor .* cache.channel_area_old
        cache.plume_energy_old .= 30.0
        cache.ne .= 1.0 / factor
        het.apply_plume_energy_change!(cache, landmark)
        gamma = landmark ? 2.0 : 5 / 3
        # For an isolated adiabatic change, W/n^gamma is invariant and Te
        # decreases on expansion. Pure dilution alone would keep Te constant.
        @test cache.nϵ ./ cache.ne .^ gamma ≈ fill(30.0, 5)
        @test cache.nϵ ./ (1.5 .* cache.ne) ≈ fill(20.0 * factor^(1 - gamma), 5)
        cache.plume_energy_old .= cache.nϵ
        cache.channel_area_old .= cache.channel_area
        cache.channel_area ./= factor
        het.apply_plume_energy_change!(cache, landmark)
        @test cache.nϵ ≈ fill(30.0, 5)
    end
end

@testset "Heavy integration applies electron expansion and compression work" begin
    for acceleration in (-1.0e11, 1.0e11)
        config = het.Config(
            thruster = het.SPT_100, domain = (0.0, 0.08), discharge_voltage = 300.0,
            anode_mass_flow_rate = 5.0e-6, solve_plume = true,
        )
        sim = het.SimParams(
            grid = het.EvenGrid(16), dt = 1.0e-8, duration = 1.0e-8,
            num_save = 2, verbose = false,
        )
        params = het.setup_simulation(config, sim)
        (; cache, grid, fluid_array, fluid_containers) = params
        empty!(params.reaction_groups)
        cache.∇ϕ .= 0.0
        cache.Tev .= 40.0
        for fluid in fluid_array
            fluid.density .= 1.0e18 * fluid.species.element.m
            if fluid.species.Z != 0
                fluid.momentum .= 1.0e4 .* fluid.density
            end
        end
        het.update_heavy_species!(params)
        @. cache.nϵ = 1.5 * cache.ne * cache.Tev
        het.update_plume_geometry!(params)
        area_before = copy(cache.channel_area)
        ne_before = copy(cache.ne)
        energy_before = copy(cache.nϵ)
        source = function (fluids, _)
            # Isolate temporal geometry change by canceling axial transport.
            for fluid in fluid_array
                fill!(fluid.dens_ddt, 0.0)
            end
            for fluid in fluids.isothermal
                @. fluid.mom_ddt = acceleration * fluid.density
            end
            return nothing
        end
        @test !het.integrate_heavy_species!(fluid_containers, params, source, sim.dt)
        cells = het.interior_cells(grid.cell_centers)
        ratio = area_before ./ cache.channel_area
        @test maximum(abs, ratio .- 1.0) > 1.0e-3
        @test (cache.ne .* cache.channel_area)[cells] ≈ (ne_before .* area_before)[cells]
        @test cache.Tev[cells] ≈ 40.0 .* ratio[cells] .^ (2 / 3)
        @test cache.nϵ[cells] ≈ energy_before[cells] .* ratio[cells] .^ (5 / 3)
        @test cache.pe[cells] ≈ cache.ne[cells] .* cache.Tev[cells]
        ion = only(fluid_containers.isothermal)
        @test ion.vel_prim[cells] ≈ fill(1.0e4 + acceleration * sim.dt, length(cells))
    end
end

@testset "Saved plume states share one geometry" begin
    config = het.Config(
        thruster = het.SPT_100, domain = (0.0, 0.08), discharge_voltage = 300.0,
        anode_mass_flow_rate = 5.0e-6, solve_plume = true,
    )
    sim = het.SimParams(
        grid = het.EvenGrid(24), dt = 1.0e-9, duration = 1.0e-8,
        num_save = 4, adaptive = false, verbose = false,
    )
    sol = het.run_simulation(config, sim)
    @test sol.retcode == :success
    for frame in sol.frames
        @test maximum(frame.channel_area) > first(frame.channel_area)
        # The two endpoint values are independently interpolated to faces in
        # Frame, so products of those endpoint averages need not obey Ohm's law.
        cells = 2:(length(sol.grid) - 1)
        @test ((frame.ji .- het.e .* frame.ne .* frame.ue) .* frame.channel_area)[cells] ≈
            fill(frame.discharge_current[], length(cells))
        @test all(isfinite, frame.Tev)
        @test all(>=(config.min_Te), frame.Tev)
    end
end
