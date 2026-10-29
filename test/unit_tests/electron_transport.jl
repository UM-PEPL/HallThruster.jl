using Test
using HallThruster: HallThruster as het

@testset "Conductive flux uses temperature gradients" begin
    grid = het.Grid1D([0.0, 0.1, 0.3, 0.6, 1.0])
    cache = het.allocate_arrays(grid.num_cells, 0)
    cache.ne .= [2.0, 4.0, 3.0, 8.0, 5.0, 7.0]
    temperature = @. 10 + 2 * grid.cell_centers
    @. cache.nϵ = 1.5 * cache.ne * temperature
    cache.κ .= 3.0
    het.energy_flux_coefficients!(grid, cache, :dirichlet)
    flux = @. cache.energy_flux_left * cache.nϵ[1:(end - 1)] +
        cache.energy_flux_right * cache.nϵ[2:end]
    @test flux ≈ fill(-6.0, length(grid.edges))

    # The LANDMARK mean-energy law remains q = -10/9 μneϵ ∇ϵ.
    cache.μ .= 0.7
    cache.Tev .= temperature
    het.LANDMARK_conductivity()(cache.κ, (; cache))
    het.energy_flux_coefficients!(grid, cache, :dirichlet)
    for face in eachindex(grid.edges)
        w = (grid.edges[face] - grid.cell_centers[face]) / grid.dz_edge[face]
        mean_ne_temperature = (1 - w) * cache.ne[face] * temperature[face] +
            w * cache.ne[face + 1] * temperature[face + 1]
        expected = -(10 / 9) * 0.7 * 1.5 * mean_ne_temperature * 3.0
        actual = cache.energy_flux_left[face] * cache.nϵ[face] +
            cache.energy_flux_right[face] * cache.nϵ[face + 1]
        @test actual ≈ expected
    end
end

@testset "Anode sheath energy transmission" begin
    grid = het.Grid1D([0.0, 0.5, 1.0])
    cache = het.allocate_arrays(grid.num_cells, 0)
    cache.ne .= 10.0
    cache.nϵ .= 1.5 * 10.0 * 6.0
    cache.ue .= -0.5
    for voltage in (0.0, 3.0, 6.0, 30.0)
        cache.Vs[] = voltage
        het.energy_flux_coefficients!(grid, cache, :sheath)
        flux = cache.energy_flux_left[1] * cache.nϵ[1] +
            cache.energy_flux_right[1] * cache.nϵ[2]
        @test flux ≈ -5.0 * (12.0 + voltage)
    end
end

@testset "Closed tapered domain conserves electron energy at flow reversal" begin
    for implicit in (0.0, 0.5, 1.0), variable_area in (false, true)
        grid = het.Grid1D([0.0, 0.1, 0.3, 0.6, 1.0])
        cache = het.allocate_arrays(grid.num_cells, 0)
        cache.ne .= 1.0
        cache.nϵ .= [10.0, 10.0, 30.0, 15.0, 20.0, 20.0]
        cache.ue .= [0.0, 0.0, 4.0, -4.0, 0.0, 0.0]
        cache.κ .= [0.0, 0.0, 3.0, 4.0, 0.0, 0.0]
        cache.channel_area .= variable_area ? [1.0, 1.1, 1.7, 2.5, 4.0, 4.3] : ones(6)
        volume = cache.channel_area .* grid.dz_cell
        before = sum((volume .* cache.nϵ)[2:(end - 1)])
        # Insulated endpoints have zero velocity and conductivity. Ghost values
        # are extrapolated, but the interior solve determines the energy budget.
        A = cache.Aϵ
        A.d[1] = A.d[end] = 1.0
        A.du[1] = A.dl[end] = -1.0
        cache.bϵ[1] = cache.bϵ[end] = 0.0
        het.setup_energy_system!(A, cache.bϵ, grid, cache, :dirichlet, implicit, 0.002)
        het.tridiagonal_solve!(cache.nϵ, A, cache.bϵ)
        after = sum((volume .* cache.nϵ)[2:(end - 1)])
        @test after ≈ before rtol = 1.0e-13
    end
end

@testset "Constant conductive power in a tapered duct" begin
    grid = het.Grid1D(range(0.0, 1.0; length = 101))
    cache = het.allocate_arrays(grid.num_cells, 0)
    cache.channel_area .= 1 .+ grid.cell_centers
    cache.ne .= 1.0
    cache.κ .= 3.0
    @. cache.nϵ = 1.5 * (10 + log(cache.channel_area))
    before = copy(cache.nϵ)
    het.setup_energy_system!(cache.Aϵ, cache.bϵ, grid, cache, :dirichlet, 0.0, 1.0)
    @test maximum(abs, cache.bϵ[2:(end - 1)] .- before[2:(end - 1)]) < 1.0e-3
end

@testset "Electron-ion scattering adds charge-squared populations" begin
    grid = het.Grid1D([0.0, 0.5, 1.0])
    propellants = [
        het.Propellant(het.Xenon, 0.0; allowed_charges = [1, 2], excited_ion_levels = Dict(1 => [1])),
        het.Propellant(het.Hydrogen, 0.0; allowed_charges = [-1]),
    ]
    groups = [het.allocate_fluids(p, grid) for p in propellants]
    ions = [fluid for group in groups for fluid in group.isothermal]
    cache = het.allocate_arrays(grid.num_cells, 0)
    cache.Tev .= 8.0
    # Start with pure Xe2+: ne = 2 ni, not ni.
    ions[3].density .= 1.0e17 * het.Xenon.m
    fluids = het.FluidContainerSet([group.continuity[1] for group in groups], ions)
    het.update_heavy_species_cache!(fluids, cache, false, 1.0)
    het.freq_electron_ion!(cache.νei, cache.ne, cache.Tev, ions)
    expected = 2.9e-12 * 4.0e17 * het.coulomb_logarithm(2.0e17, 8.0, 2) / 8.0^1.5
    @test cache.νei[2] ≈ expected
    @test cache.Z_eff[2] ≈ 2.0

    # Keep many scattering ions while positive/negative charge nearly cancels.
    ions[3].density .= 0.0
    ions[1].density .= 2.0e17 * het.Xenon.m
    ions[4].density .= 1.9e17 * het.Hydrogen.m
    het.update_heavy_species_cache!(fluids, cache, false, 1.0)
    het.freq_electron_ion!(cache.νei, cache.ne, cache.Tev, ions)
    expected = 2.9e-12 * 3.9e17 * het.coulomb_logarithm(cache.ne[2], 8.0, 1) / 8.0^1.5
    @test cache.νei[2] ≈ expected
    @test cache.Z_eff[2] ≈ 39.0
    # Splitting a population among electronic levels cannot alter scattering.
    ions[1].density .*= 0.5
    ions[2].density .= ions[1].density
    het.freq_electron_ion!(cache.νei, cache.ne, cache.Tev, ions)
    @test cache.νei[2] ≈ expected
end

@testset "Repelled negative ions cannot subtract wall losses" begin
    grid = het.Grid1D([0.0, 0.5, 1.0])
    xe = het.allocate_fluids(het.Propellant(het.Xenon, 0.0), grid)
    hydrogen = het.allocate_fluids(het.Propellant(het.Hydrogen, 0.0; allowed_charges = [-1]), grid)
    xe.isothermal[1].density .= 2.0e17 * het.Xenon.m
    hydrogen.isothermal[1].density .= 1.0e17 * het.Hydrogen.m
    cache = het.allocate_arrays(grid.num_cells, 0)
    cache.ne .= 1.0e17
    cache.Tev .= 4.0
    cache.m_eff .= het.Xenon.m
    cache.wall_transition .= 1.0
    params = (;
        grid, cache, thruster = het.SPT_100, plume_loss_scale = 1.0,
        fluid_containers = het.FluidContainerSet([xe.continuity; hydrogen.continuity], [xe.isothermal; hydrogen.isothermal]),
    )
    model = het.WallSheath(het.BoronNitride, 0.5)
    het.freq_electron_wall!(cache.radial_loss_frequency, model, params)
    @test all(>(0), cache.radial_loss_frequency)
    reference = copy(cache.radial_loss_frequency)
    hydrogen.isothermal[1].density .= 0.0
    het.freq_electron_wall!(cache.radial_loss_frequency, model, params)
    @test cache.radial_loss_frequency ≈ reference
    het.wall_power_loss!(cache.wall_losses, model, params)
    @test all(>(0), cache.wall_losses[2:(end - 1)])
end

@testset "Negative ions below the anode barrier reflect without mass leakage" begin
    grid = het.Grid1D([0.0, 0.5, 1.0])
    propellant = het.Propellant(het.Xenon, 0.0; allowed_charges = [-1, 1])
    fluids = het.allocate_fluids(propellant, grid)
    negative, positive = fluids.isothermal
    negative.density .= 1.0e17 * het.Xenon.m
    positive.density .= 2.0e17 * het.Xenon.m
    cache = het.allocate_arrays(grid.num_cells, 0)
    cache.Tev .= 2.0
    cache.channel_area .= 1.0
    cache.Vs[] = 10.0
    for velocity in (-100.0, -1000.0), reconstruct in (false, true)
        negative.momentum .= velocity .* negative.density
        het.apply_left_boundary!(fluids, propellant, cache, :sheath, 0.0, false, 1.0)
        het.update_primitive_velocity!(negative)
        het.compute_edge_states_isothermal!(negative, reconstruct)
        het.compute_fluxes_isothermal!(negative, grid)
        @test negative.momentum[1] == -negative.momentum[2]
        @test negative.flux_dens[1] == 0.0
    end

    # Both RK evaluations must use a reflected state even as reactions and
    # transport change the adjacent interior population in the predictor.
    params = (;
        cache, grid, ion_wall_losses = false, reconstruct = true,
        landmark = false, min_number_density = 1.0,
        reaction_groups = het.ElectronImpactGroup[],
        ion_production_cost_multiplier = Dict{Symbol, Float64}(),
        fluid_array = [fluids.continuity; fluids.isothermal],
        fluids_by_propellant = [fluids], simulation = (; CFL = 0.5),
    )
    boundary_fluxes = Float64[]
    source = (containers, _) -> push!(boundary_fluxes, containers.isothermal[1].flux_dens[1])
    @test !het.step_heavy_species!(fluids, params, source, 1.0e-7)
    @test boundary_fluxes == [0.0, 0.0]
end

@testset "Electrical boundary conditions with an output current limiter" begin
    for type in (:NoCircuit, :LoadResistor, :ParallelRL, :ParallelRC, :ParallelRLSeriesC),
            limit in (Inf, 0.05, 0.0)
        grid = het.Grid1D([0.0, 0.02, 0.05, 0.08])
        cache = het.allocate_arrays(grid.num_cells, 0)
        cache.ne .= 1.0e18
        cache.Tev .= 4.0
        cache.channel_area .= 1.0e-3
        cache.B .= 0.01
        cache.νc .= 1.0e7
        cache.νan .= 1.0e7
        circuit = het.CircuitModel(type; R = 20.0, L = 0.01, C = 1.0e-3, limit_current = limit)
        for element in circuit.elements
            element.type == 'C' && (element.V = 300.0)
        end
        params = (;
            cache, grid, anom_smoothing_iters = 0, landmark = false,
            anode_bc = :dirichlet, iteration = [10], min_number_density = 1.0,
            cathode_coupling_voltage = 20.0, discharge_voltage = 300.0,
            filter_circuit = circuit, dt = [1.0e-9],
            simulation = (; current_control = het.NoController()),
        )
        het.update_electrical_vars!(params)
        @test (cache.ϕ[1] + cache.ϕ[2]) / 2 ≈ cache.Vd[]
        @test (cache.ϕ[end - 1] + cache.ϕ[end]) / 2 ≈ 20.0 atol = 1.0e-10
        @test cache.Id[] <= limit
        if isfinite(limit)
            @test cache.Id[] ≈ limit atol = 1.0e-12
            @test cache.Vd[] < 300.0
        end
    end
end
