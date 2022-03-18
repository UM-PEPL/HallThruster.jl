get_species(sim) = [Species(sim.propellant, i) for i in 0:(sim.ncharge)]

function configure_simulation(sim)
    fluids = sim.fluids
    species = [fluids[i].species for i in 1:length(fluids)]
    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(Symbol(fluid.species) => fluid_range
                              for (fluid, fluid_range) in zip(fluids, fluid_ranges))

    return species, fluids, fluid_ranges, species_range_dict
end

function allocate_arrays(grid, fluids) #rewrite allocate arrays as function of set of equations, either 1, 2 or 3
    # Number of variables in the state vector U
    nvariables = 0
    for i in 1:length(fluids)
        if fluids[i].conservation_laws.type == _ContinuityOnly
            nvariables += 1
        elseif fluids[i].conservation_laws.type == _IsothermalEuler
            nvariables += 2
        elseif fluids[i].conservation_laws.type == _EulerEquations
            nvariables += 3
        end
    end

    ncells = grid.ncells + 2
    nedges = grid.ncells + 1

    U = zeros(nvariables + 1, ncells)
    A = Tridiagonal(ones(nedges-1), ones(nedges), ones(nedges-1)) #for potential
    b = zeros(nedges) #for potential equation
    Aϵ = Tridiagonal(ones(ncells-1), ones(ncells), ones(ncells-1)) #for energy
    bϵ = zeros(ncells) #for energy
    B = zeros(ncells)
    νan = zeros(ncells)
    νc = zeros(ncells)
    νei = zeros(ncells)
    νen = zeros(ncells)
    νw = zeros(ncells)
    μ = zeros(ncells)
    ϕ = zeros(nedges)
    ϕ_cell = zeros(ncells)
    ∇ϕ = zeros(ncells)
    ne = zeros(ncells)
    Tev = zeros(ncells)
    pe = zeros(ncells)
    ∇pe = zeros(ncells)
    ue = zeros(ncells)

    cache = (; A, b, Aϵ, bϵ, B, νan, νc, μ, ϕ, ϕ_cell, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, νw)
    return U, cache
end

function update!(dU, U, p, t)
    update_heavy_species!(dU, U, p, t)
    #update_electron_energy_explicit!(dU, U, p, t)
end

left_edge(i) = i - 1
right_edge(i) = i

function electron_density(U, params, i)
    ne = 0.0
    index = params.index
    @inbounds for Z in 1:params.config.ncharge
        ne += Z * U[index.ρi[Z], i] / params.config.propellant.m
    end
    return ne
end

function inlet_neutral_density(sim)
    un = sim.neutral_velocity
    A = channel_area(sim.geometry)
    m_atom = sim.propellant.m
    nn = sim.inlet_mdot / un / A / m_atom
    return nn
end

function precompute_bfield!(B, zs)
    B_max = 0.015
    L_ch = 0.025
    for (i, z) in enumerate(zs)
        B[i] = B_field(B_max, z, L_ch)
    end
end

function make_keys(fluid_range, subscript)
    len = length(fluid_range)
    if len == 1
        return (Symbol("ρ$(subscript)"))
    elseif len == 2
        return (
            Symbol("ρ$(subscript)"),
            Symbol("ρ$(subscript)u$(subscript)")
        )
    elseif len == 3
        return (
            Symbol("ρ$(subscript)"),
            Symbol("ρ$(subscript)u$(subscript)"),
            Symbol("ρ$(subscript)E$(subscript)")
        )
    else
        throw(ArgumentError("Too many equations on fluid (this should be unreachable)"))
    end
end

function configure_index(fluid_ranges)
    lf = fluid_ranges[end][end]

    ncharge = length(fluid_ranges)-1
    solve_ion_temp = length(fluid_ranges[2]) == 3

    keys_neutrals = (:ρn, )
    values_neutrals = (1, )

    if solve_ion_temp
        keys_ions = (:ρi, :ρiui, :ρiuiEi)
        values_ions = (
            [f[1] for f in fluid_ranges[2:end]]...,
            [f[2] for f in fluid_ranges[2:end]]...,
            [f[3] for f in fluid_ranges[2:end]]...,
        )
    else
        keys_ions = (:ρi, :ρiui)
        values_ions = (
            [f[1] for f in fluid_ranges[2:end]],
            [f[2] for f in fluid_ranges[2:end]],
        )
    end

    keys_fluids = (keys_neutrals..., keys_ions...)
    values_fluids = (values_neutrals..., values_ions...)
    keys_electrons = (:nϵ, :Tev, :ne, :pe, :ϕ, :grad_ϕ, :ue)
    values_electrons = lf .+ collect(1:7)
    index_keys = (keys_fluids..., keys_electrons..., :lf)
    index_values = (values_fluids..., values_electrons..., lf)
    index = NamedTuple{index_keys}(index_values)
    @show index
    return index
end

function configure_fluids(config)
    propellant = config.propellant
    species = [propellant(i) for i in 0:config.ncharge]
    neutral_fluid = Fluid(species[1], ContinuityOnly(u = config.neutral_velocity, T = config.neutral_temperature))
    ion_eqns = if config.solve_ion_energy
        EulerEquations()
    else
        IsothermalEuler(T = config.ion_temperature)
    end
    ion_fluids = [Fluid(species[i+1], ion_eqns) for i in 1:config.ncharge]
    fluids = [neutral_fluid; ion_fluids]
    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(Symbol(fluid.species) => fluid_range
                              for (fluid, fluid_range) in zip(fluids, fluid_ranges))
    return fluids, fluid_ranges, species, species_range_dict
end

function initial_condition!(U, z_cell, IC!, fluids)
    for (i, z) in enumerate(z_cell)
        @views IC!(U[:, i], z, fluids, z_cell[end])
    end
end

function run_simulation(sim, config, alg) #put source and Bcs potential in params

    fluids, fluid_ranges, species, species_range_dict = configure_fluids(config)

    index = configure_index(fluid_ranges)
    landmark = load_landmark()

    use_restart = config.restart_file !== nothing

    if use_restart
        U, grid, B = read_restart(config.restart_file)
        _, cache = allocate_arrays(grid, fluids)
        cache.B .= B
    else
        grid = sim.grid
        U, cache = allocate_arrays(grid, fluids)
        initial_condition!(@views(U[1:index.nϵ, :]), grid.cell_centers, sim.initial_condition, fluids)
        precompute_bfield!(cache.B, grid.cell_centers)
    end

    scheme = sim.scheme
    source_term! = sim.source_term!
    timestep = sim.timestepcontrol[1]
    adaptive = sim.timestepcontrol[2]
    tspan = (0.0, sim.end_time)

    # Load ionization reactions fro file
    if config.ionization_coeffs == :LANDMARK
        if config.ncharge > 1
            throw(ArgumentError("LANDMARK ionization table does not support multiply-charged ions. Please use :BOLSIG or reduce ncharge to 1."))
        else
            ionization_reactions = [IonizationReaction(species[1], species[2], landmark.rate_coeff)]
        end
    elseif config.ionization_coeffs == :BOLSIG
        ionization_reactions = load_ionization_reactions(species)
    elseif config.ionization_coeffs == :BOLSIG_FIT
		ionization_reactions = ionization_fits_Xe(config.ncharge)
		loss_coeff = loss_coeff_fit
    else
        throw(ArgumentError("Invalid ionization reactions selected. Please choose either :LANDMARK or :BOLSIG"))
    end

    BCs = sim.boundary_conditions

    precompute_bfield!(cache.B, grid.cell_centers)

    loss_coeff = config.ionization_coeffs == :BOLSIG_FIT ? loss_coeff_fit : landmark.loss_coeff

    # callback for calling the update_values! function at each timestep
    update_callback = DiscreteCallback(Returns(true), update_values!, save_positions=(false,false))

    # Choose which cache variables to save and set up saving callback
    fields_to_save = (:μ, :Tev, :ϕ, :∇ϕ, :ne, :pe, :ue, :∇pe, :νan, :νc, :νen, :νei, :νw, :ϕ_cell)

    function save_func(u, t, integrator)
        (; μ, Tev, ϕ, ∇ϕ, ne, pe, ue, ∇pe, νan, νc, νen, νei, νw, ϕ_cell) = integrator.p.cache
        return deepcopy((; μ, Tev, ϕ, ∇ϕ, ne, pe, ue, ∇pe, νan, νc, νen, νei, νw, ϕ_cell))
    end

    saved_values = SavedValues(
        Float64, NamedTuple{fields_to_save, NTuple{length(fields_to_save), Vector{Float64}}}
    )
    saving_callback = SavingCallback(save_func, saved_values, saveat = sim.saveat)


    niters = round(Int, tspan[2] / timestep)

    # Progress callback
    if config.progress_interval > 0
        progress_bar = ProgressBar(N = niters ÷ config.progress_interval, description = "Running simulation")
    else
        progress_bar = nothing
    end

    # Assemble callback set
    if sim.callback !== nothing
        callbacks = CallbackSet(update_callback, saving_callback, sim.callback)
    else
        callbacks = CallbackSet(update_callback, saving_callback)
    end

    # Simulation parameters
    params = (;
        config = config,
        propellant = config.propellant,
        ϕ_L = config.anode_potential,
        ϕ_R = config.cathode_potential,
        Te_L = config.anode_Te,
        Te_R = config.cathode_Te,
        L_ch = config.geometry.channel_length,
        A_ch = channel_area(config.geometry.outer_radius, config.geometry.inner_radius),
        αϵ = config.radial_loss_coeffs,
        αw = config.wall_collision_coeff,
        δ = config.ion_diffusion_coeff,
        un = config.neutral_velocity,
        Tn = config.neutral_temperature,
        mdot_a = config.anode_mass_flow_rate,
        OVS = config.verification,
        anom_model = config.anom_model,
        loss_coeff = loss_coeff,
        reactions = ionization_reactions,
        implicit_energy = config.implicit_energy,
        index, cache, fluids, fluid_ranges, species_range_dict, z_cell=grid.cell_centers,
        z_edge=grid.edges, cell_volume=grid.cell_volume, source_term!,
        scheme, BCs, dt=timestep, progress_bar,
        iteration = [-1]
    )

    # Compute maximum allowed iterations
    maxiters = Int(ceil(1000 * tspan[2] / timestep))

    # Choose which function to use for ODE
    # If implicit, we update electron energy in the callback, not using diffeq
    if config.implicit_energy > 0
        f = ODEFunction(update_heavy_species!)
    else
	    f = ODEFunction(update!)
    end

    #make values in params available for first timestep
    update_values!(U, params)

    # Set up ODE problem and solve
    prob = ODEProblem{true}(f, U, tspan, params)
	sol = try
        solve(
            prob, alg; saveat=sim.saveat, callback=callbacks,
            adaptive=adaptive, dt=timestep, dtmax=10*timestep, dtmin = timestep/10, maxiters = maxiters,
	    )
    catch e
        if config.progress_interval > 0
            stop(progress_bar)
        end
        println("There was an error")
        throw(e)
    end

    # Print some diagnostic information
    if sol.retcode == :NaNDetected
        println("Simulation failed with NaN detected at t = $(sol.t[end])")
    elseif sol.retcode == :InfDetected
        println("Simulation failed with Inf detected at t = $(sol.t[end])")
    end

    if config.progress_interval > 0
        stop(progress_bar)
    end

    # Return the solution
    return HallThrusterSolution(sol, params, saved_values.saveval)
end