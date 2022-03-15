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

    ncells = grid.ncells
    nedges = grid.ncells + 1

    U = zeros(nvariables + 1, ncells + 2) # need to allocate room for ghost cells
    A = Tridiagonal(ones(nedges-1), ones(nedges), ones(nedges-1)) #for potential
    b = zeros(nedges) #for potential equation
    Aϵ = Tridiagonal(ones(ncells+1), ones(ncells+2), ones(ncells+1)) #for energy
    bϵ = zeros(ncells+2) #for energy
    B = zeros(ncells + 2)
    νan = zeros(ncells + 2)
    νc = zeros(ncells + 2)
    μ = zeros(ncells + 2)
    ϕ = zeros(nedges)
    ∇ϕ = zeros(ncells + 2)
    ne = zeros(ncells + 2)
    Tev = zeros(ncells + 2)
    pe = zeros(ncells + 2)
    ue = zeros(ncells + 2)
    BC_L = zeros(nvariables+1)
    BC_R = zeros(nvariables+1)

    cache = (; BC_L, BC_R,  A, b, Aϵ, bϵ, B, νan, νc, μ, ϕ, ∇ϕ, ne, Tev, pe, ue)
    return U, cache
end

function update_heavy_species!(dU, U, params, t)
    ####################################################################
    #extract some useful stuff from params

    (;index, z_edge, propellant, scheme, fluids, species_range_dict, reactions) = params
    (;ue, μ, BC_L, BC_R) = params.cache

    ncells = size(U, 2) - 2

    mi = propellant.m

    ##############################################################
    #FLUID MODULE

    #fluid BCs now in update_values struct

    ncharge = params.config.ncharge

    un, Tn, γn, Rn = fluids[1].conservation_laws.u, fluids[1].conservation_laws.T, γ(fluids[1]), R(fluids[1])
    Ti = fluids[2].conservation_laws.T

    @views ρn = U[index.ρn, :]
    @views nϵ = U[index.nϵ, :]

    coupled = params.config.electron_pressure_coupled


    first_ind = 2
    last_ind = ncells+1

    # Compute heavy species source terms
    @inbounds for i in first_ind:last_ind

        #Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        ρn_L = i == first_ind ? SA[BC_L[index.ρn]] : SA[U[index.ρn, i-1]]
        ρn_R = i == last_ind ? SA[BC_R[index.ρn]] : SA[U[index.ρn, i+1]]
        ρn_0 = SA[U[index.ρn, i]]

        Fn_L = scheme.flux_function(ρn_L, ρn_0, fluids[1], coupled)
        Fn_R = scheme.flux_function(ρn_0, ρn_R, fluids[1], coupled)

        dU[index.ρn, i] = -(Fn_R[1] - Fn_L[1])/Δz
        neL = 0.0
        neR = 0.0
        ne0 = 0.0
        for Z in 1:ncharge
            if i == first_ind
                neL += Z * BC_L[index.ρi[Z]]/mi
            else
                neL += Z * U[index.ρi[Z], i-1] / mi
            end

            ne0 += Z * U[index.ρi[Z], i] / mi

            if i == last_ind
                neR += Z * BC_R[index.ρi[Z]]/mi
            else
                neR += Z * U[index.ρi[Z], i+1] / mi
            end
        end

        ϵL = i == first_ind ? params.Te_L : max(0.0, nϵ[i-1]/neL)
        ϵ0 = max(0.0, nϵ[i]/ne0)
        ϵR = i == last_ind ? params.Te_R : max(0.0, nϵ[i+1]/neR)

        for Z in 1:ncharge

            UL = if i == first_ind
                SA[BC_L[index.ρi[Z]], BC_L[index.ρiui[Z]]]
            else
                SA[U[index.ρi[Z], i-1], U[index.ρiui[Z], i-1]]
            end

            U0 = SA[U[index.ρi[Z], i], U[index.ρiui[Z], i]]

            UR = if i == last_ind
                SA[BC_R[index.ρi[Z]], BC_R[index.ρiui[Z]]]
            else
                SA[U[index.ρi[Z], i+1], U[index.ρiui[Z], i+1]]
            end

            fluid = fluids[Z+1]

            FL = scheme.flux_function(UL, U0, fluid, coupled, ϵL, ϵ0, neL, ne0)
            FR = scheme.flux_function(U0, UR, fluid, coupled, ϵ0, ϵR, ne0, neR)

            F_mass, F_momentum = FR[1] - FL[1], FR[2] - FL[2]

            Q_accel = -Z * e * U[index.ρi[Z], i] * ue[i] / μ[i] / mi

            dU[index.ρi[Z], i] = -F_mass / Δz
            dU[index.ρiui[Z], i] = -F_momentum / Δz + Q_accel
        end

        # ionization reactions
        for r in reactions
            reactant_index = species_range_dict[r.reactant.symbol][1]
            product_index = species_range_dict[r.product.symbol][1]
            ρ_reactant = U[reactant_index, i]
            k = r.rate_coeff
            ρdot = k(ϵ0) * ρ_reactant * ne0
            dU[reactant_index, i] -= ρdot
            dU[product_index, i] += ρdot
        end
    end
end

function update!(dU, U, p, t)
    update_heavy_species!(dU, U, p, t)
    update_electron_energy_explicit!(dU, U, p, t)
end

left_edge(i) = i - 1
right_edge(i) = i

function electron_density(U, params)
    ne = 0.0
    index = params.index
    @inbounds for i in 1:params.config.ncharge
        ne += i * U[index.ρi[i]]
    end
    return ne
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
    #ϕ_hallis, grad_ϕ_hallis = load_hallis_for_input()

    BCs = sim.boundary_conditions

    precompute_bfield!(cache.B, grid.cell_centers)

    loss_coeff = config.ionization_coeffs == :BOLSIG_FIT ? loss_coeff_fit : landmark.loss_coeff

    params = (;
        config = config,
        propellant = config.propellant,
        ϕ_L = config.anode_potential,
        ϕ_R = config.cathode_potential,
        Te_L = config.anode_Te,
        Te_R = config.cathode_Te,
        L_ch = config.geometry.channel_length,
        A_ch = channel_area(config.geometry.outer_radius, config.geometry.inner_radius),
        νϵ = config.radial_loss_coefficients,
        νw = config.wall_collision_frequencies,
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
        scheme, BCs, dt=timestep,
    )

    U = deepcopy(U)
    params = deepcopy(params)

    #PREPROCESS
    #make values in params available for first timestep
    update_values!(U, params)

    function save_func(u, t, integrator)
        (; μ, Tev, ϕ, ∇ϕ, ne, pe, ue) = integrator.p.cache
        return deepcopy((; μ, Tev, ϕ, ∇ϕ, ne, pe, ue))
    end

    saved_values = SavedValues(Float64, NamedTuple{(:μ, :Tev, :ϕ, :∇ϕ, :ne, :pe, :ue), NTuple{7, Vector{Float64}}})

    discrete_callback = DiscreteCallback(Returns(true), update_values!, save_positions=(false,false))
    saving_callback = SavingCallback(save_func, saved_values, saveat = sim.saveat)

    if sim.callback !== nothing
        cb = CallbackSet(discrete_callback, saving_callback, sim.callback)
    else
        cb = CallbackSet(discrete_callback, saving_callback)
    end

    implicit_energy = config.implicit_energy

    maxiters = Int(ceil(1000 * tspan[2] / timestep))

    if implicit_energy > 0
        f = ODEFunction(update_heavy_species!)
    else
	    f = ODEFunction(update!)
    end

    prob = ODEProblem{true}(f, U, tspan, params)
	sol = solve(
		prob, alg; saveat=sim.saveat, callback=cb,
		adaptive=adaptive, dt=timestep, dtmax=10*timestep, dtmin = timestep/10, maxiters = maxiters,
	)

    if sol.retcode == :NaNDetected
        println("Simulation failed with NaN detected at t = $(sol.t[end])")
    elseif sol.retcode == :InfDetected
        println("Simulation failed with Inf detected at t = $(sol.t[end])")
    end

    return HallThrusterSolution(sol, params, saved_values.saveval)
end

function inlet_neutral_density(sim)
    un = sim.neutral_velocity
    A = channel_area(sim.geometry)
    m_atom = sim.propellant.m
    nn = sim.inlet_mdot / un / A / m_atom
    return nn
end

function initial_condition!(U, z_cell, IC!, fluids)
    for (i, z) in enumerate(z_cell)
        @views IC!(U[:, i], z, fluids, z_cell[end])
    end
end

############################################################################################
using DelimitedFiles

function load_hallis_output(output_path)
    output_headers = [
        :z, :ne, :ϕ, :Te, :Ez, :Br, :nn, :ndot, :μe, :μen, :μbohm, :μwall, :μei,
    ]
    output = DataFrame(readdlm(output_path, Float64), output_headers)
    output.ωce = output.Br * 1.6e-19 / 9.1e-31
    replace!(output.nn, 0.0 => 1e12)
    return output[1:end-1, :]
end


function load_hallis_for_input()
    hallis = load_hallis_output("landmark/Av_PLOT_HALLIS_1D_01.out")
    ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, hallis.ϕ)
    grad_ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, -hallis.Ez)
    return ϕ_hallis, grad_ϕ_hallis
end

