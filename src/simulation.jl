struct HyperbolicScheme{F,L}
    flux_function::F  # in-place flux function
    limiter::L # limiter
    reconstruct::Bool
end

Base.@kwdef mutable struct MultiFluidSimulation{IC,B1,B2,B3,B4,S,F,L,CB,SP,BP,VF} #could add callback, or autoselect callback when in MMS mode
    grid::Grid1D
    fluids::Vector{Fluid}     # An array of user-defined fluids.
    # This will give us the capacity to more easily do shock tubes (and other problems)
    # without Hall thruster baggage
    initial_condition::IC
    boundary_conditions::Tuple{B1,B2,B3,B4}   # Tuple of left and right boundary conditions, subject to the approval of PR #10
    end_time::Float64    # How long to simulate
    scheme::HyperbolicScheme{F,L} # Flux, Limiter
    source_term!::S  # Source term function. This can include reactons, electric field, and MMS terms
    source_potential!::SP #potential source term
    boundary_potential!::BP #boundary conditions potential
    saveat::Vector{Float64} #when to save
    timestepcontrol::Tuple{Float64,Bool} #sets timestep (first argument) if second argument false. if second argument (adaptive) true, given dt is ignored.
    callback::CB
    solve_energy::Bool
    verification::VF
end

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
    F = zeros(nvariables + 1, nedges)
    UL = zeros(nvariables + 1, nedges)
    UR = zeros(nvariables + 1, nedges)
    Q = zeros(nvariables + 1, ncells+2)
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

    cache = (; F, UL, UR, Q, A, b, Aϵ, bϵ, B, νan, νc, μ, ϕ, ∇ϕ, ne, Tev, pe, ue)
    return U, cache
end

function update_heavy_species!(dU, U, params, t)
    ####################################################################
    #extract some useful stuff from params

    (;index, z_edge, propellant, fluid_ranges, fluids, species_range_dict, reactions) = params
    (;Q, ue, μ) = params.cache

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

    # Compute heavy species source terms
    @inbounds for i in 2:(ncells + 1)

        #Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]
        # Neutral fluxes and source
        ## OLD:
        #dU[index.ρn, i] = (F[index.ρn, left] - F[index.ρn, right]) / Δz + Q[index.ρn, i]
        ## NEW:
        fL = ρn[i-1] * un
        fR = ρn[i+1] * un
        a = sqrt(γn*Rn*Tn)
        sn = max(un + a, un - a)

        Fn = (fR - fL - sn * (ρn[i-1] - 2ρn[i] + ρn[i+1])) / 2 / Δz
        ## Currently don't worry about source term, that comes next
        #Qn = -ne * ρn[i] * sum(rxn.rate_coeff(ϵ) for rxn in reactions if rxn.reactant.Z == 0)
        dU[index.ρn, i] = -Fn# + Q[index.ρn, i]

        ## Remove diffusion term for now
        #η = δ * sqrt(2*e*Tev[i]/(3*mi))

        ## OLD
        #=for j in first_ion_index:last_ion_index
            # Ion fluxes and source
            dU[j, i] = (F[j, left] - F[j, right]) / Δz + Q[j, i]
            # Compute current charge state to make sure sound speed in diffusion term is correct
            Z = ((j - first_ion_index) ÷ num_ion_equations) + 1
            # Add optional diffusion term for the ions in interior cells
            dU[j, i] += sqrt(Z) * η*(U[j, i-1] - 2U[j, i] + U[j, i+1])/Δz²
        end=#

        ## NEW
        neL = 0.0
        neR = 0.0
        ne0 = 0.0
        for Z in 1:ncharge
            neL += Z * U[index.ρi[Z], i-1]
            neR += Z * U[index.ρi[Z], i+1]
            ne0 += Z * U[index.ρi[Z], i]
        end

        neL, neR, ne0 = neL/mi, neR/mi, ne0/mi
        ϵL = max(0.0, nϵ[i-1]/neL)
        ϵ0 = max(0.0, nϵ[i]/ne0)
        ϵR = max(0.0, nϵ[i+1]/neR)

        for Z in 1:ncharge

            aL = sqrt((γn * kB * Ti + Z * e * ϵL)/mi)
            aR = sqrt((γn * kB * Ti + Z * e * ϵR)/mi)

            uL = U[index.ρiui[Z], i-1] / U[index.ρi[Z], i-1]
            uR = U[index.ρiui[Z], i+1] / U[index.ρi[Z], i+1]

            sL = max(abs(uL + aL), abs(uL - aL))
            sR = max(abs(uR + aR), abs(uR - aR))

            f_mass_L = U[index.ρiui[Z], i-1]
            f_mass_R = U[index.ρiui[Z], i+1]

            F_mass = (
                f_mass_R - f_mass_L -
                sL * U[index.ρi[Z], i-1] +
                (sL + sR) * U[index.ρi[Z], i] -
                sR * U[index.ρi[Z], i+1]
            ) / 2

            f_momentum_L = U[index.ρiui[Z], i-1]^2 / U[index.ρi[Z], i-1] + Z * e * nϵ[i-1] + U[index.ρi[Z], i-1] * Rn * Ti
            f_momentum_R = U[index.ρiui[Z], i+1]^2 / U[index.ρi[Z], i+1] + Z * e * nϵ[i+1] + U[index.ρi[Z], i+1] * Rn * Ti

            F_momentum = (
                f_momentum_R - f_momentum_L -
                sL * U[index.ρiui[Z], i-1] +
                (sL + sR) * U[index.ρiui[Z], i] -
                sR * U[index.ρiui[Z], i+1]
            ) / 2

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

        # Electron source term
        dU[index.nϵ, i] = source_electron_energy_landmark(U, params, i)
    end
end

function update_electron_energy!(dU, U, params, t)
    #########################################################
    #ELECTRON SOLVE

    ncells = size(U, 2) - 2

    (;index, z_cell, z_edge) = params
    (;μ, Tev, ue, ne) = params.cache
    nϵ = @views U[index.nϵ, :]
    implicit_energy = params.config.implicit_energy

    mi = params.propellant.m

    @inbounds for i in 2:ncells+1

        # Upwinded first derivatives
        zL = z_cell[i-1]
        z0 = z_cell[i]
        zR = z_cell[i+1]

        neL = sum(U[index.ρi[Z], i-1] for Z in 1:params.config.ncharge) / mi
        ne0 = sum(U[index.ρi[Z], i] for Z in 1:params.config.ncharge) / mi
        neR = sum(U[index.ρi[Z], i+1] for Z in 1:params.config.ncharge) / mi

        ϵL = nϵ[i-1] / neL
        ϵ0 = nϵ[i] / ne0
        ϵR = nϵ[i+1] / neR

        dϵ_dz⁺ = diff(ϵL, ϵ0, zL, z0)
        dϵ_dz⁻ = diff(ϵ0, ϵR, zL, z0)
        d²ϵ_dz² = uneven_second_deriv(ϵL, ϵ0, ϵR, zL, z0, zR)

        #ue⁺ = max(ue[i], 0.0) / ue[i]
        #ue⁻ = min(ue[i], 0.0) / ue[i]
        ue⁺ = smooth_if(ue[i], 0.0, 0.0, ue[i], 0.01) / ue[i]
        ue⁻ = smooth_if(ue[i], 0.0, ue[i], 0.0, 0.01) / ue[i]

        advection_term = (
            ue⁺ * diff(ue[i-1]*nϵ[i-1], ue[i]*nϵ[i], zL, z0) +
            ue⁻ * diff(ue[i]*nϵ[i], ue[i+1]*nϵ[i+1], z0, zR)
        )

        diffusion_term = (
            ue⁺ * diff(μ[i-1]*nϵ[i-1], μ[i]*nϵ[i], zL, z0) * dϵ_dz⁺ +
            ue⁻ * diff(μ[i]*nϵ[i], μ[i+1]*nϵ[i+1], z0, zR) * dϵ_dz⁻
        )

        diffusion_term += μ[i] * nϵ[i] * d²ϵ_dz²

        dU[index.nϵ, i] += - 5/3 * advection_term + 10/9 * diffusion_term
    end

    return nothing
end

function update!(dU, U, p, t)
    update_heavy_species!(dU, U, p, t)
    update_electron_energy!(dU, U, p, t)
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

update_values!(integrator) = update_values!(integrator.u, integrator.p)
function update_values!(U, params, CN = true)
    #update useful quantities relevant for potential, electron energy and fluid solve

    ncells = size(U, 2) - 2

    (;z_cell, fluids, fluid_ranges, index, scheme, source_term!, z_edge) = params
    (;F, UL, UR, Q, B, ue, Tev, ∇ϕ, ϕ, pe, ne) = params.cache
    OVS = params.OVS.energy.active

    mi = params.propellant.m

    # Edge state reconstruction and flux computation
    #@views compute_fluxes!(F, UL, UR, U, params)
    #@views compute_edge_states!(UL[1:index.lf, :], UR[1:index.lf, :], U[1:index.lf, :], scheme)
    coupled = params.config.electron_pressure_coupled
    #@views compute_fluxes!(F[1:index.lf, :], UL[1:index.lf, :], UR[1:index.lf, :], fluids, fluid_ranges, scheme, pe, coupled)

    # Apply boundary conditions
    @views apply_bc!(U[1:index.lf, :], params.BCs[1], :left, params.Te_L, mi)
    @views apply_bc!(U[1:index.lf, :], params.BCs[2], :right, params.Te_R, mi)

    # Update electron quantities
    @inbounds @views for i in 1:(ncells + 2)
        z = z_cell[i]
        OVS_ne = OVS * (params.OVS.energy.ne(z))
        OVS_Tev = OVS * (params.OVS.energy.Tev(z))

        @views ne[i] = (1 - OVS) * max(1e13, electron_density(U[:, i], params) / mi) + OVS_ne
        Tev[i] = (1 - OVS) * max(0.1, U[index.nϵ, i]/ne[i]) + OVS_Tev
        pe[i] = U[index.nϵ, i]
        params.cache.νan[i] = params.anom_model(U, params, i)
        params.cache.νc[i] = electron_collision_freq(params.cache.Tev[i], U[1, i]/mi , ne[i], mi)
        params.cache.μ[i] = (1 - params.OVS.energy.active)*electron_mobility(params.cache.νan[i], params.cache.νc[i], B[i]) #+ OVS*(params.OVS.energy.μ)
    end

    # update electrostatic potential and potential gradient on edges
    solve_potential_edge!(U, params)
    #U[index.ϕ, :] .= params.OVS.energy.active.*0.0 #avoiding abort during OVS
    #∇ϕ[1] = first_deriv_central_diff_pot(ϕ, params.z_cell, 1)
    #∇ϕ[end] = first_deriv_central_diff_pot(ϕ, params.z_cell, ncells+2)
    ∇ϕ[1] = uneven_forward_diff(ϕ[1], ϕ[2], ϕ[3], z_edge[1], z_edge[2], z_edge[3])
    ∇ϕ[end] = uneven_backward_diff(ϕ[end-2], ϕ[end-1], ϕ[end], z_edge[end-2], z_edge[end-1], z_edge[end])

    # Compute interior potential gradient and electron velocity and update source terms
    @inbounds for i in 2:(ncells + 1)
        # potential gradient
        ∇ϕ[i] = first_deriv_central_diff_pot(ϕ, params.z_cell, i)
        #∇ϕ[i] = uneven_central_diff(ϕ[i-1], ϕ[i], ϕ[i+1], z_cell[i-1], z_cell[i], z_cell[i+1])
        ue[i] = (1 - OVS) * electron_velocity(U, params, i) + OVS * (params.OVS.energy.ue)

        #source term (includes ionization and acceleration as well as energy source temrs)
        #@views source_term!(Q[:, i], U, params, i)
    end

    # Neumann condition for electron velocity
    ue[1] = ue[2]
    ue[end] = ue[end-1]

    # Dirchlet BCs for electron energy
    apply_bc_electron!(U, params.BCs[3], :left, index)
    apply_bc_electron!(U, params.BCs[4], :right, index)

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

function run_simulation(sim, config) #put source and Bcs potential in params

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
    else
        throw(ArgumentError("Invalid ionization reactions selected. Please choose either :LANDMARK or :BOLSIG"))
    end
    #ϕ_hallis, grad_ϕ_hallis = load_hallis_for_input()

    BCs = sim.boundary_conditions

    precompute_bfield!(cache.B, grid.cell_centers)

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
        loss_coeff = landmark.loss_coeff,
        reactions = ionization_reactions,
        implicit_energy = config.implicit_energy,
        index, cache, fluids, fluid_ranges, species_range_dict, z_cell=grid.cell_centers,
        z_edge=grid.edges, cell_volume=grid.cell_volume, source_term!,
        scheme, BCs, dt=timestep,
    )

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

    #=if implicit_energy > 0
        #alg = AutoTsit5(Rosenbrock23())
        alg = SSPRK22()
        prob = ODEProblem{true}(update_heavy_species!, U, tspan, params)
        sol = solve(prob, alg; saveat=sim.saveat, callback=cb,
        adaptive=false, dt=timestep, dtmax=timestep, maxiters = maxiters)
        #=splitprob = SplitODEProblem{true}(update_electron_energy!, update_heavy_species!, U, tspan, params)
        sol = solve(splitprob, KenCarp47(); saveat=sim.saveat, callback=cb,
        adaptive=adaptive, dtmin = 1e-12, dt=timestep, dtmax = timestep, maxiters = maxiters)=#
    else
        #alg = SSPRK22()
        alg = AutoTsit5(Rosenbrock23())
        prob = ODEProblem{true}(update!, U, tspan, params)
        sol = solve(prob, alg; saveat=sim.saveat, callback=cb,
        adaptive=adaptive, dt=timestep, dtmax=timestep, maxiters = maxiters)
    end=#
	alg = AutoTsit5(Rosenbrock23())
    #alg = Rosenbrock23()
	f = ODEFunction(update!)
    dU = copy(U)
    j_func = (dU, U) -> f(dU, U, params, 0.0)
    J = ForwardDiff.jacobian(j_func, dU, U) |> sparse
    #jac_sparsity = Symbolics.jacobian_sparsity((dU, u) -> f(dU0, U, params, 0.0), dU0, U)
    prob = ODEProblem{true}(f, U, tspan, params, jac_prototype=J)
	#modelingtoolkitize(prob)
	sol = solve(
		prob, alg; saveat=sim.saveat, callback=cb,
		adaptive=adaptive, dt=timestep, dtmax=10*timestep, dtmin = timestep/10, maxiters = maxiters,

	)

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

