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

    #Dual = ForwardDiff.Dual

    U = zeros(nvariables + 7, ncells + 2) # need to allocate room for ghost cells
    F = zeros(nvariables + 1, nedges)
    UL = zeros(nvariables + 1, nedges)
    UR = zeros(nvariables + 1, nedges)
    Q = zeros(nvariables + 1)
    A = Tridiagonal(ones(ncells), ones(ncells+1), ones(ncells)) #for potential
    b = zeros(ncells + 1) #for potential equation
    B = zeros(ncells + 2)
    νan = zeros(ncells + 2)
    νc = zeros(ncells + 2)
    μ = zeros(ncells + 2)

    cache = (; F, UL, UR, Q, A, b, B, νan, νc, μ)
    return U, cache
end

function update_heavy_species!(dU, U, params, t) #get source and BCs for potential from params
    ####################################################################
    #extract some useful stuff from params

    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    index = params.index

    F, UL, UR, Q = params.cache.F, params.cache.UL, params.cache.UR, params.cache.Q
    B = params.cache.B

    z_cell, z_edge, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    scheme = params.scheme
    source_term! = params.source_term!

    ncells = size(U, 2) - 2

    mi = m(fluids[1])

    ##############################################################
    #FLUID MODULE

    #fluid BCs
    @views apply_bc!(U[1:index.lf, :], params.BCs[1], :left, params.Te_L, mi)
    @views apply_bc!(U[1:index.lf, :], params.BCs[2], :right, params.Te_R, mi)

    #fluid computations, electron in implicit
    @views compute_edge_states!(UL[1:index.lf, :], UR[1:index.lf, :], U[1:index.lf, :], scheme)
    @views compute_fluxes!(
        F[1:index.lf, :],
        UL[1:index.lf, :],
        UR[1:index.lf, :],
        fluids,
        fluid_ranges,
        scheme,
        U[index.pe, :],
        params.config.electron_pressure_coupled
    )

    # Compute heavy species source terms
    @inbounds  for i in 2:(ncells + 1)
        @turbo Q .= 0.0

        #fluid source term (includes ionization and acceleration and energy)
        source_term!(Q, U, params, i)

        #Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        # Neutral fluxes and source
        dU[index.ρn, i] = (F[index.ρn, left] - F[index.ρn, right]) / Δz + Q[index.ρn]
        η = params.config.ion_diffusion_coeff * sqrt(2*e*U[index.Tev, i]/(3*mi))

        @turbo for j in index.ρi[1]:index.lf
            # Ion fluxes and source
            dU[j, i] = (F[j, left] - F[j, right]) / Δz + Q[j]
            # Add optional diffusion term for the ions in interior cells
            dU[j, i] += η*(U[j, i-1] - 2U[j, i] + U[j, i+1])/(Δz)^2
        end

        # Electron source term
        dU[index.nϵ, i] = Q[index.nϵ]
    end

    if !params.implicit_energy
        update_electron_energy!(dU, U, params, t)
    end
end

function update_electron_energy!(dU, U, params, t)
    #########################################################
    #ELECTRON SOLVE

    index = params.index
    μ = params.cache.μ
    z_cell = params.z_cell
    ncells = size(U, 2) - 2
    F = params.cache.F
    z_edge = params.z_edge

    #apply_bc_electron!(U, params.BCs[3], :left, index)
    #apply_bc_electron!(U, params.BCs[4], :right, index)

    # Dirchlet BCs for electron energy
    U[index.nϵ, 1] = params.Te_L * U[index.ne, 1]
    U[index.nϵ, end] = params.Te_R * U[index.ne, end]

    @inbounds for i in 2:ncells+1

        left = left_edge(i)
        right = right_edge(i)

        # Upwinded first derivatives
        ue = U[index.ue, i]
        ne = U[index.ne, i]

        Δz = z_edge[right] - z_edge[left] #boundary cells same size with upwind

        ue⁺ = max(ue, 0.0) / ue
        ue⁻ = min(ue, 0.0) / ue
        #ue⁺ = smooth_if(ue, 0.0, 0.0, ue, 0.01) / ue
        #ue⁻ = smooth_if(ue, 0.0, ue, 0.0, 0.01) / ue


        advection_term = (
            ue⁺ * (ue * U[index.nϵ, i] - U[index.ue, i-1] * U[index.nϵ, i-1]) -
            ue⁻ * (ue * U[index.nϵ, i] - U[index.ue, i+1] * U[index.nϵ, i+1])
        ) / Δz

        diffusion_term_1 = -(
            ue⁺ * (-(μ[i-1] * U[index.nϵ, i-1] - μ[i] * U[index.nϵ, i]) * (U[index.Tev, i-1] - U[index.Tev, i])) -
            ue⁻ * (μ[i+1] * U[index.nϵ, i+1] - μ[i] * U[index.nϵ, i]) * (U[index.Tev, i+1] - U[index.Tev, i])
        ) / Δz^2

        if i == 2 #2
            diffusion_term_2 = μ[i] * U[index.nϵ, i] * (8U[index.Tev, i-1] - 12U[index.Tev, i] + 4U[index.Tev, i+1])
            diffusion_term_2 /= 3*Δz^2
        elseif i == ncells + 1 #ncells + 1
            diffusion_term_2 = μ[i] * U[index.nϵ, i] * (4U[index.Tev, i-1] - 12U[index.Tev, i] + 8U[index.Tev, i+1])
            diffusion_term_2 /= 3*Δz^2
        else
            # Central second derivatives
            diffusion_term_2 = μ[i] * U[index.nϵ, i] * (U[index.Tev, i-1] - 2U[index.Tev, i] + U[index.Tev, i+1])
            diffusion_term_2 /= Δz^2
        end

        dU[index.nϵ, i] += -5/3 * advection_term + 10/9 * (diffusion_term_1 + diffusion_term_2)
    end

    return nothing
end

left_edge(i) = i - 1
right_edge(i) = i

function electron_density(U, fluid_ranges)
    ne = 0.0
    @inbounds for (i, f) in enumerate(fluid_ranges)
        if i == 1
            continue # neutrals do not contribute to electron density
        end
        charge_state = i - 1
        ne += charge_state * U[f[1]]
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

function update_values!(U, params)

    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    index = params.index

    B = params.cache.B

    z_cell, z_edge, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    ncells = size(U, 2) - 2

    #update useful quantities relevant for potential, electron energy and fluid solve
    L_ch = params.L_ch
    mi = m(fluids[1])

    OVS = params.OVS.energy.active

    @inbounds @views for i in 1:(ncells + 2)
        z = z_cell[i]
        OVS_ne = OVS * (params.OVS.energy.ne(z))
        OVS_Tev = OVS * (params.OVS.energy.Tev(z))

        U[index.ne, i] = (1 - OVS) * max(1e13, electron_density(U[:, i], fluid_ranges) / mi) + OVS_ne
        U[index.Tev, i] = (1 - OVS) * max(0.1, U[index.nϵ, i]/U[index.ne, i]) + OVS_Tev
        U[index.pe, i] = U[index.nϵ, i]
        params.cache.νan[i] = params.anom_model(i, U, params)
        params.cache.νc[i] = electron_collision_freq(U[index.Tev, i], U[1, i]/mi , U[index.ne, i], mi)
        params.cache.μ[i] = (1 - params.OVS.energy.active)*electron_mobility(params.cache.νan[i], params.cache.νc[i], B[i]) #+ OVS*(params.OVS.energy.μ)
    end

    # update electrostatic potential
    solve_potential_edge!(U, params)

    @inbounds @views for i in 1:(ncells + 2)
        @views U[index.grad_ϕ, i] = first_deriv_central_diff_pot(U[index.ϕ, :], params.z_cell, i)
        U[index.ue, i] = (1 - params.OVS.energy.active)*electron_velocity(U, params, i) + params.OVS.energy.active*(params.OVS.energy.ue)
    end

    U[index.ue, 1] = U[index.ue, 2]
    U[index.ue, end] = U[index.ue, end-1]



end

#=Config = @NamedTuple begin
    propellant::Gas
end=#

condition(u,t,integrator) = t < 1
function affect!(integrator)
    update_values!(integrator.u, integrator.p)
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

    if sim.callback !== nothing
        cb = CallbackSet(DiscreteCallback(condition, affect!, save_positions=(false,false)), sim.callback)
    else
        cb = DiscreteCallback(condition, affect!, save_positions=(false,false))
    end

    implicit_energy = sim.solve_energy

    maxiters = Int(ceil(1000 * tspan[2] / timestep))

    if implicit_energy
        splitprob = SplitODEProblem{true}(update_electron_energy!, update_heavy_species!, U, tspan, params)
        sol = solve(splitprob, KenCarp47(); saveat=sim.saveat, callback=cb,
        adaptive=adaptive, dt=timestep, dtmax = timestep)
    else
        #alg = SSPRK22()
        alg = AutoTsit5(Rosenbrock23())
        prob = ODEProblem{true}(update_heavy_species!, U, tspan, params)
        sol = solve(prob, alg; saveat=sim.saveat, callback=cb,
        adaptive=adaptive, dt=timestep, dtmax=timestep, maxiters = maxiters)
    end

    return HallThrusterSolution(sol, params)
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

