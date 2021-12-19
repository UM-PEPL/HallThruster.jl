struct HyperbolicScheme{F,L}
    flux_function::F  # in-place flux function
    limiter::L # limiter
    reconstruct::Bool
end

Base.@kwdef mutable struct MultiFluidSimulation{IC,IC_E,B1,B2,S,F,L,CB,SP,BP} #could add callback, or autoselect callback when in MMS mode
    grid::Grid1D
    fluids::Vector{Fluid}     # An array of user-defined fluids.
    # This will give us the capacity to more easily do shock tubes (and other problems)
    # without Hall thruster baggage
    initial_condition::IC
    initial_condition_E::IC_E
    boundary_conditions::Tuple{B1,B2}   # Tuple of left and right boundary conditions, subject to the approval of PR #10
    end_time::Float64    # How long to simulate
    scheme::HyperbolicScheme{F,L} # Flux, Limiter
    source_term!::S  # Source term function. This can include reactons, electric field, and MMS terms
    source_potential!::SP #potential source term
    boundary_potential!::BP #boundary conditions potential
    saveat::Vector{Float64} #when to save
    timestepcontrol::Tuple{Float64,Bool} #sets timestep (first argument) if second argument false. if second argument (adaptive) true, given dt is ignored.
    callback::CB
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

function allocate_arrays(sim) #rewrite allocate arrays as function of set of equations, either 1, 2 or 3
    # Number of variables in the state vector U
    nvariables = 1 #for electron energy
    for i in 1:length(sim.fluids)
        if sim.fluids[i].conservation_laws.type == :ContinuityOnly
            nvariables += 1
        elseif sim.fluids[i].conservation_laws.type == :IsothermalEuler
            nvariables += 2
        elseif sim.fluids[i].conservation_laws.type == :EulerEquations
            nvariables += 3
        end
    end

    ncells = sim.grid.ncells
    nedges = sim.grid.ncells + 1

    #Dual = ForwardDiff.Dual

    U = zeros(nvariables, ncells + 2) # need to allocate room for ghost cells
    F = zeros(nvariables, nedges)
    UL = zeros(nvariables, nedges)
    UR = zeros(nvariables, nedges)
    Q = zeros(nvariables)
    A = Tridiagonal(ones(ncells - 1), ones(ncells), ones(ncells - 1)) #for potential
    b = zeros(ncells) #for potential equation
    ϕ = zeros(ncells + 2) #for potential equation, need to add ghost cells to potential as well
    pe = zeros(ncells + 2)
    B = zeros(ncells + 2)
    ne = zeros(ncells + 2)
    νan = zeros(ncells + 2)
    νc = zeros(ncells + 2)
    μ = zeros(ncells + 2)

    L_ch = 0.025
    Tev = map(x -> Te_func(x, L_ch), sim.grid.cell_centers)

    cache = (; F, UL, UR, Q, A, b, ϕ, Tev, pe, ne, B, νan, νc, μ)
    return U, cache
end

function update_exp!(dU, U, params, t) #get source and BCs for potential from params
    ####################################################################
    #extract some useful stuff from params
    fluids, fluid_ranges = params.fluids, params.fluid_ranges

    F, UL, UR, Q = params.cache.F, params.cache.UL, params.cache.UR, params.cache.Q
    ϕ, Tev, B = params.cache.ϕ, params.cache.Tev, params.cache.B

    z_cell, z_edge, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    scheme = params.scheme
    source_term! = params.source_term!

    ncells = size(U, 2) - 2

    ####################################################################
    #PREPROCESS
    #calculate useful quantities relevant for potential, electron energy and fluid solve
    L_ch = 0.025
    fluid = fluids[1].species.element
    @inbounds for i in 1:(ncells + 2)
        #update electron temperature from energy using old density
        params.cache.Tev[i] = max(0, U[4, i]/3*2/params.cache.ne[i])
        params.cache.ne[i] = max(0, electron_density(@view(U[:, i]), fluid_ranges) / fluid.m)
        params.cache.pe[i] = electron_pressure(params.cache.ne[i], Tev[i])
        params.cache.νan[i] = get_v_an(z_cell[i], B[i], L_ch)
        params.cache.νc[i] = get_v_c(Tev[i], U[1, i]/fluid.m , params.cache.ne[i], fluid.m)
        params.cache.μ[i] = cf_electron_transport(params.cache.νan[i], params.cache.νc[i], B[i])
    end

    #=
    println("ne: ", params.cache.ne)
    println("pe: ", params.cache.pe)
    println("νan: ", params.cache.νan)
    println("νc: ", params.cache.νc)
    println("Tev K: ", params.cache.Tev .*e ./kB)
    =#

    #update internal electron energy according to new electron density
    U[4, :] .= 3/2 .* params.cache.ne .* params.cache.Tev

    ####################################################################
    #POTENTIAL MODULE
    #println("ϕ before solve: ", params.cache.ϕ)

    solve_potential!(ϕ, U, params)

    #println("ϕ after solve: ", params.cache.ϕ)


    ##############################################################
    #FLUID MODULE

    #fluid BCs
    apply_bc!(@views(U[1:3, :]), params.BCs[1], :left)
    apply_bc!(@views(U[1:3, :]), params.BCs[2], :right)

    #electron BCs
    Tev_anode = 3 #eV 
    Tev_cathode = 3 #eV
    left_state = [3/2*params.cache.ne[1]*Tev_anode]
    right_state = [3/2*params.cache.ne[end]*Tev_cathode]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    U[4, begin] = left_state[1]
    U[4, end] = right_state[1]
    #apply_bc!(@views(U[4, :]), BCs[1], :left)
    #apply_bc!(@views(U[4, :]), BCs[2], :right)

    #fluid computations, electron in implicit
    compute_edge_states!(@views(UL[1:3, :]), @views(UR[1:3, :]), @views(U[1:3, :]), scheme)
    compute_fluxes!(@views(F[1:3, :]), @views(UL[1:3, :]), @views(UR[1:3, :]), fluids, fluid_ranges, scheme)

    # Compute heavy species source terms
    @inbounds for i in 2:(ncells + 1)
        @turbo Q .= 0.0

        #fluid source term
        source_term!(@views(Q[1:3]), @views(U[1:3, :]), params, ϕ, Tev, i)

        #electron source term
        Q[4] = source_electron_energy!(@views(Q[4]), @views(U), params, i)

        #@show Q[4] 

        # Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        @tturbo @views @. dU[:, i] = (F[:, left] - F[:, right]) / Δz + Q #should be fine, F[4, :] should be 0
    end

    return nothing
end

function update_imp!(dU, U, params, t)
    F, UL, UR, Q = params.cache.F, params.cache.UL, params.cache.UR, params.cache.Q

    #electron BCs
    Tev_anode = 3 #eV 
    Tev_cathode = 3 #eV
    left_state = [3/2*params.cache.ne[1]*Tev_anode]
    right_state = [3/2*params.cache.ne[end]*Tev_cathode]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    U[4, begin] = left_state[1]
    U[4, end] = right_state[1]
    #println("U after BCs applied: ", U[4, :])

    
    #electron computations, fluid in explicit
    scheme = HallThruster.HyperbolicScheme(HallThruster.upwind_electron!, identity, false)
    compute_edge_states!(@views(UL[4, :]), @views(UR[4, :]), @views(U[4, :]), scheme)
    #println("UL and UR after comp edge: ", UL[4, :], UR[4, :])
    compute_fluxes_electron!(@views(F[4, :]), @views(UL[4, :]), @views(UR[4, :]), [HallThruster.Electron], [1:1], scheme, params)
    #println("Flux after compute fluxes: ", F[4, :])
    

    ncells = size(U, 2) - 2
    z_edge = params.z_edge

    @inbounds for i in 2:(ncells + 1)
        # Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        dU[4, i] = (F[4, left] - F[4, right]) / Δz
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

function run_simulation(sim) #put source and Bcs potential in params
    species, fluids, fluid_ranges, species_range_dict = configure_simulation(sim)
    grid = sim.grid

    U, cache = allocate_arrays(sim)

    println("E before initial cond applied", U[4, :])

    initial_condition!(@views(U[1:3, :]), @views(U[4, :]), grid.cell_centers, sim.initial_condition,
    sim.initial_condition_E, fluid_ranges, fluids)

    println("E after initial cond applied", U[4, :])

    scheme = sim.scheme
    source_term! = sim.source_term!
    timestep = sim.timestepcontrol[1]
    adaptive = sim.timestepcontrol[2]
    tspan = (0.0, sim.end_time)

    reactions = load_ionization_reactions(species)
    landmark = load_landmark()

    BCs = sim.boundary_conditions

    precompute_bfield!(cache.B, grid.cell_centers)

    params = (; cache, fluids, fluid_ranges, species_range_dict, z_cell=grid.cell_centers,
              z_edge=grid.edges, cell_volume=grid.cell_volume, source_term!, reactions,
              scheme, BCs, dt=timestep, source_potential! = sim.source_potential!, 
              boundary_potential! = sim.boundary_potential!, landmark)

    #PREPROCESS
    #make values in params available for first implicit timestep
    ncells = size(U, 2) - 2
    L_ch = 0.025
    fluid = fluids[1].species.element
    @inbounds for i in 1:(ncells + 2)
        params.cache.ne[i] = electron_density(@view(U[:, i]), fluid_ranges) / fluid.m
        params.cache.pe[i] = electron_pressure(params.cache.ne[i], params.cache.Tev[i])
        params.cache.νan[i] = get_v_an(params.z_cell[i], params.cache.B[i], L_ch)
        params.cache.νc[i] = get_v_c(params.cache.Tev[i], U[1, i]/fluid.m , params.cache.ne[i], fluid.m)
        params.cache.μ[i] = cf_electron_transport(params.cache.νan[i], params.cache.νc[i], params.cache.B[i])
    end

    prob = SplitODEProblem{true}(update_imp!, update_exp!, U, tspan, params)
    #tmp_prob = remake(prob, u0=convert.(eltype(params),prob.u0), p=params)
    sol = solve(prob, KenCarp3(autodiff = false); saveat=sim.saveat, callback=sim.callback,
    adaptive=adaptive, dt=timestep)
    #=sol = solve(prob, KenCarp3(); saveat=sim.saveat, callback=sim.callback,
                adaptive=adaptive, dt=timestep)=#
    return sol
end

function inlet_neutral_density(sim)
    un = sim.neutral_velocity
    A = channel_area(sim.geometry)
    m_atom = sim.propellant.m
    nn = sim.inlet_mdot / un / A / m_atom
    return nn
end

function initial_condition!(U, E, z_cell, IC!, IC_E!, fluid_ranges, fluids)
    for (i, z) in enumerate(z_cell)
        @views IC!(U[:, i], z, fluids, z_cell[end])
        IC_E!(E, U[:, i], z, z_cell[end], fluid_ranges, fluids, i)
    end
end