using Test

struct HyperbolicScheme{F,L}
    flux_function::F  # in-place flux function
    limiter::L # limiter
    reconstruct::Bool
end

Base.@kwdef mutable struct MultiFluidSimulation{IC,B1,B2,S,F,L,CB,SP,BP} #could add callback, or autoselect callback when in MMS mode
    grid::Grid1D
    fluids::Vector{Fluid}     # An array of user-defined fluids.
    # This will give us the capacity to more easily do shock tubes (and other problems)
    # without Hall thruster baggage
    initial_condition::IC
    boundary_conditions::Tuple{B1,B2}   # Tuple of left and right boundary conditions, subject to the approval of PR #10
    end_time::Float64    # How long to simulate
    scheme::HyperbolicScheme{F,L} # Flux, Limiter
    source_term!::S  # Source term function. This can include reactons, electric field, and MMS terms
    source_potential!::SP #potential source term
    boundary_potential!::BP #boundary conditions potential
    saveat::Vector{Float64} #when to save
    timestepcontrol::Tuple{Float64,Bool} #sets timestep (first argument) if second argument false. if second argument (adaptive) true, given dt is ignored.
    callback::CB
    solve_energy::Bool
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
    nvariables = 0
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

    U = zeros(nvariables + 5, ncells + 2) # need to allocate room for ghost cells
    F = zeros(nvariables + 1, nedges)
    UL = zeros(nvariables + 1, nedges)
    UR = zeros(nvariables + 1, nedges)
    Q = zeros(nvariables + 1)
    A = Tridiagonal(ones(ncells - 1), ones(ncells), ones(ncells - 1)) #for potential
    b = zeros(ncells) #for potential equation
    B = zeros(ncells + 2)
    νan = zeros(ncells + 2)
    νc = zeros(ncells + 2)
    μ = zeros(ncells + 2)
   
    L_ch = 0.025

    cache = (; F, UL, UR, Q, A, b, B, νan, νc, μ)
    return U, cache
end

function update_exp!(dU, U, params, t) #get source and BCs for potential from params
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

    ####################################################################
    #PREPROCESS
    #calculate useful quantities relevant for potential, electron energy and fluid solve
    L_ch = 0.025
    fluid = fluids[1].species.element

    #=
    @inbounds for i in 1:(ncells + 2)
        #update electron temperature from energy using old density
        if params.solve_energy
            U[index.Tev, i] = max(1, U[index.nϵ, i]/3*2/U[index.ne, i])
        end
        U[index.ne, i] = max(1e-10, electron_density(@view(U[:, i]), fluid_ranges) / fluid.m)
        U[index.pe, i] = electron_pressure(U[index.ne, i], U[index.Tev, i])
        params.cache.νan[i] = get_v_an(z_cell[i], B[i], L_ch)
        params.cache.νc[i] = get_v_c(U[index.Tev, i], U[1, i]/fluid.m , U[index.ne, i], fluid.m)
        params.cache.μ[i] = cf_electron_transport(params.cache.νan[i], params.cache.νc[i], B[i])
        
    end
    =#

    ####################################################################
    #POTENTIAL MODULE
    #solve_potential!(U, params)

    ##############################################################
    #FLUID MODULE

    #fluid BCs
    apply_bc!(@views(U[1:index.lf, :]), params.BCs[1], :left)
    apply_bc!(@views(U[1:index.lf, :]), params.BCs[2], :right)

    #=
    #electron BCs
    Tev_anode = 3 #eV 
    Tev_cathode = 3 #eV
    left_state = [3/2*1e18*Tev_anode]
    right_state = [3/2*1e18*Tev_cathode]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    U[4, begin] = left_state[1]
    U[4, end] = right_state[1]
    #apply_bc!(@views(U[4, :]), BCs[1], :left)
    #apply_bc!(@views(U[4, :]), BCs[2], :right)=#

    #fluid computations, electron in implicit
    compute_edge_states!(@views(UL[1:index.lf, :]), @views(UR[1:index.lf, :]), @views(U[1:index.lf, :]), scheme)
    compute_fluxes!(@views(F[1:index.lf, :]), @views(UL[1:index.lf, :]), @views(UR[1:index.lf, :]), fluids, fluid_ranges, scheme)

    # Compute heavy species source terms
    @inbounds for i in 2:(ncells + 1)
        @turbo Q .= 0.0

        #fluid source term
        source_term!(Q, U, params, i)

        #@show Q[4]

        # Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        @tturbo @views @. dU[1:index.nϵ, i] = (F[:, left] - F[:, right]) / Δz + Q
    end

    return nothing
end

function update_imp!(dU, U, params, t)
    
    F, UL, UR, Q = params.cache.F, params.cache.UL, params.cache.UR, params.cache.Q
    index = params.index

    #electron BCs
    #
    Tev_anode = 3 #eV
    Tev_cathode = 3 #eV
    left_state = [3/2*1e18*Tev_anode]
    right_state = [3/2*1e18*Tev_cathode]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Dirichlet(right_state))

    U[index.nϵ, begin] = left_state[1]
    U[index.nϵ, end] = right_state[1]
    
    #electron computations, fluid in explicit
    scheme = HallThruster.HyperbolicScheme(HallThruster.upwind_electron!, identity, false)
    compute_edge_states!(@views(UL[index.nϵ, :]), @views(UR[index.nϵ, :]), @views(U[index.nϵ, :]), scheme)
    #println("UL and UR after comp edge: ", UL[4, :], UR[4, :])
    println("BEFORE CALC below: #########################")
    a =  UL[index.nϵ, :]
    b =  UR[index.nϵ, :]
    c =  U[index.nϵ, :]
    @show size(F)
    @show F
    compute_fluxes_electron!(@views(F[index.nϵ, :]), @views(UL[index.nϵ, :]), @views(UR[index.nϵ, :]), U, [HallThruster.Electron], [1:1], scheme, params)
    println("AFTER CALC below: #########################")
    @show F
    @show @test UL[index.nϵ, :] == a
    @show @test UR[index.nϵ, :] == b
    @show @test U[index.nϵ, :] == c
    #println("Flux after compute fluxes: ", F[4, :])
    

    ncells = size(U, 2) - 2
    z_edge = params.z_edge

    @inbounds for i in 2:(ncells + 1)
        # Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        #could save additional quantities of interest
        #=
        U[7, i] = electron_velocity(params, i)
        U[8, i] = - first_deriv_central_diff(params.cache.ϕ, params.z_cell, i)
        grad_pe = first_deriv_central_diff(params.cache.pe, params.z_cell, i)
        U[9, i] = grad_pe/e/params.cache.ne[i]=#

        #dU[index.nϵ, i] = (F[index.nϵ, left] - F[index.nϵ, right]) / Δz
        #@show F[index.nϵ, left]
        #@show F[index.nϵ, right]
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

    lf = fluid_ranges[end][end]
    index = (;lf = lf, nϵ = lf+1, Tev = lf+2, ne = lf+3, pe = lf+4, ϕ = lf+5)

    initial_condition!(@views(U[1:index.nϵ, :]), grid.cell_centers, sim.initial_condition, fluids)

    scheme = sim.scheme
    source_term! = sim.source_term!
    timestep = sim.timestepcontrol[1]
    adaptive = sim.timestepcontrol[2]
    tspan = (0.0, sim.end_time)

    reactions = load_ionization_reactions(species)
    landmark = load_landmark()

    BCs = sim.boundary_conditions

    precompute_bfield!(cache.B, grid.cell_centers)

    params = (; index, cache, fluids, fluid_ranges, species_range_dict, z_cell=grid.cell_centers,
              z_edge=grid.edges, cell_volume=grid.cell_volume, source_term!, reactions,
              scheme, BCs, dt=timestep, source_potential! = sim.source_potential!, 
              boundary_potential! = sim.boundary_potential!, landmark, solve_energy = sim.solve_energy)

    #PREPROCESS
    #make values in params available for first implicit timestep
    ncells = size(U, 2) - 2
    L_ch = 0.025
    fluid = fluids[1].species.element
    @inbounds for i in 1:(ncells + 2)
        U[index.ne, i] = max(1e-10, electron_density(@view(U[:, i]), fluid_ranges) / fluid.m)
        U[index.Tev, i] = max(1, U[index.nϵ, i]/3*2/U[index.ne, i])
        U[index.pe, i] = electron_pressure(U[index.ne, i], U[index.Tev, i])
        params.cache.νan[i] = get_v_an(params.z_cell[i], params.cache.B[i], L_ch)
        params.cache.νc[i] = get_v_c(U[index.Tev, i], U[1, i]/fluid.m , U[index.ne, i], fluid.m)
        params.cache.μ[i] = cf_electron_transport(params.cache.νan[i], params.cache.νc[i], params.cache.B[i])
    end

    solve_potential!(U, params)

    #tmp_prob = remake(prob, u0=convert.(eltype(params),prob.u0), p=params)
    if sim.solve_energy
        splitprob = SplitODEProblem{true}(update_imp!, update_exp!, U, tspan, params)
        sol = solve(splitprob, KenCarp3(autodiff = false); saveat=sim.saveat, callback=sim.callback,
        adaptive=adaptive, dt=timestep)
    else 
        prob = ODEProblem{true}(update_exp!, U, tspan, params)
        sol = solve(prob, Tsit5(); saveat=sim.saveat, callback=sim.callback,
        adaptive=adaptive, dt=timestep)
    end

    return sol
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