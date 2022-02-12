struct HyperbolicScheme{F,L}
    flux_function::F  # in-place flux function
    limiter::L # limiter
    reconstruct::Bool
end

Base.@kwdef mutable struct MultiFluidSimulation{IC,B1,B2,B3,B4,S,F,L,CB,SP,BP} #could add callback, or autoselect callback when in MMS mode
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

    U = zeros(nvariables + 7, ncells + 2) # need to allocate room for ghost cells
    F = zeros(nvariables + 1, nedges)
    UL = zeros(nvariables + 1, nedges)
    UR = zeros(nvariables + 1, nedges)
    Q = zeros(nvariables + 1)
    A = Tridiagonal(ones(ncells+1), ones(ncells+2), ones(ncells+1)) #for potential
    b = zeros(ncells+2) #for potential equation
    B = zeros(ncells + 2)
    νan = zeros(ncells + 2)
    νc = zeros(ncells + 2)
    μ = zeros(ncells + 2)

    L_ch = 0.025

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
    ε₀ = 3.0

    ##############################################################
    #FLUID MODULE

    #fluid BCs
    @views apply_bc!(U[1:index.lf, :], params.BCs[1], :left, ε₀, mi)
    @views apply_bc!(U[1:index.lf, :], params.BCs[2], :right, ε₀, mi)

    #fluid computations, electron in implicit
    @views compute_edge_states!(UL[1:index.lf, :], UR[1:index.lf, :], U[1:index.lf, :], scheme)
    @views compute_fluxes!(F[1:index.lf, :], UL[1:index.lf, :], UR[1:index.lf, :], fluids, fluid_ranges, scheme, U[index.pe, :])

    # Compute heavy species source terms
    @inbounds  for i in 2:(ncells + 1)
        @turbo Q .= 0.0

        #fluid source term (includes ionization and acceleration)
        source_term!(Q, U, params, i)

        #Compute dU/dt
        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        # Ion and neutral fluxes
        @turbo @views @. dU[1:index.lf, i] = (F[1:index.lf, left] - F[1:index.lf, right]) / Δz + Q[1:index.lf]

        # Ion diffusion term
        η = 0.0 * sqrt(2*e*U[index.Tev, i]/(3*mi))
        dU[index.lf, i] += η*(U[index.lf, i-1] - 2U[index.lf, i] + U[index.lf, i+1])/(Δz)^2

    end

    if !params.implicit_energy
        update_electron_energy!(dU, U, params, t)
    end

end

function update_electron_energy!(dU, U, params, t)
    #########################################################
    #ELECTRON SOLVE

    index = params.index
    Q = params.cache.Q
    μ = params.cache.μ
    z_edge = params.z_edge
    ncells = size(U, 2) - 2
    F = params.cache.F

    apply_bc_electron!(U, params.BCs[3], :left, index)
    apply_bc_electron!(U, params.BCs[4], :right, index)

    for i in 2:ncells+1

        source_electron_energy_landmark!(Q, U, params, i)

        dU[index.nϵ, i] = Q[index.nϵ]

        left = left_edge(i)
        right = right_edge(i)

        Δz = z_edge[right] - z_edge[left]

        # Upwinded first derivatives
        ue = U[index.ue, i]
        ne = U[index.ne, i]
        ue⁺ = max(ue, 0.0) / ue
        ue⁻ = min(ue, 0.0) / ue


        advection_term = (
            ue⁺ * (ue * U[index.nϵ, i] - U[index.ue, i-1] * U[index.nϵ, i-1]) -
            ue⁻ * (ue * U[index.nϵ, i] - U[index.ue, i+1] * U[index.nϵ, i+1])
        ) / Δz

        # this was attempting to use ∇⋅(niui) = ∇⋅(neue)
        #=advection_term = (
            U[index.Tev] * (F[2, right] - F[2, left]) +
            ne * ue⁺ * ue * (U[index.Tev, i] - U[index.Tev, i-1]) -
            ne * ue⁻ * ue * (U[index.Tev, i] - U[index.Tev, i+1])
        ) / Δz=#

        diffusion_term = -(
            ue⁺ * (-(μ[i-1] * U[index.nϵ, i-1] - μ[i] * U[index.nϵ, i]) * (U[index.Tev, i-1] - U[index.Tev, i])) -
            ue⁻ * (μ[i+1] * U[index.nϵ, i+1] - μ[i] * U[index.nϵ, i]) * (U[index.Tev, i+1] - U[index.Tev, i])
        )

        # Central second derivatives
        diffusion_term += μ[i] * U[index.nϵ, i] * (U[index.Tev, i-1] - 2U[index.Tev, i] + U[index.Tev, i+1])
        diffusion_term /= Δz^2

        dU[index.nϵ, i] += -5/3 * advection_term + 10/9 * diffusion_term
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

condition(u,t,integrator) = t < 1
function affect!(integrator)
    U, params = integrator.u, integrator.p

    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    index = params.index

    B = params.cache.B

    z_cell, z_edge, cell_volume = params.z_cell, params.z_edge, params.cell_volume
    ncells = size(U, 2) - 2

    #update useful quantities relevant for potential, electron energy and fluid solve
    L_ch = 0.025
    fluid = fluids[1].species.element

    @inbounds @views for i in 1:(ncells + 2) #pay attention as to whether J or eV in electron energy equ. 
        @views U[index.ne, i] = max(1e13, electron_density(U[:, i], fluid_ranges) / fluid.m)
        U[index.Tev, i] = max(0.1, U[index.nϵ, i]/U[index.ne, i])
        U[index.pe, i] = U[index.nϵ, i]
        @views U[index.grad_ϕ, i] = first_deriv_central_diff_pot(U[index.ϕ, :], params.z_cell, i)
        U[index.ue, i] = electron_velocity(U, params, i)
        params.cache.νan[i] = get_v_an(z_cell[i], B[i], L_ch)
        params.cache.νc[i] = electron_collision_freq(U[index.Tev, i], U[1, i]/fluid.m , U[index.ne, i], fluid.m)
        params.cache.μ[i] = electron_mobility(params.cache.νan[i], params.cache.νc[i], B[i])
    end

    # update electrostatic potential
    solve_potential!(U, params)
end

function run_simulation(sim) #put source and Bcs potential in params
    species, fluids, fluid_ranges, species_range_dict = configure_simulation(sim)
    grid = sim.grid

    U, cache = allocate_arrays(sim)

    lf = fluid_ranges[end][end]
    index = (;lf = lf, nϵ = lf+1, Tev = lf+2, ne = lf+3, pe = lf+4, ϕ = lf+5, grad_ϕ = lf+6, ue = lf+7)

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

    OVS = Array{Union{Nothing, Bool}}(nothing, 1)
    OVS[1] = false

    params = (; OVS, index, cache, fluids, fluid_ranges, species_range_dict, z_cell=grid.cell_centers,
              z_edge=grid.edges, cell_volume=grid.cell_volume, source_term!, reactions,
              scheme, BCs, dt=timestep, source_potential! = sim.source_potential!,
              boundary_potential! = sim.boundary_potential!, landmark, implicit_energy = sim.solve_energy)

    #PREPROCESS
    #make values in params available for first implicit timestep
    ncells = size(U, 2) - 2
    L_ch = 0.025
    fluid = fluids[1].species.element
    @inbounds for i in 1:(ncells + 2)
        @views U[index.ne, i] = max(1e-10, electron_density(U[:, i], fluid_ranges) / fluid.m)
        U[index.Tev, i] = max(0.1, U[index.nϵ, i]/U[index.ne, i])
        U[index.pe, i] = U[index.nϵ, i]/3*2
        @views U[index.grad_ϕ, i] = first_deriv_central_diff(U[index.ϕ, :], params.z_cell, i)
        params.cache.νan[i] = get_v_an(params.z_cell[i], params.cache.B[i], L_ch)
        params.cache.νc[i] = electron_collision_freq(U[index.Tev, i], U[1, i]/fluid.m , U[index.ne, i], fluid.m)
        params.cache.μ[i] = electron_mobility(params.cache.νan[i], params.cache.νc[i], params.cache.B[i])
        U[index.ue, i] = electron_velocity(U, params, i)
    end

    solve_potential!(U, params)
    cb = DiscreteCallback(condition, affect!, save_positions=(false,false))

    implicit_energy = sim.solve_energy

    if implicit_energy
        splitprob = SplitODEProblem{true}(update_electron_energy!, update_heavy_species!, U, tspan, params)
        sol = solve(splitprob, KenCarp4(autodiff=false); saveat=sim.saveat, callback=cb,
        adaptive=adaptive, dt=timestep, dtmax = timestep)
    else
        #alg = SSPRK22()
        alg = AutoTsit5(Rosenbrock23())
        prob = ODEProblem{true}(update_heavy_species!, U, tspan, params)
        sol = solve(prob, alg; saveat=sim.saveat, callback=cb,
        adaptive=adaptive, dt=timestep, dtmax=timestep)
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