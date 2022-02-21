#for OVS
Base.@kwdef mutable struct Verification
    potential ::Bool
    fluid :: Bool
    energy ::Bool
end

Base.@kwdef struct MagneticField
    z::Vector{Float64}
    B::Vector{Float64}
end

function no_source!(args...) end

Base.@kwdef mutable struct HallThrusterConfig{F, L, A<:AnomalousTransportModel, S}
    name::String
    propellant::Gas
    ncells::Int  = 50
    ncharge::Int = 1
    restart_file::Union{String, Nothing} = nothing
    magnetic_field::Union{MagneticField, Nothing} = nothing
    simulation_time::Float64
    nsave::Int                    = 100
    adaptive::Bool                = true
    dt::Float64                   = 2e-9
    dtmax::Float64                = 2e-9
    implicit_energy::Bool         = false
    cathode_Te::Float64           = 2.0
    anode_Te::Float64             = 2.0
    cathode_potential::Float64    = 0.0
    anode_potential::Float64      = 300
    anode_mass_flow_rate::Float64 = 5e-3
    geometry::Geometry1D
    anom_model::A = TwoZoneBohm(1/160, 1/16)
    radial_loss_coefficients::Tuple{Float64, Float64} = (0.1, 0.1)
    wall_collision_frequencies:: Tuple{Float64, Float64} = (1e7, 0.0)
    energy_equation::Symbol = :LANDMARK
    flux::F = HLLE!
    limiter::L = no_limiter
    reconstruct::Bool = false
    ionization_coeffs::Symbol = :BOLSIG
    floor_Te::Float64 = 2.0
    floor_ne::Float64 = 1e12
    solve_ion_energy::Bool = false
    ion_temperature::Float64
    neutral_velocity::Float64
    neutral_temperature::Float64
    OVS::Verification = Verification(false, false, false)
    source_term!::S = no_source!
    ion_diffusion_coeff::Float64 = 0.0
    coupled_method::Bool = true
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

function allocate_arrays(grid, fluids) #rewrite allocate arrays as function of set of equations, either 1, 2 or 3
    # Number of variables in the state vector U
    nvariables = 0
    for i in 1:length(fluids)
        if fluids[i].conservation_laws.type == :ContinuityOnly
            nvariables += 1
        elseif fluids[i].conservation_laws.type == :IsothermalEuler
            nvariables += 2
        elseif fluids[i].conservation_laws.type == :EulerEquations
            nvariables += 3
        end
    end

    ncells = grid.ncells
    nedges = grid.ncells + 1

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

function inlet_neutral_density(config)
    un = config.neutral_velocity
    A = channel_area(config.geometry)
    m_atom = config.propellant.m
    nn = config.anode_mass_flow_rate / un / A / m_atom
    return nn
end

function initialize_solution!(U, index, grid, fluids, fluid_ranges, config)
    mi = config.propellant.m

    function Te_func(z, L_ch) #will be gone soon
        return 30 * exp(-(2 * (z - L_ch) / (0.033/0.025 * L_ch))^2)
    end

    for i in 1:length(grid.cell_centers)
        U[index.ρn, i] = inlet_neutral_density(config) * mi
        for (f, r) in zip(fluids[2:end], fluid_ranges[2:end])
            U[r[1]] = config.floor_ne * mi
            U[r[2]] = config.floor_ne * config.neutral_velocity * mi
            if length(r) == 3
                U[r[3]] = mi * config.floor_ne * (0.5 * config.neutral_velocity^2 + cv(f) * config.ion_temperature)
            end
        end
        z = grid.cell_centers[i]
        L_ch = config.geometry.channel_length
        domain = config.geometry.domain
        U[index.ne, i] = electron_density(U[:, i], fluid_ranges)
        U[index.Tev, i] = Te_func(z, L_ch)
        U[index.nϵ, i] = config.anode_Te * U[index.ne, i]
        U[index.ϕ, i] = 300 - (300 / (domain[end] - domain[1])) * (z - domain[1])
        U[index.grad_ϕ, i] = 0.0
        U[index.ue, i] = 0.0
    end
end

"""
    compute_bfield!(B, field::MagneticField, z_cell)
interpolate user-provided B-field onto grid
"""
function compute_bfield!(B, z_cell, field::MagneticField)
    itp = LinearInterpolation(field.z, field.B)
    B .= itp.(z_cell)
end

function configure_simulation(config)

    fluids, fluid_ranges, species, species_range_dict = configure_fluids(config)
    landmark = load_landmark()
    loss_coeff = landmark.loss_coeff

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

    use_restart = config.restart_file !== nothing

    index = configure_index(fluid_ranges)

    grid = generate_grid(config.geometry, config.ncells)
    U, cache = allocate_arrays(grid, fluids)
    initialize_solution!(U, index, grid, fluids, fluid_ranges, config)

    if use_restart
       initialize_from_restart!(U, cache, config.restart_file, grid, fluid_ranges)
    else
        if isnothing(config.magnetic_field)
            throw(ArgumentError("Magnetic field must be provided for cold starts."))
        else
            compute_bfield!(cache.B, grid.cell_centers, config.magnetic_field)
        end
    end

    scheme = HyperbolicScheme(config.flux, config.limiter, config.reconstruct)

    apply_acceleration! = if config.coupled_method
        apply_ion_acceleration_coupled!
    else
        apply_ion_acceleration!
    end

    function ion_source_term!(Q, U, params, i)
        apply_acceleration!(Q, U, params, i)
        apply_reactions!(Q, U, params, i)
        config.source_term!(Q, U, params, i)
    end

    params = (
        mi = config.propellant.m,
        νw = config.wall_collision_frequencies,
        νϵ = config.radial_loss_coefficients,
        ncharge = config.ncharge,
        mdot_a = config.anode_mass_flow_rate,
        L_ch = config.geometry.channel_length,
        A_ch = channel_area(config.geometry.outer_radius, config.geometry.inner_radius),
        ϕ_L = config.anode_potential,
        ϕ_R = config.cathode_potential,
        Te_L = config.anode_Te,
        Te_R = config.cathode_Te,
        OVS = config.OVS,
        z_cell = grid.cell_centers,
        z_edge = grid.edges,
        cell_volume = grid.cell_volume,
        ion_source_term! = ion_source_term!,
        implicit_energy = config.implicit_energy,
        anom_model = config.anom_model,
        coupled = config.coupled_method,
        δ = config.ion_diffusion_coeff,
        reactions = ionization_reactions,
        index, cache, fluids, fluid_ranges, species_range_dict, scheme, loss_coeff,
    )

    @show species_range_dict
    @show fluid_ranges

    println("Simulation configured.")

    return U, params
end