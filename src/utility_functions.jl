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

function precompute_bfield!(B, zs)
    B_max = 0.015
    L_ch = 0.025
    for (i, z) in enumerate(zs)
        B[i] = B_field(B_max, z, L_ch)
    end
end