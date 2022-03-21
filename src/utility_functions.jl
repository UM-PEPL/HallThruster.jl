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