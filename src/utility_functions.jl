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

"""
    anode_sheath(Te_a, j_d, un_ia, n_ea)
calculate anode sheath following model presented in Sahu 2020, 
Full fluid moment model for low temperature magnetized plasmas
"""
function anode_sheath(Te_a, j_ed, n_ea)
    if Te_a < 0
        Te_a = 0.1
    end
    U_a = -kB*Te_a/e*log(4(j_ed)/(n_ea*sqrt(8*kB*Te_a/pi/mₑ)))
    return U_a
end

function compute_ecurrent_density_anode(sol)
    index = sol.params.index
    current = zeros(length(sol.t))
    mi = sol.params.propellant.m
    for i in 1:length(sol.t)
        (;ue, ne) = sol.savevals[i]
        current[i] = ne[2] * ue[2]*HallThruster.e
    end
    return current
end

#=
#should try this with cell 1 and cell 2
for i in 1:length(anode_ϕ)
    anode_sheath[i] = HallThruster.anode_sheath(sol.savevals[i].Tev[2], current_je_anode[i], sol.savevals[i].ne[2])
end=#
