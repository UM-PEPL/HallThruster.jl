left_edge(i) = i - 1
right_edge(i) = i

@inline electron_density(U, p, i) = sum(Z * U[p.index.ρi[Z], i] for Z in 1:p.config.ncharge) / p.config.propellant.m

@inline inlet_neutral_density(config) = config.anode_mass_flow_rate / config.neutral_velocity / channel_area(config.geometry)
