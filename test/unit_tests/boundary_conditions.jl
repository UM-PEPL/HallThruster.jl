using HallThruster: HallThruster as het

e = het.e

config = het.Config(
    ncharge = 2,
    domain = (0, 1),
    discharge_voltage = 300.0,
    propellant = het.Xenon,
    anode_mass_flow_rate = 5.0e-4,
    neutral_velocity = 100,
    thruster = het.SPT_100,
    LANDMARK = true,
    conductivity_model = het.LANDMARK_conductivity(),
    anode_boundary_condition = :dirichlet,
    neutral_temperature_K = 300.0,
    ion_temperature_K = 300.0,
    background_pressure_Torr = 6.0e-4,
    background_temperature_K = 150.0,
    neutral_ingestion_multiplier = 1.5,
)

fluids, fluid_ranges, species, species_range_dict, is_velocity_index = het.configure_fluids(config)
index = het.configure_index(fluids, fluid_ranges)

mi = config.propellant.m

background_neutral_velocity = 0.25 * sqrt(
    8 * het.kB *
        config.background_temperature_K / π / mi
)
background_neutral_density = mi * config.background_pressure_Torr / het.kB /
    config.background_temperature_K
background_neutral_flux = background_neutral_density * background_neutral_velocity

params = (;
    het.params_from_config(config)...,
    Te_L = 3.0,
    Te_R = 3.0,
    A_ch = config.thruster.geometry.channel_area,
    index,
    z_cell = [0.0, 1.0, 2.0],
    cache = (
        Tev = [3.0, 3.0],
        ϕ = [300.0, 300.0],
        channel_area = ones(2) * config.thruster.geometry.channel_area,
    ),
    fluids,
    fluid_ranges,
    species,
    species_range_dict,
    is_velocity_index,
    background_neutral_density,
    background_neutral_velocity,
    background_neutral_flux,
)

u_bohm_1 = sqrt((het.kB * config.ion_temperature_K + e * params.Te_L) / mi)
u_bohm_2 = sqrt(
    (het.kB * config.ion_temperature_K + 2 * e * params.Te_L) /
        mi
)

ni_1 = 1.0e17
ni_2 = 1.0e16
Te_1 = 5.0

ne = ni_1 + 2 * ni_2

# State where anode ion velocity is less than bohm velocity
U_1 = [
    mi * ni_1 * 2,
    mi * ni_1,
    -mi * ni_1 * u_bohm_1 / 2,
    mi * ni_2,
    -mi * ni_2 * u_bohm_2 / 2,
    ne * Te_1,
]

# State where anode ion velocity is greater than bohm velocity
U_2 = [
    mi * ni_1,
    mi * ni_1,
    -mi * ni_1 * u_bohm_1 * 2,
    mi * ni_2,
    -mi * ni_2 * u_bohm_2 * 2,
    ne * Te_1,
]

U_b = zeros(length(U_1))

U1 = [U_b U_1 U_1 U_b]
U2 = [U_b U_2 U_2 U_b]

un = config.neutral_velocity

# when ion velocity at left boundary is less than bohm speed, it should be accelerated
# to reach bohm speed
het.left_boundary_state!(U_b, U1, params)
@test U_b[index.ρn] ≈
    het.inlet_neutral_density(config) -
    (U_b[index.ρiui[1]] + U_b[index.ρiui[2]]) / config.neutral_velocity +
    background_neutral_density * background_neutral_velocity / un *
    config.neutral_ingestion_multiplier
@test U_b[index.ρiui[1]] / U_b[index.ρi[1]] == -u_bohm_1
@test U_b[index.ρiui[2]] / U_b[index.ρi[2]] == -u_bohm_2

# when ion velocity at left boundary is greater than bohm speed, ions have Neumann BC
het.left_boundary_state!(U_b, U2, params)
@test U_b[index.ρn] ≈
    het.inlet_neutral_density(config) -
    (U_b[index.ρiui[1]] + U_b[index.ρiui[2]]) / config.neutral_velocity +
    background_neutral_density * background_neutral_velocity / un *
    config.neutral_ingestion_multiplier
@test U_b[index.ρiui[1]] == U_2[index.ρiui[1]]
@test U_b[index.ρiui[2]] == U_2[index.ρiui[2]]
@test U_b[index.ρiui[1]] / U_b[index.ρi[1]] == -2 * u_bohm_1
@test U_b[index.ρiui[2]] / U_b[index.ρi[2]] == -2 * u_bohm_2

# Right boundary condition should be Neumann for all species
het.right_boundary_state!(U_b, U1, params)
@test U_b[index.ρn] ≈ U_1[index.ρn]
@test U_b[index.ρi[1]] ≈ U_1[index.ρi[1]]
@test U_b[index.ρiui[1]] ≈ U_1[index.ρiui[1]]
@test U_b[index.ρi[2]] ≈ U_1[index.ρi[2]]
@test U_b[index.ρiui[2]] ≈ U_1[index.ρiui[2]]

het.right_boundary_state!(U_b, U2, params)
@test U_b[index.ρn] ≈ U_2[index.ρn]
@test U_b[index.ρi[1]] ≈ U_2[index.ρi[1]]
@test U_b[index.ρiui[1]] ≈ U_2[index.ρiui[1]]
@test U_b[index.ρi[2]] ≈ U_2[index.ρi[2]]
@test U_b[index.ρiui[2]] ≈ U_2[index.ρiui[2]]
