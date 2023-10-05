@testset "Boundary conditions" begin
    e = HallThruster.e
    me = HallThruster.me

    config = (
        ncharge = 2,
        propellant = HallThruster.Xenon,
        anode_mass_flow_rate = 5e-4,
        neutral_velocity = 100,
        thruster = HallThruster.SPT_100,
        min_number_density = 1e6,
        min_electron_temperature = 1.0,
        LANDMARK = true,
        anode_boundary_condition = :dirichlet,
        solve_background_neutrals = :true,
        neutral_temperature = 300.0,
        ion_temperature = 300.0,
        background_pressure = 6e-4,
        background_neutral_temperature = 150.0,
    )

    fluids, fluid_ranges, species, species_range_dict = HallThruster.configure_fluids(config)
    index = HallThruster.configure_index(fluids, fluid_ranges)

    mi = config.propellant.m

    background_neutral_velocity = 0.25 * sqrt(8 * HallThruster.kB * config.background_neutral_temperature / π / mi)
    background_neutral_density = mi * config.background_pressure / HallThruster.kB / config.background_neutral_temperature
    background_neutral_flux = background_neutral_density * background_neutral_velocity

    params = (;
        Te_L = 3.0,
        Te_R = 3.0,
        A_ch = config.thruster.geometry.channel_area,
        config,
        index,
        z_cell = [0., 1., 2.],
        cache = (
            Tev = [3.0, 3.0],
            ϕ = [300.0, 300.0],
            channel_area = ones(2) * config.thruster.geometry.channel_area
        ),
        num_neutral_fluids = 1,
        fluids,
        fluid_ranges,
        species,
        species_range_dict,
        background_neutral_density,
        background_neutral_velocity,
        background_neutral_flux,
    )

    u_bohm_1 = sqrt((HallThruster.kB * config.ion_temperature + e * params.Te_L) / mi)
    u_bohm_2 = sqrt((HallThruster.kB * config.ion_temperature + 2 * e * params.Te_L) / mi)

    ni_1 = 1e17
    ni_2 = 1e16
    ui_1_1 = u_bohm_1 / 2
    ui_1_2 = u_bohm_1 * 2
    ui_2_1 = u_bohm_2 / 2
    ui_2_2 = u_bohm_2 * 2
    Te_1 = 5.0
    Te_2 = 1.0

    ne = ni_1 + 2 * ni_2

    # State where anode ion velocity is less than bohm velocity
    U_1 = [
        mi * ni_1 * 2,
        mi * ni_1,
        -mi * ni_1 * u_bohm_1 / 2,
        mi * ni_2,
        -mi * ni_2 * u_bohm_2 / 2,
        ne * Te_1
    ]

    # State where anode ion velocity is greater than bohm velocity
    U_2 = [
        mi * ni_1,
        mi * ni_1,
        -mi * ni_1 * u_bohm_1 * 2,
        mi * ni_2,
        -mi * ni_2 * u_bohm_2 * 2,
        ne * Te_1
    ]

    U_b = zeros(length(U_1))

    U1 = [U_b U_1 U_1 U_b]
    U2 = [U_b U_2 U_2 U_b]

    Vs = 0.0

    un = config.neutral_velocity

    # when ion velocity at left boundary is less than bohm speed, it should be accelerated
    # to reach bohm speed
    HallThruster.left_boundary_state!(U_b, U1, params)
    @test U_b[index.ρn[1]] ≈ HallThruster.inlet_neutral_density(config) - (U_b[index.ρiui[1]] + U_b[index.ρiui[2]]) / config.neutral_velocity + background_neutral_density * background_neutral_velocity / un
    @test U_b[index.ρiui[1]] / U_b[index.ρi[1]] == -u_bohm_1
    @test U_b[index.ρiui[2]] / U_b[index.ρi[2]] == -u_bohm_2

    # when ion velocity at left boundary is greater than bohm speed, ions have Neumann BC
    HallThruster.left_boundary_state!(U_b, U2, params)
    @test U_b[index.ρn[1]] ≈ HallThruster.inlet_neutral_density(config) - (U_b[index.ρiui[1]] + U_b[index.ρiui[2]]) / config.neutral_velocity + background_neutral_density * background_neutral_velocity / un
    @test U_b[index.ρiui[1]] == U_2[index.ρiui[1]]
    @test U_b[index.ρiui[2]] == U_2[index.ρiui[2]]
    @test U_b[index.ρiui[1]] / U_b[index.ρi[1]] == -2 * u_bohm_1
    @test U_b[index.ρiui[2]] / U_b[index.ρi[2]] == -2 * u_bohm_2

    # Right boundary condition should be Neumann for all species
    HallThruster.right_boundary_state!(U_b, U1, params)
    @test U_b[index.ρn[1]] ≈ U_1[index.ρn[1]]
    @test U_b[index.ρi[1]] ≈ U_1[index.ρi[1]]
    @test U_b[index.ρiui[1]] ≈ U_1[index.ρiui[1]]
    @test U_b[index.ρi[2]] ≈ U_1[index.ρi[2]]
    @test U_b[index.ρiui[2]] ≈ U_1[index.ρiui[2]]

    HallThruster.right_boundary_state!(U_b, U2, params)
    @test U_b[index.ρn[1]] ≈ U_2[index.ρn[1]]
    @test U_b[index.ρi[1]] ≈ U_2[index.ρi[1]]
    @test U_b[index.ρiui[1]] ≈ U_2[index.ρiui[1]]
    @test U_b[index.ρi[2]] ≈ U_2[index.ρi[2]]
    @test U_b[index.ρiui[2]] ≈ U_2[index.ρiui[2]]
end
