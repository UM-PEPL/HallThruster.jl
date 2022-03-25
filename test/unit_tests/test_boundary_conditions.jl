@testset "Boundary conditions" begin
    e = HallThruster.e
    me = HallThruster.me

    index = (
        ρn = 1,
        ρi = [2, 4],
        ρiui = [3, 5],
        nϵ = 6,
    )

    config = (
        ncharge = 2,
        propellant = HallThruster.Xenon,
        anode_mass_flow_rate = 5e-4,
        neutral_velocity = 100,
        thruster = HallThruster.SPT_100,
        min_number_density = 1e6,
        min_electron_temperature = 1.0
    )

    params = (;
        Te_L = 3.0,
        Te_R = 3.0,
        A_ch = config.thruster.geometry.channel_area,
        config,
        index,
    )

    mi = config.propellant.m

    u_bohm_1 = sqrt(e * 2/3 * params.Te_L / mi)
    u_bohm_2 = sqrt(2) * u_bohm_1

    ni_1 = 1e17
    ni_2 = 1e16
    ui_1_1 = u_bohm_1 / 2
    ui_1_2 = u_bohm_1 * 2
    ui_2_1 = u_bohm_2 / 2
    ui_2_2 = u_bohm_2 * 2
    Te_1 = 5.0
    Te_2 = 1.0

    ne = ni_1 + 2 * ni_2

    U_1 = [
        mi * ni_1 * 2,
        mi * ni_1,
        -mi * ni_1 * ui_1_1,
        mi * ni_2,
        -mi * ni_2 * ui_2_1,
        ne * Te_1
    ]

    U_2 = [
        mi * ni_1,
        mi * ni_1,
        -mi * ni_1 * ui_1_2,
        mi * ni_2,
        -mi * ni_2 * ui_2_2,
        ne * Te_1
    ]

    U_b = zeros(length(U_1))

    HallThruster.left_boundary_state!(U_b, U_1, params)

    @test U_b[index.ρn] ≈ HallThruster.inlet_neutral_density(config) - (U_1[index.ρiui[1]] + U_1[index.ρiui[2]]) / config.neutral_velocity
    @test U_b[index.ρiui[1]] == U_1[index.ρiui[1]]
    @test U_b[index.ρiui[2]] == U_1[index.ρiui[2]]
    @test U_b[index.ρiui[1]] / U_b[index.ρi[1]] == -u_bohm_1
    @test U_b[index.ρiui[2]] / U_b[index.ρi[2]] == -u_bohm_2
    @test U_b[index.nϵ] == HallThruster.electron_density([U_b;;], params, 1) * params.Te_L

    HallThruster.left_boundary_state!(U_b, U_2, params)
    @test U_b[index.ρn] ≈ HallThruster.inlet_neutral_density(config) - (U_2[index.ρiui[1]] + U_2[index.ρiui[2]]) / config.neutral_velocity
    @test U_b[index.ρiui[1]] == U_2[index.ρiui[1]]
    @test U_b[index.ρiui[2]] == U_2[index.ρiui[2]]
    @test U_b[index.ρiui[1]] / U_b[index.ρi[1]] == -2 * u_bohm_1
    @test U_b[index.ρiui[2]] / U_b[index.ρi[2]] == -2 * u_bohm_2
    @test U_b[index.nϵ] == HallThruster.electron_density([U_b;;], params, 1) * params.Te_L

    HallThruster.right_boundary_state!(U_b, U_1, params)
    @test all(U_b[1:end-1] .≈ U_1[1:end-1])
    @test U_b[index.nϵ] == HallThruster.electron_density([U_b;;], params, 1) * params.Te_R

    HallThruster.right_boundary_state!(U_b, U_2, params)
    @test all(U_b[1:end-1] .≈ U_2[1:end-1])
    @test U_b[index.nϵ] == HallThruster.electron_density([U_b;;], params, 1) * params.Te_R
end