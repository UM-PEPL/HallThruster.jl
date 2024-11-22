@testset "Density calculation functions" begin
    mi = HallThruster.Xenon.m

    params = (
        index = (; ρi = [1, 2, 3]),
        config = (;
            ncharge = 3, propellant = HallThruster.Xenon, neutral_velocity = 100,
            thruster = HallThruster.SPT_100, anode_mass_flow_rate = 5e-6,
        ),
    )

    U = mi * [1e16; 2e16; 3e16;;]

    @test HallThruster.electron_density(U, params, 1) == 1e16 + 4e16 + 9e16

    A_ch = params.config.thruster.geometry.channel_area
    ρn = params.config.anode_mass_flow_rate / params.config.neutral_velocity / A_ch
    @test HallThruster.inlet_neutral_density(params.config) == ρn
end
