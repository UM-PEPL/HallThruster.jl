using HallThruster: HallThruster as het

function test_density_calculations()
    @testset "Density calculation functions" begin
        mi = HallThruster.Xenon.m

        config = het.Config(;
            ncharge = 3, propellant = HallThruster.Xenon, neutral_velocity = 100,
            thruster = HallThruster.SPT_100, anode_mass_flow_rate = 5e-6, discharge_voltage = 300.0,
            domain = (0.0, 1.0),
        )

        params = (;
            het.params_from_config(config)...,
            index = (; ρi = [1, 2, 3]),
        )

        U = mi * [1e16; 2e16; 3e16;;]

        @test HallThruster.electron_density(U, params, 1) == 1e16 + 4e16 + 9e16

        A_ch = config.thruster.geometry.channel_area
        ρn = config.anode_mass_flow_rate / config.neutral_velocity / A_ch
        @test HallThruster.inlet_neutral_density(config) == ρn
    end
end
