using HallThruster: HallThruster as het

function test_errors()
    return @testset "Errors and Warnings" begin
        Landmark_config = het.Config(;
            thruster = het.SPT_100,
            domain = (0.0, 0.08),
            discharge_voltage = 300.0,
            anode_mass_flow_rate = 5.0e-6,
            LANDMARK = true,
        )

        @test_throws ErrorException het.run_simulation(
            Landmark_config; dt = 5.0e-9, duration = 4.0e-9,
            grid = het.EvenGrid(2), nsave = 10,
        )

        config = het.Config(;
            thruster = het.SPT_100,
            domain = (0.0, 0.08),
            discharge_voltage = 300.0,
            anode_mass_flow_rate = 5.0e-6,
        )

        @test_logs (
            :warn,
            "CFL for adaptive timestepping set higher than stability limit of 0.8. Setting CFL to 0.799.",
        ) match_mode = :any het.run_simulation(
            config; dt = 5.0e-9, duration = 0.0e-9, grid = het.EvenGrid(2),
            nsave = 10, adaptive = true, CFL = 0.9,
        )
    end
end

test_errors()
