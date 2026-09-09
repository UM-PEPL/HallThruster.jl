using HallThruster: HallThruster as het

Landmark_config = het.Config(;
    thruster = het.SPT_100,
    domain = (0.0, 0.08),
    discharge_voltage = 300.0,
    anode_mass_flow_rate = 5.0e-6,
    LANDMARK = true,
)

@test_throws ErrorException het.run_simulation(
    Landmark_config,
    het.SimParams(;
        dt = 5.0e-9, duration = 4.0e-9,
        grid = het.EvenGrid(2), num_save = 10,
    ),
)

config = het.Config(;
    thruster = het.SPT_100,
    domain = (0.0, 0.08),
    discharge_voltage = 300.0,
    anode_mass_flow_rate = 5.0e-6,
)

params = het.setup_simulation(
    config,
    het.SimParams(;
        dt = 5.0e-9, duration = 0.0e-9, grid = het.EvenGrid(9),
        num_save = 10, adaptive = true, CFL = 0.9,
    ),
)
@test params.simulation.CFL == 0.9
