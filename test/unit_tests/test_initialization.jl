@testset "Initialization" begin

    domain = (0.0, 0.08)
    thruster = HallThruster.SPT_100
    config = (;
        ncharge = 3,
        propellant = HallThruster.Xenon,
        thruster,
        domain,
        cathode_Te = 5.0,
        anode_Te = 3.0,
        neutral_velocity = 300.0,
        neutral_temperature = 100.0,
        ion_temperature = 300.0,
        initial_condition = HallThruster.DefaultInitialization(),
        anode_mass_flow_rate = 3e-6,
        discharge_voltage = 500.0
    )

    mi = config.propellant.m

    ncells = 100
    fluids, fluid_ranges, species, species_range_dict = HallThruster.configure_fluids(config)
    grid = HallThruster.generate_grid(config.thruster.geometry, ncells, domain)
    U, cache = HallThruster.allocate_arrays(grid, fluids)
    index = HallThruster.configure_index(fluid_ranges)

    params = (;
        index,
        cache,
        z_cell = grid.cell_centers,
        config,
    )

    HallThruster.initialize!(U, params)

    @test abs((U[index.ρn, 1] / U[index.ρn,end]) - 100) < 1

    ne = [HallThruster.electron_density(U, params, i) for i in 1:ncells+2]
    ϵ = U[index.nϵ, :] ./ ne

    ui = [
        U[index.ρiui[Z], :] ./ U[index.ρi[Z], :] for Z in 1:config.ncharge
    ]

    @test ui[1][1] ≈ -sqrt(HallThruster.e * config.anode_Te / mi)
    @test ui[2][1] ≈ -sqrt(2 * HallThruster.e * config.anode_Te / mi)
    @test ui[3][1] ≈ -sqrt(3 * HallThruster.e * config.anode_Te / mi)

    @test abs(2/3 * ϵ[1] - config.anode_Te) < 0.1
    @test abs(2/3 * ϵ[end] - config.cathode_Te) < 0.1

    struct NewInitialization <: HallThruster.InitialCondition end
    @test_throws ArgumentError HallThruster.initialize!(U, params, NewInitialization())
end

