using HallThruster: HallThruster as het

function test_boundaries()
    e = het.e

    config = het.Config(
        ncharge = 2,
        domain = (0, 2),
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

    _, params = het.setup_simulation(config, het.SimParams(duration = 1.0e-3, grid = het.EvenGrid(2)))

    # Anode ion velocity less than bohm velocity
    prop = config.propellants[1]
    mi = prop.gas.m
    Ti = prop.ion_temperature_K
    mdot_a = prop.flow_rate_kg_s
    un = prop.velocity_m_s

    ni_1 = 1.0e17
    ni_2 = 1.0e16

    u_bohm_1 = sqrt((het.kB * Ti + e * params.Te_L) / mi)
    u_bohm_2 = sqrt((het.kB * Ti + 2 * e * params.Te_L) / mi)

    (; continuity, isothermal) = params.fluid_containers
    @. continuity[1].density = ni_1 * mi * 2
    @. isothermal[1].density = ni_1 * mi
    @. isothermal[1].momentum = -mi * ni_1 * u_bohm_1 / 2
    @. isothermal[2].density = ni_2 * mi
    @. isothermal[2].momentum = -mi * ni_2 * u_bohm_2 / 2

    ingestion_density = params.ingestion_density
    anode_bc = params.anode_bc

    het.apply_left_boundary!(params.fluid_containers, params.cache, Ti, mdot_a, ingestion_density, anode_bc)
    het.apply_right_boundary!(params.fluid_containers)

    # When ion velocity at left boundary is less than then bohm speed, it should be accelerated to reach the bohm speed
    nn_B = het.background_neutral_density(config)
    un_B = het.background_neutral_velocity(config)
    boundary_ion_flux = [ion.momentum[1] for ion in isothermal]
    @test continuity[1].density[1] ≈ het.inlet_neutral_density(config) -
        sum(boundary_ion_flux) / un + nn_B * un_B / un * config.neutral_ingestion_multiplier
    @test boundary_ion_flux[1] / isothermal[1].density[1] ≈ -u_bohm_1
    @test boundary_ion_flux[2] / isothermal[2].density[1] ≈ -u_bohm_2

    # Neumann BC for all species
    for fluid in [continuity..., isothermal...]
        @test fluid.density[end] == fluid.density[end - 1]
        @test fluid.momentum[end] == fluid.momentum[end - 1]
    end

    # Anode ion velocity greater than Bohm speed
    @. continuity[1].density = ni_1 * mi * 2
    @. isothermal[1].density = ni_1 * mi
    @. isothermal[1].momentum = -mi * ni_1 * u_bohm_1 * 2
    @. isothermal[2].density = ni_2 * mi
    @. isothermal[2].momentum = -mi * ni_2 * u_bohm_2 * 2

    # when ion velocity at left boundary is greater than bohm speed, ions have Neumann BC
    het.apply_left_boundary!(params.fluid_containers, params.cache, Ti, mdot_a, ingestion_density, anode_bc)
    het.apply_right_boundary!(params.fluid_containers)
    boundary_ion_flux = [ion.momentum[1] for ion in isothermal]
    @test continuity[1].density[1] ≈ het.inlet_neutral_density(config) -
        sum(boundary_ion_flux) / un +
        nn_B * un_B / un * config.neutral_ingestion_multiplier
    @test boundary_ion_flux[1] / isothermal[1].density[1] ≈ -2 * u_bohm_1
    @test boundary_ion_flux[2] / isothermal[2].density[1] ≈ -2 * u_bohm_2

    # Neumann BC for all species
    for fluid in [continuity..., isothermal...]
        @test fluid.density[end] == fluid.density[end - 1]
        @test fluid.momentum[end] == fluid.momentum[end - 1]
    end
    return
end

test_boundaries()
