using HallThruster: HallThruster as het

function test_boundaries()
    e = het.e

    V_d = 300.0
    V_cc = 20.0
    MIN_NUMBER_DENSITY = 1.0

    config = het.Config(
        ncharge = 2,
        domain = (0, 0.08),
        discharge_voltage = V_d,
        propellant = het.Xenon,
        anode_mass_flow_rate = 5.0e-6,
        neutral_velocity = 150,
        thruster = het.SPT_100,
        LANDMARK = true,
        conductivity_model = het.LANDMARK_conductivity(),
        anode_boundary_condition = :dirichlet,
        neutral_temperature_K = 300.0,
        ion_temperature_K = 300.0,
        background_pressure_Torr = 5.0e-5,
        background_temperature_K = 150.0,
        neutral_ingestion_multiplier = 1.5,
        cathode_coupling_voltage = V_cc,
    )

    params = het.setup_simulation(config, het.SimParams(duration = 1.0e-3, grid = het.EvenGrid(20)))

    # Check boundary handling of electron temperature solver
    # (electron temperature is solved, including boundaries, during setup)
    (; Tev, ϕ) = params.cache
    @test 0.5 * (Tev[1] + Tev[2]) ≈ params.Te_L
    @test 0.5 * (Tev[end] + Tev[end - 1]) ≈ params.Te_R

    # Check potential (also solved during setup)
    @test 0.5 * (ϕ[1] + ϕ[2]) ≈ V_d
    @test 0.5 * (ϕ[end - 1] + ϕ[end]) ≈ V_cc

    # Anode ion velocity less than bohm velocity
    prop = config.propellants[1]
    mi = prop.gas.m
    Ti = prop.ion_temperature_K
    un = prop.velocity_m_s
    mdot_a = prop.flow_rate_kg_s

    ni_1 = 1.0e17
    ni_2 = 1.0e16

    u_bohm_1 = sqrt((het.kB * Ti + e * params.Te_L) / mi)
    u_bohm_2 = sqrt((het.kB * Ti + 2 * e * params.Te_L) / mi)

    (; continuity, isothermal) = params.fluid_containers
    ingestion_flow_rate = params.ingestion_flow_rates[1]
    anode_bc = params.anode_bc
    prop = config.propellants[1]

    inlet_density = het.inlet_neutral_density(prop, config.thruster.geometry.channel_area)
    nn_B = het.background_neutral_density(prop, config)
    un_B = het.background_neutral_velocity(prop, config)

    # ====================================================================
    # Case 1:
    # When ion velocity at left boundary is less than then bohm speed,
    # it should be accelerated to reach the bohm speed
    # ====================================================================
    @. continuity[1].density = ni_1 * mi * 2
    @. isothermal[1].density = ni_1 * mi
    @. isothermal[1].momentum = -mi * ni_1 * u_bohm_1 / 2
    @. isothermal[2].density = ni_2 * mi
    @. isothermal[2].momentum = -mi * ni_2 * u_bohm_2 / 2

    het.apply_left_boundary!(params.fluid_containers, prop, params.cache, anode_bc, ingestion_flow_rate)
    het.apply_right_boundary!(params.fluid_containers)

    # Edge boundary state should equal average of ghost cell and first interior cell
    boundary_ion_flux = [(ion.momentum[1] + ion.momentum[2]) / 2 for ion in isothermal]
    boundary_ion_dens = [(ion.density[1] + ion.density[2]) / 2 for ion in isothermal]
    boundary_ion_vel = [flux / dens for (flux, dens) in zip(boundary_ion_flux, boundary_ion_dens)]

    neutral_dens_edge = 0.5 * (continuity[1].density[1] + continuity[1].density[2])
    recombined_ion_flux = -sum(boundary_ion_flux) / un
    background_flux = nn_B * un_B / un * config.neutral_ingestion_multiplier

    @test neutral_dens_edge ≈ inlet_density + recombined_ion_flux + background_flux
    @test boundary_ion_vel[1] ≈ -u_bohm_1
    @test boundary_ion_vel[2] ≈ -u_bohm_2

    # Right BC
    for fluid in isothermal
        interior_density = fluid.density[end - 1]
        interior_flux = fluid.momentum[end - 1]
        interior_velocity = interior_flux / interior_density

        if interior_velocity > 0
            fluid.density[end] = interior_density
            fluid.momentum[end] = interior_flux
        else
            fluid.density[end] = MIN_NUMBER_DENSITY * fluid.species.element.m
            fluid.momentum[end] = MIN_NUMBER_DENSITY * fluid.species.element.m * interior_velocity
        end
    end

    # ====================================================================
    # Case 2:
    # Anode ion velocity greater than Bohm speed
    # Ions should have a Neumann BC
    # ====================================================================
    @. continuity[1].density = ni_1 * mi * 2
    @. isothermal[1].density = ni_1 * mi
    @. isothermal[1].momentum = -mi * ni_1 * u_bohm_1 * 2
    @. isothermal[2].density = ni_2 * mi
    @. isothermal[2].momentum = -mi * ni_2 * u_bohm_2 * 2

    het.apply_left_boundary!(params.fluid_containers, prop, params.cache, anode_bc, ingestion_flow_rate)
    het.apply_right_boundary!(params.fluid_containers)

    # Edge boundary state should equal average of ghost cell and first interior cell
    boundary_ion_flux = [(ion.momentum[1] + ion.momentum[2]) / 2 for ion in isothermal]
    boundary_ion_dens = [(ion.density[1] + ion.density[2]) / 2 for ion in isothermal]
    boundary_ion_vel = [flux / dens for (flux, dens) in zip(boundary_ion_flux, boundary_ion_dens)]
    neutral_dens_edge = 0.5 * (continuity[1].density[1] + continuity[1].density[2])

    recombined_ion_flux = -sum(boundary_ion_flux) / un
    background_flux = nn_B * un_B / un * config.neutral_ingestion_multiplier

    @test neutral_dens_edge ≈ inlet_density + recombined_ion_flux + background_flux
    @test boundary_ion_vel[1] ≈ -2 * u_bohm_1
    @test boundary_ion_vel[2] ≈ -2 * u_bohm_2


    for fluid in isothermal
        interior_density = fluid.density[end - 1]
        interior_flux = fluid.momentum[end - 1]
        interior_velocity = interior_flux / interior_density

        if interior_velocity > 0
            fluid.density[end] = interior_density
            fluid.momentum[end] = interior_flux
        else
            fluid.density[end] = MIN_NUMBER_DENSITY * fluid.species.element.m
            fluid.momentum[end] = MIN_NUMBER_DENSITY * fluid.species.element.m * interior_velocity
        end
    end

    return
end

test_boundaries()
