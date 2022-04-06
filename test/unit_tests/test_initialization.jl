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

    @test ui[1][1] ≈ -sqrt(2/3 * HallThruster.e * config.anode_Te / mi)
    @test ui[2][1] ≈ -sqrt(2/3 * 2 * HallThruster.e * config.anode_Te / mi)
    @test ui[3][1] ≈ -sqrt(2/3 * 3 * HallThruster.e * config.anode_Te / mi)

    @test abs(ϵ[1] - config.anode_Te) < 0.1
    @test abs(ϵ[end] - config.cathode_Te) < 0.1
    
    z_cell = params.z_cell
    xlabel = "z (cm)"
    p1 = plot(z_cell * 100, U[index.ρn, :] ./ mi; xlabel, label = "", title = "Neutral number density (m⁻³)")
    p2 = plot(;xlabel, yaxis = :log, title = "Ion and electron density (m⁻³)")
    p3 = plot(;xlabel, title = "Ion velocity (km/s)")
    p4 = plot(z_cell * 100, ϵ; xlabel, title = "Electron energy (3/2 Te, eV)", label = "")

    for Z in 1:config.ncharge
        plot!(p2, z_cell * 100, U[index.ρi[Z], :] / mi, label = "Z = $Z")
        plot!(p3, z_cell * 100, ui[Z] / 1000, label = "")
    end

    plot!(p2, z_cell * 100, ne, label = "Electrons", lw = 2)

    plots = (p1, p2, p3, p4)

    foreach(plots) do p
        vline!(p, [config.thruster.geometry.channel_length * 100], lw = 2, lc = :red, ls = :dash, label = "")
    end

    p = plot(p1, p2, p3, p4, layout = (2, 2), size = (1000, 600))
    savefig(p, "docs/src/assets/intialization.svg")
    display(p)
end

