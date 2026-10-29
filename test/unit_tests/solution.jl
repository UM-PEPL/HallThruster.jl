using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

function test_solution_serialization()

    ncells = 50

    config = het.Config(;
        thruster = het.SPT_100,
        discharge_voltage = 300.0,
        domain = (0.0, 0.8),
        anode_mass_flow_rate = 5.0e-6,
    )

    simparams = het.SimParams(
        grid = het.UnevenGrid(ncells),
        duration = 1.0e-6,
        num_save = 1000,
        verbose = false
    )

    cfg_ser = het.serialize(config)

    sol = het.run_simulation(config, simparams)
    avg = het.time_average(sol)
    @test sol[:nu_anom] == sol[:nu_an]
    @test sol[:B] === sol.frames[1].B

    # The deprecated single-propellant accessors and Config-only entry point are gone.
    @test all(field -> field ∉ het.valid_fields(), (:nn, :ni, :ui, :niui))
    @test_throws ArgumentError sol[:nn]
    @test_throws MethodError sol[:ni, 1]
    @test !hasmethod(het.run_simulation, Tuple{het.Config})

    frame = avg.frames[1]
    neutral_state = frame.neutrals[:Xe]

    neu = het.serialize(neutral_state)
    test_roundtrip(het.SpeciesState, neu)

    ions = frame.ions[:Xe]
    @test typeof(ions) == Vector{het.SpeciesState}
    test_roundtrip(typeof(ions), ions)

    ion_dict = frame.ions
    @test typeof(ion_dict) == het.OrderedDict{Symbol, Vector{het.SpeciesState}}

    ion_dict_ser = het.serialize(frame.ions)
    @test typeof(ion_dict_ser) == het.OrderedDict{String, Vector{het.OrderedDict{String, Any}}}

    ion_dict_2 = het.deserialize(typeof(ion_dict), ion_dict_ser)
    for k in keys(ion_dict)
        for i in eachindex(ion_dict[k])
            @test struct_eq(ion_dict[k][i], ion_dict_2[k][i])
        end
    end

    test_roundtrip(het.Frame, frame)
    test_roundtrip(het.Solution, avg)

    @testset "Excited ions in performance metrics" begin
        metrics_sol = deepcopy(avg)
        metrics_frame = only(metrics_sol.frames)
        ground_ion = metrics_frame.ions[:Xe][1]
        metrics = (het.thrust, het.mass_eff, het.anode_eff, het.voltage_eff)
        baseline = [metric(metrics_sol, 1) for metric in metrics]

        excited_ion = het.SpeciesState(;
            n = 0.5 .* ground_ion.n,
            nu = 0.5 .* ground_ion.nu,
            u = copy(ground_ion.u),
            m = ground_ion.m,
            Z = ground_ion.Z,
            excited_level = 1,
            energy_eV = 13.4,
        )
        metrics_frame.excited_states[het.Xenon(1, 1).symbol] = excited_ion
        ground_ion.n .*= 0.5
        ground_ion.nu .*= 0.5

        # Excited neutrals must not contribute to ion performance metrics,
        # even when they have a substantial density and flux.
        metrics_frame.excited_states[het.Xenon(0, 1).symbol] = het.SpeciesState(;
            n = 100 .* excited_ion.n,
            nu = 100 .* excited_ion.nu,
            u = copy(excited_ion.u),
            m = excited_ion.m,
            Z = 0,
            excited_level = 1,
            energy_eV = 8.3,
        )
        @test [metric(metrics_sol, 1) for metric in metrics] ≈ baseline

        # The same total ion beam can reside entirely in an excited state.
        excited_ion.n .+= ground_ion.n
        excited_ion.nu .+= ground_ion.nu
        ground_ion.n .= 0.0
        ground_ion.nu .= 0.0
        ground_ion.u .= 0.0
        @test [metric(metrics_sol, 1) for metric in metrics] ≈ baseline
    end

    # Explicitly tracked neutral and ion levels are emitted separately from the
    # ground-state collections and retain their derived energies.
    excited_propellant = het.Propellant(
        het.Xenon, 1.0e-6;
        max_charge = 1,
        excited_levels = [1],
        excited_ion_levels = Dict(1 => [2]),
    )
    excited_grid = het.Grid1D(range(0.0, 1.0; length = 3))
    excited_fluids = [het.allocate_fluids(excited_propellant, excited_grid)]
    energies = het.OrderedDict(
        het.Xenon(0).symbol => 0.0,
        het.Xenon(0, 1).symbol => 8.3,
        het.Xenon(1).symbol => 12.1,
        het.Xenon(1, 2).symbol => 13.4,
    )
    neutrals, ions, excited_states = het._get_species_states(
        excited_fluids, energies,
    )

    @test collect(keys(neutrals)) == [:Xe]
    @test length(ions[:Xe]) == 1
    @test collect(keys(excited_states)) == [
        het.Xenon(0, 1).symbol, het.Xenon(1, 2).symbol,
    ]
    @test excited_states[het.Xenon(0, 1).symbol].energy_eV == 8.3
    @test excited_states[het.Xenon(1, 2).symbol].energy_eV == 13.4
    @test excited_states[het.Xenon(1, 2).symbol].excited_level == 2

    profile_grid = het.Grid1D([0.0, 0.04, 0.08])
    velocity_profile = het.LinearInterpolation([0.0, 0.08], [150.0, 250.0])
    profiled_propellant = het.Propellant(
        het.Xenon, 1.0e-6; velocity_m_s = velocity_profile,
    )
    profiled_fluids = [het.allocate_fluids(profiled_propellant, profile_grid)]
    profiled_neutrals, _, _ = het._get_species_states(
        profiled_fluids, het.OrderedDict(het.Xenon(0).symbol => 0.0, het.Xenon(1).symbol => 12.1),
    )
    output_grid = het.copy_and_remove_ghosts(profile_grid.cell_centers)
    @test profiled_neutrals[:Xe].u ≈ velocity_profile.(output_grid)
    @test profiled_neutrals[:Xe].nu ≈
        profiled_neutrals[:Xe].n .* profiled_neutrals[:Xe].u

    emission = het.PhotonEmission(;
        upper = het.Xenon(0, 1).symbol,
        lower = het.Xenon(0).symbol,
        frequency = 2.0e7,
        energy_eV = 8.3,
        emission_rate = fill(4.0e18, 4),
    )
    test_roundtrip(het.PhotonEmission, het.serialize(emission))

    emission_sol = deepcopy(sol[[1, end]])
    for (i, frame) in enumerate(emission_sol.frames)
        push!(frame.photon_emissions, deepcopy(emission))
        frame.photon_emissions[1].emission_rate .= i == 1 ? 0.0 : 4.0e18
    end
    @test only(het.time_average(emission_sol).frames).photon_emissions[1].emission_rate == fill(4.0e18, 4)
    @test only(het.time_average(emission_sol[2]).frames).photon_emissions[1].emission_rate == fill(4.0e18, 4)

    return
end

test_solution_serialization()
