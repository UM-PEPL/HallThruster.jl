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

    # Explicitly tracked neutral and ion levels are emitted separately from the
    # ground-state collections and retain their derived energies.
    excited_propellant = het.Propellant(
        het.Xenon, 1.0e-6;
        max_charge = 1,
        excited_levels = [1],
        excited_ion_levels = Dict(1 => [2]),
    )
    excited_fluids = [het.allocate_fluids(excited_propellant, 2)]
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

    emission = het.PhotonEmission(;
        upper = het.Xenon(0, 1).symbol,
        lower = het.Xenon(0).symbol,
        frequency = 2.0e7,
        energy_eV = 8.3,
        emission_rate = fill(4.0e18, 4),
    )
    test_roundtrip(het.PhotonEmission, het.serialize(emission))

    return
end

test_solution_serialization()
