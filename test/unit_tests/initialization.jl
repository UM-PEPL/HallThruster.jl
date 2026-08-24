using HallThruster: HallThruster as het

function restart_fixture()
    ncells = 3
    z = collect(0.0:(ncells + 1))
    xenon = het.Propellant(
        het.Xenon, 1.0e-6;
        max_charge = 1,
        excited_levels = [1],
        excited_ion_levels = Dict(1 => [1]),
    )
    krypton = het.Propellant(
        het.Krypton, 1.0e-6;
        max_charge = 1,
        excited_levels = [2],
    )
    propellants = [xenon, krypton]
    fluids_by_propellant = [het.allocate_fluids(prop, ncells) for prop in propellants]
    energies = het.OrderedDict(
        het.Xenon(0).symbol => 0.0,
        het.Xenon(0, 1).symbol => 8.3,
        het.Xenon(1).symbol => 12.1,
        het.Xenon(1, 1).symbol => 13.4,
        het.Krypton(0).symbol => 0.0,
        het.Krypton(0, 2).symbol => 10.2,
        het.Krypton(1).symbol => 14.0,
    )
    cache = (;
        ne = zeros(length(z)),
        nϵ = zeros(length(z)),
        Tev = zeros(length(z)),
        ∇ϕ = zeros(length(z)),
        ϕ = zeros(length(z)),
    )
    params = (;
        grid = (; cell_centers = z),
        cache,
        propellants,
        fluids_by_propellant,
        species_energies_eV = energies,
    )

    state(n, nu = n; kwargs...) = Dict(
        "n" => fill(Float64(n), length(z)),
        "nu" => fill(Float64(nu), length(z)),
        (string(key) => value for (key, value) in kwargs)...,
    )
    frame = Dict(
        "z" => z,
        "neutrals" => Dict(
            "Xe" => state(2.0),
            "Kr" => state(3.0),
        ),
        "ions" => Dict(
            "Xe" => [state(4.0, 5.0; Z = 1)],
            "Kr" => [state(6.0, 7.0; Z = 1)],
        ),
        # Xe(*) is present and validated. Xe(1+,*) and Kr(2*) are intentionally
        # absent to verify that newly configured levels start empty on restart.
        "excited_states" => Dict(
            string(het.Xenon(0, 1).symbol) => state(
                8.0, 9.0; Z = 0, excited_level = 1, energy_eV = 8.3,
            ),
        ),
        "ne" => fill(10.0, length(z)),
        "Tev" => fill(11.0, length(z)),
        "potential" => fill(12.0, length(z)),
        "E" => fill(13.0, length(z)),
    )
    return params, frame
end

@testset "Explicit-state restart handling" begin
    params, frame = restart_fixture()
    het.initialize_from_restart!(params, frame)

    xenon_fluids, krypton_fluids = params.fluids_by_propellant
    @test all(≈(2.0), xenon_fluids.continuity[1].density ./ het.Xenon.m)
    @test all(≈(8.0), xenon_fluids.continuity[2].density ./ het.Xenon.m)
    @test all(≈(4.0), xenon_fluids.isothermal[1].density ./ het.Xenon.m)
    @test all(≈(5.0), xenon_fluids.isothermal[1].momentum ./ het.Xenon.m)
    @test all(iszero, xenon_fluids.isothermal[2].density)
    @test all(iszero, xenon_fluids.isothermal[2].momentum)

    @test all(≈(3.0), krypton_fluids.continuity[1].density ./ het.Krypton.m)
    @test all(iszero, krypton_fluids.continuity[2].density)
    @test all(≈(6.0), krypton_fluids.isothermal[1].density ./ het.Krypton.m)
    @test all(≈(7.0), krypton_fluids.isothermal[1].momentum ./ het.Krypton.m)
    @test params.cache.ne == frame["ne"]
    @test params.cache.Tev == frame["Tev"]
    @test params.cache.ϕ == frame["potential"]
    @test params.cache.∇ϕ == .-frame["E"]

    # Missing required ground states and incompatible excited-state metadata
    # should identify the offending species rather than surface a KeyError.
    missing_ground = deepcopy(frame)
    delete!(missing_ground["neutrals"], "Kr")
    err = try
        het.initialize_from_restart!(first(restart_fixture()), missing_ground)
        nothing
    catch exception
        exception
    end
    @test err isa ArgumentError
    @test occursin("ground-state Kr neutral", sprint(showerror, err))

    wrong_energy = deepcopy(frame)
    wrong_energy["excited_states"][string(het.Xenon(0, 1).symbol)]["energy_eV"] = 9.0
    err = try
        het.initialize_from_restart!(first(restart_fixture()), wrong_energy)
        nothing
    catch exception
        exception
    end
    @test err isa ArgumentError
    @test occursin("active chemistry uses 8.3 eV", sprint(showerror, err))

    malformed = deepcopy(frame)
    malformed["ions"]["Xe"][1]["nu"] = [1.0, 2.0]
    err = try
        het.initialize_from_restart!(first(restart_fixture()), malformed)
        nothing
    catch exception
        exception
    end
    @test err isa ArgumentError
    @test occursin("expected 5", sprint(showerror, err))

    # Restart documents carry the same schema marker as ordinary JSON inputs;
    # future versions must be rejected before their frame layout is interpreted.
    future_restart = tempname() * ".json"
    open(future_restart, "w") do io
        het.JSON.write_json(
            io, Dict(
                "serialization_version" => het.SERIALIZATION_VERSION + 1,
                "output" => Dict("average" => frame),
            )
        )
    end
    err = try
        het.initialize_from_restart!(first(restart_fixture()), future_restart)
        nothing
    catch exception
        exception
    end
    @test err isa ArgumentError
    @test occursin("unsupported serialization version", sprint(showerror, err))
    rm(future_restart)
end
