using HallThruster: HallThruster as het

@testset "Competing destruction channels share a timestep budget" begin
    limits = Float64[]
    for positive_states in (1, 3)
        grid = het.Grid1D([0.0, 0.04, 0.08])
        propellant = het.Propellant(
            het.Xenon, 5.0e-6; allowed_charges = [-1, 1],
            excited_ion_levels = Dict(1 => collect(1:(positive_states - 1))),
        )
        group = het.allocate_fluids(propellant, grid)
        fluids = [group.continuity; group.isothermal]
        neutral = het.ground_neutral(group)
        negative = first(group.isothermal)
        positives = group.isothermal[2:end]
        neutral.density .= propellant.gas.m * 4.0e17
        negative.density .= propellant.gas.m * 1.0e15
        for positive in positives
            positive.density .= propellant.gas.m * 3.0e17 / positive_states
        end
        cache = het.allocate_arrays(grid.num_cells, 0)
        ne = 3.0e17 - 1.0e15
        cache.nϵ .= 10.0 * ne
        k_ei = 2.0e-13
        reactions = [
            het.ElectronImpactReaction(
                0.7, propellant.gas(-1), [propellant.gas(0)], fill(k_ei, 256),
            ),
        ]
        reaction_groups = het.build_electron_impact_groups(
            reactions, het.reactant_indices(reactions, fluids),
            het.product_indices(reactions, fluids), fluids,
        )
        cache.reaction_rate_index_limit[] = het.common_rate_index_limit(reaction_groups)
        params = (; cache, fluids_by_propellant = [group])
        ei_loss = k_ei * ne
        mn_loss = het._K_MUTUAL_NEUTRALIZATION * 3.0e17
        ad_loss = het._K_NEUTRAL_DETACHMENT * 4.0e17

        # Repeating the stage must reset its budget rather than accumulate it.
        for _ in 1:2
            for fluid in fluids
                fill!(fluid.dens_ddt, 0.0)
                fill!(fluid.mom_ddt, 0.0)
            end
            het.apply_reaction_groups!(fluids, reaction_groups, cache, false, 1.0, Dict())
            @test cache.dt_iz[] ≈ inv(ei_loss)
            het.apply_mutual_neutralization!(params)
            @test cache.dt_iz[] ≈ inv(ei_loss + mn_loss)
            het.apply_associative_detachment!(params)
            @test cache.dt_iz[] ≈ inv(ei_loss + mn_loss + ad_loss)
            @test negative.loss_frequency[2] ≈ ei_loss + mn_loss + ad_loss
            @test negative.loss_frequency[[1, end]] == [0.0, 0.0]
        end
        push!(limits, cache.dt_iz[])

        # Wall losses also compete with volume reactions for the same ions.
        wall_frequency = 1.0e6
        wall_base = fill(wall_frequency * sqrt(propellant.gas.m), grid.num_cells)
        het.apply_ion_wall_losses!(group, wall_base, het.interior_cells(neutral.density))
        het.update_reaction_timestep!(cache, fluids)
        @test cache.dt_iz[] ≈ inv(wall_frequency + het._K_MUTUAL_NEUTRALIZATION * 1.0e15)
    end
    # Partitioning a fixed ion population into excited levels changes no rates.
    @test limits[1] ≈ limits[2]
end

@testset "Mutual neutralization charge conservation" begin
    grid = het.Grid1D([0.0, 0.04, 0.08])
    cases = (
        ([-1], [1]),
        ([-1], [2]),
        ([-2], [1]),
        ([-2], [2]),
        ([-2, -1], [1, 2, 3]),
    )

    for (negative_charges, positive_charges) in cases,
            same_gas in (false, true), excited in (false, true)
        @testset "$(negative_charges), $(positive_charges), same gas=$(same_gas), excited=$(excited)" begin
            make_propellant(gas, charges) = het.Propellant(
                gas, 0.0; allowed_charges = charges,
                excited_ion_levels = Dict(Z => [1] for Z in charges if excited),
            )
            propellants = if same_gas
                [make_propellant(het.Xenon, [negative_charges; positive_charges])]
            else
                [
                    make_propellant(het.Xenon, negative_charges),
                    make_propellant(het.Krypton, positive_charges),
                ]
            end
            groups = [het.allocate_fluids(propellant, grid) for propellant in propellants]
            ions = [fluid for group in groups for fluid in group.isothermal]
            fluids = [fluid for group in groups for fluid in [group.continuity; group.isothermal]]

            for ion in ions
                n = 1.0e16 * (1 + abs(ion.species.Z) + ion.species.excited_level)
                ion.density .= n * ion.species.element.m
                ion.momentum .= ion.density .* (ion.species.Z < 0 ? -200.0 : 300.0)
            end
            negative_density = sum(
                ion.density[2] / ion.species.element.m for ion in ions
                    if ion.species.Z == -1; init = 0.0,
            )
            positive_density = sum(
                ion.density[2] / ion.species.element.m for ion in ions
                    if ion.species.Z == 1; init = 0.0,
            )
            cache = (; dt_iz = fill(Inf))
            params = (; fluids_by_propellant = groups, cache)
            het.apply_mutual_neutralization!(params)

            for ion in ions
                partner_density = if ion.species.Z == -1
                    positive_density
                elseif ion.species.Z == 1
                    negative_density
                else
                    0.0
                end
                expected_loss = het._K_MUTUAL_NEUTRALIZATION .* partner_density .* ion.density
                expected_loss[[1, end]] .= 0.0
                @test ion.dens_ddt ≈ -expected_loss
                if abs(ion.species.Z) > 1
                    @test all(iszero, ion.mom_ddt)
                end
            end

            for cell in eachindex(grid.cell_centers)
                charge_sources = [
                    ion.species.Z * ion.dens_ddt[cell] / ion.species.element.m
                        for ion in ions
                ]
                mass_sources = [fluid.dens_ddt[cell] for fluid in fluids]
                @test abs(sum(charge_sources)) <= 1.0e-12 * sum(abs, charge_sources)
                @test abs(sum(mass_sources)) <= 1.0e-12 * sum(abs, mass_sources)
            end

            # Each parent gas gains exactly the mass removed from its ions.
            for group in groups
                expected_neutral_source = zeros(grid.num_cells)
                for ion in group.isothermal
                    expected_neutral_source .-= ion.dens_ddt
                end
                @test het.ground_neutral(group).dens_ddt ≈ expected_neutral_source
            end
            if negative_density == 0 || positive_density == 0
                @test cache.dt_iz[] == Inf
            else
                @test isfinite(cache.dt_iz[])
            end
        end
    end
end

@testset "Primitive velocity cache" begin
    # Primitive velocity is derived once per RK stage and stored with its fluid;
    # density and momentum remain the authoritative conservative state.
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 1)
    grid = het.Grid1D(range(0.0, 1.0; length = 4))
    ion = only(het.allocate_fluids(propellant, grid).isothermal)
    ion.density .= [0.0, 2.0, 4.0, 5.0, 0.0]
    ion.momentum .= [1.0, 6.0, -8.0, 0.0, -1.0]

    @test het.update_primitive_velocity!(ion) === nothing
    @test ion.vel_prim == [0.0, 3.0, -2.0, 0.0, 0.0]
end

@testset "Profiled neutral flux" begin
    grid = het.Grid1D([0.0, 0.04, 0.08])
    propellant = het.Propellant(
        het.Xenon, 5.0e-6;
        velocity_m_s = het.LinearInterpolation([0.0, 0.08], [100.0, 300.0]),
        temperature_K = het.LinearInterpolation([0.0, 0.08], [400.0, 800.0]),
    )
    neutral = only(het.allocate_fluids(propellant, grid).continuity)
    neutral.dens_L .= [1.0, 2.0, 3.0]
    neutral.dens_R .= [4.0, 5.0, 6.0]

    het.compute_fluxes_continuity!(neutral, grid)

    expected_flux = @. 0.5 * (
        neutral.vel_L * (neutral.dens_L + neutral.dens_R) -
            neutral.wave_speed * (neutral.dens_R - neutral.dens_L)
    )
    @test neutral.flux_dens ≈ expected_flux
    @test neutral.max_timestep[] ≈ minimum(grid.dz_edge ./ neutral.wave_speed)
end

@testset "Isothermal ion flux" begin
    grid = het.Grid1D([0.0, 0.4, 1.0])
    temperature = 900.0
    propellant = het.Propellant(
        het.Xenon, 0.0; max_charge = 1, ion_temperature_K = temperature,
    )
    ion = only(het.allocate_fluids(propellant, grid).isothermal)
    RT = het.R0 / het.Xenon.M * temperature

    # For p = ρRT, both the pressure flux and the eigenspeeds must use the
    # isothermal sound speed sqrt(RT), without an adiabatic gamma factor.
    @test ion.sound_speed ≈ sqrt(RT)

    ion.dens_L .= [1.0, 2.0, 3.0]
    ion.dens_R .= [1.5, 2.5, 3.5]
    ion.vel_L .= [-2.0, 1.0, 4.0]
    ion.vel_R .= [-1.0, 3.0, 2.0]
    het.compute_fluxes_isothermal!(ion, grid)

    smax = @. max(abs(ion.vel_L), abs(ion.vel_R)) + sqrt(RT)
    expected_density_flux = @. 0.5 * (
        ion.dens_L * ion.vel_L + ion.dens_R * ion.vel_R -
            smax * (ion.dens_R - ion.dens_L)
    )
    expected_momentum_flux = @. 0.5 * (
        ion.dens_L * (ion.vel_L^2 + RT) +
            ion.dens_R * (ion.vel_R^2 + RT) -
            smax * (ion.dens_R * ion.vel_R - ion.dens_L * ion.vel_L)
    )

    @test ion.flux_dens ≈ expected_density_flux
    @test ion.flux_mom ≈ expected_momentum_flux
    @test ion.max_timestep[] ≈ minimum(grid.dz_edge ./ smax)
end

@testset "Zero-density heavy species" begin
    @test het.primitive_velocity(0.0, 0.0) == 0.0
    @test isfinite(het.primitive_velocity(0.0, 0.0))
    @test het.primitive_velocity(6.0, 2.0) == 3.0

    ncells = 3
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 1)
    grid = het.Grid1D(range(0.0, 1.0; length = ncells + 1))
    fluids = het.allocate_fluids(propellant, grid)
    fluid_arr = [fluids.continuity[1], fluids.isothermal[1]]

    rxn = het.ElectronImpactReaction(0.0, het.Xenon(0), [het.Xenon(1)], ones(256))
    groups = het.build_electron_impact_groups([rxn], [1], [[2]], fluid_arr)
    num_grid_cells = ncells + 2
    fluid_arr[2].density .= propellant.gas.m * 1.0e18
    cache = (;
        inelastic_losses = zeros(num_grid_cells),
        νiz = zeros(num_grid_cells),
        νex_explicit = zeros(num_grid_cells),
        ϵ = zeros(num_grid_cells),
        ne = zeros(num_grid_cells),
        K = zeros(num_grid_cells),
        nϵ = fill(1.0e19, num_grid_cells),
        cell_cache_1 = zeros(num_grid_cells),
        reaction_loss_frequency = zeros(num_grid_cells),
        reaction_rate_indices = zeros(Int, num_grid_cells),
        reaction_rate_fractions = zeros(num_grid_cells),
        reaction_rate_index_limit = [254],
        dt_iz = [Inf],
    )

    het.apply_reaction_groups!(fluid_arr, groups, cache, false, 1.0, Dict())

    @test cache.dt_iz[] == Inf
    @test all(isfinite, cache.νiz)
    @test all(isfinite, cache.inelastic_losses)
    @test all(iszero, cache.νiz)
    @test all(iszero, cache.νex_explicit)
    @test all(iszero, cache.inelastic_losses)
    @test all(iszero, fluid_arr[1].dens_ddt)
    @test all(iszero, fluid_arr[2].dens_ddt)

    # A populated reactant produces the expected grouped chemistry timestep.
    fluid_arr[1].density[2:(end - 1)] .= propellant.gas.m
    het.apply_reaction_groups!(fluid_arr, groups, cache, false, 1.0, Dict())
    @test cache.dt_iz[] ≈ 1.0e-18
end

@testset "Ion production cost multiplier" begin
    ncells = 3
    propellant = het.Propellant(
        het.Xenon, 0.0; max_charge = 2, excited_levels = [1],
    )
    grid = het.Grid1D(range(0.0, 1.0; length = ncells + 1))
    fluid_set = het.allocate_fluids(propellant, grid)
    fluids = [fluid_set.continuity; fluid_set.isothermal]
    mass = propellant.gas.m

    ground_neutral = findfirst(fluid -> fluid.species == het.Xenon(0), fluids)
    excited_neutral = findfirst(fluid -> fluid.species == het.Xenon(0, 1), fluids)
    singly_charged = findfirst(fluid -> fluid.species == het.Xenon(1), fluids)
    doubly_charged = findfirst(fluid -> fluid.species == het.Xenon(2), fluids)

    fluids[ground_neutral].density .= 2.0 * mass
    fluids[singly_charged].density .= 3.0 * mass
    fluids[doubly_charged].density .= mass

    reactions = [
        het.ElectronImpactReaction(
            2.0, het.Xenon(0), [het.Xenon(1)], ones(256),
        ),
        het.ElectronImpactReaction(
            3.0, het.Xenon(0), [het.Xenon(0, 1)], ones(256),
        ),
        het.ElectronImpactReaction(
            5.0, het.Xenon(1), [het.Xenon(2)], ones(256),
        ),
    ]
    groups = het.build_electron_impact_groups(
        reactions, het.reactant_indices(reactions, fluids),
        het.product_indices(reactions, fluids), fluids,
    )

    num_grid_cells = ncells + 2
    cache = (;
        inelastic_losses = zeros(num_grid_cells),
        νiz = zeros(num_grid_cells),
        νex_explicit = zeros(num_grid_cells),
        ϵ = zeros(num_grid_cells),
        ne = zeros(num_grid_cells),
        K = zeros(num_grid_cells),
        nϵ = fill(50.0, num_grid_cells),
        cell_cache_1 = zeros(num_grid_cells),
        reaction_loss_frequency = zeros(num_grid_cells),
        reaction_rate_indices = zeros(Int, num_grid_cells),
        reaction_rate_fractions = zeros(num_grid_cells),
        reaction_rate_index_limit = [254],
        dt_iz = [Inf],
    )

    # ne = 5: ground -> Xe+ costs 5 * 2 * 2, excitation costs 5 * 2 * 3,
    # and Xe+ -> Xe2+ costs 5 * 3 * 5. Only the first term is multiplied.
    interior = 2:(num_grid_cells - 1)
    multipliers = Dict(:Xe => 4.0, :Kr => 100.0)
    het.apply_reaction_groups!(fluids, groups, cache, false, 1.0, multipliers)
    @test all(cache.inelastic_losses[interior] .== 185.0)

    # An absent Xe entry falls back to one; another species' value is not used.
    het.apply_reaction_groups!(
        fluids, groups, cache, false, 1.0, Dict(:Kr => 100.0),
    )
    @test all(cache.inelastic_losses[interior] .== 125.0)
end

@testset "Empty heavy-species populations" begin
    ncells = 3
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 1)
    grid = het.Grid1D(range(0.0, 1.0; length = ncells + 1))
    fluids = het.allocate_fluids(propellant, grid)
    cache = het.allocate_arrays(length(first(fluids.continuity).density), 0)

    # User-provided initial conditions and restarts may contain an empty cell
    # before the normal density limiter runs; derived plasma fields must stay finite.
    het.update_heavy_species_cache!(fluids, cache, false, 3.0)
    @test all(isfinite, cache.avg_neutral_vel)
    @test all(isfinite, cache.avg_ion_vel)
    @test all(isfinite, cache.m_eff)
    @test all(isfinite, cache.Z_eff)
    @test all(iszero, cache.avg_neutral_vel)
    @test all(iszero, cache.avg_ion_vel)
    @test all(==(propellant.gas.m), cache.m_eff)
    @test all(==(1.0), cache.Z_eff)
    @test all(==(3.0), cache.ne)

    # A configuration with no charged fluids cannot provide the ion properties
    # required by the electron and wall models, so fail with a useful error.
    neutral_only = (; continuity = fluids.continuity, isothermal = typeof(fluids.isothermal)())
    @test_throws ArgumentError het.update_heavy_species_cache!(
        neutral_only, cache, false, 3.0,
    )
end

@testset "Ion acceleration timestep" begin
    ncells = 3
    propellant = het.Propellant(het.Xenon, 0.0, max_charge = 3)
    fluid_grid = het.Grid1D(range(0.0, 1.0; length = ncells + 1))
    fluids = het.allocate_fluids(propellant, fluid_grid).isothermal
    grid = (; dz_cell = [1.0, 0.4, 0.7, 0.5, 1.0])
    cache = (; ∇ϕ = [0.0, -2.0, 4.0, -8.0, 0.0], dt_E = fill(0.0))

    for (j, fluid) in enumerate(fluids)
        fluid.density .= j
    end

    het.apply_ion_acceleration!(fluids, grid, cache)

    expected_dt = Inf
    for fluid in fluids
        qe_m = fluid.species.Z * het.e / fluid.species.element.m
        for i in 2:(ncells + 1)
            qE_m = -qe_m * cache.∇ϕ[i]
            expected_dt = min(expected_dt, abs(grid.dz_cell[i] / qE_m))
            @test fluid.mom_ddt[i] == qE_m * fluid.density[i]
        end
    end

    @test cache.dt_E[] == sqrt(expected_dt)
end
