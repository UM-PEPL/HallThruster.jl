using Test, Documenter, HallThruster, StaticArrays, BenchmarkTools, Symbolics, Statistics, LinearAlgebra

doctest(HallThruster)

@testset "Gas and species tests" begin
    include("unit_tests/gas_and_species.jl")
end

@testset "Conservation law systems, fluids, and fluxes" begin
    include("unit_tests/conservation_laws.jl")
end

@testset "Flux limiters" begin
    include("unit_tests/flux_limiters.jl")
end

@testset "Ionization" begin
    include("unit_tests/ionization.jl")
end

@testset "Utility functions" begin
    include("unit_tests/utility_funcs.jl")
end

@testset "Electron transport" begin
    include("unit_tests/electron_transport.jl")
end

@testset "Order verification (potential and gradients)" begin
    include("order_verification/ovs_funcs.jl")
    include("order_verification/ovs_potential.jl")
    refinements = refines(4, 10, 2)

    # check that first-order convergence is maintained for L1, L2, and L∞ norms
    for p in (1, 2, Inf)
        (slope_ϕ,), norms_ϕ =  test_refinements(OVS_Potential.verify_potential, refinements, p)
        (slope_∇ϕ, slope_∇pe, slope_ue), norms_grad =  test_refinements(OVS_Potential.verify_gradients, refinements, p)

        tol = 0.15

        # Check that potential and gradients are first order or better
        @test abs(slope_ϕ - 1.0) < tol || slope_ϕ > 1.0
        # Check that potential gradient and ue are first order or better
        @test abs(slope_∇ϕ - 1.0) < tol || slope_∇ϕ > 1.0
        @test abs(slope_ue - 1.0) < tol || slope_ue > 1.0
        # Check that pressure gradient is second order or better
        @test abs(slope_∇pe - 2.0) < tol || slope_∇pe > 2.0
    end
end

@testset "Order verification (electron energy)" begin
    include("order_verification/ovs_funcs.jl")
    include("order_verification/ovs_energy.jl")
    refinements = refines(6, 20, 2)

    # Test spatial order of accuracy of implicit solver and crank-nicholson on L1, L2, and L∞ norms
    slopes_nϵ, norms_nϵ = test_refinements(OVS_Energy.verify_energy, refinements, [1, 2, Inf])

    # Check that electron energy is solved to at least first order
    for slope in slopes_nϵ
        @test abs(slope - 2.0) < 0.2 || slope > 2.0
    end
end

@testset "Order verification (neutrals and ions)" begin
    include("order_verification/ovs_funcs.jl")
    include("order_verification/ovs_ions.jl")
    refinements = refines(5, 40, 2)

    limiter = HallThruster.minmod
    flux_names = ("HLLE", "Rusanov")
    fluxes = (HallThruster.HLLE, HallThruster.rusanov)

    for (flux, flux_name) in zip(fluxes, flux_names)
        for reconstruct in (false, )
            scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct)

            slopes, norms = test_refinements(ncells -> OVS_Ions.solve_ions(ncells, scheme, false), refinements, [1, 2, Inf])

            theoretical_order = 1 + reconstruct

            for slope in slopes
                @test abs(slope - theoretical_order) < 0.2
            end
        end
    end
end

#=
@testset "Simulation setup tests" begin
    @test SPT_100 isa HallThruster.Geometry1D
    @test HallThruster.channel_area(SPT_100) == π * (0.05^2 - 0.0345^2)

    species = [
        HallThruster.Species(HallThruster.Xenon, 0),
        HallThruster.Species(HallThruster.Xenon, 1),
        HallThruster.Species(HallThruster.Xenon, 2),
        HallThruster.Species(HallThruster.Xenon, 3),
    ]

    @test HallThruster.get_species(simulation) == species

    _, fluids, fluid_ranges, species_range_dict = HallThruster.configure_simulation(simulation)

    @test fluids == [
        HallThruster.Fluid(species[1], HallThruster.ContinuityOnly(u = 300.0, T = 500.0)),
        HallThruster.Fluid(species[2], HallThruster.IsothermalEuler(T = 500.0)),
        HallThruster.Fluid(species[3], HallThruster.IsothermalEuler(T = 500.0)),
        HallThruster.Fluid(species[4], HallThruster.IsothermalEuler(T = 500.0)),
    ]

    @test fluid_ranges == [1:1, 2:3, 4:5, 6:7]

    @test species_range_dict == Dict{HallThruster.Species, UnitRange{Int64}}(
        species[1] => fluid_ranges[1],
        species[2] => fluid_ranges[2],
        species[3] => fluid_ranges[3],
        species[4] => fluid_ranges[4]
    )

    z_cell, z_edge = HallThruster.generate_grid(SPT_100, simulation.ncells)
    @test z_cell[1] == z_edge[1] && z_cell[end] == z_edge[end]
    @test z_cell[2] == 0.5 * (z_edge[2] + z_edge[1])
    @test z_edge[2] - z_edge[1] == (SPT_100.domain[2] - SPT_100.domain[1]) / simulation.ncells
    @test z_cell[3] - z_cell[2] == (SPT_100.domain[2] - SPT_100.domain[1]) / simulation.ncells
    
    U, (F, UL, UR, Q) = HallThruster.allocate_arrays(simulation)
    @test size(U, 1) == size(F, 1) == size(UL, 1) == size(UR, 1) == size(Q, 1)
    nvariables = size(U, 1)
    @test nvariables == 1 + 6 + 3
    
    @test size(U, 2) == simulation.ncells+2
    @test size(UL, 2) == size(UR, 2) == size(F, 2) == simulation.ncells+1
    
    mdot = 5e-6 # kg/s
    un = 300 # m/s
    A = π * (0.05^2 - 0.0345^2) # m^2
    m_atom = HallThruster.Xenon.M / HallThruster.NA

    @test m_atom == HallThruster.Xenon.m
    nn = mdot / un / A / m_atom

    @test nn == HallThruster.inlet_neutral_density(simulation)
    
    HallThruster.initial_condition!(U, z_cell, simulation, fluid_ranges)

    @test U[end, :] == ϕ_func.(z_cell)
    @test U[end-1, :] == 6 .* ni_func.(z_cell)
    @test U[end-2, :] == Te_func.(z_cell)
    
    @show maximum(U[end-2, :])

    @test all(U[1, :] .== nn)
    @test U[2, :] == ni_func.(z_cell)
    @test U[3, :] == un .* ni_func.(z_cell)
    @test U[4, :] == ni_func.(z_cell)
    @test U[5, :] == un .* ni_func.(z_cell)
    @test U[6, :] == ni_func.(z_cell)
    @test U[7, :] == un .* ni_func.(z_cell)

    cache = (F, UL, UR, Q)

    scheme = simulation.scheme

    reactions = HallThruster.load_ionization_reactions(species)

    params = (;
        cache,
        fluids,
        fluid_ranges,
        species_range_dict,
        z_cell,
        z_edge,
        reactions,
        scheme
    )
    #dU = zeros(size(U))
    #@time HallThruster.update!(dU, U, params, 0.0)
end

# TODO: using any of the SSP methods, this fails sometimes and succeeds others, in a way that seems independent of CFL number
@testset "Freestream preservation" begin
    include("freestream_preservation.jl")
    test_preservation(0.9)
end
=#
######################################
#computations for MMS OVS