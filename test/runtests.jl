using HallThruster
using Test
using Documenter
using Symbolics
using Statistics
using DelimitedFiles
using LinearAlgebra
using Unitful
using SparseArrays
using Alert

doctest(HallThruster)

HallThruster.example_simulation(;ncells=20, duration=1e-7, dt=1e-8, nsave=2)

include("unit_tests/test_gas.jl")
include("unit_tests/test_conservation_laws.jl")
include("unit_tests/test_limiters.jl")
include("unit_tests/test_reactions.jl")
include("unit_tests/test_misc.jl")
include("unit_tests/test_electrons.jl")
include("unit_tests/test_boundary_conditions.jl")
include("unit_tests/test_walls.jl")
include("unit_tests/test_initialization.jl")
include("unit_tests/test_restarts.jl")
include("unit_tests/test_json.jl")

@testset "Order verification (electron energy)" begin
    include("order_verification/ovs_funcs.jl")
    include("order_verification/ovs_energy.jl")

    vfunc = x -> OVS_Energy.verify_energy(x; make_plots = false)
    refinements = refines(6, 20, 2)

    # Test spatial order of accuracy of implicit solver and crank-nicholson on L1, L2, and L∞ norms
    slopes_nϵ, norms_nϵ = test_refinements(vfunc, refinements, [1, 2, Inf])

    # Check that electron energy equation is solved to at least first order in space
    for slope in slopes_nϵ
        @show slope
        @test slope > 0.8
    end
end

@testset "Order verification (neutrals and ions)" begin
    include("order_verification/ovs_funcs.jl")
    include("order_verification/ovs_ions.jl")
    refinements = refines(8, 10, 1.43)

    limiter = HallThruster.van_leer
    flux_names = ["HLLE", "Rusanov", "Global Lax-Friedrichs"]
    fluxes = [HallThruster.HLLE, HallThruster.rusanov, HallThruster.global_lax_friedrichs,]

    for (flux, flux_name) in zip(fluxes, flux_names)
        for reconstruct in [false, true]
            if reconstruct && flux_name == "HLLE"
                continue
            end

            scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct)

            # Which L-P norms to check
            norms_to_check = [1, 2, Inf]
            orders, norms = test_refinements(ncells -> OVS_Ions.solve_ions(ncells, scheme; make_plots = false), refinements, norms_to_check)

            for order in orders
                if !reconstruct
                    println("No reconstruction, $flux_name: ", order)
                else
                    println("With reconstruction, $flux_name: ", order)
                end

                # Check that we achieve the desired order of accuracy
                if (reconstruct)
                    @test order > 1.5
                else
                    @test order > 0.85
                end
            end
        end
    end
end
