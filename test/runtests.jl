using HallThruster
using Test
using Documenter
using Statistics
using DelimitedFiles
using LinearAlgebra
using Unitful
using SparseArrays

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

using Symbolics
include("order_verification/ovs_funcs.jl")
include("order_verification/ovs_energy.jl")
include("order_verification/ovs_ions.jl")

@testset "Order verification (electron energy)" begin
    vfunc = x -> OVS_Energy.verify_energy(x)
    refinements = refines(6, 20, 2)

    cases = ["implicit", "Crank-Nicholson"]
    norms_to_test = [1, 2, Inf]
    num_norms = length(norms_to_test)

    # Test spatial order of accuracy of implicit solver and crank-nicholson on L1, L2, and L∞ norms
    slopes_nϵ, norms_nϵ = test_refinements(vfunc, refinements, norms_to_test)

    # Check that electron energy equation is solved to at least first order in space
    for (i, slope) in enumerate(slopes_nϵ)
        norm_ind = mod1(i, num_norms)
        case_ind = ((i-1) ÷ num_norms) + 1
        println("Electron energy ($(cases[case_ind]), $(norms_to_test[norm_ind])-norm): ", slope)
        @test slope > 0.8
    end
end

@testset "Order verification (neutrals and ions)" begin
    refinements = refines(6, 10, 2)

    limiter = HallThruster.van_leer
    flux_names = ["HLLE", "Rusanov", "Global Lax-Friedrichs"]
    #flux_names = ["Global Lax-Friedrichs"]
    fluxes = [HallThruster.HLLE, HallThruster.rusanov, HallThruster.global_lax_friedrichs,]

    # Which L-P norms to check
    cases = ["ρn", "ρi", "ρiui"]
    norms_to_test = [1, 2, Inf]
    num_norms = length(norms_to_test)

    for (flux, flux_name) in zip(fluxes, flux_names)
        for reconstruct in [false, true]
            if reconstruct && flux_name == "HLLE"
                continue
            end
            scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct)
            orders, norms = test_refinements(ncells -> OVS_Ions.solve_ions(ncells, scheme), refinements, norms_to_test)
            for (i, (order, norm)) in enumerate(zip(orders, norms))

                norm_ind = mod1(i, num_norms)
                case_ind = ((i-1) ÷ num_norms) + 1
                case_str = "($(cases[case_ind]), $(norms_to_test[norm_ind])-norm)"

                if !reconstruct
                    println("No reconstruction, $flux_name ($case_str): ", order)
                else
                    println("With reconstruction, $flux_name ($case_str):  ", order)
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
