using HallThruster
using Test
using Documenter
using Symbolics
using Statistics
using DelimitedFiles
using LinearAlgebra
using Unitful
using SparseArrays

doctest(HallThruster)

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

@testset "Order verification (electron energy)" begin
    include("order_verification/ovs_funcs.jl")
    include("order_verification/ovs_energy.jl")

    vfunc = x -> OVS_Energy.verify_energy(x)
    refinements = refines(6, 20, 2)

    # Test spatial order of accuracy of implicit solver and crank-nicholson on L1, L2, and L∞ norms
    slopes_nϵ, norms_nϵ = test_refinements(vfunc, refinements, [1, 2, Inf])

    tol = 0.15
    order = 1.0

    # Check that electron energy is solved to at least first order
    for slope in slopes_nϵ
        @show slope
        @test abs(slope - order) < tol || slope ≥ order
    end
end

@testset "Order verification (neutrals and ions)" begin
    include("order_verification/ovs_funcs.jl")
    include("order_verification/ovs_ions.jl")
    refinements = refines(5, 40, 2)

    limiter = HallThruster.van_leer
    flux_names = ("HLLE", "Rusanov", "Global Lax-Friedrichs")
    fluxes = (HallThruster.HLLE, HallThruster.rusanov, HallThruster.global_lax_friedrichs,)

    for (flux, flux_name) in zip(fluxes, flux_names)
        for reconstruct in (false, true)

            if reconstruct && flux_name == "HLLE"
                continue
            end

            scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct)

            slopes, norms = test_refinements(ncells -> OVS_Ions.solve_ions(ncells, scheme, false), refinements, [1, 2, Inf])

            # Check that gradient reconstruction gets us to at least ~1.75-order accuracy
            theoretical_order = 1 + reconstruct
            for slope in slopes
                if !reconstruct
                    println("No reconstruction, $flux_name: ", slope)
                else
                    println("With reconstruction, $flux_name: ", slope)
                end

                @test abs(slope - theoretical_order) < 0.25 || slope > theoretical_order
            end
        end
    end
end
