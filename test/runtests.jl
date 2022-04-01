using Test, Documenter, HallThruster, StaticArrays, BenchmarkTools, Symbolics, Statistics, LinearAlgebra

doctest(HallThruster)

include("unit_tests/test_gas.jl")
include("unit_tests/test_conservation_laws.jl")
include("unit_tests/test_limiters.jl")
include("unit_tests/test_ionization.jl")
include("unit_tests/test_misc.jl")
include("unit_tests/test_electrons.jl")
include("unit_tests/test_boundary_conditions.jl")
include("unit_tests/test_walls.jl")

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
    WENO_names = ("WENO off", "WENO on")
    WENOs = (false, true)

    for (flux, flux_name) in zip(fluxes, flux_names)
        for (WENO, WENO_name) in zip(WENOs, WENO_names)
            for reconstruct in (false, )
                scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct, WENO)

                slopes, norms = test_refinements(ncells -> OVS_Ions.solve_ions(ncells, scheme, false), refinements, [1, 2, Inf])

                theoretical_order = 1 + reconstruct

                for slope in slopes
                    if WENO
                        println("WENO slope: ", slope)
                    else
                        @test abs(slope - theoretical_order) < 0.2
                    end
                end
            end
        end
    end
end
