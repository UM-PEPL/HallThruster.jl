using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/order_verification/ovs_funcs.jl")
#=
@testset "Order verification (electron energy)" begin
    include("$(het.TEST_DIR)/order_verification/ovs_energy.jl")
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
        case_ind = ((i - 1) ÷ num_norms) + 1
        println("Electron energy ($(cases[case_ind]), $(norms_to_test[norm_ind])-norm): ",
            slope,)
        @test slope > 0.8
    end
end
=#

@testset "Order verification (neutrals and ions)" begin
    include("$(het.TEST_DIR)/order_verification/ovs_ions.jl")
    refinements = refines(5, 10, 2)

    limiter = het.van_leer
    flux_names = ["HLLE", "Rusanov", "Global Lax-Friedrichs"]
    fluxes = [het.HLLE, het.rusanov, het.global_lax_friedrichs]

    # Which L-P norms to check
    cases = ["ρn", "ρi", "ρiui"]
    norms_to_test = [1, 2, Inf]
    num_norms = length(norms_to_test)

    for (flux, flux_name) in zip(fluxes, flux_names)
        for reconstruct in [false, true]
            scheme = het.HyperbolicScheme(flux, limiter, reconstruct)
            orders, norms = test_refinements(
                ncells -> OVS_Ions.solve_ions(ncells, scheme), refinements, norms_to_test,)
            for (i, (order, norm)) in enumerate(zip(orders, norms))
                norm_ind = mod1(i, num_norms)
                case_ind = ((i - 1) ÷ num_norms) + 1
                case_str = "($(cases[case_ind]), $(norms_to_test[norm_ind])-norm)"

                if !reconstruct
                    println("No reconstruction, $flux_name ($case_str): ", order)
                else
                    println("With reconstruction, $flux_name ($case_str):  ", order)
                end

                # Check that we achieve the desired order of accuracy
                if (reconstruct)
                    @test order >= 1.5
                else
                    @test order > 0.75
                end
            end
        end
    end
end
