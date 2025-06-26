using HallThruster: HallThruster as het
using Test

include("$(het.TEST_DIR)/order_verification/ovs_energy.jl")
include("$(het.TEST_DIR)/order_verification/ovs_funcs.jl")
include("$(het.TEST_DIR)/order_verification/ovs_ions.jl")

using .OVS_Energy
using .OVS_Ions

function test_ovs_energy()
    return @testset "Order verification (electron energy)" begin
        vfunc = x -> OVS_Energy.verify_energy_all(x)
        refinements = refines(6, 20, 2)

        cases = ["Implicit", "Crank-Nicholson"]
        norms_to_test = [1, 2, Inf]
        num_norms = length(norms_to_test)

        # Test spatial order of accuracy of implicit solver and crank-nicholson on L1, L2, and L∞ norms
        slopes_nϵ, norms_nϵ = test_refinements(vfunc, refinements, norms_to_test)

        # Check that electron energy equation is solved to at least first order in space
        for (i, slope) in enumerate(slopes_nϵ)
            norm_ind = mod1(i, num_norms)
            case_ind = ((i - 1) ÷ num_norms) + 1
            println(
                "Electron energy ($(cases[case_ind]), $(norms_to_test[norm_ind])-norm): ",
                slope,
            )
            @test slope > 0.8
        end
    end
end

function test_ovs_ions()
    return @testset "Order verification (neutrals and ions)" begin
        refinements = refines(5, 10, 2)
        # Which L-P norms to check
        cases = ["ρn", "ρi", "ρiui"]
        norms_to_test = [1, 2, Inf]
        num_norms = length(norms_to_test)

        for reconstruct in [false, true]
            orders, _ = test_refinements(
                ncells -> OVS_Ions.solve_ions(ncells, reconstruct), refinements, norms_to_test,
            )
            for (i, order) in enumerate(orders)
                norm_ind = mod1(i, num_norms)
                case_ind = ((i - 1) ÷ num_norms) + 1
                case_str = "($(cases[case_ind]), $(norms_to_test[norm_ind])-norm)"

                if !reconstruct
                    println("No reconstruction: ", order)
                else
                    println("With reconstruction:  ", order)
                end

                # Check that we achieve the desired order of accuracy
                if reconstruct
                    @test order >= 1.5
                else
                    @test order > 0.5
                end
            end
        end
    end
end

test_ovs_energy()
#test_ovs_ions()
