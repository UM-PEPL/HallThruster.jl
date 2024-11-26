using HallThruster: HallThruster as het
using Test

include("regression_utils.jl")

function test_spt100_regression()
    @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 85.487,
            current = 4.616,
            ion_current = 3.935,
            max_Te = 24.033,
            max_E = 6.233e4,
            max_nn = 2.096e19,
            max_ni = 8.749e17,
        )
        check_regression_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 103.029,
            current = 4.988,
            ion_current = 4.090,
            max_Te = 27.784,
            max_E = 9.649e4,
            max_nn = 2.102e19,
            max_ni = 1.018e18,
        )
        check_regression_case(with_plume)
    end
end
