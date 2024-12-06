using HallThruster: HallThruster as het
using Test

include("regression_utils.jl")

function test_spt100_regression()
    @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 85.394,
            current = 4.572,
            ion_current = 3.876,
            max_Te = 24.306,
            max_E = 6.531e4,
            max_nn = 2.092e19,
            max_ni = 8.933e17,
        )
        check_regression_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 99.159,
            current = 4.933,
            ion_current = 3.851,
            max_Te = 27.322,
            max_E = 9.464e4,
            max_nn = 2.096e19,
            max_ni = 1.056e18,
        )
        check_regression_case(with_plume)
    end
end
