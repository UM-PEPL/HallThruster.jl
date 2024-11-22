using HallThruster: HallThruster as het
using Test

include("regression_utils.jl")

function test_spt100_regression()
    @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 85.895,
            current = 4.613,
            ion_current = 3.960,
            max_Te = 24.039,
            max_E = 62255.34,
            max_nn = 2.096e19,
            max_ni = 8.754e17,
        )
        check_regression_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 102.527,
            current = 4.990,
            ion_current = 4.072,
            max_Te = 27.709,
            max_E = 95449.049,
            max_nn = 2.1017e19,
            max_ni = 1.017e18,
        )
        check_regression_case(with_plume)
    end
end
