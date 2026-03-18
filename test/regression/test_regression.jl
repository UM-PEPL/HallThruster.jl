using HallThruster: HallThruster as het
using Test

# Regression tests
include("$(het.TEST_DIR)/regression/regression_utils.jl")

fix = false

@testset "LANDMARK regression" begin
    case1 = (;
        file = "LANDMARK case 1",
        CFL = 0.5,
        landmark_case = 1,
    )
    case2 = (;
        file = "LANDMARK case 2",
        CFL = 0.799,
        landmark_case = 2,
    )
    case3 = (;
        file = "LANDMARK case 3",
        CFL = 0.799,
        landmark_case = 3,
    )

    for case in [case1, case2, case3]
        check_regression_case(case; fix)
    end
end


@testset "SPT-100 regression" begin
    baseline = (;file = "baseline.json")
    with_plume = (;file = "with_plume.json")

    for case in [baseline, with_plume]
        check_regression_case(case; fix)
    end
end
