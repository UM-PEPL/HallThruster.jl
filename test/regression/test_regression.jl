using HallThruster: HallThruster as het
using Test

# Regression tests
include("$(het.TEST_DIR)/regression/regression_utils.jl")

fix = false

@testset "LANDMARK regression" begin
    case1 = (;
        file = "LANDMARK case 1",
        landmark_case = 1,
    )
    case2 = (;
        file = "LANDMARK case 2",
        landmark_case = 2,
    )
    case3 = (;
        file = "LANDMARK case 3",
        landmark_case = 3,
    )

    for case in [case1, case2, case3]
        check_regression_case(case; fix)
    end
end


@testset "SPT-100 regression" begin
    for file in ["spt100_tutorial.json", "spt100_baseline.json", "spt100_withplume.json"]
        check_regression_case((; file); fix)
    end
end
