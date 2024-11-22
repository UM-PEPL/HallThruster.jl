using HallThruster: HallThruster as het
using Test

function test_landmark_regression()
    case1 = (;
        file = "LANDMARK case 1",
        CFL = 0.25,
        landmark_case = 1,
        thrust = 89.876,
        current = 7.682,
        ion_current = 3.692,
        max_Te = 22.363,
        max_E = 32248.5,
        max_nn = 3.80456e19,
        max_ni = 1.0378e18,
    )

    case2 = (;
        file = "LANDMARK case 2",
        CFL = 0.799,
        landmark_case = 2,
        thrust = 100.869,
        current = 8.052,
        ion_current = 3.655,
        max_Te = 33.168,
        max_E = 44412.88,
        max_nn = 5.350e19,
        max_ni = 3.937e18,
    )
    case3 = (;
        file = "LANDMARK case 3",
        CFL = 0.799,
        landmark_case = 3,
        thrust = 99.403,
        current = 7.847,
        ion_current = 3.649,
        max_Te = 34.786,
        max_E = 4.48e4,
        max_nn = 6.249e19,
        max_ni = 4.705e18,
    )

    @testset "LANDMARK regression" begin
        check_regression_case(case1)
        check_regression_case(case2)
        check_regression_case(case3)
    end
end
