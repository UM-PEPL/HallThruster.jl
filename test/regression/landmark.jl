using HallThruster: HallThruster as het
using Test

function test_landmark_regression()
    case1 = (;
        file = "LANDMARK case 1",
        CFL = 0.5,
        landmark_case = 1,
        thrust = 97.155,
        current = 8.061,
        ion_current = 3.703,
        max_Te = 24.335,
        max_E = 3.614e+4,
        max_nn = 3.819e+19,
        max_ni = 1.051e+18,
        efficiencies = Dict(
            "Mass" => 1.012,
            "Current" => 0.471,
            "Divergence" => 1.0,
            "Voltage" => 0.873,
            "Anode" => 0.435,
        ),
    )

    case2 = (;
        file = "LANDMARK case 2",
        CFL = 0.799,
        landmark_case = 2,
        thrust = 99.44,
        current = 7.946,
        ion_current = 3.663,
        max_Te = 33.367,
        max_E = 4.081e+4,
        max_nn = 5.272e+19,
        max_ni = 3.393e+18,
        efficiencies = Dict(
            "Mass" => 0.967,
            "Current" => 0.462,
            "Divergence" => 1.0,
            "Voltage" => 0.947,
            "Anode" => 0.402,
        ),
    )
    case3 = (;
        file = "LANDMARK case 3",
        CFL = 0.799,
        landmark_case = 3,
        thrust = 97.899,
        current = 7.744,
        ion_current = 3.65,
        max_Te = 35.098,
        max_E = 4.122e+4,
        max_nn = 6.128e+19,
        max_ni = 4.114e+18,
        efficiencies = Dict(
            "Mass" => 0.95,
            "Current" => 0.472,
            "Divergence" => 1.0,
            "Voltage" => 0.95,
            "Anode" => 0.394,
        ),
    )

    return @testset "LANDMARK regression" begin
        check_regression_case(case1)
        check_regression_case(case2)
        check_regression_case(case3)
    end
end
