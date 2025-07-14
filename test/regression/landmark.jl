using HallThruster: HallThruster as het
using Test

function test_landmark_regression()
    case1 = (;
        file = "LANDMARK case 1",
        CFL = 0.5,
        landmark_case = 1,
        thrust = 100.333,
        current = 8.247,
        ion_current = 3.723,
        max_Te = 27.331,
        max_E = 3.926e+4,
        max_nn = 3.886e+19,
        max_ni = 1.182e+18,
        efficiencies = Dict(
            "Mass" => 1.009,
            "Current" => 0.46,
            "Divergence" => 1.0,
            "Voltage" => 0.902,
            "Anode" => 0.431,
        ),
    )

    case2 = (;
        file = "LANDMARK case 2",
        CFL = 0.799,
        landmark_case = 2,
        thrust = 101.167,
        current = 8.082,
        ion_current = 3.673,
        max_Te = 34.068,
        max_E = 4.145e+4,
        max_nn = 5.678e+19,
        max_ni = 4.387e+18,
        efficiencies = Dict(
            "Mass" => 0.949,
            "Current" => 0.455,
            "Divergence" => 1.0,
            "Voltage" => 0.986,
            "Anode" => 0.401,
        ),
    )
    case3 = (;
        file = "LANDMARK case 3",
        CFL = 0.799,
        landmark_case = 3,
        thrust = 100.011,
        current = 7.909,
        ion_current = 3.681,
        max_Te = 35.481,
        max_E = 4.171e+4,
        max_nn = 6.574e+19,
        max_ni = 5.131e+18,
        efficiencies = Dict(
            "Mass" => 0.936,
            "Current" => 0.466,
            "Divergence" => 1.0,
            "Voltage" => 0.988,
            "Anode" => 0.393,
        ),
    )

    return @testset "LANDMARK regression" begin
        check_regression_case(case1)
        check_regression_case(case2)
        check_regression_case(case3)
    end
end
