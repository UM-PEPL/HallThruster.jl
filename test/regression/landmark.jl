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
        thrust = 101.481,
        current = 8.128,
        ion_current = 3.674,
        max_Te = 33.068,
        max_E = 4.102e+4,
        max_nn = 5.468e+19,
        max_ni = 4.168e+18,
        efficiencies = Dict(
            "Mass" => 0.954,
            "Current" => 0.453,
            "Divergence" => 1.0,
            "Voltage" => 0.985,
            "Anode" => 0.403,
        ),
    )
    case3 = (;
        file = "LANDMARK case 3",
        CFL = 0.799,
        landmark_case = 3,
        thrust = 100.187,
        current = 7.939,
        ion_current = 3.677,
        max_Te = 34.686,
        max_E = 4.139e+4,
        max_nn = 6.398e+19,
        max_ni = 4.969e+18,
        efficiencies = Dict(
            "Mass" => 0.938,
            "Current" => 0.464,
            "Divergence" => 1.0,
            "Voltage" => 0.988,
            "Anode" => 0.395,
        ),
    )

    return @testset "LANDMARK regression" begin
        check_regression_case(case1)
        check_regression_case(case2)
        check_regression_case(case3)
    end
end
