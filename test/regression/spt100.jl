using HallThruster: HallThruster as het
using Test

function test_spt100_regression()
    return @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 87.304,
            current = 4.614,
            ion_current = 3.922,
            max_Te = 24.832,
            max_E = 6.665e4,
            max_nn = 2.088e19,
            max_ni = 9.651e17,
            efficiencies = Dict(
                "Mass" => 0.954,
                "Current" => 0.873,
                "Divergence" => 0.949,
                "Voltage" => 0.661,
                "Anode" => 0.5656,
            ),
        )
        check_regression_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 106.033,
            current = 5.022,
            ion_current = 4.058,
            max_Te = 27.418,
            max_E = 9.824e+4,
            max_nn = 2.095e+19,
            max_ni = 1.162e+18,
            efficiencies = Dict(
                "Mass" => 0.992,
                "Current" => 0.823,
                "Divergence" => 0.963,
                "Voltage" => 0.722,
                "Anode" => 0.6,
            ),
        )
        check_regression_case(with_plume)
    end
end
