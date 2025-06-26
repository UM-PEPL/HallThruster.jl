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
            thrust = 105.973,
            current = 5.028,
            ion_current = 4.068,
            max_Te = 27.094,
            max_E = 9.604e+4,
            max_nn = 2.096e+19,
            max_ni = 1.164e+18,
            efficiencies = Dict(
                "Mass" => 0.995,
                "Current" => 0.824,
                "Divergence" => 0.963,
                "Voltage" => 0.72,
                "Anode" => 0.601,
            ),
        )
        check_regression_case(with_plume)
    end
end
