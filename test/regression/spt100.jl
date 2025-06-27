using HallThruster: HallThruster as het
using Test

function test_spt100_regression()
    return @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 86.841,
            current = 4.614,
            ion_current = 3.896,
            max_Te = 24.832,
            max_E = 6.62e+4,
            max_nn = 2.088e+19,
            max_ni = 9.641e+17,
            efficiencies = Dict(
                "Mass" => 0.952,
                "Current" => 0.872,
                "Divergence" => 0.949,
                "Voltage" => 0.661,
                "Anode" => 0.565,
            ),
        )
        check_regression_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 105.605,
            current = 5.023,
            ion_current = 4.056,
            max_Te = 27.233,
            max_E = 9.635e+4,
            max_nn = 2.095e+19,
            max_ni = 1.166e+18,
            efficiencies = Dict(
                "Mass" => 0.995,
                "Current" => 0.827,
                "Divergence" => 0.963,
                "Voltage" => 0.72,
                "Anode" => 0.601,
            ),
        )
        check_regression_case(with_plume)
    end
end
