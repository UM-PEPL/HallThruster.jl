using HallThruster: HallThruster as het
using Test

function test_spt100_regression()
    return @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 86.979,
            current = 4.62,
            ion_current = 3.906,
            max_Te = 24.712,
            max_E = 6.548e+4,
            max_nn = 2.088e+19,
            max_ni = 9.641e+17,
            efficiencies = Dict(
                "Mass" => 0.968,
                "Current" => 0.876,
                "Divergence" => 0.949,
                "Voltage" => 0.66,
                "Anode" => 0.594,
            ),
        )
        check_regression_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 105.837,
            current = 5.012,
            ion_current = 4.057,
            max_Te = 27.566,
            max_E = 9.858e+4,
            max_nn = 2.095e+19,
            max_ni = 1.163e+18,
            efficiencies = Dict(
                "Mass" => 0.991,
                "Current" => 0.828,
                "Divergence" => 0.963,
                "Voltage" => 0.722,
                "Anode" => 0.596,
            ),
        )
        check_regression_case(with_plume)
    end
end
