using HallThruster: HallThruster as het
using Test

function test_spt100_regression()
    return @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 86.999,
            current = 4.621,
            ion_current = 3.908,
            max_Te = 24.71,
            max_E = 6.543e+4,
            max_nn = 2.088e+19,
            max_ni = 9.641e+17,
            efficiencies = Dict(
                "Mass" => 0.968,
                "Current" => 0.877,
                "Divergence" => 0.949,
                "Voltage" => 0.66,
                "Anode" => 0.594,
            ),
        )
        check_regression_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 105.119,
            current = 5.018,
            ion_current = 4.021,
            max_Te = 27.384,
            max_E = 9.806e+4,
            max_nn = 2.095e+19,
            max_ni = 1.163e+18,
            efficiencies = Dict(
                "Mass" => 0.987,
                "Current" => 0.822,
                "Divergence" => 0.963,
                "Voltage" => 0.721,
                "Anode" => 0.586,
            ),
        )
        check_regression_case(with_plume)
    end
end
