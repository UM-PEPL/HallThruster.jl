using HallThruster: HallThruster as het
using Test

include("regression_utils.jl")

function test_spt100_regression()
    @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 87.490,
            current = 4.618,
            ion_current = 3.934,
            max_Te = 24.736,
            max_E = 6.612e4,
            max_nn = 2.088e19,
            max_ni = 9.650e17,
			efficiencies = Dict(
				"Mass" => 0.957,
				"Current" => 0.877,
				"Divergence" => 0.949,
				"Voltage" => 0.660,
				"Anode" => 0.572,
			)
        )
        check_regression_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 105.644,
            current = 5.019,
            ion_current = 4.058,
            max_Te = 27.327,
            max_E = 9.731e4,
            max_nn = 2.096e19,
            max_ni = 1.167e18,
			efficiencies = Dict(
				"Mass" => 0.993,
				"Current" => 0.825,
				"Divergence" => 0.963,
				"Voltage" => 0.719,
				"Anode" => 0.597,
			)
        )
        check_regression_case(with_plume)
    end
end
