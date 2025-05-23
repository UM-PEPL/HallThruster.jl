using HallThruster: HallThruster as het
using Test

# Regression tests
include("$(het.TEST_DIR)/regression/landmark.jl")
include("$(het.TEST_DIR)/regression/spt100.jl")

function test_regression()
    @testset "Regression tests" begin
        test_landmark_regression()
        test_spt100_regression()
    end
end
