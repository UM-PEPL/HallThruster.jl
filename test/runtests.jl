# doctest(HallThruster)

using SafeTestsets

@safetestset "Doctests" begin
    using HallThruster
    using Documenter

    doctest(HallThruster)
end

@safetestset "Errors" include("unit_tests/errors.jl")
@safetestset "Gases and species" include("unit_tests/gases.jl")
@safetestset "Thrusters" include("unit_tests/thrusters.jl")
@safetestset "Linear interpolation" include("unit_tests/interpolation.jl")
@safetestset "Linear algebra" include("unit_tests/linear_algebra.jl")
@safetestset "Heavy species numerics" include("unit_tests/numerics.jl")
@safetestset "Solution objects" include("unit_tests/solution.jl")
@safetestset "Simulation setup" include("unit_tests/setup.jl")
@safetestset "Grid" include("unit_tests/grid.jl")
@safetestset "Current control" include("unit_tests/current_control.jl")
@safetestset "Reactions" include("unit_tests/reactions.jl")
@safetestset "Thermal conductivity serialization" include("unit_tests/thermal_conductivity.jl")
@safetestset "Wall losses" include("unit_tests/walls.jl")
@safetestset "Collisions and mobility" include("unit_tests/collisions.jl")
@safetestset "Boundary conditions" include("unit_tests/boundary_conditions.jl")
@safetestset "JSON serialization" include("json/json.jl")
@safetestset "Prediction regression" include("regression/test_regression.jl")

# # Order verification tests
# include("$(HallThruster.TEST_DIR)/order_verification/test_ovs.jl")
# test_ovs()
