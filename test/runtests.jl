using SafeTestsets
using Documenter
using HallThruster


doctest(HallThruster)

@safetestset "Thrusters" include("unit_tests/thrusters.jl")
@safetestset "Limiters" include("unit_tests/limiters.jl")
@safetestset "Allocation" include("unit_tests/allocation.jl")
@safetestset "Interpolation.jl" include("unit_tests/interpolation.jl")
@safetestset "Configuration" include("unit_tests/configuration.jl")
@safetestset "Errors" include("unit_tests/errors.jl")
@safetestset "Linear algebra" include("unit_tests/linear_algebra.jl")
@safetestset "Collisions" include("unit_tests/collisions.jl")
@safetestset "Thermal conductivity" include("unit_tests/thermal_conductivity.jl")
@safetestset "Walls" include("unit_tests/walls.jl")
@safetestset "Reactions" include("unit_tests/reactions.jl")
@safetestset "Gases" include("unit_tests/gases.jl")
@safetestset "Boundary conditions" include("unit_tests/boundary_conditions.jl")
@safetestset "Grid" include("unit_tests/grid.jl")
@safetestset "Current control" include("unit_tests/current_control.jl")
@safetestset "Setup" include("unit_tests/setup.jl")
@safetestset "JSON" include("json/json.jl")
@safetestset "External coupling" include("external_coupling/anom.jl")
@safetestset "Output regression" include("regression/test_regression.jl")
@safetestset "Order verification" include("order_verification/test_ovs.jl")