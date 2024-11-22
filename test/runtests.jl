using HallThruster
using Test
using Documenter
using DelimitedFiles
using LinearAlgebra
using Printf
using Unitful

doctest(HallThruster)

# exercise all parts of the solver loop
HallThruster.example_simulation(; ncells = 20, duration = 1e-7, dt = 1e-8, nsave = 2)

include("$(HallThruster.TEST_DIR)/unit_tests/thrusters.jl")
test_thrusters()

include("$(HallThruster.TEST_DIR)/unit_tests/restarts.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_gas.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_conservation_laws.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/limiters.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_reactions.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/density_calculations.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/allocation.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/interpolation.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/linear_algebra.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/errors.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/collisions.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_boundary_conditions.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/walls.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/initialization.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_json.jl")

include("$(HallThruster.TEST_DIR)/unit_tests/grid.jl")
test_grid()

include("$(HallThruster.TEST_DIR)/unit_tests/current_control.jl")
test_current_control()

# Regression tests
include("$(HallThruster.TEST_DIR)/regression/test_regression.jl")
test_regression()

# Order verification tests
include("$(HallThruster.TEST_DIR)/order_verification/test_ovs.jl")
test_ovs()
