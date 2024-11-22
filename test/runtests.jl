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

include("$(HallThruster.TEST_DIR)/unit_tests/test_restarts.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_gas.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_conservation_laws.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_limiters.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_reactions.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_misc.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_geometry.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_electrons.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_boundary_conditions.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_walls.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_initialization.jl")
include("$(HallThruster.TEST_DIR)/unit_tests/test_json.jl")

include("$(HallThruster.TEST_DIR)/unit_tests/grid.jl")
test_grid()

# Regression tests
include("$(HallThruster.TEST_DIR)/test_regression.jl")
test_regression()

# Order verification tests
include("$(HallThruster.TEST_DIR)/order_verification/test_ovs.jl")
