using HallThruster
using Test
using Documenter
using DelimitedFiles
using LinearAlgebra
using Printf
using Unitful

doctest(HallThruster)

include("$(HallThruster.TEST_DIR)/unit_tests/thrusters.jl")
test_thrusters()

include("$(HallThruster.TEST_DIR)/unit_tests/thermodynamics.jl")
test_thermodynamics()

include("$(HallThruster.TEST_DIR)/unit_tests/fluxes.jl")
test_fluxes()

include("$(HallThruster.TEST_DIR)/unit_tests/schemes.jl")
test_schemes()

include("$(HallThruster.TEST_DIR)/unit_tests/limiters.jl")
test_limiters()

include("$(HallThruster.TEST_DIR)/unit_tests/allocation.jl")
test_allocation()

include("$(HallThruster.TEST_DIR)/unit_tests/interpolation.jl")
test_interpolation()

include("$(HallThruster.TEST_DIR)/unit_tests/initialization.jl")
test_initialization()

include("$(HallThruster.TEST_DIR)/unit_tests/configuration.jl")
test_configuration()

include("$(HallThruster.TEST_DIR)/unit_tests/errors.jl")
test_errors()

include("$(HallThruster.TEST_DIR)/unit_tests/linear_algebra.jl")
test_linear_algebra()

include("$(HallThruster.TEST_DIR)/unit_tests/density_calculations.jl")
test_density_calculations()

include("$(HallThruster.TEST_DIR)/unit_tests/collisions.jl")
test_collisions()

include("$(HallThruster.TEST_DIR)/unit_tests/thermal_conductivity.jl")
test_thermal_conductivity()

include("$(HallThruster.TEST_DIR)/unit_tests/walls.jl")
test_walls()

include("$(HallThruster.TEST_DIR)/unit_tests/reactions.jl")
test_reactions()

include("$(HallThruster.TEST_DIR)/unit_tests/gases.jl")
test_gases()

include("$(HallThruster.TEST_DIR)/unit_tests/boundary_conditions.jl")
test_boundaries()

include("$(HallThruster.TEST_DIR)/unit_tests/grid.jl")
test_grid()

include("$(HallThruster.TEST_DIR)/unit_tests/current_control.jl")
test_current_control()

# JSON serialization tests
include("$(HallThruster.TEST_DIR)/json/json.jl")
test_json()

# Regression tests
include("$(HallThruster.TEST_DIR)/regression/test_regression.jl")
test_regression()

# Order verification tests
include("$(HallThruster.TEST_DIR)/order_verification/test_ovs.jl")
test_ovs()
