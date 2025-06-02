module HallThruster

using LinearAlgebra: Tridiagonal
using DocStringExtensions

using DelimitedFiles: readdlm, writedlm
using PrecompileTools: @compile_workload

using JSON3
using OrderedCollections

# Packages used for making plots
using Measures: mm
using RecipesBase

# path to the HallThruster directory
const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const REACTION_FOLDER = joinpath(PACKAGE_ROOT, "reactions")
const LANDMARK_FOLDER = joinpath(PACKAGE_ROOT, "landmark")
const LANDMARK_RATES_FILE = joinpath(LANDMARK_FOLDER, "landmark_rates.csv")
const TEST_DIR = joinpath(PACKAGE_ROOT, "test")

const MIN_NUMBER_DENSITY = 1e6

include("utilities/utility_functions.jl")
include("utilities/interpolation.jl")
include("utilities/smoothing.jl")
include("utilities/statistics.jl")
include("utilities/linearalgebra.jl")
include("utilities/integration.jl")
include("utilities/macros.jl")
include("utilities/serialization.jl")
include("utilities/keywords.jl")
include("utilities/units.jl")
include("utilities/finite_differences.jl")

using .Serialization: serialize, deserialize

include("physics/physicalconstants.jl")
include("physics/gas.jl")
include("physics/fluid.jl")
include("physics/thermal_conductivity.jl")

include("walls/materials.jl")
include("walls/wall_losses.jl")
include("walls/no_wall_losses.jl")
include("walls/constant_sheath_potential.jl")
include("walls/wall_sheath.jl")

include("collisions/anomalous.jl")
include("collisions/reactions.jl")
include("collisions/ionization.jl")
include("collisions/excitation.jl")
include("collisions/elastic.jl")
include("collisions/collision_frequencies.jl")

include("thruster/geometry.jl")
include("thruster/magnetic_field.jl")
include("thruster/thruster.jl")
include("thruster/spt100.jl")

include("grid/gridspec.jl")
include("grid/grid.jl")

include("simulation/current_control.jl")
include("simulation/initialization.jl")
include("simulation/configuration.jl")
include("simulation/allocation.jl")
include("simulation/potential.jl")
include("simulation/heavy_species_fluxes.jl")
include("simulation/heavy_species_update.jl")
include("simulation/electron_energy.jl")
include("simulation/sourceterms.jl")
include("simulation/plume.jl")
include("simulation/electron_update.jl")
include("simulation/solution.jl")
include("simulation/postprocess.jl")
include("simulation/simulation.jl")
include("simulation/json.jl")
include("visualization/plotting.jl")
include("visualization/recipes.jl")

@public PYTHON_PATH

"""
$(SIGNATURES)

The absolute path to the `HallThruster` python code on your machine.

```julia
julia> using HallThruster; HallThruster.PYTHON_PATH
"/Users/archermarks/src/HallThruster/python"
```
"""
const PYTHON_PATH = joinpath(PACKAGE_ROOT, "python")

# this is an example simulatin that we can run to exercise all parts of the code. this helps to make sure most relevant
# routines are compiled at pre-compile time
function example_simulation(; ncells, duration, dt, nsave)
    config_1 = Config(;
        thruster = HallThruster.SPT_100,
        domain = (0.0, 0.08),
        discharge_voltage = 300.0,
        anode_mass_flow_rate = 5e-6,
        wall_loss_model = WallSheath(BoronNitride),
        neutral_temperature_K = 500,
    )
    sol_1 = run_simulation(
        config_1; ncells, duration, dt, nsave, verbose = false,)

    if sol_1.retcode != :success
        error()
    end

    config_2 = HallThruster.Config(;
        thruster = SPT_100,
        domain = (0.0, 0.08),
        discharge_voltage = 300.0,
        anode_mass_flow_rate = 5e-6,
        anom_model = MultiLogBohm([0.02, 0.025, 0.03], [0.0625, 0.00625, 0.0625]),
        wall_loss_model = ConstantSheathPotential(20.0, 1.0, 1.0),
        LANDMARK = true,
        conductivity_model = LANDMARK_conductivity(),
        neutral_temperature_K = 500,
        ion_wall_losses = true,
        solve_plume = true,
    )

    sol_2 = run_simulation(
        config_2; ncells, duration, dt, nsave, adaptive = true, CFL = 0.75, verbose = false,)

    if sol_2.retcode != :success
        error()
    end

    config_4 = HallThruster.Config(;
        thruster = SPT_100,
        domain = (0.0, 0.08),
        discharge_voltage = 300.0,
        anode_mass_flow_rate = 5e-6,
        anom_model = ScaledGaussianBohm(anom_scale = 0.0625, barrier_scale = 0.9, center = 1.0, width = 0.1),
        wall_loss_model = ConstantSheathPotential(20.0, 1.0, 1.0),
        conductivity_model = Braginskii(),
        neutral_temperature_K = 500,
        ion_wall_losses = false,
        solve_plume = false,
    )

    sol_4 = run_simulation(
        config_4; ncells, duration, dt, nsave, adaptive = true, CFL = 0.75, verbose = false,)

    if sol_4.retcode != :success
        error()
    end

    time_average(sol_1)
    discharge_current(sol_1)
    thrust(sol_1)
    mass_eff(sol_1)
    current_eff(sol_1)
    divergence_eff(sol_1)
    voltage_eff(sol_1)
    return sol_1
end

# Precompile statements to improve load time
@compile_workload begin
    example_simulation(; ncells = 20, duration = 1e-7, dt = 1e-8, nsave = 2)

    for file in readdir(joinpath(TEST_DIR, "precompile"), join = true)
        if splitext(file)[2] != ".json"
            continue
        end
        sol = run_simulation(file)
    end
    rm("__output.json", force = true)
end

end # module
