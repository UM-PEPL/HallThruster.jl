module HallThruster

using LinearAlgebra
using DocStringExtensions

using SparseArrays: Tridiagonal
using PartialFunctions
using QuadGK: quadgk
using DelimitedFiles: readdlm, writedlm
using Unitful: @u_str, uconvert, ustrip, Quantity

using PrecompileTools: @setup_workload, @compile_workload, @recompile_invalidations

using JSON3
using JLD2

# Packages used for making plots
using Measures: mm
using RecipesBase

# path to the HallThruster directory
const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const REACTION_FOLDER = joinpath(PACKAGE_ROOT, "reactions")
const LANDMARK_FOLDER = joinpath(PACKAGE_ROOT, "landmark")
const LANDMARK_RATES_FILE = joinpath(LANDMARK_FOLDER, "landmark_rates.csv")

include("utilities/utility_functions.jl")

include("physics/physicalconstants.jl")
include("physics/gas.jl")
include("physics/fluid.jl")
include("physics/thermal_conductivity.jl")
include("physics/thermodynamics.jl")

include("wall_loss_models/wall_losses.jl")
include("wall_loss_models/no_wall_losses.jl")
include("wall_loss_models/constant_sheath_potential.jl")
include("wall_loss_models/wall_sheath.jl")

include("numerics/finite_differences.jl")
include("numerics/limiters.jl")
include("numerics/flux.jl")

include("collisions/anomalous.jl")
include("collisions/reactions.jl")
include("collisions/ionization.jl")
include("collisions/excitation.jl")
include("collisions/elastic.jl")
include("collisions/charge_exchange.jl")
include("collisions/collision_frequencies.jl")

include("simulation/initialization.jl")
include("simulation/geometry.jl")
include("simulation/boundaryconditions.jl")
include("simulation/potential.jl")
include("simulation/update_heavy_species.jl")
include("simulation/electronenergy.jl")
include("simulation/sourceterms.jl")
include("simulation/plume.jl")
include("simulation/configuration.jl")
include("simulation/update_electrons.jl")
include("simulation/solution.jl")
include("simulation/simulation.jl")
include("simulation/json.jl")
include("simulation/restart.jl")
include("simulation/postprocess.jl")
include("visualization/plotting.jl")
include("visualization/recipes.jl")

export time_average, Xenon, Krypton

# this is an example simulatin that we can run to exercise all parts of the code. this helps to make sure most relevant
# routines are compiled at pre-compile time
function example_simulation(;ncells, duration, dt, nsave)
    config_1 = HallThruster.Config(;
        thruster = HallThruster.SPT_100,
        domain = (0.0u"cm", 8.0u"cm"),
        discharge_voltage = 300.0u"V",
        anode_mass_flow_rate = 5u"mg/s",
        wall_loss_model = WallSheath(BoronNitride),
        neutral_temperature = 500,
    )
    sol_1 = HallThruster.run_simulation(config_1; ncells, duration, dt, nsave, verbose = false)

    config_2 = HallThruster.Config(;
        thruster = HallThruster.SPT_100,
        domain = (0.0u"cm", 8.0u"cm"),
        discharge_voltage = 300.0u"V",
        anode_mass_flow_rate = 5u"mg/s",
        wall_loss_model = ConstantSheathPotential(20.0, 1.0, 1.0),
        LANDMARK = true,
        conductivity_model = LANDMARK_conductivity(),
        neutral_temperature = 500u"K",
        propellant = Krypton,
    )
    sol_2 = HallThruster.run_simulation(config_2; ncells, duration, dt, nsave, adaptive = true, CFL = 0.75, verbose = false)
    HallThruster.time_average(sol_1)
    HallThruster.discharge_current(sol_1)
    HallThruster.thrust(sol_1)
    HallThruster.mass_eff(sol_1)
    HallThruster.current_eff(sol_1)
    HallThruster.divergence_eff(sol_1)
    HallThruster.voltage_eff(sol_1)
    return sol_1, sol_2
end

# Precompile statements to improve load time
@recompile_invalidations begin
@setup_workload begin
    JSON_DIR = joinpath(PACKAGE_ROOT, "test", "json")
    JSON_FILES = readdir(JSON_DIR, join=true)
    outfile = joinpath(JSON_DIR, "_out.json")

    @compile_workload begin
        try
            sol1, sol2 = example_simulation(;ncells=20, duration=1e-7, dt=1e-8, nsave=2)
            avg1 = time_average(sol1)
            write_to_json(outfile, sol1)
            write_to_json(outfile, avg1)
            avg2 = time_average(sol2)
            write_to_json(outfile, sol2)
            write_to_json(outfile, avg2)
            for file in JSON_FILES
                sol = HallThruster.run_simulation(file, verbose = false)
                avg = time_average(sol)
                write_to_json(outfile, sol)
                write_to_json(outfile, avg)
            end
        catch err
            rm(outfile)
            error(err)
        end
        rm(outfile)
    end
end
end

end # module
