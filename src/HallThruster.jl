module HallThruster

using Accessors: @reset
using DelimitedFiles: readdlm, writedlm
using DocStringExtensions: SIGNATURES, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using JSON3: JSON3, StructTypes
using LinearAlgebra: Tridiagonal
using PrecompileTools: @compile_workload, @recompile_invalidations
using Unitful: @u_str, uconvert, ustrip, Quantity

# path to the HallThruster directory
const PACKAGE_ROOT = joinpath(splitpath(@__DIR__)[1:(end - 1)]...)
const REACTION_FOLDER = joinpath(PACKAGE_ROOT, "reactions")
const LANDMARK_FOLDER = joinpath(PACKAGE_ROOT, "landmark")
const LANDMARK_RATES_FILE = joinpath(LANDMARK_FOLDER, "landmark_rates.csv")
const TEST_DIR = joinpath(PACKAGE_ROOT, "test", "unit_tests")

include("utilities/utility_functions.jl")
include("utilities/macros.jl")
include("utilities/keywords.jl")
include("utilities/serialization.jl")

using .Serialization: serialize, deserialize, OrderedDict

include("physics/physicalconstants.jl")
include("physics/gas.jl")
include("physics/species.jl")
include("physics/fluid.jl")
include("physics/thermal_conductivity.jl")
include("physics/thermodynamics.jl")

include("wall_loss_models/wall_losses.jl")
include("wall_loss_models/no_wall_losses.jl")
include("wall_loss_models/constant_sheath_potential.jl")
include("wall_loss_models/wall_materials.jl")
include("wall_loss_models/wall_sheath.jl")

include("numerics/finite_differences.jl")
include("numerics/limiters.jl")
include("numerics/flux.jl")
include("numerics/schemes.jl")
include("numerics/edgefluxes.jl")

include("collisions/anomalous.jl")
include("collisions/reactions.jl")
include("collisions/ionization.jl")
include("collisions/excitation.jl")
include("collisions/elastic.jl")
include("collisions/collision_frequencies.jl")

include("thruster/geometry.jl")
include("thruster/magneticfield.jl")
include("thruster/thruster.jl")
include("thruster/spt100.jl")

include("grid/grid.jl")
include("grid/gridspec.jl")

include("simulation/initialization.jl")
include("simulation/configuration.jl")
include("simulation/current_control.jl")
include("simulation/allocate.jl")
include("simulation/params.jl")
include("simulation/state.jl")

include("simulation/boundaryconditions.jl")
include("simulation/potential.jl")
include("simulation/update_heavy_species.jl")
include("simulation/electronenergy.jl")
include("simulation/sourceterms.jl")
include("simulation/plume.jl")
include("simulation/update_electrons.jl")
include("simulation/solution.jl")
include("simulation/simulation.jl")
include("simulation/json.jl")

include("output/postprocess.jl")

export time_average, Xenon, Krypton
export EvenGrid, UnevenGrid

# this is an example simulatin that we can run to exercise all parts of the code. this helps to make sure most relevant
# routines are compiled at pre-compile time
function example_simulation(; ncells, duration, dt, nsave)
    config_1 = HallThruster.Config(;
        thruster = HallThruster.SPT_100,
        domain = (0.0, 8.0),
        discharge_voltage = 300,
        anode_mass_flow_rate = 5e-6,
        wall_loss_model = WallSheath(BoronNitride, 1.0),
        neutral_temperature_K = 500.0,
    )
    sol_1 = run_simulation(config_1; ncells, duration, dt, nsave, verbose = false)

    config_2 = HallThruster.Config(;
        thruster = HallThruster.SPT_100,
        domain = (0.0, 8.0),
        discharge_voltage = 300.0,
        anode_mass_flow_rate = 5e-6,
        wall_loss_model = ConstantSheathPotential(20.0, 1.0, 1.0),
        LANDMARK = true,
        conductivity_model = LANDMARK_conductivity(),
        neutral_temperature_K = 500.0,
    )

    sol_2 = run_simulation(
        config_2; ncells, duration, dt, nsave, adaptive = true, CFL = 0.75, verbose = false,
    )

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
    outfile = joinpath(TEST_DIR, "_out.json")

    rm(outfile; force = true)
    for file in readdir(TEST_DIR; join = true)
        if isdir(file) || splitext(file)[end] != ".json"
            continue
        end
        sim = run_simulation(file; verbose = false)
        avg = time_average(sim)
        write_to_json(outfile, sim)
        write_to_json(outfile, avg)
    end
    rm(outfile; force = true)
end

# Entrypoint for static compilation
Base.@ccallable function main()::Cint
    return 0
end

end # module
