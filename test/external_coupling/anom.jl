using HallThruster: HallThruster as het
using Libdl

#==============================================================================
    Setup - compile and load libraries
==============================================================================#
function compile(src, obj; compiler = "gfortran")
    compilation_command = `$compiler $src -o $obj -shared -fPIC`
    run(compilation_command)
    return
end

# Library names must be constants
const src_f90 = "$(@__DIR__)/anom.f90"
const obj_f90 = "$(@__DIR__)/anom_f90.so"

const src_cpp = "$(@__DIR__)/anom.cpp"
const obj_cpp = "$(@__DIR__)/anom_cpp.so"

# Compile the source code for each case and load the necessary symbols
compile(src_f90, obj_f90; compiler = "gfortran")
lib_f90 = Libdl.dlopen(obj_f90)
two_zone_bohm_f90 = Libdl.dlsym(lib_f90, :two_zone_bohm)

compile(src_cpp, obj_cpp; compiler = "gcc")
lib_cpp = Libdl.dlopen(obj_cpp)
two_zone_bohm_cpp = Libdl.dlsym(lib_cpp, :two_zone_bohm)

#==============================================================================
    Standalone tests - common parameters
==============================================================================#

num_cells = 10
z = LinRange(0, 0.08, num_cells) |> collect
B = ones(num_cells)
nu_an = zeros(num_cells)

q_e = 1.60217663e-19
m_e = 9.1093837e-31
c = (0.5, 1.0)
channel_length = 0.025

#==============================================================================
    Fortran - standalone call
==============================================================================#

Lch_ref = Ref(channel_length)
num_cells_ref = Ref(num_cells)

@ccall $two_zone_bohm_f90(
    nu_an::Ptr{Float64},
    z::Ptr{Float64},
    B::Ptr{Float64},
    c::Ref{NTuple{2, Float64}},
    channel_length::Ref{Float64},
    num_cells::Ref{Int64}
)::Cvoid

expected = [
    (_z < channel_length ? c[1] : c[2]) * (q_e * _B / m_e)
        for (_z, _B) in zip(z, B)
]
@test all(nu_an .≈ expected)

#==============================================================================
    C/C++: standalone call
==============================================================================#
nu_an .= 0.0

@ccall $two_zone_bohm_cpp(
    nu_an::Ptr{Float64},
    z::Ptr{Float64},
    B::Ptr{Float64},
    c::Ref{NTuple{2, Float64}},
    channel_length::Float64,
    num_cells::Int32
)::Cvoid

@test all(nu_an .≈ expected)

#==============================================================================
    Fortran: integrate into anom model
==============================================================================#
struct TwoZoneBohm_F90 <: het.AnomalousTransportModel
    c::NTuple{2, Float64}
end

function (model::TwoZoneBohm_F90)(nu_an, params, _)
    (; thruster, grid, cache) = params
    channel_length = thruster.geometry.channel_length
    num_cells = length(grid.cell_centers)

    @ccall $two_zone_bohm_f90(
        nu_an::Ptr{Float64},
        grid.cell_centers::Ptr{Float64},
        cache.B::Ptr{Float64},
        model.c::Ref{NTuple{2, Float64}},
        channel_length::Ref{Float64},
        num_cells::Ref{Int64}
    )::Cvoid

    return nu_an
end

model = TwoZoneBohm_F90((0.05, 0.5))

config = het.Config(
    thruster = het.SPT_100,
    propellants = [het.Propellant(het.Xenon, 5.0e-6)],
    domain = (0.0, 0.08),
    discharge_voltage = 300.0,
    anom_model = model
)

num_cells = 10
simparams = het.SimParams(grid = het.UnevenGrid(num_cells))
params = het.setup_simulation(config, simparams)
model(params.cache.νan, params, nothing)

nu_an_sim = params.cache.νan
expected = [
    (_z < channel_length ? model.c[1] : model.c[2]) * (q_e * _B / m_e)
        for (_z, _B) in zip(params.grid.cell_centers, params.cache.B)
]

@test all(expected .≈ nu_an_sim)

#==============================================================================
    Close libraries
==============================================================================#
Libdl.dlclose(lib_f90)
Libdl.dlclose(lib_cpp)
