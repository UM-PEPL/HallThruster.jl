"""
    mutable struct EnergyOVS
Enables setting mu, ue, Tev and ne to certain values to very electron energy equation
"""
mutable struct EnergyOVS{F1, F2}
    active ::Int64
    μ::Union{Float64, Nothing}
    ue::Union{Float64, Nothing}
    Tev::F1
    ne::F2
end

"""
    mutable struct Verification
is passed to params to identify if OVS is active.
"""
mutable struct Verification{F1, F2}
    potential::Int64
    fluid::Int64
    energy::EnergyOVS{F1, F2}
end

struct HyperbolicScheme{F,L}
    flux_function::F
    limiter::L
    reconstruct::Bool
end

Base.@kwdef mutable struct MultiFluidSimulation{IC,B1,B2,B3,B4,S,F,L,CB,SP,BP,VF} #could add callback, or autoselect callback when in MMS mode
    grid::Grid1D
    fluids::Vector{Fluid}     # An array of user-defined fluids.
    # This will give us the capacity to more easily do shock tubes (and other problems)
    # without Hall thruster baggage
    initial_condition::IC
    boundary_conditions::Tuple{B1,B2,B3,B4}   # Tuple of left and right boundary conditions, subject to the approval of PR #10
    end_time::Float64    # How long to simulate
    scheme::HyperbolicScheme{F,L} # Flux, Limiter
    source_term!::S  # Source term function. This can include reactons, electric field, and MMS terms
    source_potential!::SP #potential source term
    boundary_potential!::BP #boundary conditions potential
    saveat::Vector{Float64} #when to save
    timestepcontrol::Tuple{Float64,Bool} #sets timestep (first argument) if second argument false. if second argument (adaptive) true, given dt is ignored.
    callback::CB
    solve_energy::Bool
    verification::VF
end

struct wall_material{F}
    SEE_func::F
    conducting::Bool
end