abstract type AbstractFluid end

@enum ConservationLawType begin
    _ContinuityOnly
    _IsothermalEuler
end

struct Fluid <: AbstractFluid
    # Conservative variables
    density::Vector{Float64}
    momentum::Vector{Float64}

    # Time deriviatives
    dens_ddt::Vector{Float64}
    mom_ddt::Vector{Float64}

    # Edge states
    dens_L::Vector{Float64}
    dens_R::Vector{Float64}
    mom_L::Vector{Float64}
    mom_R::Vector{Float64}

    # Fluxes
    flux_dens::Vector{Float64}
    flux_mom::Vector{Float64}

    # Data
    wave_speed::Array{Float64, 0}
    max_timestep::Array{Float64, 0}
    species::Species
    sound_speed::Float64
    const_velocity::Float64
    type::ConservationLawType

    function Fluid(type, species, num_cells; temp, vel = 0.0)
        R = species.element.R
        γ = species.element.γ
        a = sqrt(γ * R * temp)

        return new(
            # Conservative variables and time derivatives
            zeros(num_cells + 2), zeros(num_cells + 2),
            zeros(num_cells + 2), zeros(num_cells + 2),

            # Edge states
            zeros(num_cells + 1), zeros(num_cells + 1),
            zeros(num_cells + 1), zeros(num_cells + 1),

            # Fluxes
            zeros(num_cells + 1), zeros(num_cells + 1),

            # Data
            fill(max(abs(vel + a), abs(vel - a))), fill(0.0), species, a, vel, type
        )
    end
end

function temperature(f::Fluid)
    R = f.species.element.R
    γ = f.species.element.γ
    return f.sound_speed^2 / (γ * R)
end

function allocate_fluids(propellant, ncharge, ncells, u_neutral, T_neutral, T_ion)
    continuity = [
        Fluid(_ContinuityOnly, propellant(0), ncells; vel = u_neutral, temp = T_neutral),
    ]

    isothermal = [
        Fluid(_IsothermalEuler, propellant(Z), ncells; temp = T_ion) for Z in 1:ncharge
    ]
    return (; continuity, isothermal)
end

function _to_state_vector!(U, continuity, isothermal)
    @. @views U[1, :] = continuity[1].density

    for (i, fluid) in enumerate(isothermal)
        @. @views U[2 * i, :] = fluid.density
        @. @views U[2 * i + 1, :] = fluid.momentum
    end
    return
end

function _from_state_vector!(continuity, isothermal, U)
    @. @views continuity[1].density = U[1, :]

    for (i, fluid) in enumerate(isothermal)
        @. @views fluid.density = U[2 * i, :]
        @. @views fluid.momentum = U[2 * i + 1, :]
    end
    return
end
