@enum ConservationLawType begin
    _ContinuityOnly
    _IsothermalEuler
    _EulerEquations
end

struct FluidContainer
    # Conservative variables
    density::Vector{Float64}
    momentum::Vector{Float64}

    # Time derivatives
    dens_ddt::Vector{Float64}
    mom_ddt::Vector{Float64}

    # Timestepping caches
    dens_cache::Vector{Float64}
    mom_cache::Vector{Float64}

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

    function FluidContainer(type, species, num_cells; temp, vel = 0.0)
        R = species.element.R
        γ = species.element.γ
        a = sqrt(γ * R * temp)

        return new(
            # Conservative variables, caches, and time derivatives
            zeros(num_cells + 2), zeros(num_cells + 2),
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

function temperature(f::FluidContainer)
    R = f.species.element.R
    γ = f.species.element.γ
    return f.sound_speed^2 / (γ * R)
end

function allocate_fluids(p::Propellant, ncells)
    continuity = [
        FluidContainer(_ContinuityOnly, p.gas(0), ncells; vel = p.velocity_m_s, temp = p.temperature_K),
    ]

    isothermal = [
        FluidContainer(_IsothermalEuler, p.gas(Z), ncells; temp = p.ion_temperature_K) for Z in 1:p.max_charge
    ]
    return (; continuity, isothermal)
end
