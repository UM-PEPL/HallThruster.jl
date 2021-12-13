function flux(U, fluid)
    if fluid.conservation_laws.type == :ContinuityOnly
        ρ = U[1]
        u = velocity(U, fluid)
        F = (ρ * u, 0.0, 0.0)
    elseif fluid.conservation_laws.type == :IsothermalEuler
        ρ, ρu = U
        u = velocity(U, fluid)
        p = pressure(U, fluid)
        F = (ρu, ρu * u + p, 0.0)
    elseif fluid.conservation_laws.type == :EulerEquations
        ρ, ρu, ρE = U
        u = U[2] / U[1]
        p = pressure(U, fluid)
        ρH = ρE + p
        F = (ρu, ρu * u + p, ρH * u)
    end
    return F
end

function flux!(F, U, fluid)
    if fluid.conservation_laws.type == :ContinuityOnly
        ρ = U[1]
        u = velocity(U, fluid)
        F[1] = ρ * u
    elseif fluid.conservation_laws.type == :IsothermalEuler
        ρ, ρu = U
        u = velocity(U, fluid)
        p = pressure(U, fluid)
        F[1] = ρu
        F[2] = ρu * u + p
    elseif fluid.conservation_laws.type == :EulerEquations
        ρ, ρu, ρE = U
        u = U[2] / U[1]
        p = pressure(U, fluid)
        ρH = ρE + p

        F[1] = ρu
        F[2] = ρu * u + p
        F[3] = ρH * u
    end
    return F
end

function compute_primitive(U, γ)
    ρ, ρu, ρE = U
    u = ρu / ρ
    p = (γ - 1) * (ρE - 0.5 * ρu * u)
    return ρ, u, p
end

function compute_conservative(ρ, u, p, γ)
    ρE = p / (γ - 1) + 0.5 * ρ * u^2
    return ρ, ρ * u, ρE
end

# NOTE: this can be sped up significantly if we write specialized versions for each fluid type
# we're losing a lot of time (~1/4 of the run time) on the conditionals in the thermodynamics, better to do one conditional
# and then go from there. however, that would lead to about 2x more code in this section and a loss of generality. probably
# better to wait to overhaul this until the main features are in and we can think about a refactor
function HLLE!(F, UL, UR, fluid)
    γ = fluid.species.element.γ

    uL = velocity(UL, fluid)
    uR = velocity(UR, fluid)

    aL = sound_speed(UL, fluid)
    aR = sound_speed(UR, fluid)

    sL_min, sL_max = min(0, uL - aL), max(0, uL + aL)
    sR_min, sR_max = min(0, uR - aR), max(0, uR + aR)

    smin = min(sL_min, sR_min)
    smax = max(sL_max, sR_max)

    FL = flux(UL, fluid)
    FR = flux(UR, fluid)

    for i in 1:length(F)
        F[i] = 0.5 * (FL[i] + FR[i]) -
               0.5 * (smax + smin) / (smax - smin) * (FR[i] - FL[i]) +
               smax * smin / (smax - smin) * (UR[i] - UL[i])
    end
    return F
end

function upwind!(F, UL, UR, fluid::Fluid)
    uL = velocity(UL, fluid)
    uR = velocity(UR, fluid)
    avg_velocity = 0.5 * (uL + uR)

    if avg_velocity ≥ 0
        flux!(F, UL, fluid)
    else
        flux!(F, UR, fluid)
    end
    return F
end

# Generate out-of-place (allocating) versions of flux functions
for fluxfn in ["HLLE", "upwind"]
    let inplace = Symbol(fluxfn * '!'), outofplace = Symbol(fluxfn)
        eval(quote
                 $outofplace(UL, UR, fluid::Fluid) = $inplace(similar(UL), UL, UR, fluid)
             end)
    end
end

"""
	reconstruct!
Reconstruction using the MUSCL scheme. UL is the flux to the left, therefore evaluated on
the right face, UR is the flux to the right, therefore evaluated on the left face.

"""
function reconstruct!(UL, UR, U, scheme)
    nconservative, ncells = size(U)
    Ψ = scheme.limiter

    # compute left and right edge states
    for i in 2:(ncells - 1)
        for j in 1:nconservative
            u₋ = U[j, i - 1]
            uᵢ = U[j, i]
            u₊ = U[j, i + 1]
            Δu = u₊ - uᵢ
            ∇u = uᵢ - u₋
            r = Δu / ∇u
            UL[j, right_edge(i)] = uᵢ + 0.5 * Ψ(r) * ∇u
            UR[j, left_edge(i)] = uᵢ - 0.5 * Ψ(1 / r) * Δu
        end
    end

    return UL, UR
end

function compute_edge_states!(UL, UR, U, scheme)
    nconservative, ncells = size(U)

    if scheme.reconstruct
        reconstruct!(UL, UR, U, scheme)
    else
        for i in 2:(ncells - 1)
            for j in 1:nconservative
                UL[j, right_edge(i)] = U[j, i]
                UR[j, left_edge(i)] = U[j, i]
            end
        end
    end

    for j in 1:nconservative
        UL[j, 1] = U[j, 1] #2 this would be more consistent
        UR[j, end] = U[j, end] #end - 1 
    end

    return UL, UR
end

function compute_fluxes!(F, UL, UR, fluids, fluid_ranges, scheme)
    _, nedges = size(F)

    for i in 1:nedges
        for (j, (fluid, fluid_range)) in enumerate(zip(fluids, fluid_ranges))
            @views scheme.flux_function(F[fluid_range, i], UL[fluid_range, i],
                                        UR[fluid_range, i], fluid)
        end
    end
    return F
end