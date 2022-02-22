function flux(U, fluid, pe)
    if fluid.conservation_laws.type == _ContinuityOnly
        ρ = U[1]
        u = velocity(U, fluid)
        F = (ρ * u, 0.0, 0.0)
    elseif fluid.conservation_laws.type == _IsothermalEuler
        ρ, ρu = U
        u = velocity(U, fluid)
        p = pressure(U, fluid)
        F = (ρu, ρu * u + p + pe, 0.0)
    elseif fluid.conservation_laws.type == _EulerEquations
        ρ, ρu, ρE = U
        u = U[2] / U[1]
        p = pressure(U, fluid)
        ρH = ρE + p
        F = (ρu, ρu * u + p, ρH * u)
    end
    return F
end

function flux!(F, U, fluid, pe)
    if fluid.conservation_laws.type == _ContinuityOnly
        ρ = U[1]
        u = velocity(U, fluid)
        F[1] = ρ * u
    elseif fluid.conservation_laws.type == _IsothermalEuler
        ρ, ρu = U
        u = velocity(U, fluid)
        p = pressure(U, fluid)
        F[1] = ρu
        F[2] = ρu * u + p + pe
    elseif fluid.conservation_laws.type == _EulerEquations
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

# NOTE: this can be sped up significantly if we write specialized versions for each fluid type
# we're losing a lot of time (~1/4 of the run time) on the conditionals in the thermodynamics, better to do one conditional
# and then go from there. however, that would lead to about 2x more code in this section and a loss of generality. probably
# better to wait to overhaul this until the main features are in and we can think about a refactor
function HLLE!(F, UL, UR, fluid, peL, peR, coupled)
    γ = fluid.species.element.γ
    R = fluid.species.element.R
    Z = fluid.species.Z

    uL = velocity(UL, fluid)
    uR = velocity(UR, fluid)

    TL = temperature(UL, fluid)
    TR = temperature(UR, fluid)

    mi = m(fluid)

    TeL = max(0.0, peL / UL[1] * mi)
    TeR = max(0.0, peR / UR[1] * mi)

    aL = sqrt((Z * e * coupled * TeL + γ * kB * TL) / mi)
    aR = sqrt((Z * e * coupled * TeR + γ * kB * TR) / mi)

    sL_min, sL_max = min(0, uL - aL), max(0, uL + aL)
    sR_min, sR_max = min(0, uR - aR), max(0, uR + aR)

    smin = min(sL_min, sR_min)
    smax = max(sL_max, sR_max)

    FL = flux(UL, fluid, Z * coupled * e * peL)
    FR = flux(UR, fluid, Z * coupled * e * peR)

    for i in 1:length(F)
        F[i] = 0.5 * (FL[i] + FR[i]) -
               0.5 * (smax + smin) / (smax - smin) * (FR[i] - FL[i]) +
               smax * smin / (smax - smin) * (UR[i] - UL[i])
    end
    return F
end

function upwind!(F, UL, UR, fluid::Fluid, pe_L, pe_R, coupled)
    uL = velocity(UL, fluid)
    uR = velocity(UR, fluid)
    avg_velocity = 0.5 * (uL + uR)
    Z = fluid.species.Z
    if avg_velocity ≥ 0
        flux!(F, UL, fluid, Z * coupled * e * pe_L)
    else
        flux!(F, UR, fluid, Z * coupled * e * pe_R)
    end
    return F
end

# Generate out-of-place (allocating) versions of flux functions
for fluxfn in ["HLLE", "upwind"]
    let inplace = Symbol(fluxfn * '!'), outofplace = Symbol(fluxfn)
        eval(quote
                 $outofplace(UL, UR, fluid::Fluid, pe_L, pe_R) = $inplace(similar(UL), UL, UR, fluid, pe_L, pe_R)
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
    if length(U) == size(U)[1]
        ncells = length(U)
        nconservative = 1
    else
        nconservative, ncells = size(U)
    end

    if scheme.reconstruct
        reconstruct!(UL, UR, U, scheme)
    else
        for i in 2:(ncells - 1)
            if nconservative == 1
                UL[right_edge(i)] = U[i]
                UR[left_edge(i)] = U[i]
            else
                for j in 1:nconservative
                    #println("j: ", j)
                    UL[j, right_edge(i)] = U[j, i]
                    UR[j, left_edge(i)] = U[j, i]
                end
            end
        end
    end

    if nconservative == 1
        UL[1] = U[1]
        UR[end] = U[end]
    else
        for j in 1:nconservative
            UL[j, 1] = U[j, 1]
            UR[j, end] = U[j, end]
        end
    end
    return UL, UR
end

function compute_fluxes!(F, UL, UR, fluids, fluid_ranges, scheme, pe, coupled)
    _, nedges = size(F)

    for i in 1:nedges
        for (j, (fluid, fluid_range)) in enumerate(zip(fluids, fluid_ranges))
            pe_L = pe[i]
            pe_R = pe[i+1]
            @views scheme.flux_function(F[fluid_range, i], UL[fluid_range, i],
                                        UR[fluid_range, i], fluid, pe_L, pe_R, coupled)
        end
    end
    return F
end