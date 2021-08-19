function flux(U, fluid)
    F = similar(U)
    flux!(F, U, fluid)
    return F
end

function flux!(F, U, fluid)
    if fluid.conservation_laws.type == :ContinuityOnly
        n = U[1]
        u = velocity(U, fluid)
        F[1] = n * u
    elseif fluid.conservation_laws.type == :IsothermalEuler
        n, nu = U
        u = velocity(U, fluid)
        p = pressure(U, fluid)
        F[1] = nu
        F[2] = nu * u + p / m(fluid)
    elseif fluid.conservation_laws.type == :EulerEquations
        n, nu, nE = U
        u = velocity(U, fluid)
        p = pressure(U, fluid)
        nH = nE + p/m(fluid)

        F[1] = nu
        F[2] = nu * u + p/m(fluid)
        F[3] = nH * u
    end
    return F
end

function HLLE!(F, UL, UR, fluid)
    aL = sound_speed(UL, fluid)
    aR = sound_speed(UR, fluid)

    uL = velocity(UL, fluid)
    uR = velocity(UR, fluid)

    sL_min, sL_max = min(0, uL-aL), max(0, uL+aL)
    sR_min, sR_max = min(0, uR-aR), max(0, uR+aR)

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

function reconstruct!(UL, UR, U, scheme)
	nconservative, ncells = size(U)
	Ψ = scheme.limiter

    # compute left and right edge states
    for i in 2:ncells-1
        for j in 1:nconservative
            if scheme.reconstruct
                u₋ = U[j, i-1]
                uᵢ = U[j, i]
                u₊ = U[j, i+1]
                Δu = u₊ - uᵢ
                ∇u = uᵢ - u₋
                r = Δu / ∇u
                UL[j, right_edge(i)] = uᵢ + 0.5 * Ψ(r) * ∇u
                UR[j, left_edge(i)]  = uᵢ - 0.5 * Ψ(1/r) * Δu
            else
                UL[j, right_edge(i)] = U[j, i]
                UR[j, left_edge(i)]  = U[j, i]
            end
        end
    end
	@. @views UL[:, 1] = U[:, 1]
    @. @views UR[:, end] = U[:, end]

	return UL, UR
end

function compute_fluxes!(F, UL, UR, fluids, fluid_ranges, scheme)
	_, nedges = size(F)

    for i in 1:nedges
		for (fluid, fluid_range) in zip(fluids, fluid_ranges)
			@views scheme.flux_function(
                F[fluid_range, i],
				UL[fluid_range, i],
				UR[fluid_range, i],
				fluid
			)
		end
    end
    return F
end