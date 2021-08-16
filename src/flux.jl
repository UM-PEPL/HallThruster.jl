function flux(U, fluid::ContinuityFluid)
    n = U[1]
    u = velocity(U, fluid)
    return SA[n * u]
end

function flux(U, fluid::IsothermalFluid)
    n, nu = U
    u = velocity(U, fluid)
    p = pressure(U, fluid)
    return SA[nu, nu * u + p/m(fluid)]
end

function flux(U, fluid::EulerFluid)
    n, nu, nE = U
    u = velocity(U, fluid)
    p = pressure(U, fluid)
    nH = nE + p/m(fluid)

    return SA[nu, nu * u + p/m(fluid), nH * u]
end

# This loop automatically generates specialized versions of the flux functions
# for each fluid type, so that the correct size of array is returned
for fluid_type in (ContinuityFluid, IsothermalFluid, EulerFluid)
    @eval function HLLE(UL, UR, fluid::$fluid_type)
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

        f = @SVector[
            0.5 * (FL[i] + FR[i]) -
            0.5 * (smax + smin) / (smax - smin) * (FR[i] - FL[i]) + 
            smax * smin / (smax - smin) * (UR[i] - UL[i])
            for i in 1:nvars($fluid_type)
        ]
    end
end

function upwind(UL, UR, fluid::Fluid)
    uL = velocity(UL, fluid)
    uR = velocity(UR, fluid)
    avg_velocity = 0.5 * (uL + uR)

    if avg_velocity ≥ 0
        return flux(UL, fluid)
    else
        return flux(UR, fluid)
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