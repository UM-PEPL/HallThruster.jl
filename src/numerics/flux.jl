Base.@kwdef struct HyperbolicScheme{F,L}
    flux_function::F
    limiter::L
    reconstruct::Bool
    WENO::Bool
end

function flux(U::SVector{1, T}, fluid, pe = 0.0) where T
    ρ = U[1]
    u = velocity(U, fluid)
    return SA[ρ * u]
end

function flux(U::SVector{2, T}, fluid, pe = 0.0) where T
    ρ, ρu = U
    u = velocity(U, fluid)
    p = pressure(U, fluid)

    return SA[ρu, ρu * u + p + pe]
end

function flux(U::SVector{3, T}, fluid, pe = 0.0) where T
    ρ, ρu, ρE = U
    u = U[2] / U[1]
    p = pressure(U, fluid)
    ρH = ρE + p
    return SA[ρu, ρu * u + p + pe, ρH * u]
end

# use fun metaprogramming create specialized flux versions for each type of fluid
for NUM_CONSERVATIVE in 1:3
    eval(quote

    function upwind(UL::SVector{$NUM_CONSERVATIVE, T}, UR::SVector{$NUM_CONSERVATIVE, T}, fluid::Fluid, coupled = false, TeL = 0.0, TeR = 0.0, neL = 1.0, neR = 1.0, args...) where T
        uL = velocity(UL, fluid)
        uR = velocity(UR, fluid)

        peL = neL * TeL
        peR = neR * TeR

        avg_velocity = 0.5 * (uL + uR)

        Z = fluid.species.Z
        if avg_velocity ≥ 0
            return flux(UL, fluid, Z * coupled * e * peL)
        else
            return flux(UR, fluid, Z * coupled * e * peR)
        end
    end

    function rusanov(UL::SVector{$NUM_CONSERVATIVE, T}, UR::SVector{$NUM_CONSERVATIVE, T}, fluid, coupled = false, TeL = 0.0, TeR = 0.0, neL = 1.0, neR = 1.0, args...) where T
        γ = fluid.species.element.γ
        Z = fluid.species.Z

        uL = velocity(UL, fluid)
        uR = velocity(UR, fluid)

        TL = temperature(UL, fluid)
        TR = temperature(UR, fluid)

        mi = m(fluid)

        peL = neL * TeL
        peR = neR * TeR

        charge_factor = Z * e * coupled

         aL = sqrt((charge_factor * TeL + γ * kB * TL) / mi)
         aR = sqrt((charge_factor * TeR + γ * kB * TR) / mi)

        sL_max = max(abs(uL - aL), abs(uL + aL), #=abs(uL)=#)
        sR_max = max(abs(uR - aR), abs(uR + aR), #=abs(uR)=#)

        smax = max(sL_max, sR_max)

        FL = flux(UL, fluid, charge_factor * peL)
        FR = flux(UR, fluid, charge_factor * peR)

        return @SVector [0.5 * (FL[j] + FR[j]) - 0.5 * smax * (UR[j] - UL[j]) for j in 1:$(NUM_CONSERVATIVE)]
    end

    function HLLE(UL::SVector{$NUM_CONSERVATIVE, T}, UR::SVector{$NUM_CONSERVATIVE, T}, fluid, coupled = false, TeL = 0.0, TeR = 0.0, neL = 1.0, neR = 1.0, args...) where T
        γ = fluid.species.element.γ
        Z = fluid.species.Z

        uL = velocity(UL, fluid)
        uR = velocity(UR, fluid)

        peL = TeL * neL
        peR = TeR * neR

        TL = temperature(UL, fluid)
        TR = temperature(UR, fluid)

        mi = m(fluid)

        charge_factor = Z * e * coupled

        aL = sqrt((charge_factor * TeL + γ * kB * TL) / mi)
        aR = sqrt((charge_factor * TeR + γ * kB * TR) / mi)

        sL_min, sL_max = min(0, uL - aL), max(0, uL + aL)
        sR_min, sR_max = min(0, uR - aR), max(0, uR + aR)

        smin = min(sL_min, sR_min)
        smax = max(sL_max, sR_max)

        FL = flux(UL, fluid, charge_factor * peL)
        FR = flux(UR, fluid, charge_factor * peR)

        return @SVector[
            0.5 * (FL[j] + FR[j]) -
            0.5 * (smax + smin) / (smax - smin) * (FR[j] - FL[j]) +
            smax * smin / (smax - smin) * (UR[j] - UL[j])
            for j in 1:$(NUM_CONSERVATIVE)
        ]
    end

    function global_lax_friedrichs(UL::SVector{$NUM_CONSERVATIVE, T}, UR::SVector{$NUM_CONSERVATIVE, T}, fluid, coupled = false, TeL = 0.0, TeR = 0.0, neL = 1.0, neR = 1.0, λ_global = 0.0, args...) where T
        γ = fluid.species.element.γ
        Z = fluid.species.Z

        uL = velocity(UL, fluid)
        uR = velocity(UR, fluid)

        peL = TeL * neL
        peR = TeR * neR

        charge_factor = Z * e * coupled

        FL = flux(UL, fluid, charge_factor * peL)
        FR = flux(UR, fluid, charge_factor * peR)

        return @SVector[
            0.5 * (FL[j] + FR[j]) + 0.5*λ_global*UL[j] - 0.5*λ_global*UR[j]
            for j in 1:$(NUM_CONSERVATIVE)
        ]
    end

    end)
end

 #input is flux from HLLE or rusanov as first order, only on a large stencil. if not enough points available, reduce to first order
 function WENO5_compute_fluxes(f₋₂, f₋₁, f₀, f₊₁, f₊₂)
    f_hat¹ = 1/3*f₋₂ - 7/6*f₋₁ + 11/6*f₀
    f_hat² = -1/6*f₋₁ + 5/6*f₀ + 1/3*f₊₁
    f_hat³ = 1/3*f₀ + 5/6*f₊₁ - 1/6*f₊₂

    γ₁ = 1/10
    γ₂ = 3/5
    γ₃ = 3/10

    β₁ = 12/13*(f₋₂ - 2*f₋₁ + f₀).^2 + 1/4*(f₋₂ - 4*f₋₁ + 3*f₀).^2
    β₂ = 12/13*(f₋₁ - 2*f₀ + f₊₁).^2 + 1/4*(f₋₁ - f₊₁).^2
    β₃ = 12/13*(f₀ - 2*f₊₁ + f₊₂).^2 + 1/4*(3*f₀ - 4*f₊₁ + f₊₂).^2

    ϵₖ = 1e-6

    w_tilde₁ = γ₁/(ϵₖ .+ β₁).^2
    w_tilde₂ = γ₂/(ϵₖ .+ β₂).^2
    w_tilde₃ = γ₃/(ϵₖ .+ β₃).^2

    w₁ = w_tilde₁ / (w_tilde₁ + w_tilde₂ + w_tilde₃)
    w₂ = w_tilde₂ / (w_tilde₁ + w_tilde₂ + w_tilde₃)
    w₃ = w_tilde₃ / (w_tilde₁ + w_tilde₂ + w_tilde₃)

    f_hat = w₁*f_hat¹ + w₂*f_hat² + w₃*f_hat³

    return f_hat
end

function compute_edge_states!(UL, UR, U, scheme)
    (nvars,  ncells) = size(U)
    Ψ = scheme.limiter
    # compute left and right edge states
    @inbounds for i in 2:ncells-1
        #=if scheme.reconstruct
            u₋ = U[1, i-1]
            uᵢ = U[1, i]
            u₊ = U[1, i+1]
            Δu = u₊ - uᵢ
            ∇u = uᵢ - u₋
            r = Δu / ∇u

             UL[1, right_edge(i)] = uᵢ + 0.5 * Ψ(r) * ∇u
             UR[1, left_edge(i)]  = uᵢ - 0.5 * Ψ(1/r) * Δu
        else
            UL[1, right_edge(i)] = U[1, i]
            UR[1, left_edge(i)]  = U[1, i]
        end
        for j in 2:nvars
            UL[j, right_edge(i)] = U[j, i]
            UR[j, left_edge(i)]  = U[j, i]
        end=#
        for j in 1:nvars
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

    @. @views UL[:, 1]   = U[:, 1]
    @. @views UR[:, end] = U[:, end]
end

function compute_fluxes!(F, UL, UR, U, params)
    (;config, index, fluids) = params
    (;propellant, electron_pressure_coupled, scheme) = config
    nvars, ncells = size(U)

    coupled = electron_pressure_coupled

    ncharge = config.ncharge

    mi = propellant.m

    nedges = ncells-1

    compute_edge_states!(UL, UR, U, scheme)
    λ_global = zeros(ncharge+1) #max global wavespeed


    @inbounds for i in 1:nedges
        # Compute number density
        ne = 0.0
        for Z in 1:ncharge
            ni = U[index.ρi[Z], i] / mi
            ne = ne + Z * ni
        end

        # Compute electron temperature
        Te = max(params.config.min_electron_temperature, U[index.nϵ, i] / ne)

        fluid = fluids[1]
        γ = fluid.species.element.γ
        U_neutrals = SA[U[index.ρn, i]]
        u = velocity(U_neutrals, fluid)
        T = temperature(U_neutrals, fluid)
        λ_global[1] = sqrt((γ * kB * T) / mi)

        for Z in 1:ncharge

            fluid = fluids[Z+1]
            γ = fluid.species.element.γ
            U_ions = SA[U[index.ρi[Z], i], U[index.ρiui[Z], i]]

            u = velocity(U_ions, fluid)

            pe = Te* ne

            T = temperature(U_ions, fluid)
            charge_factor = Z * e * coupled

            a = sqrt((charge_factor * Te + γ * kB * T) / mi)
            s_max = max(abs(u - a), abs(u + a))
            if s_max > λ_global[ncharge+1]
                λ_global[ncharge+1] = s_max
            end
        end
    end

    @inbounds for i in 1:nedges

        # Compute number density
        neL = 0.0
        neR = 0.0
        for Z in 1:ncharge
            ni_L = UL[index.ρi[Z], i] / mi
            ni_R = UR[index.ρi[Z], i] / mi
             neL = neL + Z * ni_L
             neR = neR + Z * ni_R
        end

        # Compute electron temperature
        ϵL = max(params.config.min_electron_temperature, UL[index.nϵ, i] / neL)
        ϵR = max(params.config.min_electron_temperature, UR[index.nϵ, i] / neR)

        # Neutral flux at edge i
        left_state_n  = SA[UL[index.ρn, i]]
        right_state_n = SA[UR[index.ρn, i]]

        F[index.ρn, i] = scheme.flux_function(left_state_n, right_state_n, fluids[1], coupled, ϵL, ϵR, neL, neR, λ_global[1])[1]

        # Ion fluxes at edge i
        for Z in 1:ncharge
            left_state_i  = SA[UL[index.ρi[Z], i], UL[index.ρiui[Z], i]]
            right_state_i = SA[UR[index.ρi[Z], i], UR[index.ρiui[Z], i]]
            F_mass, F_momentum = scheme.flux_function(left_state_i, right_state_i, fluids[Z+1], coupled, ϵL, ϵR, neL, neR, λ_global[Z+1])
            F[index.ρi[Z],   i] = F_mass
            F[index.ρiui[Z], i] = F_momentum
        end
    end

    return F
end
