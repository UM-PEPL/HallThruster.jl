Base.@kwdef struct HyperbolicScheme{F, L}
    flux_function::F = rusanov
    limiter::SlopeLimiter{L} = van_leer
    reconstruct::Bool = true
end

@inline function flux(U::NTuple{1, T}, fluid) where T
    ρ = U[1]
    u = fluid.u
    return (ρ * u,)
end

@inline function flux(U::NTuple{2, T}, fluid) where T
    ρ, ρu = U
    p = pressure(U, fluid)
    return (ρu, U[2]^2 / U[1] + p)
end

@inline function flux(U::NTuple{3, T}, fluid) where T
    ρ, ρu, ρE = U
    u = ρu / ρ
    p = pressure(U, fluid)
    ρH = ρE + p
    return (ρu, U[2]^2 / U[1] + p, ρH * u)
end

# i forget why i did things this way. this seems super unnecessary, but we'll get to it
# oh yeah, it was because I removed StaticArrays
macro NTuple(ex)
    if !isa(ex, Expr)
        error("Bad input for @NTuple")
    end
    head = ex.head
    if head === :comprehension
        if length(ex.args) != 1 || !isa(ex.args[1], Expr) || ex.args[1].head != :generator
            error("Expected generator in comprehension, e.g. [f(i) for i = 1:3]")
        end
        ex = ex.args[1]
        if length(ex.args) != 2
            error("Use a one-dimensional comprehension for @$SV")
        end
        rng = eval(ex.args[2].args[2])
        exprs = (:(f($j)) for j in rng)
        return quote
            let
                f($(esc(ex.args[2].args[1]))) = $(esc(ex.args[1]))
                NTuple{$(length(rng))}($tuple($(exprs...)))
            end
        end
    elseif head === :typed_comprehension
        if length(ex.args) != 2 || !isa(ex.args[2], Expr) || ex.args[2].head != :generator
            error("Expected generator in typed comprehension, e.g. Float64[f(i) for i = 1:3]")
        end
        T = esc(ex.args[1])
        ex = ex.args[2]
        if length(ex.args) != 2
            error("Use a one-dimensional comprehension for @$SV")
        end
        rng = eval(ex.args[2].args[2])
        exprs = (:(f($j)) for j in rng)
        return quote
            let
                f($(esc(ex.args[2].args[1]))) = $(esc(ex.args[1]))
                NTuple{$(length(rng)),$T}($tuple($(exprs...)))
            end
        end
    else
        error("Expected comprehension")
    end
end

# use fun metaprogramming create specialized flux versions for each type of fluid
for NUM_CONSERVATIVE in 1:3
    eval(quote

    @inbounds @fastmath function rusanov(UL::NTuple{$NUM_CONSERVATIVE, T}, UR::NTuple{$NUM_CONSERVATIVE, T}, fluid, args...) where T
        γ = fluid.species.element.γ
        Z = fluid.species.Z

        uL = velocity(UL, fluid)
        uR = velocity(UR, fluid)
        TL = temperature(UL, fluid)
        TR = temperature(UR, fluid)
        aL = sound_speed(UL, fluid)
        aR = sound_speed(UL, fluid)

        sL_max = max(abs(uL - aL), abs(uL + aL))
        sR_max = max(abs(uR - aR), abs(uR + aR))

        smax = max(sL_max, sR_max)

        FL = flux(UL, fluid)
        FR = flux(UR, fluid)

        return @NTuple [0.5 * ((FL[j] + FR[j]) - smax * (UR[j] - UL[j])) for j in 1:$(NUM_CONSERVATIVE)]
    end

    @fastmath function HLLE(UL::NTuple{$NUM_CONSERVATIVE, T}, UR::NTuple{$NUM_CONSERVATIVE, T}, fluid, args...) where T
        γ = fluid.species.element.γ
        Z = fluid.species.Z

        uL = velocity(UL, fluid)
        uR = velocity(UR, fluid)
        TL = temperature(UL, fluid)
        TR = temperature(UR, fluid)
        aL = sound_speed(UL, fluid)
        aR = sound_speed(UL, fluid)

        sL_min, sL_max = min(0, uL - aL), max(0, uL + aL)
        sR_min, sR_max = min(0, uR - aR), max(0, uR + aR)

        smin = min(sL_min, sR_min)
        smax = max(sL_max, sR_max)

        FL = flux(UL, fluid)
        FR = flux(UR, fluid)

        return @NTuple[
            0.5 * (FL[j] + FR[j]) -
            0.5 * (smax + smin) / (smax - smin) * (FR[j] - FL[j]) +
            smax * smin / (smax - smin) * (UR[j] - UL[j])
            for j in 1:$(NUM_CONSERVATIVE)
        ]
    end

    @fastmath function global_lax_friedrichs(UL::NTuple{$NUM_CONSERVATIVE, T}, UR::NTuple{$NUM_CONSERVATIVE, T}, fluid, λ_global = 0.0, args...) where T
        γ = fluid.species.element.γ
        Z = fluid.species.Z

        uL = velocity(UL, fluid)
        uR = velocity(UR, fluid)

        FL = flux(UL, fluid)
        FR = flux(UR, fluid)

        return @NTuple[
            0.5 * (FL[j] + FR[j]) + 0.5 * λ_global * (UL[j] - UR[j])
            for j in 1:$(NUM_CONSERVATIVE)
        ]
    end

    end)
end

@inline function reconstruct(uⱼ₋₁, uⱼ, uⱼ₊₁, limiter)
    r = (uⱼ₊₁ - uⱼ) / (uⱼ - uⱼ₋₁)
    Δu = 0.25 * limiter(r) * (uⱼ₊₁ - uⱼ₋₁)
    return uⱼ - Δu, uⱼ + Δu
end

function compute_edge_states!(UL, UR, U, params; apply_boundary_conditions = false)
    (nvars,  ncells) = size(U)
    (;config, index, is_velocity_index) = params
    (;scheme) = config

    # compute left and right edge states
    if (scheme.reconstruct)
        @inbounds for j in 1:nvars
            if is_velocity_index[j] # reconstruct velocity as primitive variable instead of momentum density
                for i in 2:ncells-1
                    u₋ = U[j, i-1]/U[j-1, i-1]
                    uᵢ = U[j, i]/U[j-1, i]
                    u₊ = U[j, i+1]/U[j-1, i+1]
                    uR, uL = reconstruct(u₋, uᵢ, u₊, scheme.limiter)

                    ρL = UL[j-1, right_edge(i)] #use previously-reconstructed edge density to compute momentum
                    ρR = UR[j-1, left_edge(i)]
                    UL[j, right_edge(i)] = uL*ρL
                    UR[j, left_edge(i)] = uR*ρR
                end
            else
                for i in 2:ncells-1
                    u₋ = U[j, i-1]
                    uᵢ = U[j, i]
                    u₊ = U[j, i+1]

                    UR[j, left_edge(i)], UL[j, right_edge(i)] = reconstruct(u₋, uᵢ, u₊, scheme.limiter)
                end
            end
        end
    else
        @inbounds for i in 2:ncells-1, j in 1:nvars
            UL[j, right_edge(i)] = U[j, i]
            UR[j, left_edge(i)]  = U[j, i]
        end
    end

    if apply_boundary_conditions
        @views left_boundary_state!(UL[:, 1], U, params)
        @views right_boundary_state!(UR[:, end], U, params)
    else
        @. @views UL[:, 1] = U[:, 1]
        @. @views UR[:, end] = U[:, end]
    end
end

function compute_fluxes!(F, UL, UR, U, params; apply_boundary_conditions = false)
    (;config, index, fluids, Δz_edge, num_neutral_fluids, cache) = params
    (;λ_global) = cache
    (;propellant, scheme, ncharge) = config
    ncells = size(U, 2)

    mi = propellant.m
    nedges = ncells-1

    # Reconstruct the states at the left and right edges using MUSCL scheme
    compute_edge_states!(UL, UR, U, params; apply_boundary_conditions)

    # Compute maximum wave speed in domain and use this to update the max allowable timestep, if using adaptive timestepping
    @inbounds for i in 1:nedges
        # Compute wave speeds for each component of the state vector.
        # The only wave speed for neutrals is the neutral convection velocity
        for j in 1:num_neutral_fluids
            neutral_fluid = fluids[j]
            U_neutrals = (U[index.ρn[j], i],)
            u = velocity(U_neutrals, neutral_fluid)
            λ_global[j] = abs(u)
        end

        # Ion wave speeds
        for Z in 1:ncharge
            fluid_ind = Z + num_neutral_fluids
            fluid = fluids[fluid_ind]
            γ = fluid.species.element.γ
            UL_ions = (UL[index.ρi[Z], i], UL[index.ρiui[Z], i],)
            UR_ions = (UR[index.ρi[Z], i], UR[index.ρiui[Z], i],)

            uL = velocity(UL_ions, fluid)
            TL = temperature(UL_ions, fluid)
            uR = velocity(UR_ions, fluid)
            TR = temperature(UR_ions, fluid)
            aL = sound_speed(UL_ions, fluid)
            aR = sound_speed(UL_ions, fluid)

            # Maximum wave speed
            s_max = max(abs(uL + aL), abs(uL - aL), abs(uR + aR), abs(uR - aR))

            # a Δt / Δx = 1 for CFL condition, user-supplied CFL number restriction applied later, in update_values
            dt_max = Δz_edge[i] / s_max
            params.cache.dt_u[i] = dt_max

            # Update maximum wavespeeds and maximum allowable timestep
            λ_global[fluid_ind] = max(s_max, λ_global[fluid_ind])
        end
    end

    @inbounds for i in 1:nedges
        # Neutral fluxes at edge i
        for j in 1:params.num_neutral_fluids
            left_state_n  = (UL[index.ρn[j], i],)
            right_state_n = (UR[index.ρn[j], i],)

            F[index.ρn[j], i] = scheme.flux_function(left_state_n, right_state_n, fluids[j], λ_global[j])[1]
        end

        # Ion fluxes at edge i
        for Z in 1:ncharge
            left_state_i  = (UL[index.ρi[Z], i], UL[index.ρiui[Z], i],)
            right_state_i = (UR[index.ρi[Z], i], UR[index.ρiui[Z], i],)
            fluid_ind = Z + num_neutral_fluids
            F_mass, F_momentum = scheme.flux_function(left_state_i, right_state_i, fluids[fluid_ind], λ_global[fluid_ind])
            F[index.ρi[Z],   i] = F_mass
            F[index.ρiui[Z], i] = F_momentum
        end
    end

    return F
end
