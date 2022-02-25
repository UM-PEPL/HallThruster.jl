#should be able to use global params variables now
function electron_velocity(U, params, i)
    (;z_cell) = params
    (;μ, ∇ϕ, ne, pe) = params.cache

    ∇pe = first_deriv_central_diff(pe, z_cell, i)
    #∇pe = uneven_central_diff(pe[i-1], pe[i], pe[i+1], z_cell[i-1], z_cell[i], z_cell[i+1])
    uₑ = μ[i] * (∇ϕ[i] - ∇pe/ne[i])
    return uₑ
end

"""
    first_deriv_central_diff(u::Vector{Float64}, z_cell::Vector{Float64}, i::Int64)

returns the first derivative second order central difference approximation at location i. 
if i == 1, returns right one sided second order approx, elseif i == length(array), 
returns left one sided second order approx. 
"""

function first_deriv_central_diff(u, z_cell, i) #central second order approx of first derivative
    if i == 1
        grad = (-3*u[i] + 4*u[i+1] - u[i+2])/(abs(z_cell[i+1]-z_cell[i+3])) #second order one sided for boundary, or adapt for non constant stencil
    elseif i == length(u)
        grad = (u[i-2] - 4*u[i-1] + 3*u[i])/(abs(z_cell[i-3]-z_cell[i-1])) #second order one sided for boundary
    else
        grad = (u[i+1] - u[i-1])/(abs(z_cell[4]-z_cell[2])) #centered difference
    end

    return grad
end

function first_deriv_central_diff_pot(u, z_cell, i) #central difference of first deriv
    if i == 1
        grad = (u[2] - u[1])/(z_cell[3] - z_cell[2])
    elseif i == length(u)
        grad = (u[i-1] - u[i-2])/(z_cell[3] - z_cell[2])
    else 
        grad = (u[i] - u[i-1])/(z_cell[3] - z_cell[2])
    end
    return grad
end

function second_deriv_central_diff_energy(U, z_cell, params, i)
    index = params.index
    #do once with 1/2 e^2 and once with e only
    μ⁺ = 10.0 #*3/2*HallThruster.kB/HallThruster.e*(U[index.Tev, i] + U[index.Tev, i+1])/2
    μ⁻ = 10.0 #*3/2*HallThruster.kB/HallThruster.e*(U[index.Tev, i] + U[index.Tev, i-1])/2
    nϵ⁺ = 1.0
    nϵ⁻ = 1.0

    BC_left = 3/2*4/3*1e18*50000*HallThruster.kB
    BC_right = 3/2*4/3*1e18*50000*HallThruster.kB
    u = zeros(length(z_cell))
    u = U[index.Tev, :]*3/2*HallThruster.kB/HallThruster.e

    #μ⁺ = (params.cache.μ[i] + params.cache.μ[i+1])/2
    #μ⁻ = (params.cache.μ[i] + params.cache.μ[i-1])/2
    #nϵ⁺ = (u[i] + u[i+1])/2 
    #nϵ⁻ = (u[i] + u[i-1])/2
    #grad = 10/9/(z_cell[3]-z_cell[2])^2*(μ⁻*nϵ⁻*u[i-1] - (μ⁻*nϵ⁻ + μ⁺*nϵ⁺)*u[i] + μ⁺*nϵ⁺*u[i+1])
    #grad = 10/9*((μ⁻*nϵ⁻ * (-u[i] + u[i-1])/(z_cell[i]-z_cell[i-1])^2) + (μ⁺*nϵ⁺ * (-u[i] + u[i+1])/(z_cell[i+1]-z_cell[i])^2))
    #grad = 10/9*((μ⁻*nϵ⁻ * (-u[i] + u[i-1])/((z_cell[i]-z_cell[i-1])*((z_cell[i+1]-z_cell[i])))) + (μ⁺*nϵ⁺ * (-u[i] + u[i+1])/((z_cell[i+1]-z_cell[i])*(z_cell[i]-z_cell[i-1]))))
    if i == 2 
        grad = 10/9*((μ⁻*nϵ⁻ * (-u[i] + u[i-1])/((z_cell[i]-z_cell[i-1])*((z_cell[i]-z_cell[i-1])))) + (μ⁺*nϵ⁺ * (-u[i] + u[i+1])/((z_cell[i+1]-z_cell[i])*(z_cell[i]-z_cell[i-1]))))
    elseif i == length(z_cell)-1
        grad = 10/9*((μ⁻*nϵ⁻ * (-u[i] + u[i-1])/((z_cell[i]-z_cell[i-1])*((z_cell[i+1]-z_cell[i])))) + (μ⁺*nϵ⁺ * (-u[i] + u[i+1])/((z_cell[i+1]-z_cell[i])*(z_cell[i+1]-z_cell[i]))))
    else
        grad = 10/9*((μ⁻*nϵ⁻ * (-u[i] + u[i-1])/(z_cell[i]-z_cell[i-1])^2) + (μ⁺*nϵ⁺ * (-u[i] + u[i+1])/(z_cell[i+1]-z_cell[i])^2))
    end
    return grad
end

function second_deriv_central_diff_gen(u, z_cell, i)
    if i == 2 
        grad = (-u[i] + u[i-1])/((z_cell[i]-z_cell[i-1])*((z_cell[i]-z_cell[i-1]))) + (-u[i] + u[i+1])/((z_cell[i+1]-z_cell[i])*(z_cell[i]-z_cell[i-1]))
    elseif i == length(z_cell)-1
        grad = (-u[i] + u[i-1])/((z_cell[i]-z_cell[i-1])*((z_cell[i+1]-z_cell[i]))) + (-u[i] + u[i+1])/((z_cell[i+1]-z_cell[i])*(z_cell[i+1]-z_cell[i]))
    else
        grad = (-u[i] + u[i-1])/(z_cell[i]-z_cell[i-1])^2 + (-u[i] + u[i+1])/(z_cell[i+1]-z_cell[i])^2
    end
    return grad
end 


"""
    first_deriv_facereconstr_2order(u::Vector{Float64}, z_cell::Vector{Float64}, i::Int64)

returns the first derivative second order face reconstruction at i. stencil limits, l = -1/2,
r = 1/2
if i == 1, returns right one sided second order approx, elseif i == length(array), 
returns left one sided second order approx. 
"""

function first_deriv_facereconstr_2order(u, z_cell, i)
    grad = (u[i+1]-u[i])/abs(z_cell[i+1] - z_cell[i]) #centered difference
    #grad = (-u[i] + u[i-1])/abs(z_cell[i] - z_cell[i-1])
    return grad
end

"""
    S_wall(params)

wall heat loss. electron density multiplied by electron wall collision frequency and mean electron energy loss due to wall collision, 
which is assumed 2 Tev + sheath potential, using from Eq. ..., from Fundamentals of Electric Propulsion, Goebel and Katz, 2008.
"""

function S_wall_Bohm(params, i) #hara mikellides 2018
    σ = 1e-10 #electron collision area
    νₑ_w = 1 #needs to be added
    fluid = params.fluids[1].species.element
    Δϵ_w = 2*params.Tev[i] + params.Tev[i]*log(1-σ)/sqrt(2*pi*mₑ/fluid.m)
    return params.ne[i]*νₑ_w*Δϵ_w
end

function S_wall_simple(ϵ, i) #landmark and Hara non-oscillatory
    return 1e7*exp(-20/ϵ[i])*ϵ[i] #also anomalous energy loss, #different cases for ν\_ϵ    
end

function S_coll(U, params, i) #landmark table
    index = params.index
    fluid = params.fluids[1].species.element
    (; ne, Tev) = params.cache.Tev
    neutral_density = U[1, i]/fluid.m
    W = params.loss_coeff(Tev[i])
    return neutral_density*ne[i]*W
end

function energy_crank_nicholson!(U, params)
    (;Aϵ, bϵ, μ, ue, ne, Tev) = params.cache
    (;z_cell, dt, index) = params
    implicit = params.config.implicit_energy
    explicit = 1 - implicit
    ncells = length(z_cell)

    nϵ = @views U[index.nϵ, :]
    Aϵ.d[1] = 1.0
    Aϵ.du[1] = 0.0
    Aϵ.d[end] = 1.0
    Aϵ.dl[end-1] = 0.0

    bϵ[1] = params.Te_L * ne[1]
    bϵ[end] = params.Te_R * ne[end]

    for i in 2:ncells-1
        Q = source_electron_energy_landmark(U, params, i)

        # Upwinded first derivatives
        z0 = z_cell[i-1]
        z1 = z_cell[i]
        z2 = z_cell[i+1]
        Δz⁻ = z1 - z0
        Δz⁺ = z2 - z1

        ue⁺ = smooth_if(ue[i], 0.0, 0.0, ue[i], 0.01) / ue[i]
        ue⁻ = smooth_if(ue[i], 0.0, ue[i], 0.0, 0.01) / ue[i]

        dϵ_dz⁺ = diff(Tev[i-1], Tev[i], z0, z1)
        dϵ_dz⁻ = diff(Tev[i], Tev[i+1], z1, z2)

        d²ϵ_dz² = uneven_second_deriv(Tev[i-1], Tev[i], Tev[i+1], z0, z1, z2)

        advection_term = (
            ue⁺ * diff(ue[i-1]*nϵ[i-1], ue[i]*nϵ[i], z0, z1) +
            ue⁻ * diff(ue[i]*nϵ[i], ue[i+1]*nϵ[i+1], z1, z2)
        )

        diffusion_term = (
            ue⁺ * diff(μ[i-1]*nϵ[i-1], μ[i]*nϵ[i], z0, z1) * dϵ_dz⁺ +
            ue⁻ * diff(μ[i]*nϵ[i], μ[i+1]*nϵ[i+1], z1, z2) * dϵ_dz⁻
        )

        diffusion_term += μ[i] * nϵ[i] * d²ϵ_dz²

        F = 5/3 * advection_term - 10/9 * (diffusion_term)

        # Explicit part
        bϵ[i] = nϵ[i] + dt * (Q - explicit * F)

        # Implicit contribution of advection term
        Aϵ.d[i] = 5/3 * ue[i] * (ue⁺ / Δz⁻ - ue⁻ / Δz⁺)
        Aϵ.dl[i-1] = -5/3 * -ue[i-1] * ue⁺ / Δz⁻
        Aϵ.du[i] = 5/3 * ue[i+1] * ue⁻ / Δz⁺

        # Implicit contribution of first part of diffusion term
        Aϵ.d[i] -= 10/9 * μ[i] * (ue⁺ / Δz⁻ * dϵ_dz⁺ - ue⁻ / Δz⁺ * dϵ_dz⁻)
        Aϵ.dl[i-1] -= 10/9 * μ[i-1] * (-ue⁺ / Δz⁻ * dϵ_dz⁺)
        Aϵ.du[i] -= 10/9 * μ[i+1] * (ue⁻ / Δz⁺ * dϵ_dz⁻)

        # Implicit contribution of second part of diffusion term
        Aϵ.d[i] -= 10/9 * μ[i] * d²ϵ_dz²

        # Time stepping component
        Aϵ.d[i] = 1.0 + implicit * dt * Aϵ.d[i]
        Aϵ.dl[i-1] *= implicit * dt
        Aϵ.du[i] *= implicit * dt
    end

    tridiagonal_solve!(nϵ, Aϵ, bϵ)
end