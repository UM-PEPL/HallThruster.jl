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


function update_electron_energy!(dU, U, params, t)
    #########################################################
    #ELECTRON SOLVE

    ncells = size(U, 2) - 2

    (;index, z_cell, z_edge) = params
    (;Tev, ue, ne, B, μ) = params.cache
    nϵ = @views U[index.nϵ, :]
    implicit_energy = params.config.implicit_energy

    mi = params.propellant.m

    neL = 0.0
    ne0 = sum(U[index.ρi[Z], 1] for Z in 1:params.config.ncharge) / mi
    neR = sum(U[index.ρi[Z], 2] for Z in 1:params.config.ncharge) / mi

    ϵL = 0.0
    ϵ0 = nϵ[1] / ne0
    ϵR = nϵ[2] / neR

    νc_L = 0.0
    νc_0 = electron_collision_freq(ϵ0, U[index.ρn, 1] / mi, ne0, mi)
    νc_R = electron_collision_freq(ϵR, U[index.ρn, 2] / mi, neR, mi)

    νan_L = 0.0
    νan_0 = params.anom_model(U, params, 1)
    νan_R = params.anom_model(U, params, 2)

    μL = 0.0
    μ0 = electron_mobility(νan_0, νc_0, params.cache.B[1])
    μR = electron_mobility(νan_R, νc_R, params.cache.B[2])

    ueL = 0.0
    ue0 = electron_velocity(U, params, 1)
    ueR = electron_velocity(U, params, 2)

    @inbounds for i in 2:ncells+1

        zL, z0, zR = z_cell[i-1], z_cell[i], z_cell[i+1]
        neL, ne0, neR = ne0, neR, sum(U[index.ρi[Z], i+1] for Z in 1:params.config.ncharge) / mi
        ϵL, ϵ0, ϵR = ϵ0, ϵR, nϵ[i+1] / neR

        νc_L, νc_0, νc_R = νc_0, νc_R, electron_collision_freq(ϵR, U[index.ρn, i+1] / mi, neR, mi)
        νan_L, νan_0, νan_R = νan_0, νan_R, params.anom_model(U, params, i+1)
        μL, μ0, μR = μ0, μR, electron_mobility(νan_R, νc_R, B[i+1])
        ueL, ue0, ueR = ue0, ueR, electron_velocity(U, params, i+1)

        advection_term = uneven_central_diff(ue[i-1] * nϵ[i-1], ue[i] * nϵ[i], ue[i+1] * nϵ[i+1], zL, z0, zR)
        diffusion_term = uneven_central_diff(μL * nϵ[i-1], μ0 * nϵ[i], μR * nϵ[i+1], zL, z0, zR)

        d²ϵ_dz² = uneven_second_deriv(ϵL, ϵ0, ϵR, zL, z0, zR)

        diffusion_term += μ0 * nϵ[i] * d²ϵ_dz²

        source_term = source_electron_energy_landmark(U, params, i)

        dU[index.nϵ, i] = - 5/3 * advection_term + 10/9 * diffusion_term + source_term
    end

    return nothing
end
