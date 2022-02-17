#should be able to use global params variables now
function electron_velocity(U, params, i)
    index = params.index
    if i == 1
        inew = 1
    else
        inew = i
    end
    @views grad_nϵ = first_deriv_central_diff(U[index.nϵ, :], params.z_cell, i)
    uₑ = params.cache.μ[i]*(U[index.grad_ϕ, inew] - grad_nϵ/U[index.ne, i])
    return uₑ
end

function e_heat_conductivity(params, i) #without nϵ term, that is in state vector
    ν = params.cache.νan[i] + params.cache.νc[i]
    ω = e*params.cache.B[i]/mₑ
    #return 4.7/(mₑ*(1/(ν))*ω^2) #Braginskii closure, should be Tev
    return 10/9*params.cache.μ[i]
    #return 10/9*10.0
end

function flux_electron!(F, US, fluid, params, U, i)
    index = params.index
    nϵ = US
    uₑ = U[index.ue, i]
    ϵ = U[index.Tev, i]
    @views grad_Tev = first_deriv_facereconstr_2order(U[index.nϵ, :], params.z_cell, i)
    κₑ = e_heat_conductivity(params, i)
    F = 5/3*nϵ*uₑ
    return F
end

function compute_fluxes_electron!(F, UL, UR, U, fluid, fluid_ranges, scheme, params)
    nedges = length(F)

    for i in 1:nedges
        @views F[i] = scheme.flux_function(F[i], UL[i],
                                        UR[i], fluid, params, U, i)
    end
    return F
end

function upwind_electron!(F, UL, UR, fluid, params, U, i)
    index = params.index
    uL = U[index.ue, right_edge(i+1)] 
    uR = U[index.ue, left_edge(i+1)]
    avg_velocity = 0.5 * (uL + uR)

    if avg_velocity ≥ 0
        F = flux_electron!(F, UL, fluid, params, U, i)
    else
        F = flux_electron!(F, UR, fluid, params, U, i)
    end
    return F
end

function HLLE_electron!(F, UL, UR, fluid, params, U, i)
    index = params.index
    uL = U[index.ue, right_edge(i+1)]
    uR = U[index.ue, left_edge(i+1)]
    
    aL = sqrt(5/3*1.526662276445587e7*U[index.Tev, right_edge(i+1)])
    aR = sqrt(5/3*1.526662276445587e7*U[index.Tev, left_edge(i+1)])

    sL_min, sL_max = min(0, uL - aL), max(0, uL + aL)
    sR_min, sR_max = min(0, uR - aR), max(0, uR + aR)

    smin = min(sL_min, sR_min)
    smax = max(sL_max, sR_max)

    FL = flux_electron!(F, UL, fluid, params, U, i)
    FR = flux_electron!(F, UR, fluid, params, U, i)

    F = 0.5 * (FL + FR) -
    0.5 * (smax + smin) / (smax - smin) * (FR - FL) +
    smax * smin / (smax - smin) * (UR - UL)

    return F
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
    neutral_density = U[1, i]/fluid.m
    W = params.landmark.loss_coeff(U[index.Tev, i])
    return neutral_density*U[index.ne, i]*W
end
