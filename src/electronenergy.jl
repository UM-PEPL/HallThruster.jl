#should be able to use global params variables now
function electron_velocity(U, params, i)
    index = params.index
    (;μ, ∇ϕ, ne) = params.cache
    if i == 1
        inew = 1
    else
        inew = i
    end
    @views grad_nϵ = first_deriv_central_diff(U[index.nϵ, :], params.z_cell, i)
    uₑ = μ[i]*(∇ϕ[inew]- grad_nϵ/ne[i])
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
