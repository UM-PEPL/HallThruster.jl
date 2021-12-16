#should be able to use global params variables now
function electron_velocity(params, i)
    grad_ϕ = first_deriv_central_diff(params.cache.ϕ, params.z_cell, i)
    grad_pe = first_deriv_central_diff(params.cache.pe, params.z_cell, i)
    uₑ = -params.cache.μ[i]*(-grad_ϕ + grad_pe/e/params.cache.ne[i])
    return uₑ
end

function e_heat_conductivity(params, i)
    ν = params.cache.νan[i] + params.cache.νc[i]
    ω = e*params.cache.B[i]/mₑ
    return 4.7*params.cache.ne[i]*params.cache.Tev[i]*e/kB/(mₑ*(1/(ν))*ω^2)
end

function flux_electron(U, fluid, params, i)
    ϵ = U #3/2 ne kB Te
    grad_Te = first_deriv_central_diff(params.cache.Tev, params.z_cell, i)*e/kB#Te is inside ϵ, so could rewrite this
    uₑ = electron_velocity(params, i) #use Ohms law
    κₑ = e_heat_conductivity(params, i) #from the braginskii closure, need Te and cyclotron frequency, and tau e, which come from anomalous transport
    F = (ϵ*5/3*uₑ + κₑ*grad_Te, 0, 0)
    return F
end

#have to parallel updates for now, ϵ and Tev makes both messy

function flux_electron!(F, U, fluid, params, i)
    ϵ = U #3/2 ne kB Te
    grad_Te = first_deriv_central_diff(params.cache.Tev, params.z_cell, i)*e/kB#Te is inside ϵ, so could rewrite this
    uₑ = electron_velocity(params, i) #use Ohms law
    κₑ = e_heat_conductivity(params, i) #from the braginskii closure, need Te and cyclotron frequency, and tau e, which come from anomalous transport
    F = ϵ*5/3*uₑ + κₑ*grad_Te
    return F
end

function compute_fluxes_electron!(F, UL, UR, fluids, fluid_ranges, scheme, params)
    nedges = length(F)

    for i in 1:nedges
        for (j, (fluid, fluid_range)) in enumerate(zip(fluids, fluid_ranges))
            @views scheme.flux_function(F[i], UL[i],
                                        UR[i], fluid, params, i)
        end
    end
    return F
end

function upwind_electron!(F, UL, UR, fluid, params, i)
    #redo velocity calculation
    uL = electron_velocity(params, left_edge(i+1))
    uR = electron_velocity(params, right_edge(i+1)) 
    avg_velocity = 0.5 * (uL + uR)

    if avg_velocity ≥ 0
        flux_electron!(F, UL, fluid, params, i)
    else
        flux_electron!(F, UR, fluid, params, i)
    end
    return F
end


function source_electron_energy!(QE, U, params, i)
    uₑ = electron_velocity(params, i)
    grad_pe = first_deriv_central_diff(params.cache.pe, params.z_cell, i)
    ν = params.cache.νan[i] + params.cache.νc[i]
    QE = grad_pe*uₑ + mₑ*params.cache.ne[i]*ν*uₑ^2 + S_wall_simple(U[4, :], i) + S_coll(U, params, i) #resistive heating collisions, u has to be total u not just z, azimuthal component dominating
    return QE
end

"""
    first_deriv_central_diff(u::Vector{Float64}, z_cell::Vector{Float64}, i::Int64)

returns the first derivative second order central difference approximation at location i. 
if i == 1, returns right one sided second order approx, elseif i == length(array), 
returns left one sided second order approx. 
"""

function first_deriv_central_diff(u::Vector{Float64}, z_cell::Vector{Float64}, i::Int64) #central second order approx of first derivative
    if i == 1
        grad = (-3*u[i] + 4*u[i+1] - u[i+2])/(abs(z_cell[i]-z_cell[i+2])) #second order one sided for boundary
    elseif i == length(u)
        grad = (u[i-2] - 4*u[i-1] + 3*u[i])/(abs(z_cell[i-2]-z_cell[i])) #second order one sided for boundary
    else
        grad = (u[i+1] - u[i-1])/(abs(z_cell[i+1]-z_cell[i-1])) #centered difference
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

function S_wall_simple(E, i) #landmark and Hara non-oscillatory
    return 0.1*10e7*exp(-20/E[i])*E[i]
end

function S_coll(U, params, i) #landmark table
    fluid = params.fluids[1].species.element
    neutral_density = U[1, i]/fluid.m
    W = params.landmark.loss_coeff(params.cache.Tev[i])
    return neutral_density*params.cache.ne[i]*W
end
